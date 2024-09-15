require(MASS)
require(splines)

ols.function <- function(X, y, vcov = FALSE) {
  res <- list()
  if(vcov) {
    res$vcov <- solve(crossprod(X)) 
    res$coeff <- res$vcov%*%crossprod(X,y)
  } else {
    res$coeff <- try(solve(crossprod(X), crossprod(X,y)), silent = TRUE) 
  }
  res 
}

bddp <- function(y, X, prior, mcmc, standardise = TRUE) {
    multinom <- function(prob) {
      probs <- t(apply(prob,1,cumsum)) 
      res <- rowSums(probs - runif(nrow(probs)) < 0) + 1 
      return(res)  
    }
    
    yt <- y
    if(standardise == TRUE) {
      yt <- (y-mean(y))/sd(y)
    }
    n <- length(y)
    k <- ncol(X)
    
    m <- prior$m0
    S <- prior$S0
    nu <- prior$nu
    psi <- prior$Psi
    a <- prior$a
    b <- prior$b
    alpha <- prior$alpha
    L <- prior$L
    
    nburn <- mcmc$nburn
    nsave <- mcmc$nsave
    nskip <- mcmc$nskip
    nsim <- nburn + nsave*nskip
    
    p <- ns <- rep(0, L)
    v <- rep(1/L,L)
    v[L] <- 1
    
    z <- matrix(NA_real_, nrow = nsim, ncol = n, dimnames = list(1:nsim, 1:n))
    z_tmp <- vector(length = n)
    
    z[1,] <- rep(1,n)
    
    beta <- matrix(0, nrow = L, ncol = k)
    aux <- ols.function(X, yt)$coeff
    if(!inherits(aux, "try-error")) {
      for(l in 1:L) {
        beta[l,] <- aux
      }
    }
    
    tau <- rep(1/var(yt),L)
    prop <- prob <- matrix(NA_real_, nrow = n, ncol = L)
    
    P <- Tau <- matrix(NA_real_, nrow = nsim, ncol = L, dimnames = list(1:nsim, 1:L))
    Beta <- array(NA_real_,c(nsim,L,k), dimnames = list(1:nsim, 1:L, colnames(X)))
    
    Beta[1,,] <- beta
    Tau[1,] <- tau
    
    mu <- mvrnorm(1, mu = m, Sigma = S)
    Sigmainv <- rWishart(1, df = nu, solve(nu*psi))[,,1] 
    
    for(i in 2:nsim) {
      cumv <- cumprod(1-v)
      p[1] <- v[1]
      p[2:L] <- v[2:L]*cumv[1:(L-1)]
      
      for(l in 1:L) {
        prop[,l] <- p[l]*dnorm(yt, mean = X%*%beta[l,], sd = sqrt(1/tau[l]))
      }
      
      prob <- prop/rowSums(prop)
      
      z_tmp <- multinom(prob)
      
      ns <- sapply(1:L, function(x, v) sum(v == x), v = z_tmp)
      
     v[1:(L-1)] <- rbeta(L-1, 1 + ns[1:(L-1)], alpha + rev(cumsum(rev(ns[-1]))))
      
      Sigmainv_mu <- Sigmainv%*%mu
      
      for(l in 1:L) {
        X_l  <- matrix(X[z_tmp == l, ], ncol = k, nrow = ns[l])
        yt_l <- yt[z_tmp == l]
        V <- solve(Sigmainv + tau[l]*crossprod(X_l))
        mu1 <- V%*%(Sigmainv_mu + tau[l]*crossprod(X_l, yt_l))
        beta[l,] <- mvrnorm(1, mu = mu1, Sigma = V)
        
        aux <- yt_l - X_l%*%beta[l,]
        tau[l] <- rgamma(1, shape = a + (ns[l]/2), rate = b + 0.5*crossprod(aux))
      }
      
      S_inv <- solve(S)
      Vaux <- solve(S_inv + L*Sigmainv)
      if(k == 1) {
        meanmu <- Vaux%*%(S_inv%*%m + Sigmainv%*%sum(beta))
      } else {
        meanmu <- Vaux%*%(S_inv%*%m + Sigmainv%*%colSums(beta))
      }
      mu <- mvrnorm(1, mu = meanmu, Sigma = Vaux)
      
      Vaux1 <- 0
      for(l in 1:L) {
        Vaux1 <- Vaux1 + tcrossprod(beta[l,] - mu)
      }
      
      Sigmainv <- rWishart(1, nu + L, solve(nu*psi + Vaux1))[,,1]
      
      P[i,] <- p
      z[i,] <- z_tmp
      Beta[i,,] <- beta
      Tau[i,] <- tau
    }
    
    if (standardise == TRUE) {
      Beta[,,1] <- sd(y)*Beta[,,1] + mean(y)
      if(k > 1) {
        Beta[,,2:k] <- sd(y)*Beta[,,2:k]
      }
      Sigma2 <- var(y)*(1/Tau)
    } else {
      Sigma2 <- 1/Tau
    }
    
    res <- list()
    #latent component indicator
    res$z <- z[seq(nburn+1, nsim, by = nskip),]
    #weights
    res$P <- P[seq(nburn+1, nsim, by = nskip),]
    #regression coefficients
    res$Beta <- Beta[seq(nburn+1, nsim, by = nskip),,,drop = FALSE]
    #Normal components variances
    res$Sigma2 <- Sigma2[seq(nburn+1, nsim, by = nskip),]
    res
  }

inf_criteria <-
  function(y, X, res){
    n <- length(y)
    
    if(is.null(res$P)){
      L <- 1
    } else{
      L <- ncol(res$Beta)
    }
    
    if(L > 1){
      p <- res$P
    }
    if(L == 1){
      p <- NULL
    }
    
    beta <- res$Beta
    sigma2 <- res$Sigma2
    niter <- nrow(beta)
    
    if(L == 1){
      term <- matrix(0, nrow = niter, ncol = n)
      for(k in 1:niter) {
        term[k,] <- dnorm(y, mean = X%*%beta[k,], sd = sqrt(sigma2[k]))
      }        
    }
    
    if(L > 1){
      term_1 <- array(0, c(niter, L, n))
      term <- matrix(0, nrow = niter, ncol = n)
      
      for(i in 1:n) {
        for(l in 1:L) {
          term_1[,l,i] <- p[,l]*dnorm(y[i], mean = c(X[i,]%*%t(beta[,l,])), sd = sqrt(sigma2[,l]))
        }
        term[,i] <- apply(term_1[,,i], 1, function(x) sum(x))
      }
    }
    
    term
  }


waicnp <-
  function(y, X, res, L, termsum = NULL) {
    n <- length(y)  
    p <- res$P
    beta <- res$Beta
    sigma2 <- res$Sigma2
    niter <- nrow(p)
    
    if(is.null(termsum)) {
      termsum <- inf_criteria(y, X, res)
    }
    logtermsum <- log(termsum)
    
    lpd <- sum(log(apply(exp(logtermsum),2,mean)))
    p2 <- sum(apply(logtermsum,2,var))
    waic <- -2*(lpd-p2)
    
    res <- list()
    res$pW <- p2
    res$WAIC <- waic
    res
  }

lpml <- function(y, X, res, L, termsum = NULL) {
    n <- length(y)
    p <- res$P
    beta <- res$Beta
    sigma2 <- res$Sigma2
    niter <- nrow(p)
    
    if(is.null(termsum)) {
      termsum <- inf_criteria(y, X, res)
    }
    
    aux <- 1/termsum
    omegabari <- apply(aux, 2, mean)
    omegabari_1 <- sqrt(niter) * omegabari
    omegatilde <- matrix(0, nrow = niter, ncol = n)
    
    for(i in 1:n) {
      omegatilde[,i] <- pmin(aux[,i], omegabari_1[i])  
    }
    
    sum_omegatilde <- apply(omegatilde,2,sum)
    sum_term_omegatilde <- apply(termsum*omegatilde, 2, sum)
    cpo <- sum_term_omegatilde/sum_omegatilde
    
    lpml <- sum(log(cpo))
    
    res <- list()
    res$cpo <- cpo
    res$lpml <- lpml
    res 
  }


