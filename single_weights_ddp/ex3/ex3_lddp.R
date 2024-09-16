require(ggplot2)
require(mcclust.ext)
require(nor1mix)

setwd("/Users/vandainacio/Dropbox/BNPregression/Code/single_weights_DDP")
source("main_function_mcmc.R")

# Generate data
set.seed(1010101)
n <- 600
p <- 2
x <- matrix(0, n, p)

x[,1] <- runif(n, -2, 2)
x[,2] <- runif(n, -2, 2)

ind <- sin(x[, 1]*x[, 2]*pi/2) <= 0

m_true <- (-1)*(1 - ind) + ind
y <- m_true +rnorm(n, 0, 0.1)

# Prediction
xgrid <- seq(-2, 2, .1)
x_new <- cbind(rep(xgrid, length(xgrid)), rep(xgrid, each = length(xgrid)))
n_new <- dim(x_new)[1]

ind_new <- sin(x_new[, 1]*x_new[, 2]*pi/2) <= 0
m_true_new <- (-1)*(1 - ind_new) + ind_new

y_grid <- seq(-1.4, 1.4, .02)
m2 <- length(y_grid)
f_true_new <- dnorm(matrix(y_grid, nrow = m2, ncol = n_new), t(matrix(m_true_new, nrow = n_new, ncol = m2)), 0.1)

# LDDP model
# Design and Prediction matrices
X <- cbind(rep(1, n),
           x[, 1],
           x[, 2]
           )
k <- ncol(X)

X_predict <-cbind(rep(1, n_new),
                  x_new[, 1],
                  x_new[, 2]
                  )

# Data-driven prior
fit_lm <- lm(y ~ X - 1)
prior <- list(m0 = fit_lm$coefficients,
              S0 = solve(t(X)%*%X)*(sigma(fit_lm)^2), 
              nu = k + 2, 
              Psi = 30*(solve(t(X)%*%X)*(sigma(fit_lm)^2)),
              a = 2,
              b = (sigma(fit_lm)^2)/2,
              alpha = 1, 
              L = 20)

# MCMC configuration
mcmc <- list(nburn = 5000, nsave = 50000, nskip = 1)

# Run MCMC
set.seed(123)
fit_lddp <- bddp(y = y,
                 X = X,
                 prior = prior, 
                 mcmc = mcmc,
                 standardise = FALSE)

# Estimated clustering and plot
psm_lddp <- comp.psm(fit_lddp$z)
output_vi_lddp <- minVI(psm_lddp, fit_lddp$z)

ggplot() +
  geom_point(aes(x = x[,1], y = x[,2], color = as.factor(output_vi_lddp$cl))) +
  theme_bw() +
  labs( x = "x_1", y = "x_2", color = "Cluster")

# Predictive regression function and plot
m_pred_lddp <- matrix(0, nrow = n_new, ncol = mcmc$nsave)
for(l in 1:mcmc$nsave){
  for(j in 1:n_new){  
    m_pred_lddp[j, l] <- sum(fit_lddp$P[l,]*X_predict[j,]%*%t(fit_lddp$Beta[l,,]))  
  }
}

m_pred_lddp_m <- apply(m_pred_lddp, 1, mean)
m_pred_lddp_l <- apply(m_pred_lddp, 1, quantile, prob = 0.025)
m_pred_lddp_h <- apply(m_pred_lddp, 1, quantile, prob = 0.975)

df <- data.frame(x1 = x_new[,1], x2 = x_new[,2], y = m_pred_lddp_m)
ggplot(df) +
  geom_raster(aes(x1, x2, fill =  y),interpolate = TRUE) +
  theme_classic()

xnew2 <- 1.5
ggplot() +
  geom_point(aes(x = x_new[x_new[, 2] == xnew2, 1], y = m_pred_lddp_m[x_new[, 2] == xnew2]), col = "black") +
  geom_point(aes(x = x_new[x_new[, 2] == xnew2, 1], y = m_true_new[x_new[, 2] == xnew2]), col = "red") +
  geom_ribbon(aes(x = x_new[x_new[, 2] == xnew2, 1],ymin=m_pred_lddp_l[x_new[, 2] == xnew2], 
                  ymax=m_pred_lddp_h[x_new[, 2] == xnew2]), alpha = 0.2) +
  theme_bw() +
  labs( x = "x_1", y = "y") 
ylab("y")

# Empirical l2 prediction error
l2_err_lddp <- sum(((m_true_new - m_pred_lddp_m)^2)/n_new)^.5
l2_err_lddp

# Empirical coverage for prediction
ec_pred_lddp <- sum((m_true_new >= m_pred_lddp_l)&(m_true_new <= m_pred_lddp_h))/n_new
ec_pred_lddp

# Credible interval length for prediction
ci_length_lddp <- m_pred_lddp_h - m_pred_lddp_l
mean(ci_length_lddp)

# Predictive density function and plot
index <- seq(1, mcmc$nsave, by = 25)
dpred_lddp <- array(0,c(length(index), m2, n_new))

for(k in 1:length(index)){
  mu <- tcrossprod(X_predict, fit_lddp$Beta[(k*25) - 24, , ])
    
  for(l in 1:n_new) {
    aux <- norMix(mu = c(mu[l, ]), sigma = sqrt(fit_lddp$Sigma2[(k*25) - 24,]), w = fit_lddp$P[(k*25) - 24,])
    dpred_lddp[k, , l] <- dnorMix(y_grid, aux)
    }
}

dpred_lddp_m <- apply(dpred_lddp, c(2,3), mean)
dpred_lddp_l <- apply(dpred_lddp, c(2,3), quantile, prob = 0.025)
dpred_lddp_h <- apply(dpred_lddp, c(2,3), quantile, prob = 0.975)

# Estimated l1 distance for density
l1_dist_lddp <- colSums(abs(f_true_new - dpred_lddp_m))*(y_grid[2] - y_grid[1])
mean(l1_dist_lddp)


