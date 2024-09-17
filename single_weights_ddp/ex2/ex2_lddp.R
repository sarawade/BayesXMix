require(ggplot2)
require(mcclust.ext)
require(nor1mix)

setwd("/Users/vandainacio/Dropbox/BNPregression/Code/single_weights_DDP")
source("main_function_mcmc.R")

# Generate data
set.seed(1010101)
n <- 400
p <- 1
x <- matrix(0, n, p)
y <- rep(0, n)

x[, 1] <- runif(n, -2, 10)

m  <-  rep(0, n)
m[x <= 2] <- 0
m[x > 2 & x <= 5] <- 2*(x[x > 2 & x <=5 ] -2)
m[x > 5] <- 6
eps <- rep(0, n)
eps[x <= 2] <- rnorm(sum(x <= 2), 0, 0.2)
eps[x > 2 & x <= 5] <- rnorm(sum(x > 2 & x <= 5), 0, 0.05)
eps[x > 5] <- rnorm(sum(x > 5), 0, (x[x > 5] - 5)^2/15 + 0.01)
y <- m + eps

# Prediction
n_new <- 1201
x_new <- matrix(0, n_new, p)

x_new[, 1] <- seq(-2, 10, 12/(n_new - 1))

m_true_new <- rep(0, n_new)
m_true_new[x_new > 2 & x_new <= 5] <- 2*(x_new[x_new > 2 & x_new <= 5] - 2)
m_true_new[x_new > 5] <- 6

y_grid <- seq(-1, 10, .02)
m2 <- length(y_grid)
f_true_new <- matrix(0, m2, n_new)
f_true_new[, x_new <= 2] <- dnorm(matrix(y_grid, nrow = m2, ncol = sum(x_new <= 2)),
                                  t(matrix(m_true_new[x_new <= 2],
                                           nrow = sum(x_new <= 2), ncol = m2)), 0.2)
f_true_new[, x_new > 2 & x_new <= 5] <- dnorm(matrix(y_grid, nrow = m2, 
                                                     ncol = sum(x_new > 2 & x_new <= 5)),
                                              t(matrix(m_true_new[x_new > 2 & x_new <= 5],
                                                       nrow = sum(x_new > 2 & x_new <= 5),
                                                       ncol = m2)), 0.05)
f_true_new[, x_new > 5] <- dnorm(matrix(y_grid, nrow = m2, ncol = sum(x_new > 5)),
                                 t(matrix(m_true_new[x_new > 5],
                                          nrow = sum(x_new > 5), ncol = m2)),
                                 t(matrix((x_new[x_new > 5] - 5)^2/15 + 0.01,
                                          nrow = sum(x_new > 5), ncol = m2)))

# LDDP model
# Design and Prediction matrices
X <- cbind(rep(1, n), x[, 1])
k <- ncol(X)
X_predict <-cbind(rep(1, n_new), x_new[, 1])

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
mcmc <- list(nburn = 5000, nsave = 5000, nskip = 1)

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
  geom_point(aes(x = x[, 1], y = y, color = as.factor(output_vi_lddp$cl))) +
  theme_bw() +
  labs( x = "x_1", y = "y", color = "Cluster") +
  geom_line(aes(x = x[, 1], y = m), col ="red")

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

ggplot() +
  geom_line(aes(x = x_new[,1], y = m_pred_lddp_m), col = "black") +
  geom_line(aes(x = x_new[,1], y = m_true_new), col = "red") +
  geom_ribbon(aes(x = x_new[,1], ymin = m_pred_lddp_l, ymax = m_pred_lddp_h),
              alpha = 0.2) +
  theme_bw() +
  labs( x = "x_1", y = "y") 

# Empirical l2 prediction error
l2_err_lddp <- sum(((m_true_new - m_pred_lddp_m)^2)/n_new)^.5
l2_err_lddp

# Empirical coverage for prediction
ec_pred_lddp <- sum((m_true_new >= m_pred_lddp_l)&(m_true_new <= m_pred_lddp_h))/n_new
ec_pred_lddp

# Credible interval length for prediction
ci_length_lddp <- m_pred_lddp_h - m_pred_lddp_l
mean(ci_length_lddp)

# Predictive densities and plots
dpred_lddp <- array(0, c(mcmc$nsave, m2, n_new))

for(k in 1:mcmc$nsave){
  mu <- tcrossprod(X_predict, fit_lddp$Beta[k, , ])
  
  for(l in 1:n_new) {
    aux <- norMix(mu = c(mu[l, ]), sigma = sqrt(fit_lddp$Sigma2[k,]), w = fit_lddp$P[k,])
    dpred_lddp[k, , l] <- dnorMix(y_grid, aux)
  }
}

dpred_lddp_m <- apply(dpred_lddp, c(2,3), mean)
dpred_lddp_l <- apply(dpred_lddp, c(2,3), quantile, prob = 0.025)
dpred_lddp_h <- apply(dpred_lddp, c(2,3), quantile, prob = 0.975)

# Estimated l1 distance for density
l1_dist_lddp <- colSums(abs(f_true_new - dpred_lddp_m))*(y_grid[2] - y_grid[1])
mean(l1_dist_lddp)

f_pred_lddp <- data.frame(x = rep(x_new, each = m2), 
                          y = rep(y_grid, n_new), 
                          density = c(dpred_lddp_m))

ggplot(f_pred_lddp) +
  geom_raster(aes(x, y, fill = density), interpolate = TRUE) +
  scale_fill_gradient(low = "white", high ="red") +
  theme_classic() +
  ylim(-1.15, 10.25)
