require(ggplot2)
require(mcclust.ext)
require(nor1mix)

setwd("/Users/vandainacio/Dropbox/BNPregression/Code/single_weights_DDP")
source("main_function_mcmc.R")

# Generate data
set.seed(1010101)
n <- 800
p <- 1
x <- matrix(0,n,p)

t <- runif(n, min = 0, max = 2*pi) 
x[, 1] <- cos(t)

y <- sin(t) + rnorm(n, 0, 0.05)
m_true <- sin(t)

# Prediction
x_new <- seq(-1, 1, pi/(2^8))
n_new <- length(x_new)
x_new <- matrix(x_new, n_new, p)

m_true_new <- rep(0, n_new)

# True density
y_grid <- seq(-1.2, 1.2, .02)
m2 <- length(y_grid)
f_true_new <- 0.5*dnorm(matrix(y_grid, nrow = m2, ncol = n_new),t(matrix(sin(acos(x_new)), nrow = n_new, ncol = m2)), 0.05) + 
  0.5*dnorm(matrix(y_grid, nrow = m2, ncol = n_new),t(matrix(sin(acos(x_new) + pi), nrow = n_new, ncol = m2)), 0.05)

# LDDP model
# Design and Prediction matrices
X <- cbind(rep(1, n), x[ , 1])
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
mcmc <- list(nburn = 5000, nsave = 10000, nskip = 1)

# Run MCMC
fit_lddp <- bddp(y = y,
                 X = X,
                 prior = prior, 
                 mcmc = mcmc,
                 standardise = FALSE)

# Estimated clustering
psm_lddp <- comp.psm(fit_lddp$z)
output_vi_lddp <- minVI(psm_lddp, fit_lddp$z)

ggplot() +
  geom_point(aes(x = x[,1], y = y, color = as.factor(output_vi_lddp$cl))) +
  theme_bw() +
  labs( x = "x_1", y = "y", color = "Cluster") +
  geom_point(aes(x = x[,1], y = m_true), col ="red")

# MCMC realisations predictive density
dpred_lddp <- array(0, c(mcmc$nsave, m2, n_new))

for(k in 1:mcmc$nsave){
  mu <- tcrossprod(X_predict, fit_lddp$Beta[k, , ])
  
  for(l in 1:n_new) {
    aux <- norMix(mu = c(mu[l,]), sigma = sqrt(fit_lddp$Sigma2[k,]), w = fit_lddp$P[k,])
    dpred_lddp[k, , l] <- dnorMix(y_grid, aux)
  }
}

dpred_lddp_m <- apply(dpred_lddp, c(2,3), mean)
dpred_lddp_l <- apply(dpred_lddp, c(2,3), quantile, prob = 0.025)
dpred_lddp_h <- apply(dpred_lddp, c(2,3), quantile, prob = 0.975)

f_pred_lddp <- data.frame(x = rep(x_new, each = m2), 
                          y = rep(y_grid, n_new), 
                          density = c(dpred_lddp_m))

# Heatmap predictive densities
ggplot(f_pred_lddp) +
  geom_raster(aes(x, y, fill = density), interpolate = TRUE) +
  scale_fill_gradient(low = "white", high ="red") +
  theme_classic() 

# Plot predictive densities
ggplot() +
  geom_line(aes(x = y_grid, y = dpred_lddp_m[, inds[1]]), col = cols[1]) +
  geom_line(aes(x = y_grid, y = f_true_new[, inds[1]]), 
            col = cols[1], linetype = "dashed") +
  geom_ribbon(aes(x = y_grid, ymin = dpred_lddp_l[, inds[1]], ymax = dpred_lddp_h[, inds[1]]), 
              alpha = 0.2, fill = cols[1]) +
  
  geom_line(aes(x = y_grid, y = dpred_lddp_m[, inds[3]]), col = cols[3]) +
  geom_line(aes(x = y_grid, y = f_true_new[, inds[3]]), 
            col = cols[3], linetype = "dashed") +
  geom_ribbon(aes(x = y_grid, ymin = dpred_lddp_l[, inds[3]], ymax = dpred_lddp_h[, inds[3]]), 
              alpha = 0.2, fill = cols[3]) +
  
  geom_line(aes(x = y_grid, y = dpred_lddp_m[, inds[5]]), col = cols[5]) +
  geom_line(aes(x = y_grid, y = f_true_new[, inds[5]]), 
            col = cols[5], linetype = "dashed") +
  geom_ribbon(aes(x = y_grid, ymin = dpred_lddp_l[, inds[5]], ymax = dpred_lddp_h[, inds[5]]), 
              alpha = 0.2, fill = cols[5]) +
  theme_bw() +
  labs( x = "y", y = "Density") +
  ylim(0, 9)
