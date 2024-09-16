require(ggplot2)
require(mcclust.ext)
require(nor1mix)

setwd("/Users/vandainacio/Dropbox/BNPregression/Code/single_weights_DDP")
source("main_function_mcmc.R")

# Generate data
set.seed(1010101)
n <- 800
p <- 2
x <- matrix(0,n,p)

x[ , 1] <- runif(n, -1, 8)
x[ , 2] = (x[ , 1] - 3.5)^2/3 - 1 + rnorm(n, 0, 0.05)

comp_ind <- runif(n) < 0.5
comp_mean <- (x[ , 1] + 2)/60
eps <- rnorm(n, 0, 0.05) + comp_mean*comp_ind - comp_mean*(1 - comp_ind)
y <- 5-log(x[,1]+2) + eps
m_true <- 5 - log(x[ , 1] + 2)

# Prediction
n_new <- 800
x_new <- matrix(0, n_new, p)

x_new[ , 1] <- runif(n_new, -1, 8)
x_new[ , 2] <- (x_new[ , 1] - 3.5)^2/3 - 1 + rnorm(n_new, 0, 0.05)

m_true_new <- 5 - log(x_new[ , 1] + 2)

y_grid <- seq(2.3, 5.4,.02)
m2 <- length(y_grid)
f_true_new <- 0.5*dnorm(matrix(y_grid, nrow = m2, ncol = n_new) - t(matrix(m_true_new, nrow = n_new, ncol = m2)),
                        t(matrix((x_new[ , 1] + 2)/60, nrow = n_new, ncol = m2)), 0.05) + 
  0.5*dnorm(matrix(y_grid, nrow = m2, ncol = n_new) - t(matrix(m_true_new, nrow = n_new, ncol = m2)),
            -t(matrix((x_new[ , 1] + 2)/60, nrow = n_new, ncol = m2)), 0.05)

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
mcmc <- list(nburn = 5000, nsave = 15000, nskip = 1)

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
  geom_line(aes(x = x[, 1], y = m_true), col ="red")

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
  geom_line(aes(x = x_new[ , 1], y = m_pred_lddp_m), col = "black") +
  geom_line(aes(x = x_new[ , 1], y = m_true_new), col = "red") +
  geom_ribbon(aes(x = x_new[ , 1], ymin = m_pred_lddp_l, ymax = m_pred_lddp_h),
              alpha = 0.2) +
  theme_bw() +
  labs( x = "x_1", y = "y") +
  ylim(1.7, 5.1)

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

xs <- sort(x_new[ ,1], index.return = TRUE) 
inds <- xs$ix[c(seq(1, n_new, n_new/5), n_new)]

cols <- rainbow(6)

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
  ylim(0, 8.75) + 
  labs( x = "y", y = "Density") 