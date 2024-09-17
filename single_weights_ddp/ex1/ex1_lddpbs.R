require(nor1mix)
require(ggplot2)
require(mcclust.ext)

# Change the below to your own directory
setwd("/Users/vandainacio/Dropbox/BNPregression/Code/single_weights_DDP")
source("main_function_mcmc.R")

# Generate data
set.seed(1010101)
n <- 200
p <- 2
x <-  matrix(0, n, p)

x[,1] <- runif(n, -1, 8)
x[,2] <- (x[, 1] - 3.5)^2/3 - 1 + rnorm(n, 0, 0.05)

y <- 5 - log(x[, 1] + 2) + rnorm(n, 0, 0.05)
m_true <- 5 - log(x[, 1] + 2)

# Prediction
n_new <- 800
x_new <- matrix(0, n_new, p)

x_new[, 1] <- runif(n_new, -1, 8)
x_new[, 2] <- (x_new[, 1] - 3.5)^2/3 - 1 + rnorm(n_new, 0, 0.05)

m_true_new <-  5-log(x_new[, 1] + 2)

# True density
y_grid <- seq(2.4, 5.3, .02)
m2 <- length(y_grid)
f_true_new <- dnorm(matrix(y_grid, nrow = m2, ncol = n_new), t(matrix(m_true_new, nrow = n_new, ncol = m2)), 0.05)

mcmc <- list(nburn = 5000, nsave = 5000, nskip = 1)

# Fit LDDP-BS with no interior knots for x1 and x2
knots0 <- c()
X0 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots0, intercept = FALSE),
            bs(x[ , 2], degree = 3, knots = knots0, intercept = FALSE)
)

X0_predict <-cbind(rep(1, n_new),
                   predict(bs(x[, 1], degree = 3, knots = knots0, intercept = FALSE),
                           x_new[, 1]),
                   predict(bs(x[, 2], degree = 3, knots = knots0, intercept = FALSE), 
                           x_new[, 2])
)
k0 <- ncol(X0)

fit_lm0 <- lm(y ~ X0 - 1)
prior0 <- list(m0 = fit_lm0$coefficients,
               S0 = solve(t(X0)%*%X0)*(sigma(fit_lm0)^2), 
               nu = k0 + 2, 
               Psi = 30*(solve(t(X0)%*%X0)*(sigma(fit_lm0)^2)),
               a = 2,
               b = (sigma(fit_lm0)^2)/2,
               alpha = 1, 
               L = 20)

set.seed(123)
fit_lddp_bs_0 <- bddp(y = y,
                      X = X0,
                      prior = prior0, 
                      mcmc = mcmc,
                      standardise = FALSE)

waic0 <- waicnp(y = y, X = X0, res = fit_lddp_bs_0, L = 20, termsum = NULL)
lpml0 <- lpml(y = y, X = X0, res = fit_lddp_bs_0, L = 20, termsum = NULL) 
waic0$WAIC
lpml0$lpml

# Fit LDDP-BS with one interior knot for x1 and x2 (use the quantiles)
knots11 <- quantile(x = x[ , 1], probs = 0.5)
knots12 <- quantile(x = x[ , 2], probs = 0.5)
X1 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots11, intercept = FALSE),
            bs(x[ , 2], degree = 3, knots = knots12, intercept = FALSE)
)
k1 <- ncol(X1)
X1_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots11, intercept = FALSE),
                           x_new[ , 1]),
                   predict(bs(x[ , 2], degree = 3, knots = knots12, intercept = FALSE), 
                           x_new[ , 2])
)

fit_lm1 <- lm(y ~ X1 - 1)
prior1 <- list(m0 = fit_lm1$coefficients,
               S0 = solve(t(X1)%*%X1)*(sigma(fit_lm1)^2), 
               nu = k1 + 2, 
               Psi = 30*(solve(t(X1)%*%X1)*(sigma(fit_lm1)^2)),
               a = 2,
               b = (sigma(fit_lm1)^2)/2,
               alpha = 1, 
               L = 20)

set.seed(123)
fit_lddp_bs_1 <- bddp(y = y,
                      X = X1,
                      prior = prior1, 
                      mcmc = mcmc,
                      standardise = FALSE)

waic1 <- waicnp(y = y, X = X1, res = fit_lddp_bs_1, L = 20, termsum = NULL)
lpml1 <- lpml(y = y, X = X1, res = fit_lddp_bs_1, L = 20, termsum = NULL) 
waic1$WAIC; lpml1$lpml

# Fit LDDP-BS with two interior knots for x1 and x2 (use the quantiles)
knots21 <- quantile(x = x[ , 1], probs = c(0.33, 0.66))
knots22 <- quantile(x = x[ , 2], probs = c(0.33, 0.66))
X2 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots21, intercept = FALSE),
            bs(x[ , 2], degree = 3, knots = knots22, intercept = FALSE)
)
k2 <- ncol(X2)
X2_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots21, intercept = FALSE),
                           x_new[ , 1]),
                   predict(bs(x[ , 2], degree = 3, knots = knots22, intercept = FALSE), 
                           x_new[ , 2])
)

fit_lm2 <- lm(y ~ X2 - 1)
prior2 <- list(m0 = fit_lm2$coefficients,
               S0 = solve(t(X2)%*%X2)*(sigma(fit_lm2)^2), 
               nu = k2 + 2, 
               Psi = 30*(solve(t(X2)%*%X2)*(sigma(fit_lm2)^2)),
               a = 2,
               b = (sigma(fit_lm2)^2)/2,
               alpha = 1, 
               L = 20)

set.seed(123)
fit_lddp_bs_2 <- bddp(y = y,
                      X = X2,
                      prior = prior2, 
                      mcmc = mcmc,
                      standardise = FALSE)

waic2 <- waicnp(y = y, X = X2, res = fit_lddp_bs_2, L = 20, termsum = NULL)
lpml2 <- lpml(y = y, X = X2, res = fit_lddp_bs_2, L = 20, termsum = NULL)
waic2$WAIC;  lpml2$lpml

# Fit LDDP-BS with three interior knots for x1 and x2 (use the quantiles)
knots31 <- quantile(x = x[ , 1], probs = c(0.25, 0.5, 0.75))
knots32 <- quantile(x = x[ , 2], probs = c(0.25, 0.5, 0.75))
X3 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots31, intercept = FALSE),
            bs(x[ , 2], degree = 3, knots = knots32, intercept = FALSE)
)
k3 <- ncol(X3)
X3_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots31, intercept = FALSE),
                           x_new[ , 1]),
                   predict(bs(x[ , 2], degree = 3, knots = knots32, intercept = FALSE), 
                           x_new[ , 2])
)

fit_lm3 <- lm(y ~ X3 - 1)
prior3 <- list(m0 = fit_lm3$coefficients,
               S0 = solve(t(X3)%*%X3)*(sigma(fit_lm3)^2), 
               nu = k3 + 2, 
               Psi = 30*(solve(t(X3)%*%X3)*(sigma(fit_lm3)^2)),
               a = 2,
               b = (sigma(fit_lm3)^2)/2,
               alpha = 1, 
               L = 20)

set.seed(123)
fit_lddp_bs_3 <- bddp(y = y,
                      X = X3,
                      prior = prior3, 
                      mcmc = mcmc,  
                      standardise = FALSE)

waic3 <- waicnp(y = y, X = X3, res = fit_lddp_bs_3, L = 20, termsum = NULL)
lpml3 <- lpml(y = y, X = X3, res = fit_lddp_bs_3, L = 20, termsum = NULL) 
waic3$WAIC; lpml3$lpml

# Model with two interior knots has a LPML only one unit higher than 
# the model with no interior knots and therefore we choose the most parsimonious model (with no interior knots)

# Estimated clustering and plot
psm_lddpbs <- comp.psm(fit_lddp_bs_0$z)
output_vi_lddpbs <- minVI(psm_lddpbs, fit_lddp_bs_0$z)

ggplot() +
  geom_point(aes(x = x[, 1], y = y, color = as.factor(output_vi_lddpbs$cl))) +
  theme_bw() +
  labs( x = "x_1", y = "y", color = "Cluster") +
  geom_line(aes(x = x[, 1], y = m_true), col ="red")

# Predictive regression function and plot
m_pred_lddp_bs_0 <- matrix(0, nrow = n_new, ncol = mcmc$nsave)
for(l in 1:mcmc$nsave){
  for(j in 1:n_new){  
    m_pred_lddp_bs_0[j, l] <- sum(fit_lddp_bs_0$P[l,]*X0_predict[j,]%*%t(fit_lddp_bs_0$Beta[l,,]))  
  }
}  

m_pred_lddp_bs_0_m <- apply(m_pred_lddp_bs_0, 1, mean)
m_pred_lddp_bs_0_l <- apply(m_pred_lddp_bs_0, 1, quantile, prob = 0.025)
m_pred_lddp_bs_0_h <- apply(m_pred_lddp_bs_0, 1, quantile, prob = 0.975)

ggplot() +
  geom_line(aes(x = x_new[ , 1], y = m_pred_lddp_bs_0_m), col = "black") +
  geom_line(aes(x = x_new[ , 1], y = m_true_new), col = "red") +
  geom_ribbon(aes(x = x_new[ , 1], ymin = m_pred_lddp_bs_0_l, ymax = m_pred_lddp_bs_0_h),
              alpha = 0.2) +
  theme_bw() +
  labs( x = "x_1", y = "y")+
  ylim(1.7, 5.1)

# Empirical l2 prediction error
l2_err_lddpbs <- sum(((m_true_new - m_pred_lddp_bs_0_m)^2)/n_new)^.5
l2_err_lddpbs

# Empirical coverage for prediction
ec_pred_lddpbs <- sum((m_true_new >= m_pred_lddp_bs_0_l)&(m_true_new <= m_pred_lddp_bs_0_h))/n_new
ec_pred_lddpbs

# Credible interval length for prediction
ci_length_lddpbs <- m_pred_lddp_bs_0_h - m_pred_lddp_bs_0_l
mean(ci_length_lddpbs)

# Predictive densities and plot
dpred_lddp_bs_0 <- array(0, c(mcmc$nsave, m2, n_new))

for(k in 1:mcmc$nsave){
  mu <- tcrossprod(X0_predict, fit_lddp_bs_0$Beta[k, , ])
  
  for(l in 1:n_new) {
    aux <- norMix(mu = c(mu[l,]), sigma = sqrt(fit_lddp_bs_0$Sigma2[k,]), w = fit_lddp_bs_0$P[k,])
    dpred_lddp_bs_0[k, , l] <- dnorMix(y_grid, aux)
  }
}

dpred_lddp_bs_0_m <- apply(dpred_lddp_bs_0, c(2,3), mean)
dpred_lddp_bs_0_l <- apply(dpred_lddp_bs_0, c(2,3), quantile, prob = 0.025)
dpred_lddp_bs_0_h <- apply(dpred_lddp_bs_0, c(2,3), quantile, prob = 0.975)

# Estimated l1 distance for density
l1_dist_lddpbs <- colSums(abs(f_true_new - dpred_lddp_bs_0_m))*(y_grid[2] - y_grid[1])

cols <- rainbow(6)
ggplot() +
  geom_line(aes(x = y_grid, y = dpred_lddp_bs_0_m[, 1]), col = cols[1]) +
  geom_line(aes(x = y_grid, y = f_true_new[, 1]), 
            col = cols[1],linetype = "dashed") +
  geom_ribbon(aes(x = y_grid, ymin = dpred_lddp_bs_0_l[, 1], ymax = dpred_lddp_bs_0_h[, 1]), 
              alpha = 0.2, fill = cols[1]) +
  geom_line(aes(x = y_grid, y = dpred_lddp_bs_0_m[, 201]), col = cols[3]) +
  geom_line(aes(x = y_grid, y = f_true_new[, 201]), 
            col = cols[3],linetype = "dashed") +
  geom_ribbon(aes(x = y_grid, ymin = dpred_lddp_bs_0_l[, 201], ymax = dpred_lddp_bs_0_h[, 201]), 
              alpha = 0.2, fill = cols[3]) +
  geom_line(aes(x = y_grid, y = dpred_lddp_bs_0_m[, 401]), col = cols[5]) +
  geom_line(aes(x = y_grid, y = f_true_new[, 401]), 
            col = cols[5],linetype = "dashed") +
  geom_ribbon(aes(x = y_grid, ymin = dpred_lddp_bs_0_l[, 401], ymax = dpred_lddp_bs_0_h[, 401]),
              alpha = 0.2, fill = cols[5]) +
  theme_bw() +
  labs( x = "y", y = "Density") +
  ylim(0, 10.80) + 
  xlim(2.4,4.2)