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

x[, 1] <- runif(n, -2, 2)
x[, 2] <- runif(n, -2, 2)

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

# MCMC settings for all knots' configurations
mcmc <- list(nburn = 5000, nsave = 50000, nskip = 1)

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

# Fit LDDP-BS with four interior knots for x1 and x2 (use the quantiles)
knots41 <- quantile(x = x[ , 1], probs = c(0.2, 0.4, 0.6, 0.8))
knots42 <- quantile(x = x[ , 2], probs = c(0.2, 0.4, 0.6, 0.8))
X4 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots41, intercept = FALSE),
            bs(x[ , 2], degree = 3, knots = knots42, intercept = FALSE)
)
k4 <- ncol(X4)
X4_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots41, intercept = FALSE),
                           x_new[ , 1]),
                   predict(bs(x[ , 2], degree = 3, knots = knots42, intercept = FALSE), 
                           x_new[ , 2])
)

fit_lm4 <- lm(y ~ X4 - 1)
prior4 <- list(m0 = fit_lm4$coefficients,
               S0 = solve(t(X4)%*%X4)*(sigma(fit_lm4)^2), 
               nu = k4 + 2, 
               Psi = 30*(solve(t(X4)%*%X4)*(sigma(fit_lm4)^2)),
               a = 2,
               b = (sigma(fit_lm4)^2)/2,
               alpha = 1, 
               L = 20)

set.seed(123)
fit_lddp_bs_4 <- bddp(y = y,
                      X = X4,
                      prior = prior4, 
                      mcmc = mcmc,  
                      standardise = FALSE)

waic4 <- waicnp(y = y, X = X4, res = fit_lddp_bs_4, L = 20, termsum = NULL)
lpml4 <- lpml(y = y, X = X4, res = fit_lddp_bs_4, L = 20, termsum = NULL) 
waic4$WAIC; lpml4$lpml


# Fit LDDP-BS with five interior knots for x1 and x2 (use the quantiles)
knots51 <- quantile(x = x[ , 1], probs = c(1/6, 2/6, 3/6, 4/6, 5/6))
knots52 <- quantile(x = x[ , 2], probs = c(1/6, 2/6, 3/6, 4/6, 5/6))
X5 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots51, intercept = FALSE),
            bs(x[ , 2], degree = 3, knots = knots52, intercept = FALSE)
)
k5 <- ncol(X5)
X5_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots51, intercept = FALSE),
                           x_new[ , 1]),
                   predict(bs(x[ , 2], degree = 3, knots = knots52, intercept = FALSE), 
                           x_new[ , 2])
)

fit_lm5 <- lm(y ~ X5 - 1)
prior5 <- list(m0 = fit_lm5$coefficients,
               S0 = solve(t(X5)%*%X5)*(sigma(fit_lm5)^2), 
               nu = k5 + 2, 
               Psi = 30*(solve(t(X5)%*%X5)*(sigma(fit_lm5)^2)),
               a = 2,
               b = (sigma(fit_lm5)^2)/2,
               alpha = 1, 
               L = 20)

set.seed(123)
fit_lddp_bs_5 <- bddp(y = y,
                      X = X5,
                      prior = prior5, 
                      mcmc = mcmc,  
                      standardise = FALSE)

waic5 <- waicnp(y = y, X = X5, res = fit_lddp_bs_5, L = 20, termsum = NULL)
lpml5 <- lpml(y = y, X = X5, res = fit_lddp_bs_5, L = 20, termsum = NULL) 
waic5$WAIC; lpml5$lpml

# Model with no interior knots is preferred and so the remainder calculations are done for this model
# Estimated clustering and plot
psm_lddpbs <- comp.psm(fit_lddp_bs_0$z)
output_vi_lddpbs <- minVI(psm_lddpbs, fit_lddp_bs_0$z)

ggplot() +
  geom_point(aes(x = x[,1], y = x[,2], color = as.factor(output_vi_lddpbs$cl))) +
  theme_bw() +
  labs( x = "x_1", y = "x_2", color = "Cluster")

# Predictive regression function and plot
m_pred_lddp_bs <- matrix(0, nrow = n_new, ncol = mcmc$nsave)
for(l in 1:mcmc$nsave){
  for(j in 1:n_new){  
    m_pred_lddp_bs[j, l] <- sum(fit_lddp_bs_0$P[l,]*X0_predict[j,]%*%t(fit_lddp_bs_0$Beta[l,,]))  
  }
}

m_pred_lddp_bs_m <- apply(m_pred_lddp_bs, 1, mean)
m_pred_lddp_bs_l <- apply(m_pred_lddp_bs, 1, quantile, prob = 0.025)
m_pred_lddp_bs_h <- apply(m_pred_lddp_bs, 1, quantile, prob = 0.975)

df <- data.frame(x1 = x_new[,1], x2 = x_new[,2], y = m_pred_lddp_bs_m)
ggplot(df) +
  geom_raster(aes(x1, x2, fill =  y),interpolate = TRUE) +
  theme_classic()

xnew2 <- 1.5
ggplot() +
  geom_point(aes(x = x_new[x_new[, 2] == xnew2, 1], y = m_pred_lddp_bs_m[x_new[, 2] == xnew2]), col = "black") +
  geom_point(aes(x = x_new[x_new[, 2] == xnew2, 1], y = m_true_new[x_new[, 2] == xnew2]), col = "red") +
  geom_ribbon(aes(x = x_new[x_new[, 2] == xnew2, 1],ymin = m_pred_lddp_bs_l[x_new[, 2] == xnew2], 
                  ymax = m_pred_lddp_bs_h[x_new[, 2] == xnew2]), alpha = 0.2) +
  theme_bw() +
  labs( x = "x_1", y = "y") 
ylab("y")

# Empirical l2 prediction error
l2_err_lddpbs <- sum(((m_true_new - m_pred_lddp_bs_m)^2)/n_new)^.5
l2_err_lddpbs

# Empirical coverage for prediction
ec_pred_lddpbs <- sum((m_true_new >= m_pred_lddp_bs_l)&(m_true_new <= m_pred_lddp_bs_h))/n_new
ec_pred_lddpbs

# Credible interval length for prediction
ci_length_lddpbs <- m_pred_lddp_bs_h - m_pred_lddp_bs_l
mean(ci_length_lddpbs)

# Predictive density function and plot
index <- seq(1, mcmc$nsave, by = 25)
dpred_lddpbs <- array(0,c(length(index), m2, n_new))

for(k in 1:length(index)){
  mu <- tcrossprod(X0_predict, fit_lddp_bs_0$Beta[(k*25) - 24, , ])
  
  for(l in 1:n_new) {
    aux <- norMix(mu = c(mu[l, ]), sigma = sqrt(fit_lddp_bs_0$Sigma2[(k*25) - 24,]), w = fit_lddp_bs_0$P[(k*25) - 24,])
    dpred_lddpbs[k, , l] <- dnorMix(y_grid, aux)
  }
}

dpred_lddpbs_m <- apply(dpred_lddpbs, c(2,3), mean)
dpred_lddpbs_l <- apply(dpred_lddpbs, c(2,3), quantile, prob = 0.025)
dpred_lddpbs_h <- apply(dpred_lddpbs, c(2,3), quantile, prob = 0.975)

# Estimated l1 distance for density
l1_dist_lddpbs <- colSums(abs(f_true_new - dpred_lddpbs_m))*(y_grid[2] - y_grid[1])
mean(l1_dist_lddpbs)
