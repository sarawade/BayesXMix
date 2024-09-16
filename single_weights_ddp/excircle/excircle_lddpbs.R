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

# MCMC configuration for all knots configurations
mcmc <- list(nburn = 5000, nsave = 10000, nskip = 1)

# Design and prediction matrices for LDDPBS models for different configuration of knots
# no interior knots
knots0 <- c()
X0 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots0, intercept = FALSE)
)

X0_predict <-cbind(rep(1, n_new),
                   predict(bs(x[, 1], degree = 3, knots = knots0, intercept = FALSE),
                           x_new[, 1])
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

# Fit LDDP-BS with one interior knot (use the quantiles)
knots1 <- quantile(x = x[ , 1], probs = 0.5)
X1 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots1, intercept = FALSE)
)
k1 <- ncol(X1)
X1_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots1, intercept = FALSE),
                           x_new[ , 1])
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

# Fit LDDP-BS with two interior knots for x1 (use the quantiles)
knots2 <- quantile(x = x[ , 1], probs = c(0.33, 0.66))
X2 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots2, intercept = FALSE)
)
k2 <- ncol(X2)
X2_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots2, intercept = FALSE),
                           x_new[ , 1])
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

# Fit LDDP-BS with three interior knots for x1 (use the quantiles)
knots3 <- quantile(x = x[ , 1], probs = c(0.25, 0.5, 0.75))
X3 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots3, intercept = FALSE)
)
k3 <- ncol(X3)
X3_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots3, intercept = FALSE),
                           x_new[ , 1])
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


# Fit LDDP-BS with four interior knots for x1 (use the quantiles)
knots4 <- quantile(x = x[ , 1], probs = c(0.2, 0.4, 0.6, 0.8))
X4 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots4, intercept = FALSE)
)
k4 <- ncol(X4)
X4_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots4, intercept = FALSE),
                           x_new[ , 1])
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

# Fit LDDP-BS with five interior knots for x1 (use the quantiles)
knots5 <- quantile(x = x[ , 1], probs = c(1/6, 2/6, 3/6, 4/6, 5/6))
X5 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots5, intercept = FALSE)
)
k5 <- ncol(X5)
X5_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots5, intercept = FALSE),
                           x_new[ , 1])
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

waic5<- waicnp(y = y, X = X5, res = fit_lddp_bs_5, L = 20, termsum = NULL)
lpml5 <- lpml(y = y, X = X5, res = fit_lddp_bs_5, L = 20, termsum = NULL) 
waic5$WAIC; lpml5$lpml

# Fit LDDP-BS with six interior knots for x1 (use the quantiles)
knots6 <- quantile(x = x[ , 1], probs = c(1/7, 2/7, 3/7, 4/7, 5/7, 6/7))
X6 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots6, intercept = FALSE)
)
k6 <- ncol(X6)
X6_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots6, intercept = FALSE),
                           x_new[ , 1])
)

fit_lm6 <- lm(y ~ X6 - 1)
prior6 <- list(m0 = fit_lm6$coefficients,
               S0 = solve(t(X6)%*%X6)*(sigma(fit_lm6)^2), 
               nu = k6 + 2, 
               Psi = 30*(solve(t(X6)%*%X6)*(sigma(fit_lm6)^2)),
               a = 2,
               b = (sigma(fit_lm6)^2)/2,
               alpha = 1, 
               L = 20)

set.seed(123)
fit_lddp_bs_6 <- bddp(y = y,
                      X = X6,
                      prior = prior6, 
                      mcmc = mcmc,  
                      standardise = FALSE)

waic6 <- waicnp(y = y, X = X6, res = fit_lddp_bs_6, L = 20, termsum = NULL)
lpml6 <- lpml(y = y, X = X6, res = fit_lddp_bs_6, L = 20, termsum = NULL) 
waic6$WAIC; lpml6$lpml

# Fit LDDP-BS with seven interior knots for x1 (use the quantiles)
knots7 <- quantile(x = x[ , 1], probs = c(1/8, 2/8, 3/8, 4/8, 5/8, 6/8, 7/8))
X7 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots7, intercept = FALSE)
)
k7 <- ncol(X7)
X7_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots7, intercept = FALSE),
                           x_new[ , 1])
)

fit_lm7 <- lm(y ~ X7 - 1)
prior7 <- list(m0 = fit_lm7$coefficients,
               S0 = solve(t(X7)%*%X7)*(sigma(fit_lm7)^2), 
               nu = k7 + 2, 
               Psi = 30*(solve(t(X7)%*%X7)*(sigma(fit_lm7)^2)),
               a = 2,
               b = (sigma(fit_lm7)^2)/2,
               alpha = 1, 
               L = 20)

set.seed(123)
fit_lddp_bs_7 <- bddp(y = y,
                      X = X7,
                      prior = prior7, 
                      mcmc = mcmc,  
                      standardise = FALSE)

waic7 <- waicnp(y = y, X = X7, res = fit_lddp_bs_7, L = 20, termsum = NULL)
lpml7 <- lpml(y = y, X = X7, res = fit_lddp_bs_7, L = 20, termsum = NULL) 
waic7$WAIC; lpml7$lpml

# Fit LDDP-BS with eight interior knots for x1 (use the quantiles)
knots8 <- quantile(x = x[ , 1], probs = c(1/9, 2/9, 3/9, 4/9, 5/9, 6/9, 7/9, 8/9))
X8 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots8, intercept = FALSE)
)
k8 <- ncol(X8)
X8_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots8, intercept = FALSE),
                           x_new[ , 1])
)

fit_lm8 <- lm(y ~ X8 - 1)
prior8 <- list(m0 = fit_lm8$coefficients,
               S0 = solve(t(X8)%*%X8)*(sigma(fit_lm8)^2), 
               nu = k8 + 2, 
               Psi = 30*(solve(t(X8)%*%X8)*(sigma(fit_lm8)^2)),
               a = 2,
               b = (sigma(fit_lm8)^2)/2,
               alpha = 1, 
               L = 20)

set.seed(123)
fit_lddp_bs_8 <- bddp(y = y,
                      X = X8,
                      prior = prior8, 
                      mcmc = mcmc,  
                      standardise = FALSE)

waic8 <- waicnp(y = y, X = X8, res = fit_lddp_bs_8, L = 20, termsum = NULL)
lpml8 <- lpml(y = y, X = X8, res = fit_lddp_bs_8, L = 20, termsum = NULL) 
waic8$WAIC; lpml8$lpml

# Fit LDDP-BS with nine interior knots for x1 (use the quantiles)
knots9 <- quantile(x = x[ , 1], probs = c(1/10, 2/10, 3/10, 4/10, 5/10, 6/10, 7/10, 8/10, 9/10))
X9 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots9, intercept = FALSE)
)
k9 <- ncol(X9)
X9_predict <-cbind(rep(1, n_new),
                   predict(bs(x[ , 1], degree = 3, knots = knots9, intercept = FALSE),
                           x_new[ , 1])
)

fit_lm9 <- lm(y ~ X9 - 1)
prior9 <- list(m0 = fit_lm9$coefficients,
               S0 = solve(t(X9)%*%X9)*(sigma(fit_lm9)^2), 
               nu = k9 + 2, 
               Psi = 30*(solve(t(X9)%*%X9)*(sigma(fit_lm9)^2)),
               a = 2,
               b = (sigma(fit_lm9)^2)/2,
               alpha = 1, 
               L = 20)

set.seed(123)
fit_lddp_bs_9 <- bddp(y = y,
                      X = X9,
                      prior = prior9, 
                      mcmc = mcmc,  
                      standardise = FALSE)

waic9 <- waicnp(y = y, X = X9, res = fit_lddp_bs_9, L = 20, termsum = NULL)
lpml9 <- lpml(y = y, X = X9, res = fit_lddp_bs_9, L = 20, termsum = NULL) 
waic9$WAIC; lpml9$lpml

# Fit LDDP-BS with ten interior knots for x1  (use the quantiles)
knots10 <- quantile(x = x[ , 1], probs = c(1/11, 2/11, 3/11, 4/11, 5/11, 6/11, 7/11, 8/11, 9/11, 10/11))
X10 <- cbind(rep(1, n),
             bs(x[ , 1], degree = 3, knots = knots10, intercept = FALSE)
)
k10 <- ncol(X10)
X10_predict <-cbind(rep(1, n_new),
                    predict(bs(x[ , 1], degree = 3, knots = knots10, intercept = FALSE),
                            x_new[ , 1])
)

fit_lm10 <- lm(y ~ X10 - 1)
prior10 <- list(m0 = fit_lm10$coefficients,
                S0 = solve(t(X10)%*%X10)*(sigma(fit_lm10)^2), 
                nu = k10 + 2, 
                Psi = 30*(solve(t(X10)%*%X10)*(sigma(fit_lm10)^2)),
                a = 2,
                b = (sigma(fit_lm10)^2)/2,
                alpha = 1, 
                L = 20)

set.seed(123)
fit_lddp_bs_10 <- bddp(y = y,
                       X = X10,
                       prior = prior10, 
                       mcmc = mcmc,  
                       standardise = FALSE)

waic10 <- waicnp(y = y, X = X10, res = fit_lddp_bs_10, L = 20, termsum = NULL)
lpml10 <- lpml(y = y, X = X10, res = fit_lddp_bs_10, L = 20, termsum = NULL) 
waic10$WAIC; lpml10$lpml

# Model preferred is the model with five knots so all the remaining calculations will be done for this model

# Estimated clustering
psm_lddp_bs_5 <- comp.psm(fit_lddp_bs_5$z)
output_vi_lddp_bs_5 <- minVI(psm_lddp_bs_5, fit_lddp_bs_5$z)

ggplot() +
  geom_point(aes(x = x[,1], y = y, color = as.factor(output_vi_lddp_bs_5$cl))) +
  theme_bw() +
  labs( x = "x_1", y = "y", color = "Cluster") +
  geom_point(aes(x = x[,1], y = m_true), col ="red")

#Predictive densities
dpred_lddp_bs_5 <- array(0, c(mcmc$nsave, m2, n_new))

for(k in 1:mcmc$nsave){
  mu <- tcrossprod(X5_predict, fit_lddp_bs_5$Beta[k, , ])
  
  for(l in 1:n_new) {
    aux <- norMix(mu = c(mu[l,]), sigma = sqrt(fit_lddp_bs_5$Sigma2[k,]), w = fit_lddp_bs_5$P[k,])
    dpred_lddp_bs_5[k, , l] <- dnorMix(y_grid, aux)
  }
}

dpred_lddp_bs_5_m <- apply(dpred_lddp_bs_5, c(2,3), mean)
dpred_lddp_bs_5_l <- apply(dpred_lddp_bs_5, c(2,3), quantile, prob = 0.025)
dpred_lddp_bs_5_h <- apply(dpred_lddp_bs_5, c(2,3), quantile, prob = 0.975)

f_pred_lddp_bs_5 <- data.frame(x = rep(x_new, each = m2), 
                               y = rep(y_grid, n_new), 
                               density = c(dpred_lddp_bs_5_m))

ggplot(f_pred_lddp_bs_5) +
  geom_raster(aes(x, y, fill = density), interpolate = TRUE) +
  scale_fill_gradient(low = "white", high ="red") +
  theme_classic() 

xs = sort(x_new[,1], index.return=TRUE)
inds = xs$ix[c(seq(1,n_new,n_new/5),n_new)]
cols <- rainbow(6)

ggplot() +
  geom_line(aes(x = y_grid, y = dpred_lddp_bs_5_m[, inds[1]]), col = cols[1]) +
  geom_line(aes(x = y_grid, y = f_true_new[, inds[1]]), 
            col = cols[1], linetype = "dashed") +
  geom_ribbon(aes(x = y_grid, ymin = dpred_lddp_bs_5_l[, inds[1]], ymax = dpred_lddp_bs_5_h[, inds[1]]), 
              alpha = 0.2, fill = cols[1]) +
  
  geom_line(aes(x = y_grid, y = dpred_lddp_bs_5_m[, inds[3]]), col = cols[3]) +
  geom_line(aes(x = y_grid, y = f_true_new[, inds[3]]), 
            col = cols[3], linetype = "dashed") +
  geom_ribbon(aes(x = y_grid, ymin = dpred_lddp_bs_5_l[, inds[3]], ymax = dpred_lddp_bs_5_h[, inds[3]]), 
              alpha = 0.2, fill = cols[3]) +
  
  geom_line(aes(x = y_grid, y = dpred_lddp_bs_5_m[, inds[5]]), col = cols[5]) +
  geom_line(aes(x = y_grid, y = f_true_new[, inds[5]]), 
            col = cols[5], linetype = "dashed") +
  geom_ribbon(aes(x = y_grid, ymin = dpred_lddp_bs_5_l[, inds[5]], ymax = dpred_lddp_bs_5_h[, inds[5]]), 
              alpha = 0.2, fill = cols[5]) +
  theme_bw() +
  labs(x = "y", y = "Density") +
  ylim(0, 9)


