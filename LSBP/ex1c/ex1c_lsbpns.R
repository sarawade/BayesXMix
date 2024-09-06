##############################################
## Data Generation for StatSci paper
##############################################
library(ggplot2)
library(tidyverse)
library(mvtnorm) 
library(msm)
library(MCMCpack)
#devtools::install_github("tommasorigon/LSBP")
library(LSBP)
library(splines)

setwd("/Users/swade/Documents/GitHub/BayesXMix/LSBP/ex1c")

##############################################
## Example 1: multimodal errors
##############################################

# Load data
load("../.././data_simulation/ex1c/data_ex1_mm.RData")
ex1data  <- data.frame(x1 = x[,1],x2 = x[,2],y = y)

# MCMC parameters (iterations, burnin, truncation level)
p_splines <- 5         # Number of splines components
R         <- 20000     # Number of iterations
burn_in   <- 5000      # Burn-in period
H         <- 20        # Number of mixture components

# Set prior
# Note that empirically smaller/larger values of B_mixing seem to result in smoother/sharper boundaries
lmfit = lm(y~x1+x2, data = ex1data)
prior       <- prior_LSBP(p_kernel = p+1,p_mixing = p*p_splines+1,
                          b_mixing = rep(0,p*p_splines+1), B_mixing=diag(c(100,rep(10,p*p_splines))), 
                          b_kernel = c(lmfit$coefficients), B_kernel=summary(lmfit)$cov.unscaled, 
                          a_tau = 2, b_tau= sum(lmfit$residuals^2)/(n-p-1))

# Linear kernel and splines model for mixing weights
Basis1       <- ns(ex1data$x1,p_splines)
Basis2       <- ns(ex1data$x2,p_splines)
ex1data  <- data.frame(ex1data, BS1=Basis1, BS2=Basis2)
# the symbol '|', separates the kernel covariates and the mixing covariates, respectively
model_formula <- Formula::as.Formula(y ~ x1 +x2 | BS1.1 + BS1.2 + BS1.3 + BS1.4 + BS1.5 + BS2.1 + BS2.2 + BS2.3 + BS2.4 + BS2.5)

### Fit model
set.seed(10) # The seed is setted so that the Gibbs sampler is reproducible.
fit_Gibbs   <- LSBP_Gibbs(model_formula, data=ex1data, H=H, prior=prior, 
                          control=control_Gibbs(R=R,burn_in=burn_in,method_init="cluster"), verbose=TRUE)

save.image("ex1c_lsbpns.RData")

######## Results

### Number of clusters
ggplot()+
  geom_bar(aes(fit_Gibbs$param$nclust)) +
  theme_bw()

##### Prediction

## Compute predictive mean
newdata <- data.frame(y=0, x1=x_new[,1], x2=x_new[,2], 
                      BS1= ns(x_new[,1],knots=attr(Basis1,"knots"),Boundary.knots=attr(Basis1,"Boundary.knots")),
                      BS2= ns(x_new[,2],knots=attr(Basis2,"knots"),Boundary.knots=attr(Basis2,"Boundary.knots")))


gibbs_mean = predict(fit_Gibbs, type = "mean", newdata=newdata)

ypred_lsbp = apply(gibbs_mean,2,mean)
lpred_lsbp = apply(gibbs_mean,2,function(x) quantile(x,0.025))
upred_lsbp = apply(gibbs_mean,2,function(x) quantile(x,0.975))

## Compute predictive conditional density

# Design matrix for the kernel
X1 <- cbind(1,x_new[,1],x_new[,2])      
# Design matrix for the weights
X2  <- cbind(1,ns(x_new[,1],knots=attr(Basis1,"knots"),Boundary.knots=attr(Basis1,"Boundary.knots")),
             ns(x_new[,2], knots=attr(Basis2,"knots"), Boundary.knots=attr(Basis2,"Boundary.knots")))

# Conditional density - Gibbs sampling
est_Gibbs_all <- matrix(0,length(y_grid),n_new)

for(i in 1:length(y_grid)){  # Cycle over the y grid
  pred_Gibbs_i = matrix(0, R, n_new) 
  for(r in 1:R){      # Cycle over the iterations of the MCMC chain
    pred_Gibbs_i[r,] <- c(LSBP_density(y_grid[i],X1,X2,
                                       fit_Gibbs$param$beta_mixing[r,,],
                                       fit_Gibbs$param$beta_kernel[r,,],
                                       fit_Gibbs$param$tau[r,]))
  }
  est_Gibbs_all[i,]  = apply(pred_Gibbs_i,2,mean)
  print(paste(i/length(y_grid)*100, '%'))
}

# For a subet of points compute also pointwise credible intervals
# Design matrix for the kernel
X1<- cbind(1,x_new[inds,1],x_new[inds,2]) 
# Design matrix for the weights
X2  <- cbind(1,ns(x_new[inds,1],knots=attr(Basis1,"knots"),Boundary.knots=attr(Basis1,"Boundary.knots")),
             ns(x_new[inds,2], knots=attr(Basis2,"knots"), Boundary.knots=attr(Basis2,"Boundary.knots")))


# Conditional density - Gibbs sampling
est_Gibbs <- matrix(0,length(y_grid),length(inds))
lower_Gibbs <- matrix(0,length(y_grid),length(inds))
upper_Gibbs <- matrix(0,length(y_grid),length(inds))

for(i in 1:length(y_grid)){  # Cycle over the y grid
  pred_Gibbs_i = matrix(0, R,length(inds)) 
  for(r in 1:R){      # Cycle over the iterations of the MCMC chain
    pred_Gibbs_i[r,] <- c(LSBP_density(y_grid[i],X1,X2,
                                       fit_Gibbs$param$beta_mixing[r,,],
                                       fit_Gibbs$param$beta_kernel[r,,],
                                       fit_Gibbs$param$tau[r,]))
  }
  est_Gibbs[i,]  = apply(pred_Gibbs_i,2,mean)
  lower_Gibbs[i,]    <- apply(pred_Gibbs_i,2,function(x) quantile(x,0.025))
  upper_Gibbs[i,]    <- apply(pred_Gibbs_i,2,function(x) quantile(x,0.975))
}

save.image("ex1c_lsbpns.RData")

### Plot Prediction
#with credible intervals but without data
png("ex1c_lsbpns_pred.png",width = 500, height = 450)
ggplot() +
  geom_line(aes(x = x_new[,1], y = ypred_lsbp), col = "black") +
  geom_line(aes(x = x_new[,1], y = m_true_new), col = "red") +
  geom_ribbon(aes(x = x_new[,1], ymin=lpred_lsbp, ymax=upred_lsbp), alpha=0.2) +
  theme_bw() +
  labs( x = "x_1", y = "y") +
  ylim(1.7,5.1)
dev.off()

#with data but without credible intervals
ggplot() +
  geom_point(aes(x = x[,1], y = y)) +
  geom_line(aes(x = x_new[,1], y = ypred_lsbp), col = "black") +
  geom_line(aes(x = x_new[,1], y = m_true_new), col = "red") +
  theme_bw() +
  labs( x = "x_1", y = "y")

#Plot y_hat_true vs y_pred
ggplot() +
  geom_point(aes(x = m_true_new, y = ypred_lsbp), col = "black") +
  geom_point(aes(x = m_true_new, y = lpred_lsbp), col = "black", shape = 25) +
  geom_point(aes(x = m_true_new, y = upred_lsbp), col = "black", shape = 24) +
  geom_line(aes(x = m_true_new, y = m_true_new), col = "grey") +
  theme_bw() +
  labs( x = "True y", y = "Estimated y")

#PLot density for a single observation
i = 1
ggplot() +
  geom_line(aes(x = y_grid, y = est_Gibbs[,i]), col = "black") +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[i]]), col = "red") +
  geom_ribbon(aes(x=y_grid, ymin=lower_Gibbs[,i], ymax=upper_Gibbs[,i]), alpha=0.2) +
  theme_bw() +
  labs( x = "y", y = "Density")

#PLot density for a few observations
png("ex1c_lsbpns_fpred.png",width = 500, height = 450)
cols = rainbow(6)
ggplot() +
  geom_line(aes(x = y_grid, y = est_Gibbs[,1]), col = cols[1]) +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[1]]), col = cols[1],linetype = "dashed") +
  geom_ribbon(aes(x=y_grid, ymin=lower_Gibbs[,1], ymax=upper_Gibbs[,1]), alpha=0.2, fill = cols[1]) +
  geom_line(aes(x = y_grid, y = est_Gibbs[,3]), col = cols[3]) +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[3]]), col = cols[3],linetype = "dashed") +
  geom_ribbon(aes(x=y_grid, ymin=lower_Gibbs[,3], ymax=upper_Gibbs[,3]), alpha=0.2, fill = cols[3]) +
  geom_line(aes(x = y_grid, y = est_Gibbs[,5]), col = cols[5]) +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[5]]), col = cols[5],linetype = "dashed") +
  geom_ribbon(aes(x=y_grid, ymin=lower_Gibbs[,5], ymax=upper_Gibbs[,5]), alpha=0.2, fill = cols[5]) +
  theme_bw() +
  labs( x = "y", y = "Density")+
  ylim(0,9.5)
dev.off()

#empirical l2 prediction error
l2_err=sum(((m_true_new-ypred_lsbp)^2)/n_new)^.5

l1_err=sum((abs(m_true_new-ypred_lsbp))/n_new)

#estimated l1 distance for density
#.02 is grid with, CHANGE if grid with changes
l1_dist=colSums(abs(f_true_new-est_Gibbs_all))*(y_grid[2]-y_grid[1])

#Average l1 distance
mean(l1_dist)

#Max l1 dist
max(l1_dist)

#Min l1 dist
min(l1_dist)

# empirical coverage for prediction
ec_pred = sum((m_true_new>=lpred_lsbp)&(m_true_new<=upred_lsbp))/n_new

# credible interval length for prediction
ci_length = upred_lsbp - lpred_lsbp
mean(ci_length)
min(ci_length)
max(ci_length)

