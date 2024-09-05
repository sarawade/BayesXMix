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

setwd("/Users/swade/Documents/GitHub/BayesXMix/LSBP/ex1")

##############################################
## Example 1: normal errors
##############################################

# Load data
load("../.././data_simulation/ex1/data_ex1.RData")
ex1data  <- data.frame(x1 = x[,1],x2 = x[,2],y = y)

# MCMC parameters (iterations, burnin, truncation level)
R         <- 20000     # Number of iterations
burn_in   <- 5000      # Burn-in period
H         <- 20        # Number of mixture components

# Set prior
# Note that empirically smaller/larger values of B_mixing seem to result in smoother/sharper boundaries
lmfit = lm(y~x1+x2, data = ex1data)
prior       <- prior_LSBP(p_kernel = p+1,p_mixing = p+1, 
                          b_mixing =c(0,0,0), B_mixing=diag(c(100,10,10)), 
                          b_kernel = c(lmfit$coefficients), B_kernel=summary(lmfit)$cov.unscaled, 
                          a_tau = 2, b_tau= sum(lmfit$residuals^2)/(n-p-1))

# Local linear model
model_formula <- Formula::as.Formula(y ~ x1 +x2 )

### Fit model
set.seed(10) # The seed is setted so that the Gibbs sampler is reproducible.
fit_Gibbs   <- LSBP_Gibbs(model_formula, data=ex1data, H=H, prior=prior, 
                          control=control_Gibbs(R=R,burn_in=burn_in,method_init="cluster"), verbose=TRUE)

save.image("ex1_lsbp.RData")

######## Results

### Number of clusters
ggplot()+
  geom_bar(aes(fit_Gibbs$param$nclust)) +
  theme_bw()

##### Prediction

## Compute predictive mean
newdata     <- data.frame(y=0, x1=x_new[,1], x2=x_new[,2])

gibbs_mean = predict(fit_Gibbs, type = "mean", newdata=newdata)

ypred_lsbp = apply(gibbs_mean,2,mean)
lpred_lsbp = apply(gibbs_mean,2,function(x) quantile(x,0.025))
upred_lsbp = apply(gibbs_mean,2,function(x) quantile(x,0.975))

## Compute predictive conditional density

# Design matrix for the kernel and weights
X1           <- cbind(1,x_new[,1],x_new[,2])      
X2 = X1

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
# Design matrix for the kernel and weights for the subset
inds = c(1,2,201,202,401,402,601,602)
X1<- cbind(1,x_new[inds,1],x_new[inds,2]) 
X2 = X1

# Conditional density - Gibbs sampling
est_Gibbs <- matrix(0,length(y_grid),8)
lower_Gibbs <- matrix(0,length(y_grid),8)
upper_Gibbs <- matrix(0,length(y_grid),8)

for(i in 1:length(y_grid)){  # Cycle over the y grid
  pred_Gibbs_i = matrix(0, R, 8) 
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

save.image("ex1_lsbp.RData")

### Plot Prediction
#with credible intervals but without data
png("ex1_lsbp_pred.png",width = 500, height = 450)
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
png("ex1_lsbp_fpred.png",width = 500, height = 450)
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
  ylim(0,9.8) +
  xlim(2.4,4.2)
dev.off()


