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

setwd("/Users/swade/Documents/GitHub/BayesXMix/LSBP/excircle")

##############################################
## Example 4: Implicit/Circle
##############################################

# Load data
load("../.././data_simulation/excircle/data_ex_circle.RData")
ex1data  <- data.frame(x1 = x[,1],y = y)

# MCMC parameters (iterations, burnin, truncation level)
p_splines <- 5         # Number of splines components
R         <- 20000     # Number of iterations
burn_in   <- 5000      # Burn-in period
H         <- 20        # Number of mixture components

# Set prior
# Note that empirically smaller/larger values of B_mixing seem to result in smoother/sharper boundaries
lmfit = lm(y~x1, data = ex1data)
k0 = 10
prior       <- prior_LSBP(p_kernel = p+1,p_mixing = p*p_splines+1,
                          b_mixing = rep(0,p*p_splines+1), B_mixing=diag(c(100,rep(10,p*p_splines))), 
                          b_kernel = c(lmfit$coefficients), B_kernel=diag(c(10,1))*var(y), 
                          a_tau = 2, b_tau= sum(lmfit$residuals^2)/(n-p-1)/k0^2)

# Linear kernel and splines model for mixing weights
Basis1       <- ns(ex1data$x1,p_splines)
ex1data  <- data.frame(ex1data, BS1=Basis1)
# the symbol '|', separates the kernel covariates and the mixing covariates, respectively
model_formula <- Formula::as.Formula(y ~ x1 | BS1.1 + BS1.2 + BS1.3 + BS1.4 + BS1.5)

### Fit model
set.seed(10) # The seed is setted so that the Gibbs sampler is reproducible.
fit_Gibbs   <- LSBP_Gibbs(model_formula, data=ex1data, H=H, prior=prior, 
                          control=control_Gibbs(R=R,burn_in=burn_in,method_init="cluster"), verbose=TRUE)

save.image("excircle_lsbpns.RData")

######## Results

### Number of clusters
ggplot()+
  geom_bar(aes(fit_Gibbs$param$nclust)) +
  theme_bw()

##### Prediction

## Compute predictive mean
newdata <- data.frame(y=0, x1=x_new[,1],  
                      BS1= ns(x_new[,1],knots=attr(Basis1,"knots"),Boundary.knots=attr(Basis1,"Boundary.knots")))


gibbs_mean = predict(fit_Gibbs, type = "mean", newdata=newdata)

ypred_lsbp = apply(gibbs_mean,2,mean)
lpred_lsbp = apply(gibbs_mean,2,function(x) quantile(x,0.025))
upred_lsbp = apply(gibbs_mean,2,function(x) quantile(x,0.975))

## Compute predictive conditional density

# Design matrix for the kernel
X1 <- cbind(1,x_new[,1])      
# Design matrix for the weights
X2  <- cbind(1,ns(x_new[,1],knots=attr(Basis1,"knots"),Boundary.knots=attr(Basis1,"Boundary.knots")))

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
xs = sort(x_new[,1], index.return=TRUE)
inds = xs$ix[c(seq(1,n_new,n_new/5),n_new)]
X1<- cbind(1,x_new[inds,1]) 
# Design matrix for the weights
X2  <- cbind(1,ns(x_new[inds,1],knots=attr(Basis1,"knots"),Boundary.knots=attr(Basis1,"Boundary.knots")))


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

save.image("excircle_lsbpns.RData")

### Plot Prediction
#with credible intervals but without data
ggplot() +
  geom_line(aes(x = x_new[,1], y = ypred_lsbp), col = "black") +
  geom_line(aes(x = x_new[,1], y = m_true_new), col = "red") +
  geom_point(aes(x = x_new[,1], y = lpred_lsbp), col = "black", shape =25) +
  geom_point(aes(x = x_new[,1], y = upred_lsbp), col = "black", shape =17) +
  theme_bw() +
  labs( x = "x_1", y = "y") +
  ylim(-1.1,1.1)

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
png("excircle_lsbpns_fpred.png",width = 500, height = 450)
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
  labs( x = "y", y = "Density")
dev.off()

#PLot heatmap of density for a few observations
png("excircle_lsbpns_dens_heat",width = 500, height = 450)
df = data.frame(x = rep(x_new[,1], each = m2), y = rep(y_grid,n_new), density= c(est_Gibbs_all))
ggplot(df) +
  geom_raster(aes(x,y,fill=density),interpolate = TRUE) +
  scale_fill_gradient2(low="white", high="red") +
  theme_classic() 
dev.off()

#empirical l2 prediction error
l2_err=sum(((m_true_new-ypred_lsbp)^2)/n_new)^.5

l1_err=sum((abs(m_true_new-ypred_lsbp))/n_new)

#estimated l1 distance for density
#.02 is grid with, CHANGE if grid with changes
l1_dist=colSums(abs(f_true_new-est_Gibbs_all))*(y_grid[2]-y_grid[1])

# Compute the empirical coverage of the pointwise intervals for the densities
ec_dens = apply(f_true_new[,inds]>=lower_Gibbs & f_true_new[,inds]<=upper_Gibbs,2,mean)

coeff = max(l1_dist)
png("excircle_lsbpns_l1dist",width = 500, height = 500)
ggplot() +
  geom_line(aes(x=x_new,y=l1_dist)) +
  geom_point( aes(x=x_new[inds],y=ec_dens*coeff), color = "red",shape=8) +
  labs( x = "x", title = paste("Avg l1 distance:", round(mean(l1_dist),4)))+
  scale_y_continuous(
    name = "l1 distance",
    sec.axis = sec_axis(~./coeff, name="Empirical Coverage")
  ) +
  theme_bw()
dev.off()

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