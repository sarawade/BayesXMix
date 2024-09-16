##############################################
## Data Generation for StatSci paper
##############################################
library(ggplot2)
library(tidyverse)
library(mvtnorm) 
library(msm)
library(MCMCpack)
library(mcclust.ext)

setwd("/Users/swade/Documents/GitHub/BayesXMix/Joint/ex1")

##############################################
## Example 1: normal errors
##############################################

# Load data
load("../.././data_simulation/ex1/data_ex1.RData")

#############
#### Joint DP
#############
#Set prior parameter values

# y parameters
lmfit = lm(y~x[,1]+x[,2])
mu_theta=matrix(lmfit$coefficients,p+1,1)
a_y=2
# Note this is an overestimate of the variance within cluster (may want to divide by a smaller number)
b_y=sum(lmfit$residuals^2)/(n-p-1)
# remember! variance for beta is sigma^2*iC
iC=diag(c(10,1,1))/b_y*var(y)
C=solve(iC)

# x parameters (empirical approach described in Fraley & Raftery (2007) and Wade (2023))
mu_0=matrix(apply(x,2,mean),p,1)
a_x=matrix(2,p,1)
# guess on the number of clusters
k0 = 6
b_x=matrix(apply(x,2,var),p,1)/k0^2 #with a_x this is the prior mean of the within cluster variance
c_x=rep(1/k0^2,p)

# prior on DP concentration parameter
u_alpha=1
v_alpha=1


### Parameters for MCMC function
S=20000
burnin=5000
X = cbind(rep(1,n), x)

## Initialize chain

## Random Cluster assignments
k_init=10

config_init=apply(matrix(runif(n),nrow=n,ncol=(k_init-1))>t(matrix(c(1:(k_init-1))/k_init,nrow=(k_init-1),ncol=n)),1,sum)+1

alpha_init=u_alpha/v_alpha

#initialize parameters for y and x
phi_init=matrix(0,p+1,k_init)
sigma_y_init=matrix(0,1,k_init)
mu_x_init=matrix(0,p,k_init)
sigma_x_init=matrix(0,p,k_init)
for(c in 1:k_init){
  
  #observations in class c
  n_c=sum(config_init==c)
  X_c=matrix(X[config_init==c,],nrow=n_c)
  x_c=matrix(x[config_init==c,],nrow=n_c)
  y_c=y[config_init==c]
  
  #for intercept/slope
  C_hat=C+ t(X_c)%*%X_c
  decomp=eigen(C_hat, symmetric=TRUE)
  iC_hat=decomp$vectors%*%diag(1/decomp$values)%*%t(decomp$vectors)
  mu_hat=iC_hat%*%(C%*%mu_theta+t(X_c)%*%y_c)
  a_y_hat=a_y+n_c/2
  b_y_hat=b_y+(sum(y_c^2)+t(mu_theta)%*%C%*%mu_theta-t(mu_hat)%*%C_hat%*%mu_hat)/2
  
  #draw value from posterior
  sigma_y_init[c]=rinvgamma(1,a_y_hat, b_y_hat)
  phi_init[,c]=rmvnorm(1, mu_hat, sigma_y_init[c]*iC_hat,  method=c("chol"))
  
  #for mu_x
  c_x_hat= c_x+n_c
  mux_hat= 1/c_x_hat*(c_x*mu_0+ colSums(x_c))
  
  #for sigma_x
  a_x_hat=a_x+n_c/2
  b_x_hat=b_x+(colSums(x_c^2)+c_x*(mu_0^2)-c_x_hat*(mux_hat^2))/2
  
  #draw value from posterior
  sigma_x_init[,c]=rinvgamma(p,a_x_hat, b_x_hat)
  mu_x_init[,c]=rnorm(p,mux_hat, (sigma_x_init[,c]/c_x_hat)^.5)
}

source(".././joint/jdp_mcmc.R")

set.seed(101010)
output=jdp_mcmc(S, burnin, y, x, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, u_alpha, v_alpha, config_init, phi_init, sigma_y_init, mu_x_init, sigma_x_init,alpha_init)

save.image("ex1_jointDP.RData")

### Results

# Marginal posterior on the number clusters
config=output$config
k=rep(0,S)
for(s in 1:S){
  k[s]=length(unique(config[s,]))
}
ggplot()+
  geom_bar(aes(k)) +
  theme_bw()

# Posterior similarity matrix
psm=comp.psm(config)
plotpsm(psm)

# Clustering estimate
output_vi=minVI(psm,config)
png("ex1_jdp_clus1.png",width = 500, height = 450)
ggplot() +
  geom_point(aes(x = x[,1], y = x[,2], color = as.factor(output_vi$cl))) +
  theme_bw() +
  labs( x = "x_1", y = "x_2", color = "Cluster")
dev.off()
png("ex1_jdp_clus2.png",width = 500, height = 450)
ggplot() +
  geom_point(aes(x = x[,1], y = y, color = as.factor(output_vi$cl))) +
  theme_bw() +
  labs( x = "x_1", y = "y", color = "Cluster") +
  geom_line(aes(x = x[,1], y = m_true), col ="red")
dev.off()

### PREDICTION
#Compute prediction and predictive density estimates with no credible intervals
source(".././joint/predict_jdp.R")
output_pred=predict_jdp(S,n_new,m2,p,x_new, y_grid, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k, config, output$alpha, output$phi_y, output$sigma_y, output$mu_x, output$sigma_x )

#Compute prediction with credible intervals
source(".././joint/predict_cred_jdp.R")
output_cred_pred=predict_cred_jdp(.05,S,n_new,p,x_new, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k, output$config, output$alpha, output$phi_y, output$sigma_y, output$mu_x, output$sigma_x )

#Compute predictive density with credible bounds
source(".././joint/predict_fcred_jdp.R")
inds = c(1,2,201,202,401,402,601,602)
output_fcred_pred=predict_fcred_jdp(.05, S,length(inds),m2,p,x_new[inds,], y_grid, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k, config, output$alpha, output$phi_y, output$sigma_y, output$mu_x, output$sigma_x )

save.image("ex1_jointDP.RData")

###Plot Prediction
#with credible intervals but without data
png("ex1_jdp_pred.png",width = 500, height = 450)
ggplot() +
  geom_line(aes(x = x_new[,1], y = output_cred_pred$y_pred), col = "black") +
  geom_line(aes(x = x_new[,1], y = m_true_new), col = "red") +
  geom_ribbon(aes(x = x_new[,1], ymin=output_cred_pred$l_pred, ymax=output_cred_pred$u_pred), alpha=0.2) +
  theme_bw() +
  labs( x = "x_1", y = "y")
dev.off()

#with data but without credible intervals
ggplot() +
  geom_point(aes(x = x[,1], y = y)) +
  geom_line(aes(x = x_new[,1], y = output_pred$y_pred), col = "black") +
  geom_line(aes(x = x_new[,1], y = m_true_new), col = "red") +
  theme_bw() +
  labs( x = "x_1", y = "y")

#Plot y_hat_true vs y_pred
ggplot() +
  geom_point(aes(x = m_true_new, y = output_cred_pred$y_pred), col = "black") +
  geom_point(aes(x = m_true_new, y = output_cred_pred$l_pred), col = "black", shape = 25) +
  geom_point(aes(x = m_true_new, y = output_cred_pred$u_pred), col = "black", shape = 24) +
  geom_line(aes(x = m_true_new, y = m_true_new), col = "grey") +
  theme_bw() +
  labs( x = "True y", y = "Estimated y")

#PLot density for a single observation
i = 1
ggplot() +
  geom_line(aes(x = y_grid, y = output_fcred_pred$f_pred[,i]), col = "black") +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[i]]), col = "red") +
  geom_ribbon(aes(x=y_grid, ymin=output_fcred_pred$l_fpred[,i], ymax=output_fcred_pred$u_fpred[,i]), alpha=0.2) +
  theme_bw() +
  labs( x = "y", y = "Density")

#PLot density for a few observations
png("ex1_jdp_fpred.png",width = 500, height = 450)
cols = rainbow(6)
ggplot() +
  geom_line(aes(x = y_grid, y = output_fcred_pred$f_pred[,1]), col = cols[1]) +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[1]]), col = cols[1],linetype = "dashed") +
  geom_ribbon(aes(x=y_grid, ymin=output_fcred_pred$l_fpred[,1], ymax=output_fcred_pred$u_fpred[,1]), alpha=0.2, fill = cols[1]) +
  geom_line(aes(x = y_grid, y = output_fcred_pred$f_pred[,3]), col = cols[3]) +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[3]]), col = cols[3],linetype = "dashed") +
  geom_ribbon(aes(x=y_grid, ymin=output_fcred_pred$l_fpred[,3], ymax=output_fcred_pred$u_fpred[,3]), alpha=0.2, fill = cols[3]) +
  geom_line(aes(x = y_grid, y = output_fcred_pred$f_pred[,5]), col = cols[5]) +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[5]]), col = cols[5],linetype = "dashed") +
  geom_ribbon(aes(x=y_grid, ymin=output_fcred_pred$l_fpred[,5], ymax=output_fcred_pred$u_fpred[,5]), alpha=0.2, fill = cols[5]) +
  theme_bw() +
  labs( x = "y", y = "Density")  +
  xlim(2.4,4.2)
dev.off()

#empirical l2 prediction error
l2_err=sum(((m_true_new-output_pred$y_pred)^2)/n_new)^.5

l1_err=sum((abs(m_true_new-output_pred$y_pred))/n_new)

#estimated l1 distance for density
l1_dist=colSums(abs(f_true_new-output_pred$f_pred))*(y_grid[2]-y_grid[1])

#Average l1 distance
mean(l1_dist)

#Max l1 dist
max(l1_dist)

#Min l1 dist
min(l1_dist)

# empirical coverage for prediction
ec_pred = sum((m_true_new>=output_cred_pred$l_pred)&(m_true_new<=output_cred_pred$u_pred))/n_new

# credible interval length for prediction
ci_length = output_cred_pred$u_pred - output_cred_pred$l_pred
mean(ci_length)
min(ci_length)
max(ci_length)