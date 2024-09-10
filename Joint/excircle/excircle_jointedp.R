##############################################
## Data Generation for StatSci paper
##############################################
library(ggplot2)
library(tidyverse)
library(mvtnorm) 
library(msm)
library(MCMCpack)
library(mcclust.ext)

setwd("/Users/swade/Documents/GitHub/BayesXMix/Joint/excircle")

##############################################
## Example 1: skew-errors
##############################################

# Load data
load("../.././data_simulation/excircle/data_ex_circle.RData")

#############
#### Joint DP
#############
#Set prior parameter values

# y parameters
lmfit = lm(y~x[,1])
mu_theta=matrix(lmfit$coefficients,p+1,1)
a_y=2
# Note this is an overestimate of the variance within cluster (may want to divide by a smaller number)
k0y = 5
b_y=sum(lmfit$residuals^2)/(n-p-1)/k0y^2
# remember! variance for beta is sigma^2*iC
iC=diag(c(10,1))/b_y*var(y)
C=solve(iC)

# x parameters (empirical approach described in Fraley & Raftery (2007) and Wade (2023))
mu_0=matrix(apply(x,2,mean),p,1)
a_x=matrix(2,p,1)
# guess on the number of clusters
k0 = 6
b_x=matrix(apply(x,2,var),p,1)/k0^2 #with a_x this is the prior mean of the within cluster variance
c_x=rep(1/k0^2,p)

# prior on EDP concentration parameter
u_alpha_x=1
v_alpha_x=1
u_alpha_y=1
v_alpha_y=1


### Parameters for MCMC function
S=20000
burnin=5000
X = cbind(rep(1,n), x)

## Initial chain: Random Cluster assignments
ky_init=4
kx_init=4

config_y_init=apply(matrix(runif(n),nrow=n,ncol=(ky_init-1))>t(matrix(c(1:(ky_init-1))/ky_init,nrow=(ky_init-1),ncol=n)),1,sum)+1
config_x_init=apply(matrix(runif(n),nrow=n,ncol=(kx_init-1))>t(matrix(c(1:(kx_init-1))/kx_init,nrow=(kx_init-1),ncol=n)),1,sum)+1

alpha_y_init=u_alpha_y/v_alpha_y
alpha_x_init=rep(u_alpha_x/v_alpha_x,ky_init)

#initialize parameters for y and x
phi_init=matrix(0,p+1,ky_init)
sigma_y_init=matrix(0,1,ky_init)
mu_x_init=list()
sigma_x_init=list()
for(c in 1:ky_init){
  mu_x_init[[c]]=matrix(0, p, kx_init)
  sigma_x_init[[c]]=matrix(0,p,kx_init)
  
  #observations in class c
  n_c=sum(config_y_init==c)
  X_c=matrix(X[config_y_init==c,],nrow=n_c)
  y_c=y[config_y_init==c]
  
  #for intercept/slope
  C_hat=C+ t(X_c)%*%X_c
  decomp=eigen(C_hat, symmetric=TRUE)
  iC_hat=decomp$vectors%*%diag(1/decomp$values)%*%t(decomp$vectors)
  mu_hat=iC_hat%*%(C%*%mu_theta+t(X_c)%*%y_c)
  a_y_hat=a_y+n_c/2
  b_y_hat=b_y+(sum(y_c^2)+t(mu_theta)%*%C%*%mu_theta-t(mu_hat)%*%C_hat%*%mu_hat)/2
  
  #draw value from posterior
  sigma_y_init[c]=b_y_hat/(a_y_hat-1)
  phi_init[,c]=mu_hat
  
  
  for(h in 1:kx_init){
    
    #Draw a value for phi_y_c
    n_ch=sum(config_x_init[config_y_init==c]==h)
    x_ch=matrix(matrix(x[config_y_init==c,], nrow=n_c)[config_x_init[config_y_init==c]==h,], nrow=n_ch)
    
    #for mu_x
    c_x_hat= c_x+n_ch
    mux_hat= 1/c_x_hat*(c_x*mu_0+ colSums(x_ch))
    
    #for sigma_x
    a_x_hat=a_x+n_ch/2
    b_x_hat=b_x+(colSums(x_ch^2)+c_x*(mu_0^2)-c_x_hat*(mux_hat^2))/2
    
    #draw value from posterior
    sigma_x_init[[c]][,h]=b_x_hat/(a_x_hat-1)
    mu_x_init[[c]][,h]=mux_hat
  }
}

source(".././joint/edpyx_mcmc.R")

set.seed(101010)
output=edp_mcmc(S, burnin, y, x, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, u_alpha_x, v_alpha_x, u_alpha_y, v_alpha_y, config_y_init, config_x_init, phi_init, sigma_y_init, mu_x_init, sigma_x_init,alpha_y_init, alpha_x_init)

save.image("excircle_jointEDP.RData")

### Results

# Marginal posterior on the number clusters
config_x=output$config_x
config_y=output$config_y
k_y=rep(0,S)
k_x=list()
for(s in 1:S){
  k_y[s]=length(unique(config_y[s,]))
  k_x[[s]]=rep(0,k_y[s])
  for(i in 1:k_y[s]){
    k_x[[s]][i]=length(unique(config_x[s,config_y[s,]==i]))
  }
}
ggplot()+
  geom_bar(aes(k_y)) +
  theme_bw()

# Posterior similarity matrix for y-clustering
psm=comp.psm(config_y)
plotpsm(psm)

# y-Clustering estimate
output_vi=minVI(psm,config_y)
png("excircle_jedp_clus.png",width = 500, height = 450)
ggplot() +
  geom_point(aes(x = x[,1], y = y, color = as.factor(output_vi$cl))) +
  theme_bw() +
  labs( x = "x_1", y = "y", color = "Cluster") +
  geom_line(aes(x = x[,1], y = m_true), col ="red")
dev.off()

### PREDICTION
#Compute prediction and predictive density estimates with no credible intervals
source(".././joint/predict_edp.R")
output_pred=predict_edp(S,n_new,m2,p,x_new, y_grid, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k_y, k_x, config_y, config_x, output$alpha_x, output$alpha_y, output$phi_y, output$sigma_y, output$mu_x, output$sigma_x ) 

#Compute prediction with credible intervals
source(".././joint/predict_cred_edp.R")
output_cred_pred=predict_cred_edp(.05,S,n_new,p,x_new, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k_y, k_x, config_y, config_x, output$alpha_x, output$alpha_y, output$phi_y, output$sigma_y, output$mu_x, output$sigma_x ) 

#Compute predictive density with credible bounds
source(".././joint/predict_fcred_edp.R")
xs = sort(x_new[,1], index.return=TRUE)
inds = xs$ix[c(seq(1,n_new,n_new/5),n_new)]
output_fcred_pred=predict_fcred_edp(.05, S,length(inds),m2,p,x_new[inds,], y_grid, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k_y, k_x, config_y, config_x, output$alpha_x, output$alpha_y, output$phi_y, output$sigma_y, output$mu_x, output$sigma_x ) 

save.image("excircle_jointEDP.RData")

###Plot Prediction
#with credible intervals but without data
ggplot() +
  geom_point(aes(x = x_new[,1], y = output_cred_pred$y_pred), col = "black") +
  geom_point(aes(x = x_new[,1], y = m_true_new), col = "red") +
  geom_point(aes(x = x_new[,1], y = output_cred_pred$l_pred), col = "black", shape =25) +
  geom_point(aes(x = x_new[,1], y = output_cred_pred$u_pred), col = "black", shape =17) +
  theme_bw() +
  labs( x = "x", y = "y") +
  ylim(-1.1,1.1)

#with data but without credible intervals
ggplot() +
  geom_point(aes(x = x[,1], y = y)) +
  geom_line(aes(x = x_new[,1], y = output_pred$y_pred), col = "black") +
  geom_line(aes(x = x_new[,1], y = m_true_new), col = "red") +
  theme_bw() +
  labs( x = "x_1", y = "y")

#PLot density for a single observation
i = 1
ggplot() +
  geom_line(aes(x = y_grid, y = output_fcred_pred$f_pred[,i]), col = "black") +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[i]]), col = "red") +
  geom_ribbon(aes(x=y_grid, ymin=output_fcred_pred$l_fpred[,i], ymax=output_fcred_pred$u_fpred[,i]), alpha=0.2) +
  theme_bw() +
  labs( x = "y", y = "Density")

#PLot density for a few observations
png("excircle_jedp_dens",width = 500, height = 450)
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
  labs( x = "y", y = "Density")
dev.off()

#PLot true density heatmap.
png("excircle_jedp_dens_heat",width = 500, height = 450)
df = data.frame(x = rep(x_new[,1], each = m2), y = rep(y_grid,n_new), density= c(output_pred$f_pred))
ggplot(df) +
  geom_raster(aes(x,y,fill=density),interpolate = TRUE) +
  scale_fill_gradient2(low="white", high="red") +
  theme_classic() 
dev.off()

#empirical l2 prediction error
l2_err=sum(((m_true_new-output_pred$y_pred)^2)/n_new)^.5

l1_err=sum((abs(m_true_new-output_pred$y_pred))/n_new)

#estimated l1 distance for density
l1_dist=colSums(abs(f_true_new-output_pred$f_pred))*(y_grid[2]-y_grid[1])

png("excircle_jedp_l1dist",width = 500, height = 500)
ggplot() +
  geom_line(aes(x=x_new,y=l1_dist)) +
  labs( x = "x", y = "l1 distance", title = paste("Avg l1 distance:", round(mean(l1_dist),4)))+
  theme_bw()
dev.off()

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