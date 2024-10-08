##############################################
## Data Generation for StatSci paper
##############################################
library(ggplot2)
library(tidyverse)
library(mvtnorm) 
library(msm)
library(MCMCpack)
library(mcclust.ext)

setwd("/Users/swade/Documents/GitHub/BayesXMix/NormalizedWeights/ex3")

##############################################
## Example 3: Drawback of Conditional with dependent weights
##############################################

# Load data
load("../.././data_simulation/ex3/data_ex3.RData")

##############################################
## Results
##############################################
## See Normalized weights folder for matlab code
## to run the model with normalized weights

# Marginal posterior on the number clusters
config=as.matrix(read.csv("ex3_config_nwreg.csv", header = FALSE))
S = dim(config)[1]
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
png("ex3_nw_clus.png",width = 500, height = 450)
ggplot() +
  geom_point(aes(x = x[,1], y = x[,2], color = as.factor(output_vi$cl))) +
  theme_bw() +
  labs( x = "x_1", y = "x_2", color = "Cluster") 
dev.off()

### PREDICTION
pred=as.matrix(read.csv("ex3_pred_nwreg.csv", header = FALSE))
fpred=as.matrix(read.csv("ex3_fpred_nwreg.csv", header = FALSE))
lfpred=as.matrix(read.csv("ex3_lfpred_nwreg.csv", header = FALSE))
ufpred=as.matrix(read.csv("ex3_ufpred_nwreg.csv", header = FALSE))

inds = seq(20,1681,41*8)

###Plot Prediction

#Heat map of mean
png("ex3_nw_pred_heat.png",width = 500, height = 450)
df = data.frame(x1 =x_new[,1],x2=x_new[,2],y = pred[,1])
ggplot(df) +
  geom_raster(aes(x1,x2,fill= y),interpolate = TRUE) +
  theme_classic() 
dev.off()

#with credible intervals but without data
png("ex3_nw_pred.png",width = 500, height = 450)
xnew2 = 1.5
ggplot() +
  geom_point(aes(x = x_new[x_new[,2]==xnew2,1], y = pred[x_new[,2]==xnew2,1]), col = "black") +
  geom_point(aes(x = x_new[x_new[,2]==xnew2,1], y = m_true_new[x_new[,2]==xnew2]), col = "red") +
  geom_ribbon(aes(x = x_new[x_new[,2]==xnew2,1],ymin=pred[x_new[,2]==xnew2,2], ymax=pred[x_new[,2]==xnew2,3]), alpha=0.2) +
  theme_bw() +
  labs( x = "x_1", y = "y")
dev.off()

#Plot y_hat_true vs y_pred
ggplot() +
  geom_point(aes(x = m_true_new, y = pred[,1]), col = "black") +
  geom_point(aes(x = m_true_new, y = pred[,2]), col = "black", shape = 25) +
  geom_point(aes(x = m_true_new, y = pred[,3]), col = "black", shape = 24) +
  geom_line(aes(x = m_true_new, y = m_true_new), col = "grey") +
  theme_bw() +
  labs( x = "True y", y = "Estimated y")

#PLot density for a single observation
i = 1
ggplot() +
  geom_line(aes(x = y_grid, y = fpred[inds[i],]), col = "black") +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[i]]), col = "red") +
  geom_ribbon(aes(x=y_grid, ymin=lfpred[inds[i],], ymax=ufpred[inds[i],]), alpha=0.2) +
  theme_bw() +
  labs( x = "y", y = "Density")

#PLot density for a few observations
png("ex3_nw_fpred.png",width = 500, height = 450)
cols = rainbow(6)
ggplot() +
  geom_line(aes(x = y_grid, y = fpred[inds[1],]), col = cols[1]) +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[1]]), col = cols[1],linetype = "dashed") +
  geom_ribbon(aes(x=y_grid, ymin=lfpred[inds[1],], ymax=ufpred[inds[1],]), alpha=0.2, fill = cols[1]) +
  geom_line(aes(x = y_grid, y = fpred[inds[3],]), col = cols[3]) +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[3]]), col = cols[3],linetype = "dashed") +
  geom_ribbon(aes(x=y_grid, ymin=lfpred[inds[3],], ymax=ufpred[inds[3],]), alpha=0.2, fill = cols[3]) +
  geom_line(aes(x = y_grid, y = fpred[inds[5],]), col = cols[5]) +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[5]]), col = cols[5],linetype = "dashed") +
  geom_ribbon(aes(x=y_grid, ymin=lfpred[inds[5],], ymax=ufpred[inds[5],]), alpha=0.2, fill = cols[5]) +
  theme_bw() +
  labs( x = "y", y = "Density") 
dev.off()

#empirical l2 prediction error
l2_err=sum(((m_true_new-pred[,1])^2)/n_new)^.5

l1_err=sum((abs(m_true_new-pred[,1]))/n_new)

#estimated l1 distance for density
l1_dist=colSums(abs(f_true_new-t(fpred)))*(y_grid[2]-y_grid[1])

#Heat map of l1 distance
png("ex3_nw_l1_heat.png",width = 500, height = 450)
df = data.frame(x1 =x_new[,1],x2=x_new[,2],l1 =l1_dist)
ggplot(df) +
  geom_raster(aes(x1,x2,fill= l1),interpolate = TRUE) +
  theme_classic() 
dev.off()

#Average l1 distance
mean(l1_dist)

#Max l1 dist
max(l1_dist)

#Min l1 dist
min(l1_dist)


# empirical coverage for prediction
ec_pred =  sum((m_true_new>=pred[,2])&(m_true_new<=pred[,3]))/n_new

# credible interval length for prediction
ci_length =  pred[,3] - pred[,2]
mean(ci_length)
min(ci_length)
max(ci_length)
