##############################################
## Data Generation for StatSci paper
##############################################
library(ggplot2)
library(tidyverse)
library(mvtnorm) 
library(msm)
library(MCMCpack)
library(mcclust.ext)

setwd("/Users/swade/Documents/GitHub/BayesXMix/NormalizedWeights/ex_circle")

##############################################
## Example 4: Implicit/Circle
##############################################

# Load data
load("../.././data_simulation/excircle/data_ex_circle.RData")

##############################################
## Results
##############################################
## See Normalized weights folder for matlab code
## to run the model with normalized weights

# Marginal posterior on the number clusters
config=as.matrix(read.csv("ex_config_nwreg_circle.csv", header = FALSE))
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
png("excircle_nw_clus.png",width = 500, height = 450)
ggplot() +
  geom_point(aes(x = x[,1], y = y, color = as.factor(output_vi$cl))) +
  theme_bw() +
  labs( x = "x_1", y = "y", color = "Cluster")
dev.off()

### PREDICTION
pred=as.matrix(read.csv("ex_pred_nwreg_circle.csv", header = FALSE))
fpred=as.matrix(read.csv("ex_fpred_nwreg_circle.csv", header = FALSE))
lfpred=as.matrix(read.csv("ex_lfpred_nwreg_circle.csv", header = FALSE))
ufpred=as.matrix(read.csv("ex_ufpred_nwreg_circle.csv", header = FALSE))

###Plot Prediction
#with credible intervals but without data
ggplot() +
  geom_point(aes(x = x_new[,1], y = pred[,1]), col = "black") +
  geom_point(aes(x = x_new[,1], y = m_true_new), col = "red") +
  geom_point(aes(x = x_new[,1], y = pred[,2]), col = "black", shape =25) +
  geom_point(aes(x = x_new[,1], y = pred[,3]), col = "black", shape =17) +
  theme_bw() +
  labs( x = "x_1", y = "y")

#with data but without credible intervals
ggplot() +
  geom_point(aes(x = x[,1], y = y)) +
  geom_line(aes(x = x_new[,1], y =pred[,1]), col = "black") +
  geom_line(aes(x = x_new[,1], y = m_true_new), col = "red") +
  theme_bw() +
  labs( x = "x_1", y = "y")


#PLot density for a single observation
xs = sort(x_new[,1], index.return=TRUE)
inds = xs$ix[c(seq(1,n_new,n_new/5),n_new)]
i = 1
ggplot() +
  geom_line(aes(x = y_grid, y = fpred[inds[i],]), col = "black") +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[i]]), col = "red") +
  geom_ribbon(aes(x=y_grid, ymin=lfpred[inds[i],], ymax=ufpred[inds[i],]), alpha=0.2) +
  theme_bw() +
  labs( x = "y", y = "Density")

#PLot density for a few observations
png("excircle_nw_fpred_skewerrors.png",width = 500, height = 450)
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

#PLot true density heatmap.
png("excircle_nw_dens_heat",width = 500, height = 450)
df = data.frame(x = rep(x_new[,1], each = m2), y = rep(y_grid,n_new), density= c(t(fpred)))
ggplot(df) +
  geom_raster(aes(x,y,fill=density),interpolate = TRUE) +
  scale_fill_gradient2(low="white", high="red") +
  theme_classic() 
dev.off()

#empirical l2 prediction error
l2_err=sum(((m_true_new-pred[,1])^2)/n_new)^.5

l1_err=sum((abs(m_true_new-pred[,1]))/n_new)

#estimated l1 distance for density
l1_dist=colSums(abs(f_true_new-t(fpred)))*(y_grid[2]-y_grid[1])

# Compute the empirical coverage of the pointwise intervals for the densities
ec_dens = apply(f_true_new>=t(lfpred) & f_true_new<=t(ufpred),2,mean)

coeff = max(l1_dist)
png("excircle_nw_l1dist",width = 500, height = 500)
ggplot() +
  geom_line(aes(x=x_new,y=l1_dist)) +
  geom_line( aes(x=x_new,y=ec_dens*coeff), color ='red') +
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
ec_pred =  sum((m_true_new>=pred[,2])&(m_true_new<=pred[,3]))/n_new

# credible interval length for prediction
ci_length =  pred[,3] - pred[,2]
mean(ci_length)
min(ci_length)
max(ci_length)
