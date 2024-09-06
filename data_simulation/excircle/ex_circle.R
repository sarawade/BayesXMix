##############################################
## Data Generation for StatSci paper
##############################################
library(ggplot2)
library(tidyverse)
library(mvtnorm) 
library(msm)

setwd("/Users/swade/Documents/GitHub/BayesXMix/data_simulation/excircle")

##############################################
## Example 4: Implicit/Circle
##############################################

# Training Data
set.seed(1010101)
n = 800
p = 1
x= matrix(0,n,p)

# Covariates
t = runif(n,min=0,max=2*pi) 
x[,1] = cos(t)

ggplot() +
  geom_point(aes(x = t, y = x[,1])) +
  theme_bw() +
  xlab("t") +
  ylab("x")

# Response
y = sin(t) +rnorm(n,0,0.05)
m_true = sin(t)

ggplot() +
  geom_point(aes(x = x[,1], y = y)) +
  theme_bw() +
  xlab("x") +
  ylab("y") +
  geom_point(aes(x = x[,1], y = m_true), col ="red")

### PREDICTION

# Covariates
x_new = seq(-1,1,pi/(2^8))
n_new = length(x_new)
x_new= matrix(x_new,n_new,p)

# True regression function
m_true_new = rep(0,n_new)

# True density
y_grid=seq(-1.2,1.2,.02)
m2=length(y_grid)
f_true_new = 0.5*dnorm(matrix(y_grid,nrow=m2,ncol=n_new),t(matrix(sin(acos(x_new)),nrow=n_new,ncol=m2)), 0.05) + 0.5*dnorm(matrix(y_grid,nrow=m2,ncol=n_new),t(matrix(sin(acos(x_new)+pi),nrow=n_new,ncol=m2)), 0.05)

# Heatmap of true density
png("excircle_true_dens_heat",width = 500, height = 450)
f_pred_true = data.frame(x = rep(x_new[,1], each = m2), y = rep(y_grid,n_new), density= c(f_true_new))
ggplot(f_pred_true) +
  geom_raster(aes(x,y,fill=density),interpolate = TRUE) +
  scale_fill_gradient2(low="white", high="red") +
  theme_classic() 
dev.off()

# Save data
save.image("data_ex_circle.RData")

write.csv(cbind(x,y), file = "extrain_circle.csv", row.names = TRUE)
write.csv(x_new, file = "extest_circle.csv", row.names = TRUE)
write.csv(f_true_new, file = "extest_circle_dens.csv", row.names = TRUE)

