##############################################
## Data Generation for StatSci paper
##############################################
library(ggplot2)
library(tidyverse)
library(mvtnorm) 
library(msm)

setwd("/Users/swade/Documents/GitHub/BayesXMix/data_simulation/ex2")

##############################################
## Example 2: Drawbacks of single weights
##############################################

# Training Data
set.seed(1010101)
n = 400
p = 1
x= matrix(0,n,p)
y = rep(0,n)

# Covariates
x[,1] = runif(n,-2,10)

# Response 0 piecewise linear with covariate dependent variance
m  = rep(0,n)
m[x<=2] = 0
m[x>2 & x<=5] = 2*(x[x>2 & x<=5]-2)
m[x>5] = 6
eps = rep(0,n)
eps[x<=2] = rnorm(sum(x<=2), 0, 0.2)
eps[x>2 & x<=5] = rnorm(sum(x>2 & x<=5), 0, 0.05)
eps[x>5] = rnorm(sum(x>5), 0, (x[x>5]-5)^2/15+0.01)
y = m + eps

ggplot() +
  geom_point(aes(x = x[,1], y = y)) +
  theme_bw() +
  xlab("x") +
  ylab("y") +
  geom_point(aes(x = x[,1], y = m), col ="red")

### PREDICTION

# Covariates
n_new = 1201
x_new= matrix(0,n_new,p)
x_new[,1] = seq(-2,10,12/(n_new-1))

# True regression function
m_true_new = rep(0,n_new)
m_true_new[x_new>2 & x_new<=5] = 2*(x_new[x_new>2 & x_new<=5]-2)
m_true_new[x_new>5] = 6

ggplot() +
  geom_line(aes(x = x_new, y = m_true_new), col ="red") +
  theme_bw() +
  xlab("x_1") +
  ylab("y") 

# True density
y_grid=seq(-1,10,.02)
m2=length(y_grid)
f_true_new = matrix(0,m2,n_new)
f_true_new[,x_new<=2] = dnorm(matrix(y_grid,nrow=m2,ncol=sum(x_new<=2)),t(matrix(m_true_new[x_new<=2],nrow=sum(x_new<=2),ncol=m2)), 0.2)
f_true_new[,x_new>2 & x_new<=5] = dnorm(matrix(y_grid,nrow=m2,ncol=sum(x_new>2 & x_new<=5)),t(matrix(m_true_new[x_new>2 & x_new<=5],nrow=sum(x_new>2 & x_new<=5),ncol=m2)), 0.05)
f_true_new[,x_new>5] = dnorm(matrix(y_grid,nrow=m2,ncol=sum(x_new>5)),t(matrix(m_true_new[x_new>5],nrow=sum(x_new>5),ncol=m2)), t(matrix((x_new[x_new>5]-5)^2/15+0.01,nrow=sum(x_new>5),ncol=m2)))

# Heatmap of true density
png("excircle_true_dens_heat",width = 500, height = 450)
f_pred_true = data.frame(x = rep(x_new, each = m2), y = rep(y_grid,n_new), density= c(f_true_new))
ggplot(f_pred_true) +
  geom_raster(aes(x,y,fill=density),interpolate = TRUE) +
  scale_fill_gradient2(low="white", high="red", trans = "sqrt") +
  theme_classic() 
dev.off()

# Save data
save.image("data_ex2.RData")

write.csv(cbind(x,y), file = "ex2train.csv", row.names = TRUE)
write.csv(x_new, file = "ex2test.csv", row.names = TRUE)
write.csv(f_true_new, file = "ex2test_dens.csv", row.names = TRUE)

