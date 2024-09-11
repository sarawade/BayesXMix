##############################################
## Data Generation for StatSci paper
##############################################
library(ggplot2)
library(tidyverse)
library(mvtnorm) 
library(msm)

setwd("/Users/swade/Documents/GitHub/BayesXMix/data_simulation/ex3")

##############################################
## Example 3: Drawback of Conditional with dependent weights
##############################################

# Training Data
set.seed(1010101)
n = 600
p = 2

# Covariates
x= matrix(0,n,p)
x[,1] = runif(n,-2,2)
x[,2] = runif(n,-2,2)

ggplot() +
  geom_point(aes(x = x[,1], y = x[,2])) +
  theme_bw() +
  xlab("x_1") +
  ylab("x_2")

# Response
ind = sin(x[,1]*x[,2]*pi/2)<=0
m_true = (-1)*(1-ind) + ind
y = m_true +rnorm(n, 0, 0.1)

ggplot() +
  geom_point(aes(x = x[,1], y = x[,2], col = y)) +
  theme_bw() +
  xlab("x_1") +
  ylab("x_2") 

# Visualize on a grid
xgrid = seq(-2,2,.01)
dfgrid = data.frame(x1 = rep(xgrid,length(xgrid)), x2 = rep(xgrid,each = length(xgrid)) )
dfgrid$ind = sin(dfgrid$x1*dfgrid$x2*pi/2)<=0
dfgrid$m = (-1)*(1-dfgrid$ind)+ dfgrid$ind

# Mean function on a grid
png("ex3_true_mean_heat",width = 500, height = 450)
ggplot(dfgrid) +
  geom_raster(aes(x1,x2,fill=m),interpolate = TRUE) +
  theme_classic() +
  xlab("x_1") +
  ylab("x_2") 
dev.off()

# Mean function slice at fixed x2 value
x2value =1.5
ggplot() +
  geom_line(aes(x = dfgrid$x1[dfgrid$x2==x2value], y = dfgrid$m[dfgrid$x2==x2value]),col ="red") +
  theme_bw() +
  xlab("x_1") +
  ylab("y") 

### PREDICTION

# Covariates
xgrid = seq(-2,2,.1)
x_new= cbind(rep(xgrid,length(xgrid)),rep(xgrid,each = length(xgrid)))
n_new = dim(x_new)[1] 

# True regression function
ind_new = sin(x_new[,1]*x_new[,2]*pi/2) <=0
m_true_new = (-1)*(1-ind_new) + ind_new

# True density
y_grid=seq(-1.4,1.4,.02)
m2=length(y_grid)
f_true_new = dnorm(matrix(y_grid,nrow=m2,ncol=n_new),t(matrix(m_true_new,nrow=n_new,ncol=m2)), 0.1)

ggplot() +
  geom_point(aes(x = x_new[,1], y = x_new[,2], col = m_true_new)) +
  theme_bw() +
  xlab("x_1") +
  ylab("x_2")

df = data.frame(x1 =x_new[,1],x2=x_new[,2],y =m_true_new)
ggplot(df) +
  geom_raster(aes(x1,x2,fill= y),interpolate = TRUE) +
  theme_classic() 

# Save data
save.image("data_ex3.RData")

write.csv(cbind(x,y), file = "ex3train_terrors.csv", row.names = TRUE)
write.csv(x_new, file = "ex3test_terrors.csv", row.names = TRUE)
write.csv(f_true_new, file = "ex3test_terrors_dens.csv", row.names = TRUE)

