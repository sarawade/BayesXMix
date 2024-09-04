##############################################
## Data Generation for StatSci paper
##############################################
library(ggplot2)
library(tidyverse)
library(mvtnorm) 
library(msm)
library(sn)

setwd("/Users/swade/Documents/GitHub/BayesXMix/data_simulation/ex1b")

##############################################
## Example 1: Drawback of Joint
##############################################

# Training Data
set.seed(1010101)
n = 800
p = 2

# Covariates
x= matrix(0,n,p)
x[,1] = runif(n,-1,8)
x[,2] = (x[,1]-3.5)^2/3-1+rnorm(n,0,0.05)

ggplot() +
  geom_point(aes(x = x[,1], y = x[,2])) +
  theme_bw() +
  xlab("x_1") +
  ylab("x_2")

# Response
omega = 0.1
alpha = x[,1]
xi = - omega*sqrt(2/pi)*alpha/sqrt(1+alpha^2)
y = 5-log(x[,1]+2) +rst(n, xi=xi,omega=omega,alpha=alpha)[1:n]
m_true = 5-log(x[,1]+2)

ggplot() +
  geom_point(aes(x = x[,1], y = y)) +
  theme_bw() +
  xlab("x_1") +
  ylab("y") +
  geom_line(aes(x = x[,1], y = m_true), col ="red")

### PREDICTION
n_new = 800

# Covariates
x_new= matrix(0,n_new,p)
x_new[,1] = runif(n_new,-1,8)
x_new[,2] = (x_new[,1]-3.5)^2/3-1+rnorm(n_new,0,0.05)

# True regression function
m_true_new = 5-log(x_new[,1]+2)

# True density
y_grid=seq(2.3,5.4,.02)
m2=length(y_grid)
alpha_new = x_new[,1]
xi_new = - omega*sqrt(2/pi)*alpha_new/sqrt(1+alpha_new^2)
f_true_new = matrix(0, nrow=m2,ncol=n_new)
for(i in 1:m2){
  f_true_new[i,] = dsn(y_grid[i]-m_true_new,xi=xi_new,omega=omega,alpha=alpha_new)
}

ggplot() +
  geom_point(aes(x = x_new[,1], y = x_new[,2])) +
  theme_bw() +
  xlab("x_1") +
  ylab("x_2")

ggplot() +
  geom_line(aes(x = x_new[,1], y = m_true_new), col ="red") +
  theme_bw() +
  xlab("x_1") +
  ylab("y") 

#Plot density for a few observations
xs = sort(x_new[,1], index.return=TRUE)
inds = xs$ix[c(seq(1,n_new,n_new/5),n_new)]
cols = rainbow(6)
ggplot() +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[1]]), col = cols[1],linetype = "dashed") +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[2]]), col = cols[2],linetype = "dashed") +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[3]]), col = cols[3],linetype = "dashed") +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[4]]), col = cols[4],linetype = "dashed") +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[5]]), col = cols[5],linetype = "dashed") +
  geom_line(aes(x = y_grid, y = f_true_new[,inds[6]]), col = cols[6],linetype = "dashed") +
  theme_bw() +
  labs( x = "y", y = "Density")

# Save data
save.image("data_ex1_skew.RData")

write.csv(cbind(x,y), file = "ex1train_skewerrors.csv", row.names = TRUE)
write.csv(x_new, file = "ex1test_skewerrors.csv", row.names = TRUE)
write.csv(f_true_new, file = "ex1test_skewerrors_dens.csv", row.names = TRUE)

