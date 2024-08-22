#joint DP mixture for regression
#Simulate data and call jDP MCMC function to obtain
#posterior samples

library(mvtnorm) 
library(msm)
library(MCMCpack)


## Add increasing number of covariates

###################################################################

#SIMULATE DATA

set.seed(1010)
p=1
n=200
#probx=c(1/2,1/2)
#configx_t=(runif(n)<probx[1])+1
#mux_t=c(3,7)
#sdx_t=c(1,1)
#x=matrix(rnorm(n),n,1)*(sdx_t[1]*(configx_t==1)+sdx_t[2]*(configx_t==2))+(mux_t[1]*(configx_t==1)+mux_t[2]*(configx_t==2))
x=matrix(rnorm(n),n,1)*2+4
#two clusters
prob=c(1/2,1/2)
mu_t=c(4,6)
sd_t=c(1,1)/2
px_t=matrix(0,n,2)
px_t[,1]=dnorm(x,mu_t[1],sd_t[1])*prob[1]
px_t[,2]=dnorm(x,mu_t[2],sd_t[2])*prob[2]
std=function(x){
z=x/sum(x)
return(z)
}
px_t=t(apply(px_t,1,std))
u=runif(n)
config_true=rep(1,n)
config_true[u>px_t[,1]]=2
beta_t=matrix(c(0,1,5*.9,.1),2,2)
X=cbind(rep(1,n), x)
y=matrix(0,p,1)
sigma1=1/16
sigma2=1/8
y[config_true==1]=rnorm(sum(config_true==1),X[config_true==1,]%*%beta_t[,1],sigma1^.5)
y[config_true==2]=rnorm(sum(config_true==2),X[config_true==2,]%*%beta_t[,2],sigma2^.5)


#######################################################################
#plot data

par(mfrow=c(1,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
plot(x[,1],y,xlab="X",ylab="Y")
points(x[config_true==2,1],y[config_true==2],col=2)

par(mfrow=c(1,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
sx=sort(x[,1],index.return=T)
plot(sx$x,px_t[sx$ix,1],xlab="X",ylab="Y",type="l")
lines(sx$x,px_t[sx$ix,2],col=2)

####################################################################

#sample additional x
pmax=100-1
#configxa_t=matrix((runif(n*pmax)<probx[1])+1,nrow=n,ncol=pmax)
#xa=matrix(rnorm(n*pmax),n,pmax)*(sdx_t[1]*(configxa_t==1)+sdx_t[2]*(configxa_t==2))+(mux_t[1]*(configxa_t==1)+mux_t[2]*(configxa_t==2))
if(pmax%%2==1){
sigmax2=matrix(c(3.5,0),pmax,pmax)
}
if(pmax%%2==0){
sigmax2=matrix(c(3.5,0),pmax+1,pmax+1)
sigmax2=sigmax2[-(pmax+1),-(pmax+1)]
}
diag(sigmax2)=rep(4,pmax)
sigmax21=sigmax2-1/4*matrix(c(3.5,0),pmax,1)%*%matrix(c(3.5,0),1,pmax)
mu21=4+1/4*(x[,1]-4)%x%t(matrix(c(3.5,0),pmax,1))
csigmax21=chol(sigmax21)
xa=t(t(csigmax21)%*%matrix(rnorm(pmax*n),pmax,n))+mu21

#add to x
pa=15-1
p=p+pa
if(pa>0){
x=cbind(x,xa[,1:pa])}
X=cbind(rep(1,n), x)

#####################################################################

#####################################################################

#Set parameter values

#y parameters
#remember! variance for beta is sigma^2*iC
mu_theta=apply(beta_t,1,mean)
if(pa>0){
mu_theta=matrix(c(mu_theta,rep(0,pa)),nrow=(p+1))
}
#Option 1
#iXX=solve(t(X)%*%X)
#iC.corr=diag(1/sqrt(diag(iXX)))%*%iXX%*%diag(1/sqrt(diag(iXX)))
#iC=diag(sqrt(iC.diag))%*%iC.corr%*%diag(sqrt(iC.diag))
#Option 2
iC.diag=10*c((2.5^2)/(qnorm(.5+.5*.95^(1/(p+1)))^2),rep((.5^2)/(qnorm(.5+.5*.95^(1/(p+1)))^2),p))
iC=diag(iC.diag)
C=solve(iC)
a_y=2
b_y=1/10


# x parameters
a_x=matrix(2,p,1)
b_x=matrix(1,p,1)
mu_0=matrix(4,p,1)
c_x=rep((qnorm(.5+.5*.95^(1/p))^2)/(4^2),p)

#Precision
u_alpha=1
v_alpha=1

#######################################################################

### Parameters for MCMC function
S=20000
burnin=5000

#Initialize chain

#Option 1
#start everyone in 1 group
#config_init=rep(1,n)

#alpha_init=u_alpha/v_alpha

#compute updated parameters for y
#C_hat=C+ t(X)%*%X
#decomp=eigen(C_hat, symmetric=TRUE)
#iC_hat=decomp$vectors%*%diag(1/decomp$values)%*%t(decomp$vectors)
#mu_hat=iC_hat%*%(C%*%mu_theta+t(X)%*%y)
#a_y_hat=a_y+n/2
#b_y_hat=b_y+(sum(y^2)+t(mu_theta)%*%C%*%mu_theta-t(mu_hat)%*%C_hat%*%mu_hat)/2

#initialize parameters for y
#phi_init=mu_hat
#sigma_y_init=b_y_hat/(a_y_hat-1)


#compute updated hyperparameters for x
#c_x_hat=c_x+n
#mu_x_hat=1/c_x_hat*(c_x*mu_0+colSums(x))
#a_x_hat=a_x+n/2
#b_x_hat=b_x+1/2*(colSums(x^2)+c_x*mu_0^2-c_x_hat*mu_x_hat^2)

#initialize parameters for x
#mu_x_init=mu_x_hat
#sigma_x_init=b_x_hat/(a_x_hat-1)

#Option 2: Load previous data
#config_prev=config
#S_prev=S
#alpha_prev=alpha
#phi_prev=phi_y
#sigma_y_prev=sigma_y
#mu_x_prev=mu_x
#sigma_x_prev=sigma_x
#k_prev=k

#config_init=config_prev[S_prev,]
#phi_init=phi_prev[[S_prev]]
#sigma_y_init=sigma_y_prev[[S_prev]]
#mu_x_init=mu_x_prev[[S_prev]]
#sigma_x_init=sigma_x_prev[[S_prev]]
#alpha_init=alpha_prev[S_prev]

## Option 3: Random Cluster assignments
k_init=3

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


#####################################################################

source("./Dropbox/EDP_mix/Ex_ip/jdp_mcmc.R")

#Call MCMC edp function
set.seed(101010)
output=jdp_mcmc(S, burnin, y, x, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, u_alpha, v_alpha, config_init, phi_init, sigma_y_init, mu_x_init, sigma_x_init,alpha_init)

save.image("./Dropbox/EDP_mix/Ex_ip/jdp_ex_ip15_S20000b5000_v2_n200.RData")

#######################################################################

# Define variables
config=output$config
phi_y=output$phi_y
sigma_y=output$sigma_y
mu_x=output$mu_x
sigma_x=output$sigma_x
alpha=output$alpha


# If previous data, attach new results
#config=rbind(config_prev, config)
#phi_y_combined=list()
#phi_y_combined[1:S_prev]=phi_prev
#phi_y_combined[(S_prev+1):(S_prev+S)]=phi_y
#phi_y=phi_y_combined
#sigma_y_combined=list()
#sigma_y_combined[1:S_prev]=sigma_y_prev
#sigma_y_combined[(S_prev+1):(S_prev+S)]=sigma_y
#sigma_y=sigma_y_combined
#mu_x_combined=list()
#mu_x_combined[1:S_prev]=mu_x_prev
#mu_x_combined[(S_prev+1):(S_prev+S)]=mu_x
#mu_x=mu_x_combined
#sigma_x_combined=list()
#sigma_x_combined[1:S_prev]=sigma_x_prev
#sigma_x_combined[(S_prev+1):(S_prev+S)]=sigma_x
#sigma_x=sigma_x_combined
#alpha=c(alpha_prev, alpha)
#S=S+S_prev

k=rep(0,S)
for(s in 1:S){
	k[s]=length(unique(config[s,]))
}

#######################################################################

#Diagnostics
#choose individual
i=1
phi_hat_i=matrix(0,S,p+1)
phi_i=matrix(0,S, p+1)
phi_i[1,]=phi_y[[1]][,config[1,i]]
phi_hat_i[1,]=(phi_i[1,])
for(s in 2:S){
phi_i[s,]=phi_y[[s]][,config[s,i]]
phi_hat_i[s,]=(phi_hat_i[s-1,]*(s-1)+phi_i[s,])/s
}

par(mfrow=c(2,1),mar=c(1.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
for(j in 1:(2)){
plot(phi_i[,j], ylab=paste("beta",j-1), xlab="", main=paste("beta*",j-1,"=",round(phi_hat_i[s,j],4)) )
points(phi_hat_i[,j],col=2)
}

par(mfrow=c(1,1),mar=c(1.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
j=2
plot(phi_i[,j], ylab="beta_2",xlab="", main=paste("beta*_2=",round(phi_hat_i[s,j],4)))
points(phi_hat_i[,j],col=2)

round(phi_hat_i[s,],4)


#for sigma_y
i=5
sigmay_hat_i=matrix(0,S,1)
sigmay_i=matrix(0,S, 1)
sigmay_i[1]=sigma_y[[1]][config[1,i]]
sigmay_hat_i[1]=(sigmay_i[1])
for(s in 2:S){
sigmay_i[s]=sigma_y[[s]][config[s,i]]
sigmay_hat_i[s]=(sigmay_hat_i[s-1]*(s-1)+sigmay_i[s])/s
}

par(mfrow=c(1,1),mar=c(1.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
plot(sigmay_i, ylab="sigmay", xlab="", main=paste("sigmay*=",round(sigmay_hat_i[s],4)))
points(sigmay_hat_i,col=2)

round(sigmay_hat_i[s],4)


par(mfrow=c(3,1),mar=c(1.5,2.5,1.5,.6)+.1, mgp=c(1.5, .5, 0))
for(j in 1:(p+1)){
plot(phi_i[,j], ylab=paste("beta",j-1), xlab="", main=paste("beta*",j-1,"=",round(phi_hat_i[s,j],4)) )
points(phi_hat_i[,j],col=2)
}
plot(sigmay_i, ylab="sigmay", xlab="", main=paste("sigmay*=",round(sigmay_hat_i[s],4)))
points(sigmay_hat_i,col=2)


#### Compute subject specific parameter estimates
sigmay_hat=matrix(0,n,1)
phi_hat=matrix(0,n,p+1)
for (i in 1:n){
phi_hat[i,]=phi_y[[1]][,config[1,i]]
sigmay_hat[i]=sigma_y[[1]][config[1,i]]
for(s in 2:S){
phi_hat[i,]=phi_hat[i,]+phi_y[[s]][,config[s,i]]
sigmay_hat[i]=sigmay_hat[i]+sigma_y[[s]][config[s,i]]
}
}
phi_hat=phi_hat/S
sigmay_hat=sigmay_hat/S
beta_tf=beta_t
if(pa>0){beta_tf=rbind(beta_t,matrix(0,pa,2))}

# Compute true
phi_true=((config_true==1)%x%t(beta_tf[,1])+(config_true==2)%x%t(beta_tf[,2]))
sigmay_true=((config_true==1)*sigma1+(config_true==2)*sigma2)

# Print estimated and true for some subjects
toprint=data.frame(beta_true=phi_true[1:3,], beta_hat=phi_hat[c(1:3),])
toprint
toprint=data.frame(sigmay_true=sigmay_true[1:3], sigmay_hat=sigmay_hat[c(1:3)])
toprint
avg_diff_phi=colSums(abs(phi_hat-phi_true))/n
avg_diff_sigma=sum(abs(sigmay_hat-sigmay_true))/n
toprint=data.frame(Avg_err_beta=matrix(avg_diff_phi,1,p+1),Avg_err_sigmay=avg_diff_sigma)
toprint

## for mu
i=5
mu_hat_i=matrix(0,S,p)
mu_i=matrix(0,S, p)
mu_i[1,]=mu_x[[1]][,config[1,i]]
mu_hat_i[1,]=(mu_i[1,])
for(s in 2:S){
mu_i[s,]=mu_x[[s]][,config[s,i]]
mu_hat_i[s,]=(mu_hat_i[s-1,]*(s-1)+mu_i[s,])/s
}

par(mfrow=c(1,1),mar=c(1.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
for(j in 1:p){
plot(mu_i[,j],ylab=paste("mu",j), xlab="", main=paste("mu*",j,"=",round(mu_hat_i[s,j],4)))
points(mu_hat_i[,j],col=2)
}


round(mu_hat_i[s,],4)

## for sigma
i=5
sigmax_hat_i=matrix(0,S,p)
sigmax_i=matrix(0,S, p)
sigmax_i[1,]=sigma_x[[1]][,config[1,i]]
sigmax_hat_i[1,]=(sigmax_i[1,])
for(s in 2:S){
sigmax_i[s,]=sigma_x[[s]][,config[s,i]]
sigmax_hat_i[s,]=(sigmax_hat_i[s-1,]*(s-1)+sigmax_i[s,])/s
}


par(mfrow=c(1,1),mar=c(1.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
for(j in 1:p){
plot(sigmax_i[,j],ylab=paste("sigmax",j), xlab="", main=paste("sigmax*",j,"=",round(sigmax_hat_i[s,j],4)))
points(sigmax_hat_i[,j],col=2)
}

round(sigmax_hat_i[s,],4)



#####################################################################


#Analysis of parameters

##k
k_hat=mean(k)
k_hat
par(mfrow=c(1,1))
hist(k, freq=F)
k_hat_round=round(k_hat)
median(k)
sk=sort(k)
sk[0.025*S]
sk[0.975*S]

#alpha
alpha_hat=mean(alpha)
alpha_hat
plot(density(alpha))
sa=sort(alpha)
sa[0.025*S]
sa[0.975*S]



######################################################################

########### Inference for partitions


# Need to Reorder labels
source("./Dropbox/EDP_mix/Ex_ip/reorder_labels_jdp.R")
output_reorder=reorder_labels(S,n, config)

save.image("./Dropbox/EDP_mix/Ex_ip/jdp_ex_ip15_S20000b5000_v2_n200.RData")

#Compute estimated probability of configurations and sort by highest prob
config_p=output_reorder$config_count/S
config_sp=sort(config_p, decreasing=TRUE, index.return=TRUE)
#Total number of comfigurations with positive prob
length(config_p)

#probability of top 8 configurations
config_sp$x[1:8]
sum(config_sp$x[1:8])


#configuration of top
config_top=config[output_reorder$config_index[config_sp$ix],]

#any configurations equal to true configuration?
if(sum(apply(config_top,1,max)==2)>0){
configy_true=rowSums(config_top[apply(config_top,1,max)==2,]!=t(matrix(config_true,n,sum(apply(config_top,1,max)==2))))
sum(configy_true==n)
sum(configy_true==0)
max(configy_true)
min(configy_true)
}


########Plot top models
cl=c("black", "red", "green", "blue", "yellow", "orange", "magenta", "cyan", "gray","darkgray","blueviolet","darkmagenta","darkred", "darkorange", "darkorange4", "cornflowerblue", "darkolivegreen2", "darkolivegreen4", "brown","chartreuse","burlywood1" ) 

## TOP IN X SPACE 

#Plot top model in X space
par(mfrow=c(2,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
plot(x[,1],x[,2],xlab="x1",ylab="x2")
k_1=max(config_top[1,])
for(j in 1:k_1){
points(x[config_top[1,]==j,1], x[config_top[1,]==j,2],col=cl[j])
}
plot(x[,3],x[,4],xlab="x3",ylab="x4")
k_1=max(config_top[1,])
for(j in 1:k_1){
points(x[config_top[1,]==j,3], x[config_top[1,]==j,4],col=cl[j])
}


## TOP IN YX SPACE

#Plot top model in YX space
config_num=1
par(mfrow=c(1,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
for(i in 1:1){
plot(x[,i],y, col="white",xlab="x1",ylab="y")
k_1=k[output_reorder$config_index[config_sp$ix[config_num]]]
for(j in 1:k_1){
points(x[config_top[config_num,]==j,i], y[config_top[config_num,]==j],col=cl[j])
}
}

#Plot top model in YX space
config_num=1
par(mfrow=c(1,1), mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
for(i in 1:1){
plot(x[,i],y, col="white",xlab="x1",ylab="y",main=paste( "p( rho |y,x)=",format(config_sp$x[i],digits=4, scientific=F)))
k_1=k[config_num]
for(j in 1:k_1){
points(x[config[config_num,]==j,i], y[config[config_num,]==j],col=cl[j])
}
}
#for cambridge
par(mfrow=c(1,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
i=1
plot(x[,i],y, col="white",xlab=paste("x",i),ylab="y")
k_1=k[output_reorder$config_index[config_sp$ix[1]]]
for(j in 1:k_1){
points(x[config_top[1,]==j,i], y[config_top[1,]==j],col=cl[j])
}



## PLOT TOP 8 IN X SPACE

# Top 8 configurations in X space
par(mfrow=c(2,4), mar=c(2,2,1,1)+.1)
for(h in 1:8){
i=S-(h-1)*1000
plot(x[,1],x[,2], col="white")
k_1=k[output_reorder$config_index[config_sp$ix[i]]]
for(j in 1:k_1){
points(x[config_top[i,]==j,1], x[config_top[i,]==j,2],col=j)
}
}
par(mfrow=c(2,4), mar=c(2,2,1,1)+.1)
for(h in 1:8){
i=S-(h-1)*1000
plot(x[,3],x[,4], col="white")
k_1=k[output_reorder$config_index[config_sp$ix[i]]]
for(j in 1:k_1){
points(x[config_top[i,]==j,3], x[config_top[i,]==j,4],col=j)
}
}


##TOP 8 IN YX SPACE

# Top 8 Y models in YX space
i_covar=1
par(mfrow=c(2,4), mar=c(2,2,1,1)+.1)
for(h in 1:8){
i=S-(h-1)*1000
plot(x[,i_covar],y, col="white")
k_1=k[output_reorder$config_index[config_sp$ix[1]]]
for(j in 1:k_1){
points(x[config_top[i,]==j,i_covar], y[config_top[i,]==j],col=j)
}
}

########################################################################

### PREDICTION

## PREDICTION

#simulate new data
set.seed(101010)
m=200
x_new=matrix(rnorm(m),m,1)*2+4
#x_new=seq(-2,10,.1)
#m=length(x_new)
#x_new=matrix(x_new,m,1)

#sample additional x
mu21_new=4+1/4*(x_new[,1]-4)%x%t(matrix(c(3.5,0),pmax,1))
xa_new=t(t(csigmax21)%*%matrix(rnorm(pmax*m),pmax,m))+mu21_new

#add to x
if(pa>0){
x_new=cbind(x_new,xa_new[,1:pa])}
#if(pa>0){
#x_new=matrix(x_new,m,pa+1)}
X_new=cbind(rep(1,m), x_new)

#calculate true expected value
px_t_new=matrix(0,m,2)
px_t_new[,1]=dnorm(x_new[,1],mu_t[1],sd_t[1])*prob[1]
px_t_new[,2]=dnorm(x_new[,1],mu_t[2],sd_t[2])*prob[2]
px_t_new=t(apply(px_t_new,1,std))
y_hat_true=px_t_new[,1]*(X_new[,1:2]%*%beta_t[,1])+px_t_new[,2]*(X_new[,1:2]%*%beta_t[,2])

#Calculate true density
y_grid=seq(-2,10,.05)
m2=length(y_grid)
f_true=t(matrix(px_t_new[,1],m,m2))*dnorm(matrix(y_grid,nrow=m2,ncol=m),t(matrix(X_new[,1:2]%*%beta_t[,1],nrow=m,ncol=m2)), sigma1^.5)+t(matrix(px_t_new[,2],m,m2))*dnorm(matrix(y_grid,nrow=m2,ncol=m),t(matrix(X_new[,1:2]%*%beta_t[,2],nrow=m,ncol=m2)), sigma2^.5)


#Compute prediction and predictive density estimates with no credible intervals
source("./Dropbox/EDP_mix/Ex_ip/predict_jdp.R")
output_pred=predict_jdp(S,m,m2,p,x_new, y_grid, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k, config, alpha, phi_y, sigma_y, mu_x, sigma_x )

#Compute prediction with credible intervals
source("./Dropbox/EDP_mix/Ex_ip/predict_cred_jdp_v2.R")
output_cred_pred=predict_cred_jdp(.05,S,m,p,x_new, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k, config, alpha, phi_y, sigma_y, mu_x, sigma_x )

#Compute predictive density with credible bounds
source("./Dropbox/EDP_mix/Ex_ip/predict_fcred_jdp_v2.R")
output_fcred_pred1_10=predict_fcred_jdp(.05, S,10,m2,p,x_new[1:10,], y_grid, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k, config, alpha, phi_y, sigma_y, mu_x, sigma_x )

save.image("./Dropbox/EDP_mix/Ex_ip/jdp_ex_ip15_S20000b5000_v2_n200.RData")

################################################
###Plot Prediction
#with credible intervals but without data
par(mfrow=c(1,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
for(i in 1:1){
x_newi_s=sort(x_new[,i],index.return=TRUE)
plot(x_newi_s$x,output_cred_pred$y_pred[x_newi_s$ix],type='l',col=2, ylim=c(-3,9), xlab="x1",ylab="E[y|x]")
lines(x_newi_s$x,y_hat_true[x_newi_s$ix], col=1)
lines(x_newi_s$x,output_cred_pred$l_pred[x_newi_s$ix], col="darkgray",lty=2)
lines(x_newi_s$x,output_cred_pred$u_pred[x_newi_s$ix], col="darkgray",lty=2)
}

#with data but without credible intervals
par(mfrow=c(1,1))
for(i in 1:1){
x_newi_s=sort(x_new[,i],index.return=TRUE)
plot(x_newi_s$x,output_pred$y_pred[x_newi_s$ix],type='l',col=2, ylim=c(min(c(y, output_pred$l_pred)),max(c(y, output_pred$u_pred))), xlim=c(min(x[,i],x_new[,i]),max(x[,i],x_new[,i])))
lines(x_newi_s$x,y_hat_true[x_newi_s$ix], col=1)
points(x[,i],y)
}

#Plot y_hat_true vs y_pred
par(mfrow=c(1,1))
plot(y_hat_true, output_pred$y_pred, xlim=c(-5,10),ylim=c(-5,10))
abline(0,1, col="grey")
points(y_hat_true, output_pred$u_pred, col=4)
points(y_hat_true, output_pred$l_pred, col=4)


# plot with less observations
sx1=sort(x_new[,1],index.return=T)
sx_new=matrix(x_new[sx1$ix,],m,p)
spred=output_pred$y_pred[sx1$ix]
supred=output_cred_pred$u_pred[sx1$ix]
slpred=output_cred_pred$l_pred[sx1$ix]
strue=y_hat_true[sx1$ix]
ind=c(1,5,10,15,20,30,40,50,70,90,110,130,150,160,170,180,185,190,195,200)

library(Hmisc)

par(mfrow=c(1,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
plot(sx_new[ind,1],strue[ind],pch="*", cex=1.5,xlab="x1", ylab="E[y|x]", ylim=c(-4,8))
#segments(sx_new[ind,1],slpred[ind],x1=sx_new[ind,1],y1=supred[ind], lwd=1.5)
#points(sx_new[ind,1],spred[ind],col="red", pch=1,lwd=1.5)
d=data.frame(x=sx_new[ind,1],y=spred[ind],l=slpred[ind],u=supred[ind])
with (
  data = d
  , expr = errbar(x, y, l, u, add=T, pch=1, cap=.03,col="blue",lwd=1,cex=1)
)


#####Plot predictive Density
par(mfrow=c(1,1))
i=1
par(mfrow=c(3,3))
for(i in 1:9){
plot(y_grid, output_pred$f_pred[,i],col=2, xlim=c(-2,10), ylim=c(0,max(f_true)))
points(y_grid, f_true[,i],col=1)
ind_pred=sum(y_grid<output_pred$y_pred[i])
fy_pred=(output_pred$f_pred[ind_pred,i]*(output_pred$y_pred[i]-y_grid[ind_pred])+output_pred$f_pred[ind_pred+1,i]*(y_grid[ind_pred+1]-output_pred$y_pred[i]))/(y_grid[ind_pred+1]-y_grid[ind_pred])
points(output_pred$y_pred[i],fy_pred, pch='*',cex=2,col=2)
#ind_l=sum(y_grid<output_pred$l_pred[i])
#fy_l=(output_pred$f_pred[ind_l,i]*(output_pred$l_pred[i]-y_grid[ind_l])+output_pred$f_pred[ind_l+1,i]*(y_grid[ind_l+1]-output_pred$l_pred[i]))/(y_grid[ind_l+1]-y_grid[ind_l])
#points(output_pred$l_pred[i],fy_l, pch='*',cex=2, col='blue')
#ind_u=sum(y_grid<output_pred$u_pred[i])
#fy_u=(output_pred$f_pred[ind_u,i]*(output_pred$u_pred[i]-y_grid[ind_u])+output_pred$f_pred[ind_u+1,i]*(y_grid[ind_u+1]-output_pred$u_pred[i]))/(y_grid[ind_u+1]-y_grid[ind_u])
#points(output_pred$u_pred[i],fy_u, pch='*',cex=2, col='blue')
ind_true=sum(y_grid<y_hat_true[i])
fy_true=(f_true[ind_true,i]*(y_hat_true[i]-y_grid[ind_true])+f_true[ind_true+1,i]*(y_grid[ind_true+1]-y_hat_true[i]))/(y_grid[ind_true+1]-y_grid[ind_true])
points(y_hat_true[i],fy_true, pch='*',cex=2, col=1)
}

par(mfrow=c(1,1))
i=1
plot(y_grid, output_fcred_pred1_10$f_pred[,i],type='l', col=2, xlim=c(-2,10), ylim=c(0,max(c(f_true,output_fcred_pred1_10$u_fpred))),xlab="y", ylab="f(y|x)")
points(y_grid, f_true[,i],type='l',col=1)
points(y_grid, output_fcred_pred1_10$l_fpred[,i],type='l',col='blue',lty=2)
points(y_grid, output_fcred_pred1_10$u_fpred[,i],type='l',col='blue',lty=2)

write.table(cbind(matrix(f_true[,1:10],ncol=1),matrix(unlist(output_fcred_pred1_10),ncol=3)), file="C:\\Users\\Sara\\Documents\\EDP_mix\\Ex\\fpred_jdp.csv")

#empirical l2 prediction error
l2_err=sum(((y_hat_true-output_pred$y_pred)^2)/m)^.5
l2_err


l1_err=sum((abs(y_hat_true-output_pred$y_pred))/m)
l1_err


#print the first 10 predictions with credible intervals
round(cbind(y_hat_true[1:10], output_cred_pred$y_pred[1:10],output_cred_pred$l_pred[1:10],output_cred_pred$u_pred[1:10]),3)

write.table(cbind(y_hat_true, output_cred_pred$y_pred,output_cred_pred$l_pred,output_cred_pred$u_pred), file="C:\\Users\\Sara\\Documents\\EDP_mix\\Ex\\pred_jdp.csv")


#estimated l1 distance for density
#.05 is grid with, CHANGE if grid with changes
l1_dist=colSums(abs(f_true-output_pred$f_pred))*.05

#Average l1 distance
mean(l1_dist)

#Max l1 dist
max(l1_dist)

#Min l1 dist
min(l1_dist)


