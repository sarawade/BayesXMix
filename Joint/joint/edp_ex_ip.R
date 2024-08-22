#EDP mixture for regression
#Simulate data for second problem and call EDP MCMC function to obtain
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

par(mfrow=c(1,1))
plot(x,y)

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

#Set parameter values

#y parameters
#remember! variance for beta is sigma^2*iC
mu_theta=apply(beta_t,1,mean)
if(pa>0){
mu_theta=matrix(c(mu_theta,rep(0,pa)),nrow=(p+1))
}
#iXX=solve(t(X)%*%X)
#iC.corr=diag(1/sqrt(diag(iXX)))%*%iXX%*%diag(1/sqrt(diag(iXX)))
#iC.diag=10*c((2.5^2)/(qnorm(.5+.5*.95^(1/(p+1)))^2),rep((.5^2)/(qnorm(.5+.5*.95^(1/(p+1)))^2),p))
iC.diag=c(20,rep(1,p))
#iC=diag(sqrt(iC.diag))%*%iC.corr%*%diag(sqrt(iC.diag))
iC=diag(iC.diag)
C=solve(iC)
a_y=2
b_y=1/10

# x parameters
a_x=matrix(2,p,1)
b_x=matrix(1,p,1)
mu_0=matrix(4,p,1)
#c_x=rep((qnorm(.5+.5*.95^(1/p))^2)/(4^2),p)
c_x=rep(1/(4),p)

#Precision
u_alpha_x=1
v_alpha_x=1
u_alpha_y=1
v_alpha_y=1


### Parameters for MCMC function
S=20000
burnin=5000

#Initialize chain

#Option 1
#start everyone in 1 group
#config_x_init=rep(1,n)
#config_y_init=rep(1,n)

#alpha_y_init=u_alpha_y/v_alpha_y
#alpha_x_init=u_alpha_x/v_alpha_x

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
#mu_x_init=list()
#mu_x_init[[1]]=mu_x_hat
#sigma_x_init=list()
#sigma_x_init[[1]]=b_x_hat/(a_x_hat-1)

#Option 2: Load previous data
#configx_prev=config_x
#S_prev=S
#configy_prev=config_y
#alphay_prev=alpha_y
#alphax_prev=alpha_x
#phi_prev=phi_y
#sigma_y_prev=sigma_y
#mu_x_prev=mu_x
#sigma_x_prev=sigma_x
#ky_prev=k_y
#kx_prev=k_x

#config_y_init=configy_prev[S_prev,]
#config_x_init=configx_prev[S_prev,]
#phi_init=phi_prev[[S_prev]]
#sigma_y_init=sigma_y_prev[[S_prev]]
#mu_x_init=mu_x_prev[[S_prev]]
#sigma_x_init=sigma_x_prev[[S_prev]]
#alpha_y_init=alphay_prev[S_prev]
#alpha_x_init=alphax_prev[[S_prev]]

## Option 3: Random Cluster assignments
ky_init=3
kx_init=2

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


#####################################################################

source("edpyx_mcmc.R")

#Call MCMC edp function
set.seed(101010)
output=edp_mcmc(S, burnin, y, x, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, u_alpha_x, v_alpha_x, u_alpha_y, v_alpha_y, config_y_init, config_x_init, phi_init, sigma_y_init, mu_x_init, sigma_x_init,alpha_y_init, alpha_x_init)

save.image("edp_ex_ip15_S20000b5000_v2_n200.RData")
#######################################################################

# Define variables
config_x=output$config_x
config_y=output$config_y
phi_y=output$phi_y
sigma_y=output$sigma_y
mu_x=output$mu_x
sigma_x=output$sigma_x
alpha_y=output$alpha_y
alpha_x=output$alpha_x

k_y=rep(0,S)
k_x=list()
for(s in 1:S){
	k_y[s]=length(unique(config_y[s,]))
	k_x[[s]]=rep(0,k_y[s])
	for(i in 1:k_y[s]){
		k_x[[s]][i]=length(unique(config_x[s,config_y[s,]==i]))
	}
}

#####################################################################

#Diagnostics
# Choose an individual
i=1
phi_hat_i=matrix(0,S,p+1)
phi_i=matrix(0,S, p+1)
phi_i[1,]=phi_y[[1]][,config_y[1,i]]
phi_hat_i[1,]=(phi_i[1,])
for(s in 2:S){
phi_y_s=matrix(phi_y[[s]],p+1,k_y[s])
phi_i[s,]=phi_y_s[,config_y[s,i]]
phi_hat_i[s,]=(phi_hat_i[s-1,]*(s-1)+phi_i[s,])/s
}

par(mfrow=c(2,1),mar=c(1.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
for(j in 1:(2)){
plot(phi_i[,j], ylab=paste("beta",j-1), xlab="", main=paste("beta*",j-1,"=",round(phi_hat_i[s,j],4)) )
points(phi_hat_i[,j],col=2)
}

par(mfrow=c(1,1),mar=c(1.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
j=2
plot(phi_i[,j], ylab="beta_2", xlab="", main=paste("beta*_2=",round(phi_hat_i[s,j],4)))
points(phi_hat_i[,j],col=2)

round(phi_hat_i[s,],4)

#for sigma_y
i=1
sigmay_hat_i=matrix(0,S,1)
sigmay_i=matrix(0,S, 1)
sigmay_i[1]=sigma_y[[1]][config_y[1,i]]
sigmay_hat_i[1]=(sigmay_i[1])
for(s in 2:S){
sigmay_i[s]=sigma_y[[s]][config_y[s,i]]
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
phi_hat[i,]=phi_y[[1]][,config_y[1,i]]
sigmay_hat[i]=sigma_y[[1]][config_y[1,i]]
for(s in 2:S){
phi_y_s=matrix(phi_y[[s]],p+1,k_y[s])
phi_hat[i,]=phi_hat[i,]+phi_y_s[,config_y[s,i]]
sigmay_hat[i]=sigmay_hat[i]+sigma_y[[s]][config_y[s,i]]
}
}
phi_hat=phi_hat/S
sigmay_hat=sigmay_hat/S
beta_tf=beta_t
if(pa>0){beta_tf=rbind(beta_t,matrix(0,pa,2))}

# Compute true parameter estimates
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
i=1
mu_hat_i=matrix(0,S,p)
mu_i=matrix(0,S, p)
mu_i[1,]=mu_x[[1]][[config_y[1,i]]][,config_x[1,i]]
mu_hat_i[1,]=(mu_i[1,])
for(s in 2:S){
mu_i[s,]=mu_x[[s]][[config_y[s,i]]][,config_x[s,i]]
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
sigmax_i[1,]=sigma_x[[1]][[config_y[1,i]]][,config_x[1,i]]
sigmax_hat_i[1,]=(sigmax_i[1,])
for(s in 2:S){
sigmax_i[s,]=sigma_x[[s]][[config_y[s,i]]][,config_x[s,i]]
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


## Number of y-clusters
ky_hat=mean(k_y)
ky_hat
par(mfrow=c(1,1))
hist(k_y,freq=F)
k_y_hat_round=round(ky_hat)

# Precision parameter for y
mean(alpha_y)
plot(density(alpha_y))

#Precision parameter for x
alpha_x_vec=unlist(alpha_x[(burnin+1):(burnin+S)])
phi_vec=unlist(phi_y[(burnin+1):(burnin+S)])
sigma_y_vec=unlist(sigma_y[(burnin+1):(burnin+S)])
par(mfrow=c(2,1))
for (i in 1:(p+1)){
plot(phi_vec[seq(i,length(phi_vec),p+1)], alpha_x_vec)
}
plot(sigma_y_vec, alpha_x_vec)

##Previous plot is too confusing, average over nearby points
#number of points to use in average for the interval
n_avg=100
n_ax=length(alpha_x_vec)
par(mfrow=c(1,3),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
for (i in 1:(2)){
phi_i_s=sort(phi_vec[seq(i,length(phi_vec),p+1)],index.return=TRUE)
x_toplot=c(colSums(matrix(phi_i_s$x[1:(n_ax-n_ax%%n_avg)], nrow=n_avg))/n_avg, mean(phi_i_s$x[(n_ax-n_ax%%n_avg+1):n_ax]))
y_toplot=c(colSums(matrix(alpha_x_vec[phi_i_s$ix][1:(n_ax-n_ax%%n_avg)], nrow=n_avg))/n_avg, mean(alpha_x_vec[phi_i_s$ix][(n_ax-n_ax%%n_avg+1):n_ax]))
plot(x_toplot, y_toplot,, ylab="alpha_x",xlab=paste("beta",i-1))
}
sigmay_s=sort(sigma_y_vec,index.return=TRUE)
x_toplot=c(colSums(matrix(sigmay_s$x[1:(n_ax-n_ax%%n_avg)], nrow=n_avg))/n_avg, mean(sigmay_s$x[(n_ax-n_ax%%n_avg+1):n_ax]))
y_toplot=c(colSums(matrix(alpha_x_vec[sigmay_s$ix][1:(n_ax-n_ax%%n_avg)], nrow=n_avg))/n_avg, mean(alpha_x_vec[sigmay_s$ix][(n_ax-n_ax%%n_avg+1):n_ax]))
plot(x_toplot, y_toplot, ylab="alpha_x",xlab="sigma_y")

### k_x
sk_x=rep(0,S)
for(s  in 1:S){
sk_x[s]=sum(k_x[[s]])
}
#mean
sum(sk_x)/S
hist(sk_x)

######################################################################

########### Inference for partitions

# Need to Reorder labels
source("./Dropbox/EDP_mix/Ex_ip/reorder_labels.R")
output_reorder=reorder_labels(S,n, config_y, config_x)

save.image("./Dropbox/EDP_mix/Ex_ip/edp_ex_ip15_S20000b5000_v2_n200.RData")

#Compute estimated probability of y configurations and sort by highest prob
configy_p=output_reorder$configy_count/S
configy_sp=sort(configy_p, decreasing=TRUE, index.return=TRUE)
#Total number of y comfigurations with positive prob
length(configy_p)

#Compute estimated probability of yx configurations and sort by highest prob
configyx_p=output_reorder$configyx_count/S
configyx_sp=sort(configyx_p, decreasing=TRUE, index.return=TRUE)
#Total number of yx comfigurations with positive prob
length(configyx_p)

#probability of top 8 y configurations
configy_sp$x[1:8]
sum(configy_sp$x[1:8])

#probability of top 8 yx configurations
configyx_sp$x[1:8]
sum(configyx_sp$x[1:8])

#configuration of top
configy_top=config_y[output_reorder$configy_index[configy_sp$ix],]
configyx_topy=config_y[output_reorder$configyx_index[1,configyx_sp$ix],]
configyx_topx=config_x[output_reorder$configyx_index[2,configyx_sp$ix],]

#number in common with posterior mode and true
sum(config_true==configy_top[1,])

#any configurations equal to true configuration?
if(sum(apply(configy_top,1,max)==2)>0){
configy_true=rowSums(configy_top[apply(configy_top,1,max)==2,]!=t(matrix(config_true,n,sum(apply(configy_top,1,max)==2))))
sum(configy_true==n)
sum(configy_true==0)
max(configy_true)
min(configy_true)
}

########Plot top models

## TOP IN X SPACE 

#Plot top Y model in X space

par(mfrow=c(2,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
plot(x[,1],x[,2], xlab="x1",ylab="x2")
k_1=max(configy_top[1,])
for(j in 1:k_1){
points(x[configy_top[1,]==j,1], x[configy_top[1,]==j,2],col=j)
}
plot(x[,3],x[,4], xlab="x3",ylab="x4")
k_1=max(configy_top[1,])
for(j in 1:k_1){
points(x[configy_top[1,]==j,3], x[configy_top[1,]==j,4],col=j)
}


#Plot top Y model with X configuration in X space

par(mfrow=c(2,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
plot(x[,1],x[,2], col="white",xlab="x1",ylab="x2")
ky_1=k_y[output_reorder$configy_index[configy_sp$ix[1]]]
kx_1=k_x[[output_reorder$configy_index[configy_sp$ix[1]]]]
for(j in 1:ky_1){
for(h in 1:kx_1[j]){
points(x[configy_top[1,]==j,1][config_x[output_reorder$configy_index[configy_sp$ix[1]],configy_top[1,]==j]==h], x[configy_top[1,]==j,2][config_x[output_reorder$configy_index[configy_sp$ix[1]],configy_top[1,]==j]==h],pch=h,col=j)
}
}
plot(x[,3],x[,4], col="white",xlab="x3",ylab="x4")
ky_1=k_y[output_reorder$configy_index[configy_sp$ix[1]]]
kx_1=k_x[[output_reorder$configy_index[configy_sp$ix[1]]]]
for(j in 1:ky_1){
for(h in 1:kx_1[j]){
points(x[configy_top[1,]==j,3][config_x[output_reorder$configy_index[configy_sp$ix[1]],configy_top[1,]==j]==h], x[configy_top[1,]==j,4][config_x[output_reorder$configy_index[configy_sp$ix[1]],configy_top[1,]==j]==h],pch=h,col=j)
}
}

#Plot top YX in X space, if all the different not interesting
par(mfrow=c(1,2), mar=c(2,2,1,1)+.1)
plot(x[,1],x[,2], col="white")
ky_1=k_y[output_reorder$configyx_index[1,configyx_sp$ix[1]]]
kx_1=k_x[[output_reorder$configyx_index[2,configyx_sp$ix[1]]]]
for(j in 1:ky_1){
for(h in 1:kx_1[j]){
points(x[configyx_topy[1,]==j,1][configyx_topx[1,configyx_topy[1,]==j]==h], x[configyx_topy[1,]==j,2][configyx_topx[1,configyx_topy[1,]==j]==h],pch=h,col=j)
}
}
plot(x[,3],x[,4], col="white")
ky_1=k_y[output_reorder$configyx_index[1,configyx_sp$ix[1]]]
kx_1=k_x[[output_reorder$configyx_index[2,configyx_sp$ix[1]]]]
for(j in 1:ky_1){
for(h in 1:kx_1[j]){
points(x[configyx_topy[1,]==j,3][configyx_topx[1,configyx_topy[1,]==j]==h], x[configyx_topy[1,]==j,4][configyx_topx[1,configyx_topy[1,]==j]==h],pch=h,col=j)
}
}


## TOP IN YX SPACE

#Plot top Y model in YX space
confignum=1
par(mfrow=c(1,1), mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
for(i in 1:1){
plot(x[,i],y, col="white", xlab="x1",ylab="y", main=paste( "p( rho_y |y,x)=",format(configy_sp$x[i],digits=4, scientific=F)))
ky_1=k_y[output_reorder$configy_index[configy_sp$ix[confignum]]]
for(j in 1:ky_1){
points(x[configy_top[confignum,]==j,i], y[configy_top[confignum,]==j],col=j)
}
}

#Plot top Y model with X config in YX space
confignum=S
par(mfrow=c(1,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
for(i in 1:p){
plot(x[,i],y, col="white", xlab=paste("x",i),ylab="y")
ky_1=k_y[output_reorder$configy_index[configy_sp$ix[confignum]]]
kx_1=k_x[[output_reorder$configy_index[configy_sp$ix[confignum]]]]
for(j in 1:ky_1){
for(h in 1:kx_1[j]){
points(x[configy_top[confignum,]==j,i][config_x[output_reorder$configy_index[configy_sp$ix[confignum]],configy_top[confignum,]==j]==h], y[configy_top[confignum,]==j][config_x[output_reorder$configy_index[configy_sp$ix[confignum]],configy_top[confignum,]==j]==h],pch=h,col=j)
}
}
}

#Top YX config in YX space, if all different, not interesting
par(mfrow=c(1,1), mar=c(2,2,1,1)+.1)
for(i in 1:p){
plot(x[,i],y, col="white")
ky_1=k_y[output_reorder$configyx_index[1,configyx_sp$ix[S]]]
kx_1=k_x[[output_reorder$configyx_index[2,configyx_sp$ix[S]]]]
for(j in 1:ky_1){
for(h in 1:kx_1[j]){
points(x[configyx_topy[S,]==j,i][configyx_topx[S,configyx_topy[S,]==j]==h], y[configyx_topy[S,]==j][configyx_topx[S,configyx_topy[S,]==j]==h],pch=h,col=j)
}
}
}


## PLOT TOP 8 IN X SPACE

# Top 8 Y configurations in X space
par(mfrow=c(2,4), mar=c(2,2,1,1)+.1)
for(h in 1:8){
#i=S-(h-1)*1000
i=h
plot(x[,1],x[,2], col="white")
ky_1=k_y[output_reorder$configy_index[configy_sp$ix[i]]]
for(j in 1:ky_1){
points(x[configy_top[i,]==j,1], x[configy_top[i,]==j,2],col=j)
}
}
par(mfrow=c(2,4), mar=c(2,2,1,1)+.1)
for(h in 1:8){
i=S-(h-1)*1000
plot(x[,3],x[,4], col="white")
ky_1=k_y[output_reorder$configy_index[configy_sp$ix[i]]]
for(j in 1:ky_1){
points(x[configy_top[i,]==j,3], x[configy_top[i,]==j,4],col=j)
}
}

#Top 8 Y models with X configuration in X space
#X1 vs X2
par(mfrow=c(2,4), mar=c(2,2,1,1)+.1)
for(h in 1:8){
i=S-(h-1)*1000
plot(x[,1],x[,2], col="white",main=paste( "p( rho |y,x)=",format(configy_sp$x[i],digits=4, scientific=F)))
ky_1=k_y[output_reorder$configy_index[configy_sp$ix[i]]]
kx_1=k_x[[output_reorder$configy_index[configy_sp$ix[i]]]]
for(j in 1:ky_1){
for(h in 1:kx_1[j]){
points(x[configy_top[i,]==j,1][config_x[output_reorder$configy_index[configy_sp$ix[i]],configy_top[i,]==j]==h], x[configy_top[i,]==j,2][config_x[output_reorder$configy_index[configy_sp$ix[i]],configy_top[i,]==j]==h],pch=h,col=j)
}
}
}
#X3 vs X4
par(mfrow=c(2,4), mar=c(2,2,1,1)+.1)
for(h in 1:8){
i=S-(h-1)*1000
plot(x[,3],x[,4], col="white",main=paste( "p( rho |y,x)=",format(configy_sp$x[i],digits=4, scientific=F)))
ky_1=k_y[output_reorder$configy_index[configy_sp$ix[i]]]
kx_1=k_x[[output_reorder$configy_index[configy_sp$ix[i]]]]
for(j in 1:ky_1){
for(h in 1:kx_1[j]){
points(x[configy_top[i,]==j,3][config_x[output_reorder$configy_index[configy_sp$ix[i]],configy_top[i,]==j]==h], x[configy_top[i,]==j,4][config_x[output_reorder$configy_index[configy_sp$ix[i]],configy_top[i,]==j]==h],pch=h,col=j)
}
}
}


#Top 8 YX models in YX space, if all different, not interesting
#X1 vs X2 
par(mfrow=c(2,4), mar=c(2,2,1,1)+.1)
for(h in 1:8){
#i=S-(h-1)*1000
i=h
plot(x[,1],x[,2], col="white")
ky_1=k_y[output_reorder$configyx_index[1,configyx_sp$ix[i]]]
kx_1=k_x[[output_reorder$configyx_index[2,configyx_sp$ix[i]]]]
for(j in 1:ky_1){
for(h in 1:kx_1[j]){
points(x[configyx_topy[i,]==j,1][configyx_topx[i,configyx_topy[i,]==j]==h], x[configyx_topy[i,]==j,2][configyx_topx[i,configyx_topy[i,]==j]==h],pch=h,col=j)
}
}
}
#X3 vs X4
par(mfrow=c(2,4), mar=c(2,2,1,1)+.1)
for(h in 1:8){
i=S-(h-1)*1000
plot(x[,3],x[,4], col="white")
ky_1=k_y[output_reorder$configyx_index[1,configyx_sp$ix[i]]]
kx_1=k_x[[output_reorder$configyx_index[2,configyx_sp$ix[i]]]]
for(j in 1:ky_1){
for(h in 1:kx_1[j]){
points(x[configyx_topy[i,]==j,3][configyx_topx[i,configyx_topy[i,]==j]==h], x[configyx_topy[i,]==j,4][configyx_topx[i,configyx_topy[i,]==j]==h],pch=h,col=j)
}
}
}

##TOP 8 IN YX SPACE

# Top 8 Y models in YX space
i_covar=1
par(mfrow=c(2,4), mar=c(2,2,1,1)+.1)
for(h in 1:8){
#i=S-(h-1)*1000
i=h
plot(x[,i_covar],y, col="white")
ky_1=k_y[output_reorder$configy_index[configy_sp$ix[1]]]
for(j in 1:ky_1){
points(x[configy_top[i,]==j,i_covar], y[configy_top[i,]==j],col=j)
}
}

#Top 8 Y models with X configurations in YX space
i_chosen=1
par(mfrow=c(2,4), mar=c(2,2,1,1)+.1)
for(h in 1:8){
i=S-(h-1)*1000
plot(x[,i_chosen],y, col="white",main=paste( "p( rho |y,x)=",format(configy_sp$x[i],digits=4, scientific=F)))
ky_1=k_y[output_reorder$configy_index[configy_sp$ix[i]]]
kx_1=k_x[[output_reorder$configy_index[configy_sp$ix[i]]]]
for(j in 1:ky_1){
for(h in 1:kx_1[j]){
points(x[configy_top[i,]==j,i_chosen][config_x[output_reorder$configy_index[configy_sp$ix[i]],configy_top[i,]==j]==h], y[configy_top[i,]==j][config_x[output_reorder$configy_index[configy_sp$ix[i]],configy_top[i,]==j]==h],pch=h,col=j)
}
}
}

#Top YX configurations in YX space, if all different not interesting
i_covar=1
par(mfrow=c(2,4), mar=c(2,2,1,1)+.1)
for(h in 1:8){
i=S-(h-1)*1000
plot(x[,i_covar],y, col="white")
ky_1=k_y[output_reorder$configyx_index[1,configyx_sp$ix[i]]]
kx_1=k_x[[output_reorder$configyx_index[2,configyx_sp$ix[i]]]]
for(j in 1:ky_1){
for(h in 1:kx_1[j]){
points(x[configyx_topy[i,]==j,i_covar][configyx_topx[i,configyx_topy[i,]==j]==h], y[configyx_topy[i,]==j][configyx_topx[i,configyx_topy[i,]==j]==h],pch=h,col=j)
}
}
}


########################################################################

### PREDICTION


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

#Compute prediction and predictive density estimates without credible intervals
source("./Dropbox/EDP_mix/Ex_ip/predict_edp.R")
output_pred=predict_edp(S,m,m2,p,x_new, y_grid, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k_y, k_x, config_y, config_x, alpha_x, alpha_y, phi_y, sigma_y, mu_x, sigma_x )

#Compute credible intervals for prediction
source("./Dropbox/EDP_mix/Ex_ip/predict_cred_edp_v2.R")
output_cred_pred=predict_cred_edp(.05,S,m,p,x_new, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k_y, k_x, config_y, config_x, alpha_x, alpha_y, phi_y, sigma_y, mu_x, sigma_x )

#Compute predictive density for specific individual with credible bounds
source("./Dropbox/EDP_mix/Ex_ip/predict_fcred_edp_v2.R")
output_fcred_pred1_10=predict_fcred_edp(.05, S,10,m2,p,x_new[1:10,], y_grid, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k_y, k_x, config_y, config_x, alpha_x, alpha_y, phi_y, sigma_y, mu_x, sigma_x )

save.image("./Dropbox/EDP_mix/Ex_ip/edp_ex_ip15_S20000b5000_v2_n200.RData")


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

#no lines
#with credible intervals but without data
par(mfrow=c(1,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
for(i in 1:1){
x_newi_s=sort(x_new[,i],index.return=TRUE)
plot(x_newi_s$x,output_cred_pred$y_pred[x_newi_s$ix],pch=16, col=2, ylim=c(-3,9), xlim=c(min(x[,i],x_new[,i]),max(x[,i],x_new[,i])), xlab="x1",ylab="E[y|x]")
points(x_newi_s$x,y_hat_true[x_newi_s$ix],pch="*",col=1)
points(x_newi_s$x,output_cred_pred$l_pred[x_newi_s$ix], pch="-", cex=1, col="darkgray",lty=2)
points(x_newi_s$x,output_cred_pred$u_pred[x_newi_s$ix], pch="-", cex=1, col="darkgray",lty=2)
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
  , expr = errbar(x, y, l, u, add=T, pch=1, cap=.03,col=2,lwd=1,cex=1)
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
points(output_pred$y_pred[i],fy_pred, pch='*',cex=3,col=2)
#ind_l=sum(y_grid<output_pred2$l_y_pred[i])
#fy_l=(output_pred$f_pred[ind_l,i]*(output_pred2$l_y_pred[i]-y_grid[ind_l])+output_pred$f_pred[ind_l+1,i]*(y_grid[ind_l+1]-output_pred2$l_y_pred[i]))/(y_grid[ind_l+1]-y_grid[ind_l])
#points(output_pred2$l_y_pred[i],fy_l, pch='*',cex=3, col='blue')
#ind_u=sum(y_grid<output_pred2$u_y_pred[i])
#fy_u=(output_pred$f_pred[ind_u,i]*(output_pred2$u_y_pred[i]-y_grid[ind_u])+output_pred$f_pred[ind_u+1,i]*(y_grid[ind_u+1]-output_pred2$u_y_pred[i]))/(y_grid[ind_u+1]-y_grid[ind_u])
#points(output_pred2$u_y_pred[i],fy_u, pch='*',cex=3, col='blue')
ind_true=sum(y_grid<y_hat_true[i])
fy_true=(f_true[ind_true,i]*(y_hat_true[i]-y_grid[ind_true])+f_true[ind_true+1,i]*(y_grid[ind_true+1]-y_hat_true[i]))/(y_grid[ind_true+1]-y_grid[ind_true])
points(y_hat_true[i],fy_true, pch='*',cex=3, col=1)
}

par(mfrow=c(1,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, .5, 0))
i=1
plot(y_grid, output_fcred_pred1_10$f_pred[,i],type='l', col=2, xlim=c(-2,10), ylim=c(0,max(c(f_true,output_fcred_pred1_10$u_fpred))),xlab="y", ylab="f(y|x)")
points(y_grid, f_true[,i],type='l',col=1)
points(y_grid, output_fcred_pred1_10$l_fpred[,i],type='l',col='blue',lty=2)
points(y_grid, output_fcred_pred1_10$u_fpred[,i],type='l',col='blue',lty=2)

write.table(x_new,file="C:\\Users\\Sara\\Documents\\EDP_mix\\Ex\\x_new.csv")

write.table(cbind(matrix(f_true[,1:10],ncol=1),matrix(unlist(output_fcred_pred1_10),ncol=3)), file="C:\\Users\\Sara\\Documents\\EDP_mix\\Ex\\fpred_edp.csv")

write.table(cbind(y_hat_true, output_cred_pred$y_pred,output_cred_pred$l_pred,output_cred_pred$u_pred), file="C:\\Users\\Sara\\Documents\\EDP_mix\\Ex\\pred_edp.csv")



#empirical l2 prediction error
l2_err=sum(((y_hat_true-output_pred$y_pred)^2)/m)^.5
l2_err


l1_err=sum((abs(y_hat_true-output_pred$y_pred))/m)
l1_err

#estimated l1 distance for density
#.05 is grid with, CHANGE if grid with changes
l1_dist=colSums(abs(f_true-output_pred$f_pred))*.05

#Average l1 distance
mean(l1_dist)

#Max l1 dist
max(l1_dist)

#Min l1 dist
min(l1_dist)



