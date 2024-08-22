### Function that obtains posterior samples from  joint DP regression model

###INPUT
# S:number of iterations to save
# burnin: number of iterations to discard
# y: observed response (nx1)
# x: observed covariates (nxp)
# mu_theta, C: prior parameters for beta (p+1x1 and p+1xp+1)
# a_y, b_y: prior parameters for sigma^2_y (1x1 and 1x1)
# mu_0, c_x: prior parameters for mu  (px1 and px1)
# a_x, b_x: prior parameters for sigma^2_x (px1 and px1)
# u_alpha, v_alpha: prior parameters for alpha (1x1 and 1x1)
# config_init: initial configuration (1xn)
# phi_init: initali values of beta for each y cluster (p+1xkn_y)
# sigma_y_init: initial values of sigma_y for each y cluster (1xkn_y)
# mu_x_init, sigma_x_init: initial values of mu_x and sigma_x for each y-x cluster (list with kn_y elements, each element is pxkn_x[j])
# alpha_init: initial value of mass parameter (1x1)

###OUPUT
# config: sampled configurations (Sxn)
# phi_y, sigma_y: sampled y cluster parameters (list with S elements, each element is p+1xkn[s] and 1xkn[s])
# mu_x, sigma_x: sampled x cluster parameters (list with S elements, each element is pxkn[s] and pxkn[s])
# alpha: sampled precision parameters (Sx1)


jdp_mcmc=function(S, burnin, y, x, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, u_alpha, v_alpha, config_init, phi_init, sigma_y_init, mu_x_init, sigma_x_init,alpha_init){

#define n and p
n=dim(x)[1]
p=dim(x)[2]

#define X
X=cbind(rep(1,n), x)

#define iC
iC=solve(C)

#Empty matrices to store MCMC output

config=matrix(0,S,n)
phi_y=list()
sigma_y=list()
mu_x=list()
sigma_x=list()
alpha=matrix(0, S,1)

#Initialize
config_s=config_init
alpha_s=alpha_init
phi_s=phi_init
sigma_y_s=sigma_y_init
mu_x_s=mu_x_init
sigma_x_s=sigma_x_init

#terms that don't depend on s
muCmu= t(mu_theta)%*%C%*%mu_theta
mucmu_x=c_x*mu_0^2

for(s in 1:(S+burnin-1)){

	#Sample configurations at s+1

	#configuration to be updated
	c_update=config_s

	#unique parameters to be updated if new value is chosen
	phi_update=phi_s
	sigmay_update=sigma_y_s
	sigmax_update=sigma_x_s
	mux_update=mu_x_s

	#updated hyperparameters that dont depend on i	
	a_y_hat=a_y+1/2
	a_x_hat=a_x+1/2
	c_x_hat=c_x+1

	for(i in 1:n){
	
		#Define configuration without i
		c_ni=c_update[-i]
		c_i=c_update[i]
		uniq_c_ni=unique(c_ni)
		dn_ni=length(uniq_c_ni)
	
		#check if c_i is a class with the single element i
		singleton=(sum(c_ni==c_i)==0)


		if(singleton){
			#remove the parameters for class c_i
			phi_update[,c_i]=phi_update[,dn_ni+1]
			phi_update=matrix(phi_update[,-(dn_ni+1)], nrow=(p+1))
			sigmay_update[c_i]=sigmay_update[dn_ni+1]
			sigmay_update=sigmay_update[-(dn_ni+1)]
			sigmax_update[,c_i]=sigmax_update[,dn_ni+1]
			sigmax_update=matrix(sigmax_update[,-(dn_ni+1)], nrow=p)
			mux_update[,c_i]=mux_update[,dn_ni+1]
			mux_update=matrix(mux_update[,-(dn_ni+1)],nrow=p)
			#change c_ni, updated config
			c_ni[c_ni==(dn_ni+1)]=c_i
			c_update[c_update==(dn_ni+1)]=c_i
			c_update[i]=dn_ni+1		
		}

			

		#Set value of config for observation i
		log_prob=matrix(0,dn_ni+1,1)
		log_x=matrix(0,dn_ni+1,1)
		log_y=matrix(0,dn_ni+1,1)
		phi_update=matrix(phi_update, nrow=(p+1), ncol=dn_ni)
		mux_update=matrix(mux_update, nrow=p, ncol=dn_ni)
		sigmax_update=matrix(sigmax_update, nrow=p, ncol=dn_ni)


		#calculate probabilities for old values
		for(j in 1:dn_ni){
			n_j_ni=sum(c_ni==j)
			z_ij=sum(phi_update[,j]*X[i,])
			log_y[j]=dnorm(y[i],z_ij,sigmay_update[j]^.5, log=TRUE)
			log_x[j]=sum(dnorm(x[i,], mux_update[,j], sigmax_update[,j]^.5, log=TRUE))
			log_prob[j]=log(n_j_ni/(alpha_s+n-1))+log_y[j]+log_x[j]
		}

		#calculate probability of a new value
		#updated parameters
		C_hat=C+X[i,]%*%t(X[i,])
		decomp=eigen(C_hat, symmetric=TRUE)
		iC_hat=decomp$vectors%*%diag(1/decomp$values)%*%t(decomp$vectors)
	
		sc=1-X[i,]%*%iC_hat%*%X[i,]
		yi_hat=X[i,]%*%mu_theta
		y_sc=(y[i]-yi_hat)*(a_y*sc/b_y)^.5
		log_y[dn_ni+1]=dt(y_sc,2*a_y,log=TRUE)+.5*(log(sc)-log(b_y)+log(a_y))

		#for x
		x_sc=(x[i,]-mu_0)*(c_x*a_x/(c_x_hat*b_x))^.5
		log_x[dn_ni+1]=sum(dt(x_sc,2*a_x, log=TRUE)+.5*(log(c_x)-log(c_x+1)-log(b_x)+log(a_x)))
		log_prob[dn_ni+1]=log(alpha_s/(alpha_s+n-1))+log_y[dn_ni+1]+log_x[dn_ni+1]
	
		#need to find normalzing constant
		#f <- function (k) sum(exp(k+log_prob))-1
		#k_star=uniroot(f,lower=-10^(10),upper=10^20)$root
		#probs=exp(k_star+log_prob)
		probs=exp(log_prob)/sum(exp(log_prob))

		#set value of config according to probability
 		c_update[i]=sum(runif(1)>c(cumsum(probs[-(dn_ni+1)]),1))+1

		
		#update parameters if new chosen
		if(c_update[i]==(dn_ni+1)){
			#updated parameters
			mu_hat=iC_hat%*%(C%*%mu_theta+X[i,]*y[i])
			b_y_hat=b_y+(y[i]^2+muCmu-t(mu_hat)%*%C_hat%*%mu_hat)/2

			#updated parameters for x
			mu_x_hat=1/c_x_hat*(c_x*mu_0+x[i,])
			b_x_hat=b_x+(x[i,]^2+mucmu_x-c_x_hat*mu_x_hat^2)/2

			sigmay_new=rinvgamma(1,a_y_hat, b_y_hat)
			phi_y_new=rmvnorm(1, mu_hat, sigmay_new*iC_hat,  method=c("chol"))
			phi_update=cbind(phi_update,t(phi_y_new))
			sigmay_update=c(sigmay_update,sigmay_new)
			sigmax_new=rinvgamma(p,a_x_hat, b_x_hat)
			mu_x_new=rnorm(p, mu_x_hat, (sigmax_new/c_x_hat)^.5)
			mux_update=cbind(mux_update,mu_x_new)
			sigmax_update=cbind(sigmax_update,sigmax_new)
		}	
	}


	config_s=c_update
	dn=length(unique(c_update))

	#Draw new parameters for each cluster
	phi_s=matrix(0,p+1,dn)
	sigma_y_s=matrix(0,1,dn)
	mu_x_s=matrix(0,p,dn)
	sigma_x_s=matrix(0,p,dn)
	for(c in 1:dn){

		#observations in class c
		n_c=sum(config_s==c)
		X_c=matrix(X[config_s==c,],nrow=n_c)
		x_c=matrix(x[config_s==c],nrow=n_c)
		y_c=y[config_s==c]

		#for intercept/slope
		C_hat=C+ t(X_c)%*%X_c
		decomp=eigen(C_hat, symmetric=TRUE)
		iC_hat=decomp$vectors%*%diag(1/decomp$values)%*%t(decomp$vectors)
		mu_hat=iC_hat%*%(C%*%mu_theta+t(X_c)%*%y_c)
		a_y_hat=a_y+n_c/2
		b_y_hat=b_y+(sum(y_c^2)+muCmu-t(mu_hat)%*%C_hat%*%mu_hat)/2

		#draw value from posterior
		sigma_y_s[c]=rinvgamma(1,a_y_hat, b_y_hat)
		phi_s[,c]=rmvnorm(1, mu_hat, sigma_y_s[c]*iC_hat,  method=c("chol"))
	
		#for mu_x
		c_x_hat= c_x+n_c
		mux_hat= 1/c_x_hat*(c_x*mu_0+ colSums(x_c))
		
		#for sigma_x
		a_x_hat=a_x+n_c/2
		b_x_hat=b_x+(colSums(x_c^2)+mucmu_x-c_x_hat*mux_hat^2)/2

		#draw value from posterior
		sigma_x_s[,c]=rinvgamma(p,a_x_hat, b_x_hat)
		mu_x_s[,c]=rnorm(p,mux_hat, (sigma_x_s[,c]/c_x_hat)^.5)

	}

	###Draw other parameters

	### alpha

	#updated parameters
	nu=rbeta(1, alpha_s+1,n)
	v_alpha_post=v_alpha-log(nu)
	u_alpha_post=u_alpha+dn
	if(runif(1)<((n*v_alpha_post)/(u_alpha_post-1+n*v_alpha_post))){
		u_alpha_post=u_alpha_post-1}
	
	#draw from posterior
	alpha_s=rgamma(1,u_alpha_post, v_alpha_post)

	if (s%%10==0){
	print( paste("Number of iterations completed=", s) )
	print( paste("k=", dn) )
	#par(mfrow=c(2,3))
	#k_s=max(config_s)
	#plot(x[,1],y)
	#for(i in 1:k_s){
	#points(x[config_s==i,1],y[config_s==i],col=i)
	#}

	#plot(x[,2],y)
	#for(i in 1:k_s){
	#points(x[config_s==i,2],y[config_s==i],col=i)
	#}

	#plot(x[,3],y)
	#for(i in 1:k_s){
	#points(x[config_s==i,3],y[config_s==i],col=i)
	#}
	
	#plot(x[,1],x[,2])
	#for(i in 1:k_s){
	#points(x[config_s==i,1],x[config_s==i,2],col=i)
	#}

	#plot(x[,1],x[,3])
	#for(i in 1:k_s){
	#points(x[config_s==i,1],x[config_s==i,3],col=i)
	#}

	#plot(x[,2],x[,3])
	#for(i in 1:k_s){
	#points(x[config_s==i,2],x[config_s==i,3],col=i)
	#}
	}


	#if s is bigger than burnin, save the output
	if (s>=burnin){
		config[s+1-burnin,]=config_s
		alpha[s+1-burnin,]=alpha_s
		phi_y[[s+1-burnin]]=phi_s
		sigma_y[[s+1-burnin]]=sigma_y_s
		mu_x[[s+1-burnin]]=mu_x_s
		sigma_x[[s+1-burnin]]=sigma_x_s	
	}

}

#return output
output=list(config = config, phi_y=phi_y, sigma_y=sigma_y, mu_x=mu_x, sigma_x=sigma_x, alpha=alpha)
return( output)

}
