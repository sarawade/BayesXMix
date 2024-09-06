### Function that obtains posterior samples from EDP regression model

###INPUT
# S:number of iterations to save
# burnin: number of iterations to discard
# y: observed response (nx1)
# x: observed covariates (nxp)
# mu_theta, C: prior parameters for beta (p+1x1 and p+1xp+1)
# a_y, b_y: prior parameters for sigma^2_y (1x1 and 1x1)
# mu_0, c_x: prior parameters for mu  (px1 and px1)
# a_x, b_x: prior parameters for sigma^2_x (px1 and px1)
# u_alpha_x, v_alpha_x: prior parameters for alpha_x (1x1 and 1x1)
# u_alpha_y, v_alpha_y: prior parameters for alpha_y (1x1 and 1x1)
# config_y_init: initial configuration of y (1xn)
# config_x_init: initial configuration of x within y clusters (1xn)
# phi_init: initali values of beta for each y cluster (p+1xkn_y)
# sigma_y_init: initial values of sigma_y for each y cluster (1xkn_y)
# mu_x_init, sigma_x_init: initial values of mu_x and sigma_x for each y-x cluster (list with kn_y elements, each element is pxkn_x[j])
# alpha_y_init, alpha_x_init: initial values of mass parameters (1x1, 1xkn_y)

###OUPUT
# config_y, config_x: sampled configurations (Sxn and Sxn)
# phi_y, sigma_y: sampled y cluster parameters (list with S elements, each element is p+1xkn_y[s] and 1xkn_y[s])
# mu_x, sigma_x: sampled x cluster parameters (list with S elements, each element is a list with kn_y[s] elements, each element is pxkn_x[s,j] and pxkn_x[s,j])
# alpha_y: sampled y precision parameters (Sx1)
# alpha_x: sample x precision parameters (list with S elements, each element is 1xkn_y[s])

edp_mcmc=function(S, burnin, y, x, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, u_alpha_x, v_alpha_x, u_alpha_y, v_alpha_y, config_y_init, config_x_init, phi_init, sigma_y_init, mu_x_init, sigma_x_init,alpha_y_init, alpha_x_init){

#define n and p
n=dim(x)[1]
p=dim(x)[2]

#define X
X=cbind(rep(1,n), x)

#define iC
iC=solve(C)

#Empty matrices to store MCMC output

config_x=matrix(0,S,n)
config_y=matrix(0,S,n)
phi_y=list()
sigma_y=list()
mu_x=list()
sigma_x=list()
alpha_x=list()
alpha_y=matrix(0, S)

#Initialize
config_x_s=config_x_init
config_y_s=config_y_init
alpha_y_s=alpha_y_init
alpha_x_s=alpha_x_init
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
	c_update_x=config_x_s
	c_update_y=config_y_s

	#unique parameters to be updated if new value is chosen
	phi_update=phi_s
	sigmay_update=sigma_y_s
	sigmax_update=sigma_x_s
	mux_update=mu_x_s
	alpha_x_update=alpha_x_s

	#updated hyperparameters that dont depend on i	
	a_y_hat=a_y+1/2
	a_x_hat=a_x+1/2
	c_x_hat=c_x+1

	for(i in 1:n){
	
		#Define configuration without i
		c_ni_y=c_update_y[-i]
		c_i_y=c_update_y[i]
		uniq_c_ni_y=unique(c_ni_y)
		dn_ni_y=length(uniq_c_ni_y)

		c_ni_x=c_update_x[-i]
		c_i_x=c_update_x[i]

		#check if c_i is an x class with the single element i
		singleton_y=(sum(c_ni_y==c_i_y)==0)

		if(singleton_y){
			#remove the parameters for class c_i
			phi_update[,c_i_y]=phi_update[,dn_ni_y+1]
			phi_update=phi_update[,-(dn_ni_y+1)]
			sigmay_update[c_i_y]=sigmay_update[dn_ni_y+1]
			sigmay_update=sigmay_update[-(dn_ni_y+1)]
			sigmax_update[[c_i_y]]=sigmax_update[[dn_ni_y+1]]
			sigmax_update=sigmax_update[-(dn_ni_y+1)]
			mux_update[[c_i_y]]=mux_update[[dn_ni_y+1]]
			mux_update=mux_update[-(dn_ni_y+1)]
			#change c_ni, updated config
			c_ni_y[c_ni_y==(dn_ni_y+1)]=c_i_y
			c_update_y[c_update_y==(dn_ni_y+1)]=c_i_y
			c_update_y[i]=dn_ni_y+1	
			alpha_x_update[c_i_y]=alpha_x_update[dn_ni_y+1]
			alpha_x_update=alpha_x_update[-(dn_ni_y+1)]
			c_i_y=dn_ni_y+1
		}

		dn_ni_x=rep(0,dn_ni_y)
		for(j in 1:dn_ni_y){
			dn_ni_x[j]=length(unique(c_ni_x[c_ni_y==j]))
		}

	
		#check if c_i is a y class with the single element i
		singleton_x=(sum(c_ni_x[c_ni_y==c_i_y]==c_i_x)==0)
		
		if(singleton_x&&!singleton_y){
			mux_update[[c_i_y]][,c_i_x]=mux_update[[c_i_y]][,dn_ni_x[c_i_y]+1]
			mux_update[[c_i_y]]=mux_update[[c_i_y]][,-(dn_ni_x[c_i_y]+1)]
			sigmax_update[[c_i_y]][,c_i_x]=sigmax_update[[c_i_y]][,dn_ni_x[c_i_y]+1]
			sigmax_update[[c_i_y]]=sigmax_update[[c_i_y]][,-(dn_ni_x[c_i_y]+1)]
			#change c_ni, updated config
			c_ni_x[c_ni_y==c_i_y][c_ni_x[c_ni_y==c_i_y]==(dn_ni_x[c_i_y]+1)]=c_i_x
			c_update_x[c_update_y==c_i_y][c_update_x[c_update_y==c_i_y]==(dn_ni_x[c_i_y]+1)]=c_i_x
			c_update_x[i]=dn_ni_x[c_i_y]+1
			c_i_x=dn_ni_x[c_i_y]+1
		}
	
		#Set value of config for observation i
		log_prob_y=matrix(0,dn_ni_y+1,1)
		log_y=matrix(0,dn_ni_y+1,1)
		log_x=list()
		log_prob_x=list()
	
		#calculate probabilities for new x
		x_sc=(x[i,]-mu_0)*(c_x*a_x/(c_x_hat*b_x))^.5
		log_new_x=sum(dt(x_sc,2*a_x, log=TRUE)+.5*(log(c_x)-log(c_x+1)-log(b_x)+log(a_x)))
	
		#calculate probability of a new y
		#updated parameters
		C_hat=C+X[i,]%*%t(X[i,])
		decomp=eigen(C_hat, symmetric=TRUE)
		iC_hat=decomp$vectors%*%diag(1/decomp$values)%*%t(decomp$vectors)
	
		sc=1-X[i,]%*%iC_hat%*%X[i,]
		yi_hat=X[i,]%*%mu_theta
		y_sc=(y[i]-yi_hat)*(a_y*sc/b_y)^.5
		log_y[dn_ni_y+1]=dt(y_sc,2*a_y,log=TRUE)+.5*(log(sc)-log(b_y)+log(a_y))
		log_prob_y[dn_ni_y+1]=log(alpha_y_s/(alpha_y_s+n-1))+log_y[dn_ni_y+1]


		#calculate probabilities for old x values
		phi_update=matrix(phi_update, nrow=(p+1),ncol=dn_ni_y)
		for(j in 1:dn_ni_y){
			n_j_ni=sum(c_ni_y==j)
			z_ij=sum(phi_update[,j]*X[i,])
			log_y[j]=dnorm(y[i],z_ij,sigmay_update[j]^.5, log=TRUE)
			log_prob_y[j]= log(n_j_ni/(alpha_y_s+n-1))+log_y[j]
		
			log_x[[j]]=matrix(0,dn_ni_x[j]+1,1)
			log_prob_x[[j]]=matrix(0,dn_ni_x[j]+1,1)

			mux_update[[j]]=matrix(mux_update[[j]],nrow=(p), ncol=dn_ni_x[j])
			sigmax_update[[j]]=matrix(sigmax_update[[j]],nrow=p, ncol=dn_ni_x[j])

			for(h in 1:dn_ni_x[j]){
				n_jh_ni=sum(c_ni_x[c_ni_y==j]==h)
				log_x[[j]][h]=sum(dnorm(x[i,], mux_update[[j]][,h], sigmax_update[[j]][,h]^.5, log=TRUE))
				log_prob_x[[j]][h]=log(n_jh_ni/(alpha_x_update[j]+n_j_ni))+log_x[[j]][h]
			}
		
			log_x[[j]][dn_ni_x[j]+1]=log_new_x
			log_prob_x[[j]][dn_ni_x[j]+1]= log(alpha_x_update[j]/(alpha_x_update[j]+n_j_ni))+log_x[[j]][dn_ni_x[j]+1]
		}

		log_x[[dn_ni_y+1]]=log_new_x
		log_prob_x[[dn_ni_y+1]]=log_x[[dn_ni_y+1]]

	
		#compute log vector
		log_xy_prob=rep(0, sum(dn_ni_x+1)+1)
		cum_dn_ni_x=cumsum(c(0,dn_ni_x+1,1))
		for(j in 1:(dn_ni_y+1)){
			log_xy_prob[(cum_dn_ni_x[j]+1):(cum_dn_ni_x[j+1])]=log_prob_y[j]+log_prob_x[[j]]
		}
	
		#need to find normalizing constant
		#f <- function (k) sum(exp(k+log_xy_prob))-1
		#k_star=uniroot(f,lower=-10^(10),upper=10^20)$root
		#probs=exp(k_star+log_xy_prob)
		probs=exp(log_xy_prob)/sum(exp(log_xy_prob))

		#set value of config according to probability	
		config_chosen=sum(runif(1)>c(cumsum(probs[-length(probs)]),1))+1
		c_update_y[i]=sum(cum_dn_ni_x<config_chosen)
		c_update_x[i]=config_chosen-cum_dn_ni_x[c_update_y[i]]	

		
		#update parameters if new chosen
		#for y
		if(c_update_y[i]==(dn_ni_y+1)){
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
			mux_update[[c_update_y[i]]]=matrix(mu_x_new,nrow=p)
			sigmax_update[[c_update_y[i]]]=matrix(sigmax_new,nrow=p)
		
			#if new then n_i=1, k_i=1, so just draw from prior for alpha y
			alpha_x_new=rgamma(1,u_alpha_x, v_alpha_x)
			alpha_x_update=c(alpha_x_update, alpha_x_new)
			}	
		else{
			if(c_update_x[i]==(dn_ni_x[c_update_y[i]]+1)){
				#updated parameters
				mu_x_hat=1/c_x_hat*(c_x*mu_0+x[i,])
				b_x_hat=b_x+(x[i,]^2+mucmu_x-c_x_hat*mu_x_hat^2)/2

				sigmax_new=rinvgamma(p,a_x_hat, b_x_hat)
				mu_x_new=rnorm(p, mu_x_hat, (sigmax_new/c_x_hat)^.5)
				sigmax_update[[c_update_y[i]]]=cbind(sigmax_update[[c_update_y[i]]],sigmax_new)
				mux_update[[c_update_y[i]]]=cbind(mux_update[[c_update_y[i]]],(mu_x_new))
			}
		}
	}

	config_x_s=c_update_x
	config_y_s=c_update_y
	dn_y=length(unique(c_update_y))
	dn_x=rep(0,dn_y)
	for(j in 1:dn_y){
	dn_x[j]=length(unique(c_update_x[c_update_y==j]))
	}
	alpha_x_s=alpha_x_update

	#Draw new parameters for each cluster
	phi_s=matrix(0,p+1,dn_y)
	sigma_y_s=matrix(0,1,dn_y)
	mu_x_s=list()
	sigma_x_s=list()
	for(c in 1:dn_y){
		mu_x_s[[c]]=matrix(0, p, dn_x[c])
		sigma_x_s[[c]]=matrix(0,p,dn_x[c])

		#observations in class c
		n_c=sum(config_y_s==c)
		X_c=matrix(X[config_y_s==c,],nrow=n_c)
		y_c=y[config_y_s==c]

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
	
	
		for(h in 1:dn_x[c]){
		
		#Draw a value for phi_y_c
		n_ch=sum(config_x_s[config_y_s==c]==h)
		x_ch=matrix(matrix(x[config_y_s==c,], nrow=n_c)[config_x_s[config_y_s==c]==h,], nrow=n_ch)

		#for mu_x
		c_x_hat= c_x+n_ch
		mux_hat= 1/c_x_hat*(c_x*mu_0+ colSums(x_ch))
		
		#for sigma_x
		a_x_hat=a_x+n_ch/2
		b_x_hat=b_x+(colSums(x_ch^2)+mucmu_x-c_x_hat*mux_hat^2)/2

		#draw value from posterior
		sigma_x_s[[c]][,h]=rinvgamma(p,a_x_hat, b_x_hat)
		mu_x_s[[c]][,h]=rnorm(p,mux_hat, (sigma_x_s[[c]][,h]/c_x_hat)^.5)
		}
	}

	#Propose x-cluster switching
	##Case 1
	if(dn_y>1 & sum(dn_x>1)>0){
		# Choose x cluster among x-clusters in y-clusters with more than 1 x-clusters
		kb1=sum(dn_x[dn_x>1])
		u=runif(1)*kb1
		lj=1
		while(u>lj){lj=lj+1}
		j=sum(lj>cumsum(dn_x*(dn_x>1)))+1
		l=lj-c(0,cumsum(dn_x*(dn_x>1)))[j]
		#Choose y cluster other than j
		u=runif(1)*(dn_y-1)
		h=1
		while(u>h){h=h+1}
		h=h+(h>=j)
		#Calculate acceptance prob
		nj=sum(config_y_s==j)
		nh=sum(config_y_s==h)
		njl=sum(config_x_s[config_y_s==j]==l)
		laccep=lgamma(nj-njl)-lgamma(nj)+lgamma(nh+njl)-lgamma(nh)+lgamma(alpha_x_s[j]+nj)-lgamma(alpha_x_s[j]+nj-njl)+lgamma(alpha_x_s[h]+nh)-lgamma(alpha_x_s[h]+nh+njl)+log(alpha_x_s[h])-log(alpha_x_s[j])
		yjl=matrix(matrix(y[config_y_s==j], nrow=nj)[config_x_s[config_y_s==j]==l,], nrow=njl)
		Xjl=matrix(matrix(X[config_y_s==j,], nrow=nj)[config_x_s[config_y_s==j]==l,], nrow=njl)
		laccep=laccep+sum(dnorm(yjl,Xjl%*%phi_s[,h], sigma_y_s[h]^.5,log=T))-sum(dnorm(yjl,Xjl%*%phi_s[,j], sigma_y_s[j]^.5,log=T))
		laccep=laccep+log(kb1)-log(kb1-(dn_x[j]==2)+(dn_x[h]==1))
		accep=min(1,exp(laccep))
		if(runif(1)<accep){
			print("Switch 1 accepted")
			#Switch labels
			indjl=(config_y_s==j&config_x_s==l)
			config_y_s[indjl]=h
			config_x_s[indjl]=dn_x[h]+1
			if(dn_x[j]>l){config_x_s[(config_y_s==j&config_x_s>l)]=config_x_s[(config_y_s==j&config_x_s>l)]-1}
			#Update x parameters
			mu_x_s[[h]]=cbind(mu_x_s[[h]],mu_x_s[[j]][,l])
			sigma_x_s[[h]]=cbind(sigma_x_s[[h]],sigma_x_s[[j]][,l])
			mu_x_s[[j]]=matrix(mu_x_s[[j]][,-l],nrow=p)
			sigma_x_s[[j]]=matrix(sigma_x_s[[j]][,-l],nrow=p)
			#Update dn_x
			dn_x[j]=dn_x[j]-1
			dn_x[h]=dn_x[h]+1
		}
	}

	## Case 2 and 3
	u23=runif(1)
	ind2=sum(dn_x>1)>0
	ind3=(sum(dn_x==1)>0)&(dn_y>1)
	if(u23<(.5^ind3) & ind2){
		# Choose x cluster among x-clusters in y-clusters with more than 1 x-clusters
		kb1=sum(dn_x[dn_x>1])
		u=runif(1)*kb1
		lj=1
		while(u>lj){lj=lj+1}
		j=sum(lj>cumsum(dn_x*(dn_x>1)))+1
		l=lj-c(0,cumsum(dn_x*(dn_x>1)))[j]
		#Sample new phi
		nj=sum(config_y_s==j)
		njl=sum(config_x_s[config_y_s==j]==l)
		yjl=matrix(matrix(y[config_y_s==j], nrow=nj)[config_x_s[config_y_s==j]==l,], nrow=njl)
		Xjl=matrix(matrix(X[config_y_s==j,], nrow=nj)[config_x_s[config_y_s==j]==l,], nrow=njl)
		C_hat=C+ t(Xjl)%*%Xjl
		decomp=eigen(C_hat, symmetric=TRUE)
		iC_hat=decomp$vectors%*%diag(1/decomp$values)%*%t(decomp$vectors)
		mu_hat=iC_hat%*%(C%*%mu_theta+t(Xjl)%*%yjl)
		a_y_hat=a_y+njl/2
		b_y_hat=b_y+(sum(yjl^2)+muCmu-t(mu_hat)%*%C_hat%*%mu_hat)/2
		sigma_y_new=rinvgamma(1,a_y_hat, b_y_hat)
		phi_new=rmvnorm(1, mu_hat, sigma_y_new*iC_hat,  method=c("chol"))
		#Sample new alpha_x
		alpha_x_new=rgamma(1,u_alpha_x, v_alpha_x)
		laccep=log(alpha_y_s)+lgamma(nj-njl)-lgamma(nj)+lgamma(alpha_x_s[j]+nj)-lgamma(alpha_x_s[j]+nj-njl)+lgamma(njl)+lgamma(alpha_x_new)-lgamma(alpha_x_new+njl)+log(alpha_x_new)-log(alpha_x_s[j])
		sc=diag(njl)-Xjl%*%iC_hat%*%t(Xjl)
		decomp=eigen(sc, symmetric=TRUE)
		isc=decomp$vectors%*%diag(1/decomp$values,njl)%*%t(decomp$vectors)
		yjl_hat=Xjl%*%mu_theta
		laccep=laccep+dmvt(c(yjl-yjl_hat), sigma=(b_y/a_y*isc), df = 2*a_y, log = TRUE)-sum(dnorm(yjl,Xjl%*%phi_s[,j], sigma_y_s[j]^.5,log=T))
		laccep=laccep+log(kb1)-log(sum(dn_x==1)+1+(dn_x[j]==2))-log(dn_y)
		accep=min(1,exp(laccep))
		if(runif(1)<accep){
			print("Switch 2 accepted")
			#Switch labels
			indjl=(config_y_s==j&config_x_s==l)
			config_y_s[indjl]=dn_y+1
			config_x_s[indjl]=1
			if(dn_x[j]>l){config_x_s[(config_y_s==j&config_x_s>l)]=config_x_s[(config_y_s==j&config_x_s>l)]-1}
			#Update x parameters
			mu_x_s[[dn_y+1]]=matrix(mu_x_s[[j]][,l],nrow=p)
			sigma_x_s[[dn_y+1]]=matrix(sigma_x_s[[j]][,l],nrow=p)
			mu_x_s[[j]]=matrix(mu_x_s[[j]][,-l],nrow=p)
			sigma_x_s[[j]]=matrix(sigma_x_s[[j]][,-l],nrow=p)
			#Update y parameters
			phi_s=cbind(phi_s,t(phi_new))
			sigma_y_s=c(sigma_y_s,sigma_y_new)
			#Update alpha
			alpha_x_s=c(alpha_x_s,alpha_x_new)
			#Update dn_y
			dn_y=dn_y+1
			#Update dn_x
			dn_x[j]=dn_x[j]-1
			dn_x=c(dn_x,1)
		}
	}
	if(u23>(.5*ind2) & ind3){
	# Choose x cluster among x-clusters in y-clusters with 1 x-clusters
		k1=sum(dn_x==1)
		u=runif(1)*k1
		lj=1
		while(u>lj){lj=lj+1}
		j=sum(lj>cumsum(dn_x*(dn_x==1)))+1
		l=1
		#Choose y cluster other than j
		u=runif(1)*(dn_y-1)
		h=1
		while(u>h){h=h+1}
		h=h+(h>=j)
		#Compute acceptance probability
		nh=sum(config_y_s==h)
		njl=sum(config_y_s==j)
		laccep=-log(alpha_y_s)+lgamma(nh+njl)-lgamma(nh)+lgamma(alpha_x_s[h]+nh)-lgamma(alpha_x_s[h]+nh+njl)-lgamma(njl)-lgamma(alpha_x_s[j])+lgamma(alpha_x_s[j]+njl)+log(alpha_x_s[h])-log(alpha_x_s[j])
		yjl=matrix(y[config_y_s==j], nrow=njl)
		Xjl=matrix(X[config_y_s==j,], nrow=njl)
		C_hat=C+ t(Xjl)%*%Xjl
		decomp=eigen(C_hat, symmetric=TRUE)
		iC_hat=decomp$vectors%*%diag(1/decomp$values)%*%t(decomp$vectors)
		sc=diag(njl)-Xjl%*%iC_hat%*%t(Xjl)
		decomp=eigen(sc, symmetric=TRUE)
		isc=decomp$vectors%*%diag(1/decomp$values,njl)%*%t(decomp$vectors)
		yjl_hat=Xjl%*%mu_theta
		laccep=laccep-dmvt(c(yjl-yjl_hat), sigma=(b_y/a_y*isc), df = 2*a_y, log = TRUE)+sum(dnorm(yjl,Xjl%*%phi_s[,h], sigma_y_s[h]^.5,log=T))
		laccep=laccep+log(k1)+log(dn_y-1)-log(sum(dn_x>1)+1+(dn_x[h]==1))
		accep=min(1,exp(laccep))
		if(runif(1)<accep){
			print("Switch 3 accepted")
			#Switch labels
			indjl=(config_y_s==j)
			config_y_s[indjl]=h
			config_x_s[indjl]=dn_x[h]+1
			#Update x parameters
			mu_x_s[[h]]=cbind(mu_x_s[[h]],mu_x_s[[j]][,l])
			sigma_x_s[[h]]=cbind(sigma_x_s[[h]],sigma_x_s[[j]][,l])
			if(j<dn_y){
				config_y_s[config_y_s>j]=config_y_s[config_y_s>j]-1
				mu_x_s[j:(dn_y-1)]=mu_x_s[(j+1):(dn_y)]
				sigma_x_s[j:(dn_y-1)]=sigma_x_s[(j+1):(dn_y)]
				phi_s[j:(dn_y-1)]=phi_s[(j+1):(dn_y)]
				sigma_y_s[j:(dn_y-1)]=sigma_y_s[(j+1):(dn_y)]
				alpha_x_s[j:(dn_y-1)]=alpha_x_s[(j+1):(dn_y)]
				dn_x[j:(dn_y-1)]=dn_x[(j+1):(dn_y)]
			}
			#Remove last
			mu_x_s=mu_x_s[-(dn_y)]
			sigma_x_s=sigma_x_s[-(dn_y)]
			phi_s=phi_s[,-(dn_y)]
			sigma_y_s=sigma_y_s[-(dn_y)]
			alpha_x_s=alpha_x_s[-(dn_y)]
			dn_x=dn_x[-(dn_y)]
			#Update dn_y
			dn_y=dn_y-1
		}
	}

	###Draw other parameters

	### alpha

	#updated parameters
	nu=rbeta(1, alpha_y_s+1,n)
	v_alpha_post=v_alpha_y-log(nu)
	u_alpha_post=u_alpha_y+dn_y
	if(runif(1)<((n*v_alpha_post)/(u_alpha_post-1+n*v_alpha_post))){
		u_alpha_post=u_alpha_post-1}
	
	#draw from posterior
	alpha_y_s=rgamma(1,u_alpha_post, v_alpha_post)

	for(c in 1:dn_y){
		n_c=sum(config_y_s==c)
		nu=rbeta(1, alpha_x_s[c]+1,n_c)
		v_alpha_post=v_alpha_x-log(nu)
		u_alpha_post=u_alpha_x+dn_x[c]
		if(runif(1)<((n_c*v_alpha_post)/(u_alpha_post-1+n_c*v_alpha_post))){
			u_alpha_post=u_alpha_post-1}
	
		#draw from posterior
		alpha_x_s[c]=rgamma(1,u_alpha_post, v_alpha_post)	
	}	

	if (s%%10==0){
	print( paste("Number of iterations completed=", s) )
	print( paste("Number of y-clusters=", dn_y) )
	#n_c=rep(0,dn_y)
	#for(c in 1:dn_y){
	#	n_c[c]=sum(config_y_s==c)
	#}
	#print( paste("Cluster sizes=", n_c) )
	#par(mfrow=c(2,3))
	#k_s=max(config_y_s)
	#plot(x[,1],y)
	#for(i in 1:k_s){
	#points(x[config_y_s==i,1],y[config_y_s==i],col=i)
	#}

	#plot(x[,2],y)
	#for(i in 1:k_s){
	#points(x[config_y_s==i,2],y[config_y_s==i],col=i)
	#}

	#plot(x[,3],y)
	#for(i in 1:k_s){
	#points(x[config_y_s==i,3],y[config_y_s==i],col=i)
	#}
	
	#plot(x[,1],x[,2])
	#for(i in 1:k_s){
	#points(x[config_y_s==i,1],x[config_y_s==i,2],col=i)
	#}

	#plot(x[,1],x[,3])
	#for(i in 1:k_s){
	#points(x[config_y_s==i,1],x[config_y_s==i,3],col=i)
	#}

	#plot(x[,2],x[,3])
	#for(i in 1:k_s){
	#points(x[config_y_s==i,2],x[config_y_s==i,3],col=i)
	#}
	}

	#if s is bigger than burnin, save the output
	if (s>=burnin){
		config_x[s+1-burnin,]=config_x_s
		config_y[s+1-burnin,]=config_y_s
		alpha_y[s+1-burnin,]=alpha_y_s
		alpha_x[[s+1-burnin]]=alpha_x_s
		phi_y[[s+1-burnin]]=phi_s
		sigma_y[[s+1-burnin]]=sigma_y_s
		mu_x[[s+1-burnin]]=mu_x_s
		sigma_x[[s+1-burnin]]=sigma_x_s	
	}

}

#return output
output=list(config_y = config_y, config_x = config_x, phi_y=phi_y, sigma_y=sigma_y, mu_x=mu_x, sigma_x=sigma_x, alpha_y=alpha_y, alpha_x=alpha_x)
return( output)

} 
