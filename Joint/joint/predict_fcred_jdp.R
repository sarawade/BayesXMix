# Computes prediction and predictive density estimate for a new subject with credible intervals
# for EDP regression model

###INPUT
# a: Credible level for credible intervals of mean
# S: number of samples
# m2: length of y_grid
# p: number of covariates
# x_new: covariate of new subjects (1xp)
# y_grid: grid of y values to estimate predictive density (m2x1)
# mu_0, c_x: prior parameters for mu  (px1 and px1)
# a_x, b_x: prior parameters for sigma^2_x (px1 and px1)
# config: sampled configurations (Sxn)
# phi_y: sampled y cluster parameters (list with S elements, each element is p+1xk[s])
# mu_x, sigma_x: sampled x cluster parameters (list with S elements, each element is pxk[s] and pxk[s])
# alpha: sampled precision parameters (Sx1)

###OUTPUT
# f_pred: prediction of response for new subject (m2x1)
# l_fpred: lower bound for pointwise credible interval of predictive density (m2x1)
# u_fpred: upper bound for pointwise credible interval of predictive density (m2x1)


predict_fcred_jdp=function(a, S,m,m2,p,x_new, y_grid,mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k, config, alpha, phi_y, sigma_y, mu_x, sigma_x ){
	
#Create X matrix
x_new=matrix(x_new,m,p)
X_new=matrix(c(rep(1,m),x_new),ncol=(p+1))

#marginal for x
c_x_hat=c_x+1
x_sc=(x_new-t(matrix(mu_0,nrow=p,ncol=m)))*t(matrix((c_x*a_x/(c_x_hat*b_x))^.5,nrow=p,ncol=m))
marg_x=exp(rowSums(dt(x_sc,2*t(matrix(a_x,nrow=p,ncol=m)),log=TRUE)+t(matrix(.5*(log(c_x)-log(c_x_hat)-log(b_x)+log(a_x)),nrow=p,ncol=m))))

#marginal for y
marg_y=matrix(0,m2,m)
for(j in 1:m){
C_hat=C+X_new[j,]%*%t(X_new[j,])
decomp=eigen(C_hat, symmetric=TRUE)
iC_hat=decomp$vectors%*%diag(1/decomp$values)%*%t(decomp$vectors)
sc=1-X_new[j,]%*%iC_hat%*%X_new[j,]
yi_hat=X_new[j,]%*%mu_theta
y_sc=(y_grid-yi_hat)*(a_y*sc/b_y)^.5
marg_y[,j]=dt(y_sc,2*a_y)*(a_y*sc/b_y)^.5
}

#Initialize
f_pred_mat=matrix(0,m2,S) # for the first new data point
weight_mat=matrix(0,m,S)
weight_list = list(S)

for(s in 1:S){

	#calculate weight matrix
	weight_mat_s=matrix(0,m,k[s]+1)

	for(j in 1:k[s]){
		nj_s=sum(config[s,]==j)
		weight_mat_s[,j]=nj_s*exp(rowSums(dnorm(x_new, t(matrix(mu_x[[s]][,j], nrow=p,ncol=m)), t(matrix(sigma_x[[s]][,j]^.5,nrow=p,ncol=m)),log=T)))
	}
	weight_mat[,k[s]+1]=alpha[s]*marg_x
	
	#calculate prediction
	f_pred_s=dnorm(matrix(y_grid, nrow=m2, ncol=k[s]),matrix(1,m2,1)%*%X_new[1,]%*%phi_y[[s]],matrix(1,m2,1)%*%sigma_y[[s]]^.5 )
	f_pred_s=cbind(f_pred_s, marg_y[,1])
	
	weight_list[[s]] = weight_mat_s
	weight_mat[,s]= rowSums(weight_mat_s)
	f_pred_mat[,s]= f_pred_s%*%weight_mat_s[1,]/weight_mat[1,s]

	if(((s/S*100)%%1)==0){
		print(paste(s/S*100,"% completed"))
	}
}


###Compute Pred and Cred Intervals
f_pred=matrix(0,m2,m)
u_fpred=matrix(0,m2,m)
l_fpred=matrix(0,m2,m)

for(i in 1:m){
	f_pred[,i]=f_pred_mat%*%weight_mat[i,]/sum(weight_mat[i,])

	#calculate upper and lower 95% credible bounds for f_pred
	for(j in 1:m2){
		f_pred_i_sort=sort(f_pred_mat[j,], index.return=T)
		cweight_i_sort=cumsum(weight_mat[i,f_pred_i_sort$ix])/sum(weight_mat[i,])
		ind_l=sum(cweight_i_sort<a/2)
		l_fpred[j,i]=(f_pred_i_sort$x[ind_l]*(a/2-cweight_i_sort[ind_l])+f_pred_i_sort$x[ind_l+1]*(cweight_i_sort[ind_l+1]-a/2))/(cweight_i_sort[ind_l+1]-cweight_i_sort[ind_l])
		ind_u=sum(cweight_i_sort<(1-a/2))
		u_fpred[j,i]=(f_pred_i_sort$x[ind_u]*(1-a/2-cweight_i_sort[ind_u])+f_pred_i_sort$x[ind_u+1]*(cweight_i_sort[ind_u+1]-1+a/2))/(cweight_i_sort[ind_u+1]-cweight_i_sort[ind_u])	
	}

	if(i!=m){
		#compute f_pred_mat
		for(s in 1:S){
		  f_pred_s=dnorm(matrix(y_grid, nrow=m2, ncol=k[s]),matrix(1,m2,1)%*%X_new[i+1,]%*%phi_y[[s]],matrix(1,m2,1)%*%sigma_y[[s]]^.5 )
		  f_pred_s=cbind(f_pred_s, marg_y[,i+1])
		
			f_pred_mat[,s]=f_pred_s%*%weight_list[[s]][i+1,]/weight_mat[i+1,s]
		}
	}
}


output=list(f_pred=f_pred,l_fpred=l_fpred, u_fpred=u_fpred)
return(output)
}
