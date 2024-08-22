# Computes prediction and predictive density estimated for EDP regression model

###INPUT
# S: number of samples
# m: number of new subjects
# m2: length of y_grid
# p: number of covariates
# x_new: covariates of new subjects (mxp)
# y_grid: grid of y values to estimate predictive density (m2x1)
# mu_theta, C: prior parameters for beta (p+1x1 and p+1xp+1)
# a_y, b_y: prior parameters for sigma^2_y (1x1 and 1x1)
# mu_0, c_x: prior parameters for mu  (px1 and px1)
# a_x, b_x: prior parameters for sigma^2_x (px1 and px1)
# config: sampled configurations (Sxn)
# beta_y: sampled y cluster parameters (list with S elements, each element is p+1xk[s])
# mu_x, sigma_x: sampled x cluster parameters (list with S elements, each element is pxk[s] and pxk[s])
# alpha: sampled precision parameters (Sx1)

###OUTPUT
# y_pred: prediction of response for new subjects
# f_pred: predictive density of response for new subjects

predict_jdp=function(S,m,m2,p,x_new, y_grid, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k, config, alpha, phi_y, sigma_y, mu_x, sigma_x ){

#Create X matrix
X_new=matrix(c(rep(1,m),x_new),ncol=(p+1))

#Initialize
y_pred=matrix(0,m,1)
f_pred=matrix(0,m2,m)
norm_const=matrix(0,m,1)

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

for(s in 1:S){

	#calculate weight matrix
	weight_mat_s=matrix(0,m,k[s]+1)
	weight_mat_s[,1]=alpha[s]/(alpha[s]+n)*marg_x
	f_pred=f_pred+marg_y%*%diag(weight_mat_s[,1])

	for(j in 1:k[s]){
		nj_s=sum(config[s,]==j)
		weight_mat_s[,j+1]=nj_s/(alpha[s]+n)*exp(rowSums(dnorm(x_new, t(matrix(mu_x[[s]][,j], nrow=p,ncol=m)), t(matrix(sigma_x[[s]][,j]^.5,nrow=p,ncol=m)),log=T)))
		f_pred=f_pred+dnorm(matrix(y_grid,nrow=m2,ncol=m),t(matrix(X_new%*%phi_y[[s]][,j],nrow=m,ncol=m2)),sigma_y[[s]][j]^.5 )%*%diag(weight_mat_s[,j+1])
	}
	#calculate prediction
	y_pred_s=X_new%*%cbind(mu_theta,phi_y[[s]])
	y_pred=y_pred+rowSums(y_pred_s*weight_mat_s)	
	norm_const=norm_const+rowSums(weight_mat_s)

	if(((s/S*100)%%1)==0){
		print(paste(s/S*100,"% completed"))
	}
}

###Prediction
y_pred=y_pred/norm_const
f_pred=f_pred/t(matrix(norm_const,nrow=m,ncol=m2))


output=list(y_pred=y_pred, f_pred=f_pred)
return(output)

}
