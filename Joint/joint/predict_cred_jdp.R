# Computes prediction for a new subject with credible intervals
# for DP regression model

###INPUT
# a: Credible level for credible intervals of mean
# S: number of samples
# p: number of covariates
# x_new: covariate of new subjects (1xp)
# mu_theta, C: prior parameters for beta (p+1x1 and p+1xp+1)
# a_y, b_y: prior parameters for sigma^2_y (1x1 and 1x1)
# mu_0, c_x: prior parameters for mu  (px1 and px1)
# a_x, b_x: prior parameters for sigma^2_x (px1 and px1)
# config: sampled configurations (Sxn)
# phi_y: sampled y cluster parameters (list with S elements, each element is p+1xk[s])
# mu_x, sigma_x: sampled x cluster parameters (list with S elements, each element is pxk[s] and pxk[s])
# alpha: sampled precision parameters (Sx1)

###OUTPUT
# y_pred: prediction of response for new subject (m2x1)
# l_pred: lower bound for pointwise credible interval of prediction
# u_pred: upper bound for pointwise credible interval of prediction


predict_cred_jdp=function(a, S,m,p,x_new, mu_theta, C, a_y, b_y, mu_0, c_x, a_x, b_x, k, config, alpha, phi_y, sigma_y, mu_x, sigma_x ){

#Create X matrix
X_new=matrix(c(rep(1,m),x_new),nrow=m,,ncol=(p+1))

#marginal for x
c_x_hat=c_x+1
x_sc=(x_new-t(matrix(mu_0,nrow=p,ncol=m)))*t(matrix((c_x*a_x/(c_x_hat*b_x))^.5,nrow=p,ncol=m))
marg_x=exp(rowSums(dt(x_sc,2*t(matrix(a_x,nrow=p,ncol=m)),log=TRUE)+t(matrix(.5*(log(c_x)-log(c_x_hat)-log(b_x)+log(a_x)),nrow=p,ncol=m))))


#Initialize
#first element is for new subject
weight_mat=matrix(0,m, sum(k)+1)
weight_mat[,1]=sum(alpha/(alpha+n))*marg_x
y_pred_mat=matrix(0,m, sum(k)+1)
y_pred_mat[,1]=X_new%*%mu_theta

for(s in 1:S){

	#calculate weight matrix
	weight_mat_s=matrix(0,m,k[s])

	for(j in 1:k[s]){
		nj_s=sum(config[s,]==j)
		weight_mat_s[,j]=nj_s/(alpha[s]+n)*exp(rowSums(dnorm(x_new, t(matrix(mu_x[[s]][,j], nrow=p,ncol=m)), t(matrix(sigma_x[[s]][,j]^.5,nrow=p,ncol=m)),log=T)))
	}

	#calculate prediction
	y_pred_s=X_new%*%phi_y[[s]]
	st=2
	if(s!=1){
		st=st+sum(k[1:(s-1)])
	}
	y_pred_mat[,st:(sum(k[1:s])+1)]=y_pred_s
	weight_mat[,st:(sum(k[1:s])+1)]=weight_mat_s
	
	if(((s/S*100)%%1)==0){
		print(paste(s/S*100,"% completed"))
	}
}

###Prediction
y_pred=rowSums(y_pred_mat*weight_mat)/rowSums(weight_mat)

#calculate upper and lower 95% credible bounds for y_pred
u_pred=matrix(0,m,1)
l_pred=matrix(0,m,1)
for(i in 1:m){
	pred_i_sort=sort(y_pred_mat[i,], index.return=T)
	cweight_i_sort=cumsum(weight_mat[i,pred_i_sort$ix])/sum(weight_mat[i,])
	ind_l=sum(cweight_i_sort<a/2)
	l_pred[i]=(pred_i_sort$x[ind_l]*(a/2-cweight_i_sort[ind_l])+pred_i_sort$x[ind_l+1]*(cweight_i_sort[ind_l+1]-a/2))/(cweight_i_sort[ind_l+1]-cweight_i_sort[ind_l])
	ind_u=sum(cweight_i_sort<(1-a/2))
	u_pred[i]=(pred_i_sort$x[ind_u]*(1-a/2-cweight_i_sort[ind_u])+pred_i_sort$x[ind_u+1]*(cweight_i_sort[ind_u+1]-1+a/2))/(cweight_i_sort[ind_u+1]-cweight_i_sort[ind_u])
}


output=list(y_pred=y_pred,l_pred=l_pred, u_pred=u_pred)
return(output)
}
