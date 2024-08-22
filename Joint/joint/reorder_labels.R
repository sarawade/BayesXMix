#Function to reorder the labels

###INPUT
# S: number of samples
# n: sample size
# config_y, config_x: sample configurations (Sxn and Sxn)

###OUTPUT
# configy_reorder: sampled y configurations reorder so that the first subject is 
#   is in the first group
# configx_reorder: sampled x configurations reorered so that the subject with 
#   the smallest index is in the first group
# configy_count: number of times each y configuration indexed in configy_count
#   index is sampled (#Unique y configurationsx1)
# configyx_count: number of times each yx configuration indexed in configyx_count
#   index is sampled (#Unique yx configurationsx1)
# configy_index: indices of unique y configutations (#Unique y configurationsx1)
# configyx_index: indices of unique yx configurations (#Unique yx configurationsx2)


reorder_labels=function(S,n, config_y, config_x){

#Initialize output
configx_reorder=matrix(0,S,n)
configy_reorder=matrix(0,S,n)
configy_count=c(1)
configyx_count=c(1)
configy_index=c(1)
configyx_index=matrix(c(1,1),nrow=2)

for(s in 1:S){
	#reorder the configuration
	uniqy_s=unique(config_y[s,])
	for( h in 1:k_y[s]){
		configy_reorder[s,config_y[s,]==uniqy_s[h]]=h
		uniqx_s=unique(config_x[s,config_y[s,]==uniqy_s[h]])
		for(j in 1:k_x[[s]][uniqy_s[h]]){
			configx_reorder[s,config_y[s,]==uniqy_s[h]][config_x[s,config_y[s,]==uniqy_s[h]]==uniqx_s[j]]=j	
		}
	}
	if(s>1){
	#check if we have seen this configuration before
	beforey=colSums(matrix(t(configy_reorder[configy_index,]), nrow=n)==matrix(configy_reorder[s,],n,length(configy_index) ))==n
	if(sum(beforey)>0){
		configy_count[beforey]=configy_count[beforey]+1

		#check if we have also seen y config before
		configy_before_index=configy_index[beforey]
		configyx_before_index=configyx_index[2,configyx_index[1,]==configy_before_index]
		beforeyx=colSums(matrix(t(configx_reorder[configyx_before_index,]), nrow=n)==matrix(configx_reorder[s,],n,length(configyx_before_index) ))==n
		if(sum(beforeyx)>0){
			configyx_count[configyx_index[1,]==configy_before_index][(configyx_index[1,]==configy_before_index)[beforeyx]]=configyx_count[configyx_index[1,]==configy_before_index][(configyx_index[1,]==configy_before_index)[beforeyx]]+1
		}
		else{
			configyx_count=c(configyx_count,1)
			configyx_index=cbind(configyx_index, matrix(c(configy_before_index,s),nrow=2))
		}
	}
	else{
		configy_count=c(configy_count,1)
		configy_index=c(configy_index,s)
		configyx_count=c(configyx_count,1)
		configyx_index=cbind(configyx_index, matrix(c(s,s),nrow=2))
	}
	}	
}

output=list(configy_reorder=configy_reorder, configx_reorder=configx_reorder, configy_count=configy_count, configyx_count=configyx_count, configy_index=configy_index, configyx_index=configyx_index)
return(output)
}
