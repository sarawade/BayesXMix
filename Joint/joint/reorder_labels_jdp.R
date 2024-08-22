#Function to reorder the labels

###INPUT
# S: number of samples
# n: sample size
# config: sample configurations (Sxn and Sxn)

###OUTPUT
# config_reorder: sampled configurations reorder so that the first subject is 
#   is in the first group (Sxn)
# config_count: number of times each configuration indexed in config_index
#   is sampled (#Unique configurationsx1)
# configy_index: indices of unique configutations (#Unique configurationsx1)

reorder_labels=function(S,n, config){

#Initialize output
config_reorder=matrix(0,S,n)
config_count=c(1)
config_index=c(1)

for(s in 1:S){
	#reorder the configuration
	uniq_s=unique(config[s,])
	for( h in 1:k[s]){
		config_reorder[s,config[s,]==uniq_s[h]]=h
	}
	if(s>1){
	#check if we have seen this configuration before
	before=colSums(matrix(t(config_reorder[config_index,]), nrow=n)==matrix(config_reorder[s,],n,length(config_index) ))==n
	if(sum(before)>0){
		config_count[before]=config_count[before]+1
	}
	else{
		config_count=c(config_count,1)
		config_index=c(config_index,s)
	}
	}	
}


output=list(config_reorder=config_reorder, config_count=config_count, config_index=config_index)
return(output)
}
