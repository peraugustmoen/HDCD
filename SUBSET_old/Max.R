#####################
#######Max###########
#####################

#Function finds best location of single change (or else rejects 
#possibility of a change) using the maximum of the CUSUM statistics 
#across all variates. This is done by computing the CUSUM 
#statistic in the case that the change occurs in a variate at 
#each point in time.

Max = function(Y,beta=sqrt(2*(mad(Y[1,])^2)*log(dim(Y)[2])+2*(mad(Y[1,])^2)*log(dim(Y)[1])),alpha=sqrt(2*(mad(Y[1,])^2)*log(dim(Y)[2]))){#Y is the data matrix (each column the set of observations at a given time point, each row the set of observations across a variate). beta and alpha are the penalty values, although note that the latter is defunct. default penalty definitely too high - otherwise code is correct; penalty needed here can be slightly larger than for the summation case. 
  n=dim(Y)[2]
  p=dim(Y)[1]
  
  #print(n)
  
  S=matrix(0,nrow=p,ncol=n)
  ##first get cusum stats
  for(i in 1:p) S[i,] =cumsum(Y[i,])
  
  ##we will work with C_t=RSS_0-RSS_t for each series
  ##as this can be calculated from S
  C=matrix(0,nrow=p,ncol=n-1)
  nums <- seq(0.1,(n-1)/10,by=0.1)
  for(i in 1:p){
    C[i,]=(abs(sqrt(n/(nums*rev(nums)))*((1:(n-1))/n * S[i,n] - S[i,1:(n-1)])))/10
  }
  
  
  ##we find the maximum across all variates for each point in time
  Ctot=apply(C,2,max)
  
  if(max(Ctot)<beta) return(list(Ctot,cpt=NULL,cost=max(Ctot),affected=NULL,test.stat=max(Ctot)-beta)) ##no change
  if(max(Ctot)>beta) return(list(Ctot,cpt=which.max(Ctot),cost=max(Ctot),affected=rep(1,p),test.stat=max(Ctot)-beta)) ##output cpt location,reduction in cost relative to no change, and which series are affected.
}
