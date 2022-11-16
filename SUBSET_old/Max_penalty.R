#############################
#######Max_penalty###########
#############################

Max_penalty = function(Y,beta=sqrt(2*(mad(Y[1,])^2)*log(dim(Y)[2])+2*(mad(Y[1,])^2)*log(dim(Y)[1])),alpha=sqrt(2*(mad(Y[1,])^2)*log(dim(Y)[2]))){#Y is the data matrix - entry i,j corresponds to the observation made in variate i at time j. Note that both the beta and alpha penalty values are defunct. default penalty definitely too high - otherwise code is correct; penalty needed here can be slightly larger than for the summation case
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
  
  return(list("An Empty Slot","Another Empty Slot",max(Ctot)))
}
