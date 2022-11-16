##############
#Mean_penalty#
##############

#Function finds best location of single change (or else rejects possibility of a change) using the mean of the CUSUM statistics across all variates. This is done by computing the CUSUM statistic in the case that the change occurs in a variate at each point in time.

Mean_penalty = function(Y,beta=sqrt(2*(mad(Y[1,])^2)*log(dim(Y)[2])+2*(mad(Y[1,])^2)*log(dim(Y)[1])),alpha=sqrt(2*(mad(Y[1,])^2)*log(dim(Y)[2]))){
  n=dim(Y)[2]
  p=dim(Y)[1]
  
  S=matrix(0,nrow=p,ncol=n)
  ##first get cusum stats
  for(i in 1:p) S[i,] =cumsum(Y[i,])
  
  ##we will work with C_t=RSS_0-RSS_t for each series
  ##as this can be calculated from S
  C=matrix(0,nrow=p,ncol=n-1)
  nums <- seq(0.1,(n-1)/10,by=0.1)
  for(i in 1:p){
    C[i,]=abs(sqrt(n/(nums*rev(nums)))*((1:(n-1))/n * S[i,n] - S[i,1:(n-1)]))/10 #(abs(S[i,1:(n-1)]*sqrt((n-1):1)/sqrt(1:(n-1))-(S[i,n]-S[i,1:(n-1)])*sqrt(1:(n-1))/sqrt((n-1):1)))*sqrt(n)/abs(S[i,n])
  }
  
  
  ##total reduction in cost involves summing over series, and then normalising by the number of series in question
  Ctot=apply(C,2,sum)/p
  
  #print(Ctot) - uncomment this to see CUSUM series - a significant peak should correspond well with a true change
  
  return(list("An Empty Slot","Another Empty Slot",max(Ctot)))
}
