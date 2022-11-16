BinWeight_penalty = function(Y,alpha=sqrt(2*log(dim(Y)[2]))){#Bin-Weight method - Y is a matrix corresponding to the data: each row being the observations of a variate, each column being the observations made at a particular point in time. alpha is the `threshold penalty' above which the CUSUM statistic for a given series must exceed in order to be included in the overall test statistic.
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
    C[i,]=abs(sqrt(n/(nums*rev(nums)))*((1:(n-1))/n * S[i,n] - S[i,1:(n-1)]))/10
  }
  
  I=matrix(0,nrow=p,ncol=n-1) #updated matrix
  
  I=C*(C>alpha)
  
  ##total reduction in cost involves summing over series
  Ctot=apply(I,2,sum)
  
  #print(Ctot)

  return(list("An Empty Slot","Another Empty Slot",max(Ctot)))
}
