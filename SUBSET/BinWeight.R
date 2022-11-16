#################################
##########BinWeight##############
#################################

#BinWeight computes a binary-weighted sum of the CUSUM
#statistics in order to find the best location (including a null fit)
#of a changepoint. CUSUM statistics are weighted with a 1 if and only if
#they exceed a certain threshold, whose square is beta.

BinWeight = function(Y,beta=(dim(Y)[1])*sqrt(3*((mad(Y[1,]))^2)*log(dim(Y)[1])+((mad(Y[1,]))^2)*log(dim(Y)[2])),alpha=sqrt(2*log(dim(Y)[2]))){#Y is the data matrix - entry i,j corresponds to the jth observation made in the ith variate. beta is the penalty for flagging a change, while alpha is the threshold below which a series does not contribute to the overall test statistic. Default penalty definitely too high - otherwise code is correct; penalty needed here can be slightly larger than for the summation case
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
    C[i,]=abs(sqrt(n/(nums*rev(nums)))*((1:(n-1))/n * S[i,n] - S[i,1:(n-1)]))/10##needs checking
  }
  
  I=matrix(0,nrow=p,ncol=n-1) #updated matrix
  
  I=C*(C>alpha)
  
  ##total reduction in cost involves summing over series
  Ctot=apply(I,2,sum)
  
  #print(Ctot)
  
  if(max(Ctot)<beta) return(list(Ctot,cpt=NULL,cost=max(Ctot),affected=NULL,test.stat=max(Ctot)-beta)) ##no change
  if(max(Ctot)>beta) return(list(Ctot,cpt=which.max(Ctot),cost=max(Ctot),affected=rep(1,p),test.stat=max(Ctot)-beta)) ##output cpt location,reduction in cost relative to no change, and which series are affected.
}
