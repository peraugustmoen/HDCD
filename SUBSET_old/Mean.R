########################
#MV.single.change.CUSUM#
########################

#Function finds best location of single change (or else rejects possibility of a change) using the mean of the CUSUM statistics across all variates. This is done by computing the CUSUM statistic in the case that the change occurs in a variate at each point in time.

Mean = function(Y,beta=sqrt(2*(mad(Y[1,])^2)*log(dim(Y)[2])+2*(mad(Y[1,])^2)*log(dim(Y)[1])),alpha=sqrt(2*(mad(Y[1,])^2)*log(dim(Y)[2]))){#default penalty definitely too high - otherwise code is correct. Y should be entered as a matrix - each column the observations of a time point, and each row the observations of a variate. beta is the penalty, and alpha is a defunct second penalty argument.
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
  
  if(max(Ctot)<beta) return(list(Ctot,cpt=NULL,cost=max(Ctot),affected=NULL,test.stat=max(Ctot)-beta)) ##no change
  if(max(Ctot)>beta) return(list(Ctot,cpt=which.max(Ctot),cost=max(Ctot),affected=rep(1,p),test.stat=max(Ctot)-beta)) ##output cpt location,reduction in cost relative to no change, and which series are affected.
}
