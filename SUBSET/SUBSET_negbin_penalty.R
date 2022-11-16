SUBSET.negbin_penalty=function(Y,alpha=2*log(dim(Y)[1]),estimate.r=TRUE,r=rep(10,dim(Y)[1])){#Y is the data matrix, such that the (i,j)th entry corresponds to the jth observation of the ith variate. alpha is the penalty incurred for adding an additional variate to the affected set of variates in the sparse setting. estimate.r is a logical, such that if TRUE the method will seek to estimate the over-dispersion for each variate. If FALSE, the method will take the following argument (r) as the known values of the over-dispersions. r is a vector such that the kth entry corresponds to the (known) over-dispersion in the kth variate.
  n=dim(Y)[2]
  p=dim(Y)[1]
  
  S        =matrix(0,nrow=p,ncol=n)
  U        =matrix(0,nrow=p,ncol=n)

  #thresh = 0

  #indicators = rep(0,n-1)
  
  ##first get cusum stats
  for(i in 1:p){
    S[i,]    =cumsum(Y[i,])
    if(estimate.r){
      U[i,]    =cumsum(Y[i,]^2)
    }
  }
  
  if(estimate.r){
    for(i in 1:p){
      r[i]=abs(((S[i,n]/n)^2)/((U[i,n]-(S[i,n])^2/n)/(n-1)-S[i,n]/n +0.00001) + 0.00001)
    }
  }
  
  ##we will work with C_t=RSS_0-RSS_t for each series
  ##as this can be calculated from S
  C=matrix(0,nrow=p,ncol=n-1)
  for(i in 1:p){
    one  =2*r[i]*((1:(n-1))*log(1:(n-1))+((n-1):1)*log((n-1):1)-n*log(n))
    two  =2*r[i]*(n*log(r[i]*n+S[i,n])-(1:(n-1))*log(r[i]*(1:(n-1))+S[i,1:(n-1)])-((n-1):1)*log(((n-1):1)*r[i]+S[i,n]-S[i,1:(n-1)]))
    three=2*(S[i,1:(n-1)]*log(S[i,1:(n-1)]+0.000001)+(S[i,n]-S[i,1:(n-1)])*log(S[i,n]-S[i,1:(n-1)]+0.000001)-S[i,n]*log(S[i,n]+0.000001))
    four =2*(S[i,n]*log(r[i]*n+S[i,n])-S[i,1:(n-1)]*log(r[i]*(1:(n-1))+S[i,1:(n-1)])-(S[i,n]-S[i,1:(n-1)])*log(r[i]*((n-1):1)+(S[i,n]-S[i,1:(n-1)]))) 
    C[i,]=one+two+three+four
  }
  
  ##C is reduction in RSS if we add a changepoint. Cost of adding a changepoint per series is alpha so reduction in cost available is:
  D=pmax(C-alpha,0)
  
  Ctot=apply(C,2,sum)
  
  Dtot=apply(D,2,sum)

  beta1=max(Dtot)
  beta2= (sqrt(max(max(Ctot)-p/2,0)) -sqrt(p/2) )^2

	return(list(dense=max(Ctot),sparse=max(Dtot),beta=max(beta1,beta2),beta1=beta1,beta2=beta2))

}

