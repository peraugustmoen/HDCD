#source('modefunction.R')

####function for minimising a cost
## of sum_{over series} RSS + penalty
## for the AMOC case
## where penalty is alpha+s beta;
## alpha,beta>0 are chosen constants and s is the number of series that change.

##to do this we first calculate the cost with penalty beta for each series and each changepoint location.

## the cost is min{RSS_0,RSS_t+alpha} where RSS_0 is RSS if there are nochanges and RSS_t is RSS if there is a change at t.

## Then we sum this + beta to get the cost for a change at t.

## sum RSS_0 over series is the cost of no change.

##Y is a p by n matrix of data (p series, n time-points)
SUBSET.negbin=function(Y,beta=4*log(dim(Y)[2]),alpha=2*log(dim(Y)[1]),estimate.r=TRUE,r=rep(10,dim(Y)[1])){#Y is the inputted data matrix, such that the (i,j)th entry corresponds to the jth observation of the ith variate, beta and alpha are as defined for the SUBSET method. estimate.r is a logical such that if TRUE the method will assume that r (the over-dispersion) is unknown and will estimate it using a method of moments approach. If it is set to FALSE, then the following vector argument (r) consistes of the (known) values of (such that the kth entry of the vector corresponds to over-dispersion in the kth variate) taken by the procedure.
  n=dim(Y)[2]
  p=dim(Y)[1]
  
  S        =matrix(0,nrow=p,ncol=n)
  U        =matrix(0,nrow=p,ncol=n)

  thresh = p+sqrt(2*p*beta)

  indicators = rep(0,n-1)
  
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
  
  Ctot=apply(C,2,sum)-thresh
  
  Dtot=apply(D,2,sum)
  
  Dvalues = Dtot - beta
  Cvalues = Ctot - beta

  savings = pmax(Dvalues,Cvalues)
  
  if(max(savings)<0) return(list(Dtot,cpt=NULL,cost=0,affected=NULL,test.stat=0)) ##no change
  else{
    if(Dtot[which.max(savings)]<max(savings)){
      return(list(savings,cpt=which.max(savings),cost=max(savings),affected=rep(1,p),test.stat=max(savings)))
    }
    else{
      switch.on <- (1:p)[D[,which.max(savings)]>0]
      affected <- rep(0,p)
      affected[switch.on]=1
      return(list(savings,cpt=which.max(savings),cost=max(savings),affected,test.stat=max(savings)))
 ##output cpt location,reduction in cost relative to no change, and which series are affected.
    }
  }
}

