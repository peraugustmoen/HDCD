
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/madcorrection.R')

####function for minimising a cost
## of sum_{over series} RSS + penalty
## for the AMOC case
## where penalty is alpha+s beta;
## alpha,beta>0 are chosen constants and s is the number of series that change.

##to do this we first calculate the cost with penalty beta for each series and each changepoint location.

## the cost is min{RSS_0,RSS_t+beta} where RSS_0 is RSS if there are nochanges and RSS_t is RSS if there is a change at t.

## Then we sum this + alpha to get the cost for a change at t.

## sum RSS_0 over series is the cost of no change.

##Y is a p by n matrix of data (p series, n time-points)

SUBSET.normal_penalty=function(Y,alpha=2*log(dim(Y)[1])){#Y is the data matrix such that the (i,j)th entry corresponds to the jth entry for the ith variate; alpha is the penalty value incurred for adding another variate to the affected set in the sparse setting.
  n=dim(Y)[2]
  p=dim(Y)[1]

  #Y=t(apply(Y,1,madcorrection)) #uncomment this

  
  #indicators = rep(0,n-1)
  
  S=matrix(0,nrow=p,ncol=n)
  ##first get cusum stats
  for(i in 1:p) S[i,] =cumsum(Y[i,])
  
  ##we will work with C_t=RSS_0-RSS_t for each series
  ##as this can be calculated from S
  C=matrix(0,nrow=p,ncol=n-1)
  for(i in 1:p){
    C[i,]=S[i,1:(n-1)]^2/(1:(n-1))+(S[i,n]-S[i,1:(n-1)])^2/((n-1):1)-S[i,n]^2/n ##needs checking
  }
  
  ##C is reduction in RSS if we add a changepoint. Cost of adding a changepoint per series is beta so reduction in cost available is:
  D=pmax(C-alpha,0)
  
  Ctot=apply(C,2,sum)  
  Dtot=apply(D,2,sum)
 
  beta1=max(Dtot)
  if(max(Ctot)<p){
    beta2=0
  }
  else{
    beta2= (sqrt(max(max(Ctot)-p/2,0))-sqrt(p/2))^2
  }

	return(list(dense=max(Ctot),sparse=max(Dtot),beta=max(beta1,beta2),beta1=beta1,beta2=beta2))
  
}


the.function <- function(){
  start.time <- Sys.time()
  that.change <- MV.single.change(f[[1]],2*log(dim(f[[1]])[2]),log(log(dim(f[[1]])[2])))
  end.time <- Sys.time()
  return(end.time-start.time)
}

#the.function()
