library(HDCD)

#multiple change-points:
#multiple change-points:
p = 1000
n = 2
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
mus = matrix(0, nrow=p, ncol=n)
#set.seed(10)
etas = c(round(n/2))
randomdense =rnorm(p)
randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
randomdense = randomdense/p^(1/4)/sqrt(n)*12
k = 20
#randomsparse = rnorm(k)
#randomsparse = randomsparse/norm(randomsparse, type="2")
# randomsparse = runif(k, min=-1, max = 1)
# randomsparse = randomsparse/norm(randomsparse, type="2")
# randomsparse = randomsparse* 3/sqrt(n)*1000/k

#randomsparse = 2.5*c(2/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*rnorm(k), rep(0,p-k))
randomsparse = c(8/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))

#kk = round(p^(3/4))
#kk = 30
kk = p
#randomsparse = c(8*rep(1,kk)/sqrt(kk)*p^(1/4), rep(0, p-kk))

mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + matrix(rep(randomsparse, n-round(etas[1])) ,nrow=p)

#mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + matrix(rep(c(randomsparse, rep(0,p-k)), n-round(etas[1])) ,nrow=p)
#mus[,round(etas[2]+1):n] = mus[,round(etas[2]+1):n] + matrix(rep(randomdense, n-round(etas[2])) ,nrow=p)/0.5
#mus[,round(etas[3]+1):n] = mus[,round(etas[3]+1):n] + matrix(rep(rep(1,p)/p^(1/4)/sqrt(n)*12, n-round(etas[3])) ,nrow=p)/0.5
plot(mus[1,])
X = mus+noise
#X =mus
start = -1
stop = n-1
Xt = InspectChangepoint::cusum.transform(X)
CUSUM_X = CUSUM(X, start, stop)
CUSUM_X = CUSUM_X[1:(p*(stop-start-1))]
CUSUM_X = matrix(data=CUSUM_X, nrow=p,ncol = stop-start-1)
#CUSUM_X = CUSUM_X[, 1:(stop-start-1)]
plot(colSums(CUSUM_X^2))

CUSUM_signal = CUSUM(mus, start, stop)
CUSUM_signal = CUSUM_signal[1:(p*(stop-start-1))]
CUSUM_signal = matrix(data=CUSUM_signal, nrow=p,ncol = stop-start-1)

max(colSums(CUSUM_signal^2))
#tresholded = CUSUM_X[, ]^2
CUSUM_X = as.matrix(CUSUM_X, ncol=n-1, nrow= p)
thresholded = matrix(NA, nrow = p, ncol = n-1)
# see what s gets highest score
scores = rep(NA, floor(sqrt(p*log(n^4)))+1)
scores2 = rep(NA, floor(sqrt(p*log(n^4)))+1)
last = floor(sqrt(p*log(n^4)))+1
as = scores[]
as2 = scores2[]
nu_as = scores[]
threshes = scores[]
threshes2 = scores[]
maxes = scores[]
maxes2 = scores[]
num_inc = scores[]
maximizers = scores[]
maximizers2 = scores[]
num_inc_tmp = rep(NA, (stop-start-1))
ss = 1:(floor(sqrt(p*log(n^4))))
for (s in 1:(floor(sqrt(p*log(n^4)))+1)) {
  a = 0
  nu_a = 1
  if(s !=last){
    a = sqrt(2*log(exp(1)*p*log(n^4)/s^2))
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
  }
  num_inc[s] = 0
  as[s] = a
  nu_as[s] = nu_a
  thresholded[,] = 0
  for (i in 1:p) {
    for (j in 1:(stop-start-1)) {
      if(abs(CUSUM_X[i,j])>a){
        thresholded[i,j] = CUSUM_X[i,j]^2 - nu_a
        num_inc[s] = num_inc[s]+1
        
      }
    }
  }
  res = max(colSums(thresholded))
  maximizers[s] = which.max(colSums(thresholded))
  maxes[s] = res
  thresh = 9*sqrt(p*log(n^4))
  thresh2 = thresh
  if(s != last){
    thresh = 9*((sqrt(p*exp(-a^2/2)*log(n^4))+ log(n^4)))
    thresh2 =9*((s*log(exp(1)*p*log(n^4)/s^2)+log(n^4)))
  }
  threshes[s] = thresh
  threshes2[s] = thresh2
  scores[s] = res / thresh
  
  a = 0
  nu_a = 1
  if(s !=last){
    a = sqrt(4*log(exp(1)*p*log(n^4)/s^2))
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
  }
  #num_inc[s] = 0
  num_inc_tmp[]=0
  as2[s] = a
  #as[s] = a
  #nu_as[s] = nu_a
  thresholded[,] = 0
  for (i in 1:p) {
    for (j in 1:(stop-start-1)) {
      if(abs(CUSUM_X[i,j])>a){
        thresholded[i,j] = CUSUM_X[i,j]^2 - nu_a
        num_inc_tmp[j] = num_inc_tmp[j]+1
        
      }
    }
  }
  colsummm = colSums(thresholded)
  maxind = which.max(colsummm)
  res = colsummm[maxind]
  maxes2[s] = res
  scores2[s] = res / thresh2
  maximizers2[s] = maxind
  num_inc[s] = num_inc_tmp[maxind]
  as2[s] = a
  
}
# par(mfrow=c(1,2))
# plot(scores)
# abline(v = k, col=2)
# plot(scores2)
# abline(v = k, col=2)
par(mfrow=c(1,2))
#plot(maxes - threshes)
#abline(v = k, col=2)

diffasq = as[1:(length(as)-1)]^2 - as[2:length(as)]^2
plot(diff(maxes2),ylim = c(min(diff(maxes2), ss*2*diffasq),max(diff(maxes2), ss*2*diffasq)))
lines(ss*2*diffasq,col=2)
max(which(ss*2*diffasq <=diff(maxes2)))
max(which(ss*2*diffasq > diff(maxes2)))





colmaxes = function(x){
  p = dim(x)[1]
  n = dim(x)[2]
  
  ret = x[1,]
  
  for (j in 1:n) {
    ret[j] = max(x[,j])
  }
  return(ret)
}

colmins = function(x){
  p = dim(x)[1]
  n = dim(x)[2]
  
  ret = x[1,]
  
  for (j in 1:n) {
    ret[j] = min(x[,j])
  }
  return(ret)
}
n = 10000
p = 50
N=n
maxs = floor(sqrt(p*log(n^(1))))
ss = 1:maxs
sums = matrix(NA, nrow=N, ncol=maxs-1)
as = sums[1,]
nu_as = sums[1,]
delta_as = sums[1,]
ss = 1:maxs
for (t in 1:N) {
  y = rnorm(p)
  for (i in 1:((maxs-1))) {
    k = 2*i
    mm = c(30*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))
    mm = rep(0,p)
    x = y[] +  mm
    a = sqrt(4*log(exp(1)*p*log(n^(1))/i^2))
    as[i] = a
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
    nu_as[i] = nu_a
    #a_tilde = sqrt(2*log(p/i))
    for (j in 1:p) {
      # if(abs(x[j])>a_tilde){
      #   counter = counter+1
      # }
      if(abs(x[j])>a){
        
        #counter = counter+1
        #counter2 = counter2-1
        #if(abs(x[j])>delta_a){
        #  counter2 = counter2+2
        #}
        x[j] = x[j]^2 - nu_a 
        # if(x[j]>0){
        #   counter = counter+1
        # }
        # else{
        #   counter = counter-1
        # }
        
      }
      else{
        x[j]=0
      }
    }
    sums[t,i] = - sum(x) 
    a1 = a
    x = y[] +  mm
    a = sqrt(4*log(exp(1)*p*log(n^(1))/(i+1)^2))
    #as[i] = a
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
    #nu_as[i] = nu_a
    #a_tilde = sqrt(2*log(p/i))
    for (j in 1:p) {
      # if(abs(x[j])>a_tilde){
      #   counter = counter+1
      # }
      if(abs(x[j])>a){
        
        counter = counter+1
        counter2 = counter2-1
        if(abs(x[j])>delta_a){
          counter2 = counter2+2
        }
        x[j] = x[j]^2 - nu_a 
        # if(x[j]>0){
        #   counter = counter+1
        # }
        # else{
        #   counter = counter-1
        # }
        
      }
      else{
        x[j]=0
      }
    }
    #sums[t,i] = sums[t,i]+ sum(x) - 2*k* (a1^2 - a^2)
    sums[t,i] = sums[t,i]+ sum(x) 
    

    
  }
}
plot(colmaxes(sums) )
plot(colmins(sums))
plot(colMeans(sums))

diffasq = as[1:(length(as)-1)]^2 - as[2:length(as)]^2
ss[1:(length(ss)-2)]*diffasq
plot(colmaxes(sums) - ss[1:(length(ss)-2)]*diffasq)
plot(colmins(sums) -  ss[1:(length(ss)-2)]*diffasq)
plot(colMeans(sums))





