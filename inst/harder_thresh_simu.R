# simu with harder thresh

library(HDCD)

#multiple change-points:
#multiple change-points:
p = 1000000
n = 3
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
mus = matrix(0, nrow=p, ncol=n)
set.seed(10)
etas = c(round(n/2))
randomdense =rnorm(p)
randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
randomdense = randomdense/p^(1/4)/sqrt(n)*12
k = 1
#randomsparse = rnorm(k)
#randomsparse = randomsparse/norm(randomsparse, type="2")
# randomsparse = runif(k, min=-1, max = 1)
# randomsparse = randomsparse/norm(randomsparse, type="2")
# randomsparse = randomsparse* 3/sqrt(n)*1000/k

#randomsparse = 2.5*c(2/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*rnorm(k), rep(0,p-k))
randomsparse = c(2/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))

#kk = round(p^(3/4))
kk = 300
randomsparse = c(8*rep(1,kk)/sqrt(kk)*p^(1/4), rep(0, p-kk))

#mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + matrix(rep(randomsparse, n-round(etas[1])) ,nrow=p)


plot(mus[1,])
#X = mus+noise
#X =mus
X = noise
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

#max(colSums(CUSUM_signal^2))
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
  # num_inc[s] = 0
  # as[s] = a
  # nu_as[s] = nu_a
  # thresholded[,] = 0
  # for (i in 1:p) {
  #   for (j in 1:(stop-start-1)) {
  #     if(abs(CUSUM_X[i,j])>a){
  #       thresholded[i,j] = CUSUM_X[i,j]^2 - nu_a
  #       num_inc[s] = num_inc[s]+1
  #       
  #     }
  #   }
  # }
  #res = max(colSums(thresholded))
  #maximizers[s] = which.max(colSums(thresholded))
  #maxes[s] = res
  thresh = 2*sqrt(p*log(n^4))
  thresh2 = thresh
  if(s != last){
    thresh = 9*((sqrt(p*exp(-a^2/2)*log(n^4))+ log(n^4)))
    thresh2 =20*9*((s*log(exp(1)*p*log(n^4)/s^2)+log(n^4)))
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
        thresholded[i,j] = CUSUM_X[i,j]^2 - nu_a - 2*a^2
        num_inc_tmp[j] = num_inc_tmp[j]+1
        
      }
    }
  }
  colsummm = colSums(thresholded)
  maxind = which.min(colsummm)
  res = colsummm[maxind]
  maxes2[s] = res
  scores2[s] = res / thresh2
  maximizers2[s] = maxind
  num_inc[s] = num_inc_tmp[maxind]
  as2[s] = a
  
}

plot(maxes2, ylim = c(min(-threshes2,maxes2)-1000,0))
lines(-20*threshes2)
plot(-threshes2, ylim = c(min(-threshes2), 0))
lines(maxes2)
