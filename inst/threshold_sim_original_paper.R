library(HDCD)

#multiple change-points:
#multiple change-points:
p = 2000
n = 500
mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
#etas = c(round(n/4),round(2*n/4),round(3*n/4))
#randomdense =rnorm(p)
#randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
#randomdense = randomdense/p^(1/4)/sqrt(n)*12
k = 15
#randomsparse = runif(k, min=-1, max = 1)
#randomsparse = randomsparse/norm(randomsparse, type="2")
#randomsparse = randomsparse* 3/sqrt(n)*200/k
nn = round(n/2)
mus[,round(nn+1):n] = mus[,round(nn+1):n] + matrix(rep(c(3/sqrt(n)*75/k*sample(c(-1,1), replace=TRUE,k), rep(0,p-k)), n-nn) ,nrow=p)
#mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + matrix(rep(c(randomsparse, rep(0,p-k)), n-round(etas[1])) ,nrow=p)
#mus[,round(etas[2]+1):n] = mus[,round(etas[2]+1):n] + matrix(rep(randomdense, n-round(etas[2])) ,nrow=p)/0.5
#mus[,round(etas[3]+1):n] = mus[,round(etas[3]+1):n] + matrix(rep(rep(1,p)/p^(1/4)/sqrt(n)*12, n-round(etas[3])) ,nrow=p)/0.5
plot(mus[1,])
#X = mus+noise
X = noise

start = -1
stop = n-1
CUSUM_X = CUSUM(X, start, stop)
CUSUM_X = CUSUM_X[1:(p*(stop-start-1))]
CUSUM_X = matrix(data=CUSUM_X, nrow=p,ncol = stop-start-1)
#CUSUM_X = CUSUM_X[, 1:(stop-start-1)]
plot(colSums(CUSUM_X^2))


thresholded = CUSUM_X[,]
# see what s gets highest score
scores = rep(NA, round(sqrt(p*log(log(8*n))))+1)
last = round(sqrt(p*log(log(8*n))))+1
as = scores[]
nu_as = scores[]
threshes = scores[]
maxes = scores[]
for (s in 1:last) {
  a = 0
  nu_a = 1
  if(s !=last){
    a = 2*sqrt(log(exp(1)*p*log(log(8*n))/s^2))
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
  }
  
  as[s] = a
  nu_as[s] = nu_a
  thresholded[,] = 0
  for (i in 1:p) {
    for (j in 1:(stop-start-1)) {
      if(abs(CUSUM_X[i,j])>a){
        thresholded[i,j] = CUSUM_X[i,j]^2 - nu_a
      }
    }
  }
  res = max(colSums(thresholded))
  #res = (colSums(thresholded))[round(n/2)]
  maxes[s] = res
  if(TRUE){
    thresh = 9*max(s * log(exp(1)*p*log(log(8*n))/s^2), log(log(8*n)))
  }
  threshes[s] = thresh
  scores[s] = res / thresh
  
}
plot(scores)
plot(as)
plot(maxes)
plot(nu_as)
plot(threshes)

plot(maxes)
lines(threshes/9)
abline(v = k, col=2)

