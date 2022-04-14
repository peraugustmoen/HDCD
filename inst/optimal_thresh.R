library(HDCD)

#multiple change-points:
#multiple change-points:
p = 1000
n = 500
mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
etas = c(round(n/4),round(2*n/4),round(3*n/4))
randomdense =rnorm(p)
randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
randomdense = randomdense/p^(1/4)/sqrt(n)*12
k = 30
C0 = 9
#randomsparse = rnorm(k)
#randomsparse = randomsparse/norm(randomsparse, type="2")
# randomsparse = runif(k, min=-1, max = 1)
# randomsparse = randomsparse/norm(randomsparse, type="2")
# randomsparse = randomsparse* 3/sqrt(n)*1000/k
randomsparse = c(sample(c(-1, 1), replace=TRUE, k), rep(0,p-k))/k
randomsparse = randomsparse/sqrt(n)*4*max(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))
normm = sqrt(n)*norm(randomsparse, type="2")
#randomsparse = c(3/sqrt(n)*90/k*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))
mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + matrix(rep(randomsparse, n-round(etas[1])) ,nrow=p)
#mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + matrix(rep(c(randomsparse, rep(0,p-k)), n-round(etas[1])) ,nrow=p)
#mus[,round(etas[2]+1):n] = mus[,round(etas[2]+1):n] + matrix(rep(randomdense, n-round(etas[2])) ,nrow=p)/0.5
#mus[,round(etas[3]+1):n] = mus[,round(etas[3]+1):n] + matrix(rep(rep(1,p)/p^(1/4)/sqrt(n)*12, n-round(etas[3])) ,nrow=p)/0.5
plot(mus[1,])
X = mus+noise

start = -1
stop = etas[2]-1
CUSUM_X = CUSUM(X, start, stop)
CUSUM_X = CUSUM_X[1:(p*(stop-start-1))]
CUSUM_X = matrix(data=CUSUM_X, nrow=p,ncol = stop-start-1)
#CUSUM_X = CUSUM_X[, 1:(stop-start-1)]
plot(colSums(CUSUM_X^2))


#tresholded = CUSUM_X[, ]^2
thresholded = CUSUM_X[,]
# see what s gets highest score
scores = rep(NA, floor(sqrt(p*log(n^4)))+1)
scores2 = rep(NA, floor(sqrt(p*log(n^4)))+1)
last = floor(sqrt(p*log(n^4)))+1
as = scores[]
nu_as = scores[]
threshes = scores[]
maxes = scores[]
num_inc = scores[]
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
  maxes[s] = res
  thresh = 9*sqrt(p*log(n^4))
  if(s != last){
    thresh = 9*((sqrt(p*exp(-a^2/2)*log(n^4))+ log(n^4)))
  }
  threshes[s] = thresh
  scores[s] = res / thresh
  scores2[s] = res - thresh
  
}
plot(scores)
abline(v = k, col=2)




