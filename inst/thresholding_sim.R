library(HDCD)

#multiple change-points:
#multiple change-points:
p = 1000
n = 1000
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
mus = matrix(0, nrow=p, ncol=n)
#set.seed(10)
etas = c(round(n/2))
randomdense =rnorm(p)
randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
randomdense = randomdense/p^(1/4)/sqrt(n)*12
k = 5
#randomsparse = rnorm(k)
#randomsparse = randomsparse/norm(randomsparse, type="2")
# randomsparse = runif(k, min=-1, max = 1)
# randomsparse = randomsparse/norm(randomsparse, type="2")
# randomsparse = randomsparse* 3/sqrt(n)*1000/k

#randomsparse = 2.5*c(2/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*rnorm(k), rep(0,p-k))
randomsparse = c(15/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))

#kk = round(p^(3/4))
kk = p
randomsparse = c(rep(1,kk)/sqrt(kk)*(p*log(n))^(1/4), rep(0, p-kk))/sqrt(2)

#mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + matrix(rep(randomsparse, n-round(etas[1])) ,nrow=p)

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
plot(maxes)
# par(mfrow=c(1,2))
# plot(scores)
# abline(v = k, col=2)
# plot(scores2)
# abline(v = k, col=2)
par(mfrow=c(1,2))
#plot(maxes - threshes)
#abline(v = k, col=2)

diffasq = as[1:(length(as)-1)]^2 - as[2:length(as)]^2
plot(diff(maxes2))
lines(ss*2*diffasq,col=2)
min(which(ss*2*diffasq > diff(maxes2)))
max(which(ss*2*diffasq > diff(maxes2)))


plot(maxes2-threshes2)
abline(v = k, col=2)
plot(maxes2 - threshes2/9/9)

threshes3= nu_as[]

threshes3[1:(length(ss))] = 1*(ss[1:(length(ss))]*pmax(log(exp(1)-1 + sqrt(p*log(n^4))/ss[1:(length(ss))]), log(n^4)))



threshes3[length(nu_as)] =1* (sqrt(p*log(n^4)) + log(n^4))

plot(maxes2 - 18*threshes3)
abline(v = k, col=2)
maxs = which.max(maxes2 - 18*threshes3)
max(CUSUM_signal^2)
as2[maxs]^2
as2[round(maxs/2)]^2

maxes22 = maxes2[]
maxes22[maxes22/num_inc - 2*as2^2<0] = 0

threshes22 = threshes2[]
threshes22[maxes22/num_inc - 2*as2^2<0] = 0

plot(maxes22-threshes22)
abline(v = k, col=2)

max(CUSUM_signal^2)
maxes2[1]-threshes2[1]

# plot((maxes2-threshes2)/threshes2)
# abline(v = k, col=2)


####### test med annerledes thresholding


library(HDCD)

#multiple change-points:
#multiple change-points:
p = 1000
n = 200
mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
etas = c(round(n/2))
randomdense =rnorm(p)
randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
randomdense = randomdense/p^(1/4)/sqrt(n)*12
k = 30
#randomsparse = rnorm(k)
#randomsparse = randomsparse/norm(randomsparse, type="2")
# randomsparse = runif(k, min=-1, max = 1)
# randomsparse = randomsparse/norm(randomsparse, type="2")
# randomsparse = randomsparse* 3/sqrt(n)*1000/k

randomsparse = c(4/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*rnorm(k), rep(0,p-k))
#randomsparse = c(9/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))


#randomsparse = 2*rep(1,p)/p^(1/4)
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
thresholded = CUSUM_X[,]
# see what s gets highest score
scores = rep(NA, floor(sqrt(p*log(n^4)))+1)
scores2 = rep(NA, floor(sqrt(p*log(n^4)))+1)
last = floor(sqrt(p*log(n^4)))+1
as = scores[]
as2 = scores[]
nu_as = scores[]
threshes = scores[]
threshes2 = scores[]
maxes = scores[]
maxes2 = scores[]
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
  thresholded[,] = - nu_a*exp(-a^2/2)
  for (i in 1:p) {
    for (j in 1:(stop-start-1)) {
      if(abs(CUSUM_X[i,j])>a){
        thresholded[i,j] = CUSUM_X[i,j]^2 - nu_a*2*exp(-a^2/2)
        num_inc[s] = num_inc[s]+1
        
      }
    }
  }
  res = max(colSums(thresholded))
  maxes[s] = res
  thresh = 9*sqrt(p*log(n^4))
  thresh2 = thresh
  if(s != last){
    thresh = 9*(max(sqrt(p*exp(-a^2/2)*log(n^4)), log(n^4)))
    thresh2 = 9*((s*log(exp(1)*p*log(n^4)/s^2)+ log(n^4)))
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
  num_inc[s] = 0
  #as[s] = a
  #nu_as[s] = nu_a
  thresholded[,] = - nu_a*exp(-a^2/2)
  for (i in 1:p) {
    for (j in 1:(stop-start-1)) {
      if(abs(CUSUM_X[i,j])>a){
        thresholded[i,j] = CUSUM_X[i,j]^2 - nu_a*exp(-a^2/2)
        num_inc[s] = num_inc[s]+1
        
      }
    }
  }
  res = max(colSums(thresholded))
  maxes2[s] = res
  scores2[s] = res / thresh2
  as2[s] = a
  
}
# par(mfrow=c(1,2))
# plot(scores)
# abline(v = k, col=2)
# plot(scores2)
# abline(v = k, col=2)
par(mfrow=c(1,2))
plot(maxes - threshes)
abline(v = k, col=2)

plot(maxes2-threshes2)
abline(v = k, col=2)







plot((maxes2-threshes2)*threshes2)

ss = 1:(floor(sqrt(p*log(n^4))))


maxes[k] / (k*log(exp(1)*p*log(n^4)/k^2))


plot(scores2)
abline(v = k, col=2)
ss =  1:(floor(sqrt(p*log(n^4)))+1)
plot(ss[maxes>=9 ])
plot(as)
plot(maxes)
plot(nu_as)
plot(maxes,ylim = c(min(c(threshes, maxes)), max(threshes, maxes)))
lines(9/4*threshes)
abline(v = k, col=2)
plot(maxes-threshes)

plot(maxes[2:length(maxes)] - maxes[1:(length(maxes)-1)], type="l")
lines(maxdiffs*10,type="l",col=2)

min(which(maxes[2:length(maxes)] - maxes[1:(length(maxes)-1)] <= maxdiffs*11))

diff(nu_as)
theoreticalcuttoff = sqrt(log(n^4)*p*8*exp(-as[1:(length(as)-1)]^2/2))
lines(theoreticalcuttoff, col=2)
min(which(maxes[2:length(maxes)] - maxes[1:(length(maxes)-1)] <= theoreticalcuttoff))




# original paper:
n = etas[2]
#tresholded = CUSUM_X[, ]^2
thresholded = CUSUM_X[,]
# see what s gets highest score
scores = rep(NA, floor(sqrt(p*log(log(8*n))))+1)
last = floor(sqrt(p*log(log(8*n))))+1
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
  maxes[s] = res
  #thresh = 9*sqrt(p*log(log(8*n)))
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

plot(maxes)
lines(threshes)
abline(v = k, col=2)






# experiment:
#tresholded = CUSUM_X[, ]^2
thresholded = CUSUM_X[,]
# see what s gets highest score
scores = rep(NA, round(sqrt(p*log(n^4)))+1)
last = round(sqrt(p*log(n^4)))+1
as = scores[]
nu_as = scores[]
threshes = scores[]
maxes = scores[]
for (s in 1:(round(sqrt(p*log(n^4)))+1)) {
  a = 0
  nu_a = 1
  if(s !=last){
    a = 2*sqrt(log(exp(1)*p*log(n^4)/s^2))
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
  }
  
  as[s] = a
  nu_as[s] = nu_a
  tresholded[,] = 0
  for (i in 1:p) {
    for (j in 1:(stop-start-1)) {
      if(abs(CUSUM_X[i,j])>a){
        thresholded[i,j] = CUSUM_X[i,j]^2 - nu_a
      }
    }
  }
  res = max(colSums(thresholded))
  maxes[s] = res
  thresh = 4*sqrt(p*log(n^4))
  if(s != last){
    #thresh = max(s * log(exp(1)*p*log(n^4)/s^2), log(n^4))
    thresh = max(s * log(exp(1)*p*log(n^4)/s^2), log(n^4))
  }
  threshes[s] = thresh
  
  scores[s] = res / thresh
  if(scores[s]<1){
    scores[s] = 0
  }
  scores[s] =scores[s]/sqrt(s)
  
}
plot(scores)
plot(as)
plot(maxes)
plot(nu_as)

