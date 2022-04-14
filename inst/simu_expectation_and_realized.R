erfc = function(z){
  return(2*pnorm(z*sqrt(2), lower.tail=FALSE))
}

mu = function(c, s,aa){
  nu_aa  = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))
  EE = 1/2*(c -nu_aa +1)*erfc((aa-sqrt(c))/sqrt(2)) + exp(-(aa-sqrt(c))^2/2)*(aa+sqrt(c))/sqrt(2*pi)
  EE = EE + 1/2  *(c-nu_aa + 1)*erfc((aa+sqrt(c))/sqrt(2)) + exp(-(aa+sqrt(c))^2/2)*(aa-sqrt(c))/sqrt(2*pi)
  return(s*EE)
}
### how does threshed look?
library(HDCD)

#multiple change-points:
#multiple change-points:
p = 1000
n = 100
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
if(k<sqrt(p*log(n^4))){
  randomsparse = c(3/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))
}else{
  randomsparse = c(1.25*rep(1,k)/sqrt(k)*p^(1/4), rep(0, p-k))
}
#kk = round(p^(3/4))
#kk = p
#randomsparse = c(100*rep(1,kk)/sqrt(kk)*p^(1/4), rep(0, p-kk))

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

# max(colSums(CUSUM_signal^2))
# #tresholded = CUSUM_X[, ]^2
# thresholded = CUSUM_X[,]
# # see what s gets highest score
# scores = rep(NA, floor(sqrt(p*log(n^4)))+1)
# scores2 = rep(NA, floor(sqrt(p*log(n^4)))+1)
# last = floor(sqrt(p*log(n^4)))+1
# as = scores[]
# as2 = scores2[]
# nu_as = scores[]
# threshes = scores[]
# threshes2 = scores[]
# maxes = scores[]
# maxes2 = scores[]
# num_inc = scores[]
# maximizers = scores[]
# maximizers2 = scores[]
#for (s in 1:(floor(sqrt(p*log(n^4)))+1)) {
par(mfrow=c(2,2))
#for (s in c(1,max(2,min(round(k/2),k^(1/4))), k, p)) {
thresholded = CUSUM_X[,]
for (s in c(1,max(2,min(round(k/2))), k, p)) {
  
  
  #s = 1
  a = 0
  nu_a = 1
  if(s < sqrt(p*log(n^4))){
    a = sqrt(4*log(exp(1)*p*log(n^4)/s^2))
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
  }
  # num_inc[s] = 0
  # as2[s] = a
  #as[s] = a
  #nu_as[s] = nu_a
  #thresholded[,] = 0
  res = 1:(stop-start-1)
  
  for (j in 1:(stop-start-1)) {
    res[j] = mu(CUSUM_signal[1,j]^2, k, a)
  }
  thresh = 9*sqrt(p*log(n^4))
  if(s < sqrt(p*log(n^4))){
    thresh =9*(max(s*log(exp(1)*p*log(n^4)/s^2),log(n^4)))
  }
  
  # res = max(colSums(thresholded))
  # maxes2[s] = res
  # scores2[s] = res / thresh2
  # maximizers2[s] = which.max(colSums(thresholded))
  # as2[s] = a
  
  
  thresholded[,] = 0
  for (i in 1:p) {
    for (j in 1:(stop-start-1)) {
      if(abs(CUSUM_X[i,j])>a){
        thresholded[i,j] = CUSUM_X[i,j]^2 - nu_a
        #num_inc[s] = num_inc[s]+1

      }
    }
  }
  res_emp = colSums(thresholded)
  plot(res, xlab = sprintf("s = %d", s),ylim = c(0, max(c(res_emp, res,thresh))))
  abline(h=thresh,col=2)
  abline(h = max(res) - 3*thresh/9,col=3)
  print(sprintf("s = %d, thresh = %f, max = %f, diff = %f", s, thresh, max(res),max(res)-thresh))
  #points(res_emp,col=4, pch=2)
  lines(res_emp,col=4, lwd=2)
  #lines(res + 2*sqrt(res*thresh*9/9) + thresh/9*9, col=2, lty=2)
  #lines(res -2*sqrt(res*thresh*9/9) -thresh*9/9, col=2, lty=2)
}




# can we ever have  C^b/s <= a^2(t) and mu_t(C^b) - C1 * xi(t) > mu_s(C^b) - C1 xi(s)?
#p = 1000
logp = 50
logn4= 10
s = 10000
tos= 2*s
C1 =1
C2 = 36
t = 1
at = sqrt(4*(logp + 1 +log(logn4) -log(t^2)))
as =0
xi_t = C1*(t*at^2 +logn4)
#xi_s = C1*sqrt(log(n^4))*exp(logp/2)
if(s <sqrt(exp(logp +logn4))){
  as = sqrt(4*(logp + 1 +log(logn4) -log(s^2)))
  xi_s = C1*min(s*as^2 + logn4, sqrt(p*exp(logp +logn4)))
}else{print("FEIL")}
atos = 0
xi_tos = sqrt(exp(logp + logn4))
if(2*s < sqrt(exp(logp)*log(n^4))){
  atos = sqrt(4*(logp + 1 +log(logn4) -log((2*s)^2)))
  xi_tos = C1*(s*2*atos^2 + logn4)
}

mut = mu(at^2, s, at)
mus = mu(at^2, s, as)
mutos = mu(at^2, 2*s, atos)

cbind(mut - xi_t, mus - xi_s)
(mut - xi_t< mus - xi_s)| (mut < C2*xi_t)
(mut - xi_t ) / ( mus - xi_s)*100

cbind(mutos - xi_tos, mus - xi_s)







# what if we subtract eg 2*nu_a?

erfc = function(z){
  return(2*pnorm(z*sqrt(2), lower.tail=FALSE))
}

mu = function(c, s,aa){
  nu_aa  = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))
  EE = 1/2*(c -nu_aa +1)*erfc((aa-sqrt(c))/sqrt(2)) + exp(-(aa-sqrt(c))^2/2)*(aa+sqrt(c))/sqrt(2*pi)
  EE = EE + 1/2  *(c-nu_aa + 1)*erfc((aa+sqrt(c))/sqrt(2)) + exp(-(aa+sqrt(c))^2/2)*(aa-sqrt(c))/sqrt(2*pi)
  return(s*EE)
}
### how does threshed look?
library(HDCD)

#multiple change-points:
#multiple change-points:
p = 1000
n = 100
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
if(k<sqrt(p*log(n^4))){
  randomsparse = c(5/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))
}else{
  randomsparse = c(12*rep(1,k)/sqrt(k)*p^(1/4), rep(0, p-k))
}
#kk = round(p^(3/4))
#kk = p
#randomsparse = c(100*rep(1,kk)/sqrt(kk)*p^(1/4), rep(0, p-kk))

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

# max(colSums(CUSUM_signal^2))
# #tresholded = CUSUM_X[, ]^2
# thresholded = CUSUM_X[,]
# # see what s gets highest score
# scores = rep(NA, floor(sqrt(p*log(n^4)))+1)
# scores2 = rep(NA, floor(sqrt(p*log(n^4)))+1)
# last = floor(sqrt(p*log(n^4)))+1
# as = scores[]
# as2 = scores2[]
# nu_as = scores[]
# threshes = scores[]
# threshes2 = scores[]
# maxes = scores[]
# maxes2 = scores[]
# num_inc = scores[]
# maximizers = scores[]
# maximizers2 = scores[]
#for (s in 1:(floor(sqrt(p*log(n^4)))+1)) {
par(mfrow=c(2,3))
#for (s in c(1,max(2,min(round(k/2),k^(1/4))), k, p)) {
thresholded = CUSUM_X[,]
for (s in c(1,max(2,min(round(k/2))), k, 2*k,3*k, p)) {
  
  
  #s = 1
  a = 0
  nu_a = 1
  if(s < sqrt(p*log(n^4))){
    a = sqrt(4*log(exp(1)*p*log(n^4)/s^2))
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
  }
  # num_inc[s] = 0
  # as2[s] = a
  #as[s] = a
  #nu_as[s] = nu_a
  #thresholded[,] = 0
  res = 1:(stop-start-1)
  
  for (j in 1:(stop-start-1)) {
    res[j] = mu(CUSUM_signal[1,j]^2, k, a)
  }
  thresh = 9*sqrt(p*log(n^4))
  if(s < sqrt(p*log(n^4))){
    thresh =9*(max(s*log(exp(1)*p*log(n^4)/s^2),log(n^4)))
  }
  
  # res = max(colSums(thresholded))
  # maxes2[s] = res
  # scores2[s] = res / thresh2
  # maximizers2[s] = which.max(colSums(thresholded))
  # as2[s] = a
  
  
  thresholded[,] = 0
  for (i in 1:p) {
    for (j in 1:(stop-start-1)) {
      if(abs(CUSUM_X[i,j])>a){
        thresholded[i,j] = CUSUM_X[i,j]^2 - nu_a - a^2
        #num_inc[s] = num_inc[s]+1
        
      }
    }
  }
  res_emp = colSums(thresholded)
  plot(res_emp, xlab = sprintf("s = %d", s),ylim = c(0, max(c(res_emp, thresh))))
  abline(h=thresh,col=2)
  abline(h = max(res_emp) - 3*thresh/9,col=3)
  print(sprintf("s = %d, thresh = %f, max = %f, diff = %f", s, thresh, max(res),max(res)-thresh))
  #points(res_emp,col=4, pch=2)
  lines(res_emp,col=4, lwd=2)
  #lines(res + 2*sqrt(res*thresh*9/9) + thresh/9*9, col=2, lty=2)
  #lines(res -2*sqrt(res*thresh*9/9) -thresh*9/9, col=2, lty=2)
}











# how does it look with only noise? not too small???????

### how does threshed look?
library(HDCD)

#multiple change-points:
#multiple change-points:
p = 1000
n = 100
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
if(k<sqrt(p*log(n^4))){
  randomsparse = c(5/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))
}else{
  randomsparse = c(12*rep(1,k)/sqrt(k)*p^(1/4), rep(0, p-k))
}
#kk = round(p^(3/4))
#kk = p
#randomsparse = c(100*rep(1,kk)/sqrt(kk)*p^(1/4), rep(0, p-kk))

mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + matrix(rep(randomsparse, n-round(etas[1])) ,nrow=p)

#mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + matrix(rep(c(randomsparse, rep(0,p-k)), n-round(etas[1])) ,nrow=p)
#mus[,round(etas[2]+1):n] = mus[,round(etas[2]+1):n] + matrix(rep(randomdense, n-round(etas[2])) ,nrow=p)/0.5
#mus[,round(etas[3]+1):n] = mus[,round(etas[3]+1):n] + matrix(rep(rep(1,p)/p^(1/4)/sqrt(n)*12, n-round(etas[3])) ,nrow=p)/0.5
plot(mus[1,])
#X = mus+noise
X =noise
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

# max(colSums(CUSUM_signal^2))
# #tresholded = CUSUM_X[, ]^2
# thresholded = CUSUM_X[,]
# # see what s gets highest score
# scores = rep(NA, floor(sqrt(p*log(n^4)))+1)
# scores2 = rep(NA, floor(sqrt(p*log(n^4)))+1)
# last = floor(sqrt(p*log(n^4)))+1
# as = scores[]
# as2 = scores2[]
# nu_as = scores[]
# threshes = scores[]
# threshes2 = scores[]
# maxes = scores[]
# maxes2 = scores[]
# num_inc = scores[]
# maximizers = scores[]
# maximizers2 = scores[]
#for (s in 1:(floor(sqrt(p*log(n^4)))+1)) {
par(mfrow=c(2,3))
#for (s in c(1,max(2,min(round(k/2),k^(1/4))), k, p)) {
thresholded = CUSUM_X[,]
for (s in c(1,max(2,min(round(k/2))), k, 2*k,3*k, p)) {
  
  
  #s = 1
  a = 0
  nu_a = 1
  if(s < sqrt(p*log(n^4))){
    a = sqrt(4*log(exp(1)*p*log(n^4)/s^2))
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
  }
  # num_inc[s] = 0
  # as2[s] = a
  #as[s] = a
  #nu_as[s] = nu_a
  #thresholded[,] = 0
  res = 1:(stop-start-1)
  
  for (j in 1:(stop-start-1)) {
    res[j] = mu(CUSUM_signal[1,j]^2, k, a)
  }
  thresh = 9*sqrt(p*log(n^4))
  if(s < sqrt(p*log(n^4))){
    thresh =9*(max(s*log(exp(1)*p*log(n^4)/s^2),log(n^4)))
  }
  
  # res = max(colSums(thresholded))
  # maxes2[s] = res
  # scores2[s] = res / thresh2
  # maximizers2[s] = which.max(colSums(thresholded))
  # as2[s] = a
  
  
  thresholded[,] = 0
  for (i in 1:p) {
    for (j in 1:(stop-start-1)) {
      if(abs(CUSUM_X[i,j])>a){
        thresholded[i,j] = CUSUM_X[i,j]^2 - nu_a - a^2
        #num_inc[s] = num_inc[s]+1
        
      }
    }
  }
  res_emp = colSums(thresholded)
  plot(res_emp, xlab = sprintf("s = %d", s),ylim = c(min(res_emp), max(c(res_emp, thresh))))
  abline(h=thresh,col=2)
  abline(h = max(res_emp) - 3*thresh/9,col=3)
  print(sprintf("s = %d, thresh = %f, max = %f, diff = %f", s, thresh, max(res),max(res)-thresh))
  #points(res_emp,col=4, pch=2)
  lines(res_emp,col=4, lwd=2)
  #lines(res + 2*sqrt(res*thresh*9/9) + thresh/9*9, col=2, lty=2)
  #lines(res -2*sqrt(res*thresh*9/9) -thresh*9/9, col=2, lty=2)
}



















# simu tail bound with larger subtraction
erfc = function(z){
  return(2*pnorm(z*sqrt(2), lower.tail=FALSE))
}

mu2 = function(c, s,aa){
  nu_aa  = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))
  EE = 1/2*(c -nu_aa-aa^2 +1)*erfc((aa-sqrt(c))/sqrt(2)) + exp(-(aa-sqrt(c))^2/2)*(aa+sqrt(c))/sqrt(2*pi)
  EE = EE + 1/2  *(c-nu_aa -aa^2 + 1)*erfc((aa+sqrt(c))/sqrt(2)) + exp(-(aa+sqrt(c))^2/2)*(aa-sqrt(c))/sqrt(2*pi)
  return(s*EE)
}
p = 1000
n = 1000

a = sqrt(4*log(p*log(n)))
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))

thet = seq(0, 2.5, by=0.05)*a
sums = matrix(NA, nrow = n, ncol = length(thet))
for (t in 1:length(thet)) {
  for (i in 1:n) {
    x = rnorm(p) + thet[t]
    
    for (j in 1:p) {
      if(abs(x[j])>a){
        x[j] = x[j]^2 - nu_a - 2*a^2
      }
      else{
        x[j]=0
      }
    }
    
    sums[i,t] = sum(x)
  }
}
hist(sums[,4])
rowmaxes = rep(NA, length(thet))
rowmeans = rep(NA, length(thet))
rowmins = rep(NA, length(thet))
for (i in 1:length(thet)) {
  rowmaxes[i] = max(sums[,i])
  rowmeans[i] = mean(sums[,i])
  rowmins[i] = min(sums[,i])
}
par(mfrow=c(1,1))
EE = 1/2*(thet^2 -nu_a -2*a^2+1)*erfc((a-thet)/sqrt(2)) + exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
EE = EE + 1/2 * (thet^2-nu_a -2*a^2+ 1)*erfc((a+thet)/sqrt(2)) + exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
#EE[EE<0] = 0
plot(thet/a, rowmaxes, type="l")
lines(thet/a, rowmeans, col=2)
abline(h = 9*log(n), col=5)
lines(thet/a, p*EE + 9*log(n) + 2*sqrt(p*(abs(EE)+a^2)*9*log(n)) ,col=3)
lines(thet/a, rowmins, col=7)
lines(thet/a, p*EE - 9*log(n) -2*sqrt(p*(abs(EE)+a^2)*9*log(n)),col=4)
lines(thet/a, p*EE, col=6)
legend("topleft", legend = c("emp max", "emp mean", "thresh", 
                             "thet upper", "emp min", "thet min", "thet mean"), 
       col= c(1,2,5,3,7,4,6),lty=c(1,1,1,1,1,1,1))

p1 = pmax(1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)),0)
lines(thet/a, p*EE- 9*log(n) - 2*sqrt(p*p5*9*log(n)), col=6)
sum(rowmaxes > p*EE + 9*log(n) + 2*sqrt(p*EE*9*log(n)))
sum(rowmins < p*EE - 9*log(n) - 2*sqrt(p*EE*9*log(n)))
#lines(thet/a, p*EE, type="l", col=4)
plot(thet[1:10]/a, rowmaxes[1:10])
abline(h = 9*log(n), col=2)
lines(thet[1:10]/a, p*EE[1:10] + 9*log(n) )



# simu tail bound with zero mean


p = 100000
n = 200

a = sqrt(4*log(p*log(n)))


ss = seq(1, round(sqrt(p*log(n^4))))
sums = matrix(NA, nrow = n, ncol = length(ss))
as = rep(NA, length(ss))
for (t in 1:length(ss)) {
  s = ss[t]

  a = sqrt(4*log(exp(1)*p*log(n^4)/s^2))
  nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
  as[t] = a
  for (i in 1:n) {
    x = rnorm(p) 
    
    for (j in 1:p) {
      if(abs(x[j])>a){
        x[j] = x[j]^2 - nu_a - 2*a^2
      }
      else{
        x[j]=0
      }
    }
    
    sums[i,t] = sum(x)
  }
}
hist(sums[,4])
rowmaxes = rep(NA, length(ss))
rowmeans = rep(NA, length(ss))
rowmins = rep(NA, length(ss))
for (i in 1:length(ss)) {
  rowmaxes[i] = max(sums[,i])
  rowmeans[i] = mean(sums[,i])
  rowmins[i] = min(sums[,i])
}
par(mfrow=c(1,1))
plot(ss, rowmaxes, type="l", ylim = c(min(c(rowmins, - 9*pmax(log(n^4), ss*as^2))), 
                                      max(rowmaxes,  9*pmax(log(n^4), ss*as^2))))
lines(ss, rowmins, type="l",col=2)
lines(ss, 9*pmax(log(n^4), ss*as^2),col=3)
lines(ss,- 9*pmax(log(n^4), ss*as^2),col=4)



## hm does the previous bound really work??

logn4 = 5
logp = 50
s = round(sqrt(logn4)*exp(logp/2)/8)
a = sqrt(4*(1+logp + log(logn4) - 2*log(s)))
#a = sqrt()
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
EE = 1/2*(0 -nu_a -2*a^2/s^2+1)*erfc((a)/sqrt(2)) + exp(-(a)^2/2)*(a)/sqrt(2*pi)
EE = EE + 1/2 * (0-nu_a -2*a^2/s^2+ 1)*erfc((a)/sqrt(2)) + exp(-(a)^2/2)*(a)/sqrt(2*pi)
exp(logp)*EE
s*9*(a)
