rel = seq(0, 2, by = 0.01)

a = 4
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
maxval = rep(NA, length(rel))

K = 1000000
tmp = rep(NA, K)
for (i in 1:length(rel)) {
  tmp = rnorm(K) + a*rel[i]
  maxval[i] = max(tmp[abs(tmp)>a]^2 -nu_a)
}
plot(rel, maxval)
lines(rel, (a*rel)^2 + 15*a*rel + log(K),col=2)
abline(h = a^2)



## E(X)
erfc = function(z){
  return(2*pnorm(z*sqrt(2), lower.tail=FALSE))
}
a = 4
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
thet = seq(0.1, 2, by=0.001)* a
if(a == 0){
  thet = seq(0.1, 1)
}
EE = 1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)) + exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
EE = EE + 1/2 * (thet^2-nu_a + 1)*erfc((a+thet)/sqrt(2)) + exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
plot(thet, EE, type="l")
p1 = 1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2))
p2 = exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
p3 = 1/2 * (thet^2-nu_a + 1)*erfc((a+thet)/sqrt(2))
p4 = exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
matplot(thet, cbind(EE, p1, p2, p3, p4), type="l")
legend(x="topleft", legend = c("EE", "p1", "p2", "p3", "p4"), col=c(1,2,3,4), lty=c(1,2,3,4))

matplot(thet, cbind(EE, p1, p2, p3, p4), type="l")
lines(thet, thet^2*exp(-(a-thet)^2/2)/10)
legend(x="topleft", legend = c("EE", "p1", "p2", "p3", "p4"), col=c(1,2,3,4), lty=c(1,2,3,4))

lines(thet,thet^2)
p5 = 1/2*(thet^2)*erfc((a-thet)/sqrt(2))
matplot(cbind(EE, p1, p2, p3, p4,p5), type="l")
legend(x="topleft", legend = c("EE", "p1", "p2", "p3", "p4","p5"), col=c(1,2,3,4,5,6), lty=c(1,2,3,4,5,6))
lines(thet^2)

lines(thet, thet^2,col=2)
plot(thet^2 - nu_a - EE)
max(thet^2 - nu_a - EE)
EE/a

# diff:
a = 1
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
EE = 1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)) + exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
EE = EE + 1/2 * (thet^2-nu_a + 1)*erfc((a+thet)/sqrt(2)) + exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
plot(EE, thet^2)
#plot(diff(thet^2)/diff(EE),type="l")
plot(diff(EE)/diff(thet^2),type="l")
dd1 = 1/2*erfc((a-thet)/sqrt(2)) + 1/2*erfc((a+thet)/sqrt(2))
dd2 = 1/(2*sqrt(2*pi)*thet)*(exp(-1/2*(a-thet)^2)  -exp(-1/2*(a+thet)^2)) *
  (a^2+2-nu_a)
lines( dd1+dd2,col=2)
# tail bound
t = 1
thet = 4*a
s = 1
ls = seq(0.001, 0.24, by = 0.0001)
plot(ls, 2*exp(-ls*t-s*(a - thet - 2*ls*thet/(1-2*ls))^2*(1-2*ls)/2 + 2*s*ls^2*thet^2/(1-2*ls) +
             ls*thet^2*s - ls*nu_a*s )*sqrt(1-2*ls))
plot(ls, -ls*t-s*(a - thet - 2*ls*thet/(1-2*ls))^2*(1-2*ls)/2 + 2*s*ls^2*thet^2/(1-2*ls) +
                 ls*thet^2*s - ls*nu_a*s )




# simu tail bound

#multiple change-points:
#multiple change-points:
p = 1000
n = 1000

a = sqrt(4*log(p*log(n)))
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))

thet = seq(0, 1.2, by=0.05)*a
sums = matrix(NA, nrow = n, ncol = length(thet))
for (t in 1:length(thet)) {
  for (i in 1:n) {
    x = rnorm(p) + thet[t]

    for (j in 1:p) {
      if(abs(x[j])>a){
        x[j] = x[j]^2 - nu_a
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

EE = 1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)) + exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
EE = EE + 1/2 * (thet^2-nu_a + 1)*erfc((a+thet)/sqrt(2)) + exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
EE[EE<0] = 0
plot(thet/a, rowmaxes, type="l")
lines(thet/a, rowmeans, col=2)
abline(h = 9*log(n), col=5)
lines(thet/a, p*EE + 9*log(n) + 2*sqrt(p*EE*9*log(n)) ,col=3)
p5 = 1/2*(thet^2)*erfc((a-thet)/sqrt(2))
lines(thet/a, rowmins, col=7)
lines(thet/a, p*EE - 9*log(n) -2*sqrt(p*EE*9*log(n)),col=5)
p1 = pmax(1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)),0)
lines(thet/a, p*EE- 9*log(n) - 2*sqrt(p*p5*9*log(n)), col=6)
sum(rowmaxes > p*EE + 9*log(n) + 2*sqrt(p*EE*9*log(n)))
sum(rowmins < p*EE - 9*log(n) - 2*sqrt(p*EE*9*log(n)))
#lines(thet/a, p*EE, type="l", col=4)
plot(thet[1:10]/a, rowmaxes[1:10])
abline(h = 9*log(n), col=2)
lines(thet[1:10]/a, p*EE[1:10] + 9*log(n) )


## seems like mu is bounded above by p5
# and that max A < mu + thresh (9 log n) + 2sqrt(mu * 9 log n)




## new idea: threshold harder to guarantee a >2 (?)


#multiple change-points:
#multiple change-points:
p = 1000
n = 1000

a = sqrt(4*log(p*log(n)))
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
thet = 2*a
EE = 1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)) + exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
EE = EE + 1/2 * (thet^2-nu_a + 1)*erfc((a+thet)/sqrt(2)) + exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
EE[EE<0] = 0
thet = seq(0, 2, by=0.05)*a
sums = matrix(NA, nrow = n, ncol = length(thet))
for (t in 1:length(thet)) {
  for (i in 1:n) {
    x = rnorm(p) + thet[t]
    
    for (j in 1:p) {
      if(abs(x[j])>a){
        x[j] = (x[j]^2 - nu_a )
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

plot(thet/a, rowmaxes, type="l")
lines(thet/a, rowmeans, type="l", col=2)
lines(thet/a, rowmins, type="l", col=3)      
abline(h = 9*log(n))
lines(thet/a, p*thet^2, type="l")
# erstatt nu_a med E((Z+a)^2 | (Z+a) >a) = a*a +1






## better to decrease a?

p = 1000
n = 1000

a = sqrt(4*log(exp(1)*p*log(n)))
a2 = sqrt(4*log(exp(1)*p*log(n))/4)
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
nu_a2 = 1 + a2*exp(dnorm(a2, log=TRUE)-pnorm(a2, lower.tail = FALSE, log=TRUE))
#EE = 1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)) + exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
#EE = EE + 1/2 * (thet^2-nu_a + 1)*erfc((a+thet)/sqrt(2)) + exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
#EE[EE<0] = 0
thet = 3/4*a
sums = matrix(NA, nrow = n, ncol = 2)

for (i in 1:n) {
  y = rnorm(p) + thet
  x = y[]
  for (j in 1:p) {
    if(abs(x[j])>a){
      x[j] = (x[j]^2 - nu_a )
    }
    else{
      x[j]=0
    }
  }
  
  sums[i,1] = sum(x)
  
  x = y[]
  for (j in 1:p) {
    if(abs(x[j])>a2){
      x[j] = (x[j]^2 - nu_a2 )
    }
    else{
      x[j]=0
    }
  }
  
  sums[i,2] = sum(x)
}

hist(sums[,1])
hist(sums[,2])
hist(sums[,2]-sums[,1])
mean(sums[,2]-sums[,1])







# new simu 





n = 3
p = 10000000
s = 1000
#a = sqrt(4*log(exp(1)*p*log(n^4)))
#psi = 300
t = 1
a = sqrt(4*log(exp(1)*p*log(n^4)/t^2))
#a=1
#a=0
psi = t*log(exp(1)*p*log(n^4)/t^2)
mu = function(c, s,aa){
  nu_aa  = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))
  EE = 1/2*(c -nu_aa +1)*erfc((aa-sqrt(c))/sqrt(2)) + exp(-(aa-sqrt(c))^2/2)*(aa+sqrt(c))/sqrt(2*pi)
  EE = EE + 1/2  *(c-nu_aa + 1)*erfc((aa+sqrt(c))/sqrt(2)) + exp(-(aa+sqrt(c))^2/2)*(aa-sqrt(c))/sqrt(2*pi)
  return(s*EE)
}
#psis = c(1,2,3,4,5,6,7,8,9,10,20,30,50,60,100)
#ss = c(1,2,3,5,8,10,20,30,50,100,500,1000,10000,1000000,100000000,1e+10)
ss = seq(1, round(sqrt(p*log(n^4))))
#a=100
c_bs = rep(NA, length(ss))
c_etas = rep(NA, length(ss))
mu_inv_pen = rep(NA, length(ss))
tts = rep(NA, length(ss))
cbbetter = rep(NA, length(ss))
for (i in 1:length(ss)) {
  s = ss[i]
  minimum_mu_b =256*(psi+log(n^4))
  #minimum_mu_b =a^2
  #minimum_mu_b=0.00000001
  mudiff = function(c){
    return(minimum_mu_b - mu(c, s, a))
  }
  c_b_times_s = s*uniroot(mudiff,interval = c(0,100000000))$root
  diff = 18*(psi+log(n^4))
  mu_eta = diff + minimum_mu_b
  mudiff = function(c){
    return(mu_eta - mu(c, s, a))
  }
  
  c_eta_times_s = s*uniroot(mudiff, interval=c(0,100000000))$root
  c_bs[i] = c_b_times_s
  c_etas[i] = c_eta_times_s
  mudiff = function(c){
    return(diff - mu(c, s, a))
  }
  mu_inv_pen[i] = s*uniroot(mudiff,interval = c(0,100000000))$root
  
  c_b = c_b_times_s/s
  cbbetter[i] = 1000
  if(c_b< a^2){
    for (tt in (t+1):round(sqrt(p*log(n^4)))) {
      att = sqrt(4*log(exp(1)*p*log(n^4)/tt^2))
      if(att^2<=c_b){
        break
      }
    }
    if(TRUE){
    #if(att^2>c_b){
    
    
      tt = p
      tts[i] = tt
      cbbetter[i] = (s*c_b - minimum_mu_b) - 256*(sqrt(p*log(n^4)) - 256*(psi+log(n^4)))
    }else{
      tts[i] = tt
      cbbetter[i] = (mu(c_b, tt, att) - minimum_mu_b) - 256*(tt*att^2 - a^2)
    }
  }
  
  
  
}

plot(cbbetter)
min(cbbetter)

cbind(c_bs, c_etas)
c_etas - c_bs

plot(c_etas - c_bs)

matplot(cbind(c_etas, c_bs), type="l")
cbind(c_etas - c_bs, mu_inv_pen)
plot((c_etas - c_bs),type="l")
abline(h = diff)
plot((c_etas - c_bs),type="l")
lines(mu_inv_pen,col=2)

plot((c_etas - c_bs),type="l")
lines(c_bs, col=2)

plot((c_etas - c_bs)/c_bs)

plot(c_bs/ss - a^2)


#plot(log(((c_etas - c_bs)/c_bs)))
#plot(exp(((c_etas - c_bs)/c_bs)))







# for t < s, can we really get T(t) > (C_1 - 9) xi(t) 
# if C^b/s <a^2?
n = 20000
p = 3000
s = p

t = 1
a = sqrt(4*log(exp(1)*p*log(n^4)/t^2))
nu_a= 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))

psi = t*log(exp(1)*p*log(n^4)/t^2)
mu = function(c, s,aa){
  nu_aa  = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))
  EE = 1/2*(c -nu_aa +1)*erfc((aa-sqrt(c))/sqrt(2)) + exp(-(aa-sqrt(c))^2/2)*(aa+sqrt(c))/sqrt(2*pi)
  EE = EE + 1/2  *(c-nu_aa + 1)*erfc((aa+sqrt(c))/sqrt(2)) + exp(-(aa+sqrt(c))^2/2)*(aa-sqrt(c))/sqrt(2*pi)
  return(s*EE)
}
#psis = c(1,2,3,4,5,6,7,8,9,10,20,30,50,60,100)
#ss = c(1,2,3,5,8,10,20,30,50,100,500,1000,10000,1000000,100000000,1e+10)
ss = seq(1, round(sqrt(p*log(n^4))))

mu(a^2, s, a)
a^2*t*10

# mu(a, s, a) is never above a^2*t when s<=p!



n = 200000
p = 300
s = 10

t = 1
a = sqrt(4*log(exp(1)*p*log(n^4)/t^2))
nu_a= 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))

psi = t*log(exp(1)*p*log(n^4)/t^2)
mu = function(c, s,aa){
  nu_aa  = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))
  EE = 1/2*(c -nu_aa +1)*erfc((aa-sqrt(c))/sqrt(2)) + exp(-(aa-sqrt(c))^2/2)*(aa+sqrt(c))/sqrt(2*pi)
  EE = EE + 1/2  *(c-nu_aa + 1)*erfc((aa+sqrt(c))/sqrt(2)) + exp(-(aa+sqrt(c))^2/2)*(aa-sqrt(c))/sqrt(2*pi)
  return(s*EE)
}
#psis = c(1,2,3,4,5,6,7,8,9,10,20,30,50,60,100)
#ss = c(1,2,3,5,8,10,20,30,50,100,500,1000,10000,1000000,100000000,1e+10)
ss = seq(1, round(sqrt(p*log(n^4))))

cc = 16*s*log(exp(1)*p*log(n^4)/s^2)
ress = rep(NA, round(sqrt(p*log(n^4))))
az = rep(NA, round(sqrt(p*log(n^4))))
ss = 1:round(sqrt(p*log(n^4)))
for (t in 1:round(sqrt(p*log(n^4)))) {
  a = sqrt(4*log(exp(1)*p*log(n^4)/t^2))
  az[t] = a
  ress[t] = mu(cc/s, s, aa=a)
}
#plot(ress)
#lines(az^2*ss*9,col=2)
#plot(ress-az^2*ss*9 )
plot(ress - 9*pmax(az^2*ss, log(n^4)))
abline(v = s, col=2)







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
k = p
#randomsparse = rnorm(k)
#randomsparse = randomsparse/norm(randomsparse, type="2")
# randomsparse = runif(k, min=-1, max = 1)
# randomsparse = randomsparse/norm(randomsparse, type="2")
# randomsparse = randomsparse* 3/sqrt(n)*1000/k

#randomsparse = 2.5*c(2/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*rnorm(k), rep(0,p-k))
if(k<sqrt(p*log(n^4))){
  randomsparse = c(20/sqrt(n)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n^4)/k^2)), log(n^4)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))
}else{
  randomsparse = c(20*rep(1,k)/sqrt(k)*p^(1/4), rep(0, p-k))
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
  thresh = 128*sqrt(p*log(n^4))
  if(s < sqrt(p*log(n^4))){
    thresh =128*(max(s*log(exp(1)*p*log(n^4)/s^2),log(n^4)))
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
  plot(res, xlab = sprintf("s = %d", s),ylim = c(0, max(c(res_emp, res))))
  abline(h=thresh,col=2)
  abline(h = max(res) - 18*thresh/128,col=3)
  print(sprintf("s = %d, thresh = %f, max = %f, diff = %f", s, thresh, max(res),max(res)-thresh))
  points(res_emp,col=4, pch=2)
  lines(res + 2*sqrt(res*thresh*9/128) + thresh/128*9, col=2, lty=2)
  lines(res -2*sqrt(res*thresh*9/128) -thresh*9/128, col=2, lty=2)
}

