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


# new simu tail bound
#N = 100

# the following simu shows that for t <= log(n^4), 
#the liu statistic does not guarantee small signals to be undetected, 
# but the projection statistic does. 
p = 10
n = 1000000
N = n*1
s = floor(min(sqrt(p*log(n)), p))
k = s
C1 = 1

sums = matrix(0, nrow = N, ncol =s)
sums2 = matrix(0, nrow = N, ncol =s)
maxtots = matrix(0, nrow = N, ncol =s)
as = sums[1,]
for (iter in 1:N) {
  y = rnorm(p) 
  # diff = c()
  # if(k<s){
  #   diff = c(sqrt(16*C1)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n)/k^2)), log(n)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))
  # }
  # else{
  #   diff = rep(1, p)/sqrt(p) *sqrt(16*C1)*(p*log(n))^(1/4)
  # }
  # y = y+ diff

  #y[1:s] = y[1:s] + as[1]
  for (t in 1:s) {
    maxzz = rep(0,p)
    diff = c()
    k = t
    #if(k<s){
      #diff = c(sqrt(C1)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n)/k^2)), log(n)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))
      diff = c(sqrt(C1)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n)/k^2)), log(n)))*sample(c(-1,1), replace=TRUE,p))
    #}
    x = y + diff
    a = sqrt(4*sqrt(C1)*log(exp(1)*p*log(n)/t^2) + 8*sqrt(C1)*log(n)/t)
    as[t] = a
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
    for (j in 1:p) {
      if(abs(x[j])>a){
        
        x[j] = x[j]^2 - nu_a
        
      }
      else{
        x[j]=0
      }
    }
    sums[iter,t] = sum(x) 
    x = (y+diff)^2
    xs = sort(x, decreasing=TRUE)
    sums2[iter, t] = sum(xs[1:t])

    
    
    
  }
}
pens = 9 *(pmax(1:s*as^2, log(n)))
plot(colmaxes(sums), ylim = c(min(c(colmaxes(sums), pens)), max(c(colmaxes(sums), pens))))
points(colmaxes(sums2))
lines(pens)

fracs = rep(NA, dim(sums)[2])
# compute fraction above thresh
pens = C1 *(pmax(1:s*as^2, log(n)))
for (i in 1:dim(sums)[2]) {
  fracs[i] = mean(sums[,i] - pens[i] >0)
}

fracs








# new simu tail bound
#N = 100
p = 3
n = 10000000
N = n
s = floor(min(sqrt(p*log(n)), p))
k = min(3, s)
C1 = 9

#s = 500
#aa = sqrt(4*log(exp(1)*p*log(n)/s^2))
#nu_a = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))

#thet = seq(0, 1.2, by=0.05)*a
sums = matrix(0, nrow = N, ncol =s)
maxtots = matrix(0, nrow = N, ncol =s)
as = sums[1,]
for (iter in 1:N) {
  y = rnorm(p) 
  diff = c()
  if(k<=s){
    diff = c(sqrt(C1/4)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n)/k^2)), log(n)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))
  }
  else{
    diff = rep(1, p)/sqrt(p) *sqrt(C1/4)*(p*log(n))^(1/4)
  }
  y = y+ diff
  
  #y[1:s] = y[1:s] + as[1]
  for (t in 1:s) {
    maxzz = rep(0,p)
    x = y[]
    a = sqrt(8*log(exp(1)*p*log(n)/t^2) +8*log(n)/t)
    #a = sqrt(8*log(exp(1)*p*log(n)/t^2))
    as[t] = a
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
    for (j in 1:p) {
      if(abs(x[j])>a){
        
        x[j] = x[j]^2 - nu_a #-C1*a^2
        
      }
      else{
        x[j]=0
      }
    }
    sums[iter,t] = sum(x) #+ C1*p*2*a^2*pnorm(a, lower.tail=FALSE)
    
    
    
    
  }
}

fracs = rep(NA, dim(sums)[2])
# compute fraction above thresh
pens = C1 *((1:s*as^2/4+ log(n)))
for (i in 1:dim(sums)[2]) {
  fracs[i] = mean(sums[,i] - pens[i] >0)
}

fracs

colmaxes(sums)-pens
cbind(colmaxes(sums), pens)
k



hist(sums[,4])
rowmaxes = rep(NA, s)
rowmeans = rep(NA, s)
rowmins = rep(NA,s)
rowmaxestots = rep(NA, s)
rowmeanstots = rep(NA, s)
rowminstots = rep(NA,s)
for (i in 1:s) {
  rowmaxes[i] = max(sums[,i])
  rowmeans[i] = mean(sums[,i])
  rowmins[i] = min(sums[,i])
  rowmaxestots[i] = max(maxtots[,i])
  rowmeanstots[i] = mean(maxtots[,i])
  rowminstots[i] = min(maxtots[,i])
}

# plot(rowmaxestots)
# 
# EE = 1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)) + exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
# EE = EE + 1/2 * (thet^2-nu_a + 1)*erfc((a+thet)/sqrt(2)) + exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
# EE[EE<0] = 0

#plot(1:s, rowmaxes + p*2*pnorm(as, lower.tail=FALSE), type="l", ylim=c(min(rowmins), max(rowmaxes)))
plot(1:s, rowmaxes + 2*as^2*p*2*pnorm(as, lower.tail=FALSE),
     ylim = c(min(c(rowmins+2*as^2*p*2*pnorm(as, lower.tail=FALSE),9*(- (1:s) * as^2 - log(n)))) ,
              max(c(rowmaxes+2*as^2*p*2*pnorm(as, lower.tail=FALSE), 9*(1:s)*as^2 +9*log(n)))) , type="l")
lines(1:s, rowmeans + 2*as^2*p*2*pnorm(as, lower.tail=FALSE), type="l",col=2)
lines(1:s, rowmins+2*as^2*p*2*pnorm(as, lower.tail=FALSE), type="l",col=3)
lines(1:s, 9*(1:s * as^2 + log(n)),col=4,type="l")
lines(1:s, 9*(-(1:s) * as^2-log(n)),col=4,type="l")

lines(1:s, sqrt(s*aa^2*1:s*as^2))
lines(thet/a, rowmins, type="l",col=2)
lines(thet/a, rowmeans, type="l",col=3)
9*a^2









## counting thing


# new simu tail bound
#N = 100
p =1000
n = 1000
N = n
s = floor(min(sqrt(p*log(n)), p))
ss = 1:s
k = min(s, s)
C1 = 4

#s = 500
#aa = sqrt(4*log(exp(1)*p*log(n)/s^2))
#nu_a = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))

#thet = seq(0, 1.2, by=0.05)*a
sums = matrix(0, nrow = N, ncol =s)
counts = matrix(0, nrow = N, ncol =s)

maxtots = matrix(0, nrow = N, ncol =s)
as = sums[1,]
for (iter in ceiling(log(n)):N) {
  y = rnorm(p) 
  diff = c()
  if(k<=s){
    diff = c(sqrt(C1)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n)/k^2)), log(n)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))
  }
  else{
    diff = rep(1, p)/sqrt(p) *sqrt(C1)*(p*log(n))^(1/4)
  }
  y = y+ diff
  
  #y[1:s] = y[1:s] + as[1]
  for (t in 1:s) {
    maxzz = rep(0,p)
    x = y[]
    a = sqrt(4*log(exp(1)*p*log(n)/t^2))
    #a = sqrt(8*log(exp(1)*p*log(n)/t^2))
    as[t] = a
    count = 0
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
    for (j in 1:p) {
      if((x[j])^2> C1*(a^2 + log(n)/t))
      {
        count = count+1
      }
      if(abs(x[j])>a){
        
        x[j] = x[j]^2 - nu_a #-C1*a^2
        
      }
      else{
        x[j]=0
      }
    }
    sums[iter,t] = sum(x) #+ C1*p*2*a^2*pnorm(a, lower.tail=FALSE)
    counts[iter, t] = count
    
    
    
  }
}



colmaxes(counts)
colMeans(counts)
cbind(colmaxes(counts), ss)








# larger a

# new simu tail bound
#N = 100
p =10000
n = 100
N = n
s = floor(min(sqrt(p*log(n)), p))
ss = 1:s
k = min(100, s)
C1 = 4

#s = 500
#aa = sqrt(4*log(exp(1)*p*log(n)/s^2))
#nu_a = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))

#thet = seq(0, 1.2, by=0.05)*a
sums = matrix(0, nrow = N, ncol =s)
counts = matrix(0, nrow = N, ncol =s)

maxtots = matrix(0, nrow = N, ncol =s)
as = sums[1,]
for (iter in ceiling(log(n)):N) {
  y = rnorm(p) 
  diff = c()
  if(k<=s){
    diff = c(sqrt(C1)*1/sqrt(k)*max(c(sqrt(k*log(exp(1)*p*log(n)/k^2)), log(n)))*sample(c(-1,1), replace=TRUE,k), rep(0,p-k))
  }
  else{
    diff = rep(1, p)/sqrt(p) *sqrt(C1)*(p*log(n))^(1/4)
  }
  y = y+ diff
  
  #y[1:s] = y[1:s] + as[1]
  for (t in 1:s) {
    maxzz = rep(0,p)
    x = y[]
    a = sqrt(4*C1*log(exp(1)*p*log(n)/t^2))
    #a = sqrt(8*log(exp(1)*p*log(n)/t^2))
    as[t] = a
    count = 0
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
    for (j in 1:p) {
      if((x[j])^2> C1*(a^2 + log(n)/t))
      {
        count = count+1
      }
      if(abs(x[j])>a){
        
        x[j] = x[j]^2 - nu_a #-C1*a^2
        
      }
      else{
        x[j]=0
      }
    }
    sums[iter,t] = sum(x) #+ C1*p*2*a^2*pnorm(a, lower.tail=FALSE)
    counts[iter, t] = count
    
    
    
  }
}



colmaxes(sums)
colMeans(sums)
pens = 9*pmax(ss*log(exp(1)*p*log(n)/ss^2), log(n))
plot(ss, colmaxes(sums), ylim=c(min(c(colmaxes(sums),pens )), max(c(colmaxes(sums)),pens)))
lines(ss, pens)
abline(v = log(n))
abline(v= k)




mu = function(c, ss,aa){
  nu_aa  = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))
  EE1 = 1/2*(c -nu_aa +1)*erfc((aa-sqrt(c))/sqrt(2)) + exp(-(aa-sqrt(c))^2/2)*(aa+sqrt(c))/sqrt(2*pi)
  EE2 =  1/2  *(c-nu_aa + 1)*erfc((aa+sqrt(c))/sqrt(2)) + exp(-(aa+sqrt(c))^2/2)*(aa-sqrt(c))/sqrt(2*pi)
  EE = EE1+EE2
  return(list(ss*EE, EE1, EE2))
}
C1 = 4
tt = ceiling(log(n)):s
means = matrix(NA, ncol = length(tt), nrow = length(tt))
EE1s = matrix(NA, ncol = length(tt), nrow = length(tt))
EE2s = matrix(NA, ncol = length(tt), nrow = length(tt))

for (t in tt) {
  for (z in t:s) {
    a = sqrt(4*C1*log(exp(1)*p*log(n)/t^2))
    b = sqrt(C1*log(exp(1)*p*log(n)/z^2))
    res = mu(b^2, z, a)
    means[z - tt[1],t-tt[1]] = res[[1]]
    EE1s[z - tt[1],t-tt[1]] = res[[2]]
    EE2s[z - tt[1],t-tt[1]] = res[[3]]
    
  }
}

plot(means[,1])


p =100000
n = 100
N = n
s = floor(min(sqrt(p*log(n)), p))
#tt = ceiling(log(n)):s
tt = ceiling(log(log(n), base=2)):floor(log(s, base=2))
tt = c(tt, log(s, base=2))
k = min(100, s)
C1 = 4

means = matrix(NA, ncol = length(tt), nrow = log(p, base=2)+2)
as = tt[]


for (t in 1:length(tt)) {
  ttt = 2^(tt[t])
  already = 0
  for (hh in 1:p) {
    z = s + 2^(hh-1)
    if(z>p){
      if(already){
        z=p
        break
      }
      else{already=1;z = p}
    }
    a = sqrt(4*C1*log(exp(1)*p*log(n)/ttt^2))
    b = C1*sqrt(p*log(n))
    res = (mu(b/z, z, a))
    means[hh ,t] = res[[1]]
    as[t] = a

    
  }
}
par(mfrow=c(1,1))
plot(means[,length(tt)])
abline(h = b)
plot(means[,1])


# larger a and nu_tilde (conditional on non-zero mean)
erfc = function(z){
  return(2*pnorm(z*sqrt(2), lower.tail=FALSE))
}

# new simu tail bound
#N = 100
p =10000
n = 100
N = n
s = floor(min(sqrt(p*log(n)), p))
ss = 1:s
k = min(25, s)
C1 = 1/4

#s = 500
#aa = sqrt(4*log(exp(1)*p*log(n)/s^2))
#nu_a = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))

#thet = seq(0, 1.2, by=0.05)*a
sums = matrix(0, nrow = N, ncol =s)
counts = matrix(0, nrow = N, ncol =s)

maxtots = matrix(0, nrow = N, ncol =s)
as = sums[1,]
nu_as = sums[1,]
nu_a_tildes = sums[1,]
for (iter in ceiling(log(n)):N) {
  y = rnorm(p) 
  
  
  #y[1:s] = y[1:s] + as[1]
  for (t in 1:s) {
    
    maxzz = rep(0,p)
    x = y[]
    a_tilde= sqrt(log(exp(1)*p*log(n)/t^2))
    a = 4*sqrt(C1)*a_tilde
    #a = sqrt(8*log(exp(1)*p*log(n)/t^2))
    as[t] = a
    count = 0
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
    b = sqrt(C1)*a_tilde
    #diff = c()
    #diff = rep(b, p)
    #x = y+diff
    #y = 
    nu_a_tilde = 0.5*(b^1+1)*erfc(3*b/sqrt(2)) + 5*exp(-9*b^2/2)*b/sqrt(2*pi)
    nu_a_tilde = nu_a_tilde + 0.5*(b^1+1)*erfc(5*b/sqrt(2)) + 3*exp(-25*b^2/2)*b/sqrt(2*pi)
    nu_a_tilde = nu_a_tilde / (pnorm(3*b,lower.tail=FALSE) + pnorm(5*b, lower.tail=FALSE))
    nu_a_tildes[t] = nu_a_tilde
    for (j in 1:p) {
      if((x[j])^2> C1*(a^2 + log(n)/t))
      {
        count = count+1
      }
      if(abs(x[j])>a){
        
        x[j] = x[j]^2 - nu_a_tilde #-C1*a^2
        
      }
      else{
        x[j]=0
      }
    }
    sums[iter,t] = sum(x) #+ C1*p*2*a^2*pnorm(a, lower.tail=FALSE)
    counts[iter, t] = count
    
    
    
  }
}



colmaxes(sums)
colMeans(sums)
colmins(sums)
plot(colmins(sums))
plot(ss, colmaxes(sums), ylim=c(min(c(colmaxes(sums),9*pmax(ss*as^2, log(n)) )), max(c(colmaxes(sums)),9*pmax(ss*as^2, log(n)) )))
lines(ss, 9*pmax(ss*as^2, log(n)))
cbind(colmaxes(counts), ss)
























fracs = rep(NA, dim(sums)[2])
# compute fraction above thresh
pens = C1 *((1:s*as^2/4+ log(n)))
for (i in 1:dim(sums)[2]) {
  fracs[i] = mean(sums[,i] - pens[i] >0)
}

fracs

colmaxes(sums)-pens
cbind(colmaxes(sums), pens)
k



hist(sums[,4])
rowmaxes = rep(NA, s)
rowmeans = rep(NA, s)
rowmins = rep(NA,s)
rowmaxestots = rep(NA, s)
rowmeanstots = rep(NA, s)
rowminstots = rep(NA,s)
for (i in 1:s) {
  rowmaxes[i] = max(sums[,i])
  rowmeans[i] = mean(sums[,i])
  rowmins[i] = min(sums[,i])
  rowmaxestots[i] = max(maxtots[,i])
  rowmeanstots[i] = mean(maxtots[,i])
  rowminstots[i] = min(maxtots[,i])
}

# plot(rowmaxestots)
# 
# EE = 1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)) + exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
# EE = EE + 1/2 * (thet^2-nu_a + 1)*erfc((a+thet)/sqrt(2)) + exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
# EE[EE<0] = 0

#plot(1:s, rowmaxes + p*2*pnorm(as, lower.tail=FALSE), type="l", ylim=c(min(rowmins), max(rowmaxes)))
plot(1:s, rowmaxes + 2*as^2*p*2*pnorm(as, lower.tail=FALSE),
     ylim = c(min(c(rowmins+2*as^2*p*2*pnorm(as, lower.tail=FALSE),9*(- (1:s) * as^2 - log(n)))) ,
              max(c(rowmaxes+2*as^2*p*2*pnorm(as, lower.tail=FALSE), 9*(1:s)*as^2 +9*log(n)))) , type="l")
lines(1:s, rowmeans + 2*as^2*p*2*pnorm(as, lower.tail=FALSE), type="l",col=2)
lines(1:s, rowmins+2*as^2*p*2*pnorm(as, lower.tail=FALSE), type="l",col=3)
lines(1:s, 9*(1:s * as^2 + log(n)),col=4,type="l")
lines(1:s, 9*(-(1:s) * as^2-log(n)),col=4,type="l")

lines(1:s, sqrt(s*aa^2*1:s*as^2))
lines(thet/a, rowmins, type="l",col=2)
lines(thet/a, rowmeans, type="l",col=3)
9*a^2

