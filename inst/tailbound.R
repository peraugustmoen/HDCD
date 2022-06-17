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

## mgf

a = 2
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
thet = a
EE = 1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)) + exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
EE = EE + 1/2 * (thet^2-nu_a + 1)*erfc((a+thet)/sqrt(2)) + exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
mu = EE
ls = seq(0, 1/4, by=0.001)
plot(ls,exp(ls*(thet^2 - nu_a - mu) + 2*ls*a^2 - a^2 - thet^2)/sqrt(1-2*ls))
lines(ls, exp(ls*(thet^2 - nu_a - mu) +mu*ls/(1-2*ls) )/sqrt(1-2*ls), col=2)


## mgf 2 for thet < a

a = 4
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
thet = a
EE = 1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)) + exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
EE = EE + 1/2 * (thet^2-nu_a + 1)*erfc((a+thet)/sqrt(2)) + exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
mu = EE
ls = seq(0, 1/4, by=0.001)
aa = exp(ls*(thet^2 - nu_a - mu) + 2*ls*a^2 - a^2 - thet^2-1/2*log(1-2*ls))
## seems like, as l ->0, aa -> e^(-2a^2)
aa[1]
exp(-2*a^2)
bb = exp(-mu*(1-2*ls) - mu/2*log(1-2*ls))
cc =exp(-ls*mu - mu/2*log(1-2*ls)+ls^2*mu/(1-2*ls)) #lovende
plot(ls,aa,
     ylim = c(min(c(aa, bb,cc)), max(c(aa, bb,cc))),type="l")

lines(ls, bb,col=2)
lines(ls, cc, col=3)
sum(bb<aa)
sum(cc<aa)


# simu
n = 10000000
sims = rep(0, n)
for (i in 1:n) {
  x = rnorm(1)+ thet
  if(abs(x)>a){
    sims[i] = x^2-nu_a
  }
  # x = rnorm(p) + thet
  # 
  # for (j in 1:p) {
  #   if(abs(x[j])>a){
  #     x[j] = x[j]^2 - nu_a
  #   }
  #   else{
  #     x[j]=0
  #   }
  # }
  # 
  # sims[i] = sum(x)
}

aa_emp = rep(NA, length(ls))
for (ll in 1:length(ls)) {
  aa_emp[ll] = mean(exp(ls[ll]*(sims - EE)))
}
plot(ls, aa_emp,type="l")
lines(ls, exp(-EE*ls) +aa,col="2")

hist(sums[,4])
rowmaxes = rep(NA, length(thet))
rowmeans = rep(NA, length(thet))
rowmins = rep(NA, length(thet))
for (i in 1:length(thet)) {
  rowmaxes[i] = max(sums[,i])
  rowmeans[i] = mean(sums[,i])
  rowmins[i] = min(sums[,i])
}



#mu/thet^2

## mgf 3 for thet > a

a = 70
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
thet = max(a,0.0001)
EE = 1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)) + exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
EE = EE + 1/2 * (thet^2-nu_a + 1)*erfc((a+thet)/sqrt(2)) + exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
mu = EE
ls = seq(0, 1/4, by=0.001)
aa = exp(ls*(thet^2 - nu_a - mu) + thet^2*ls^2/(1-2*ls))/sqrt(1-2*ls)
bb = exp(mu*ls^2*(1-2*ls) - mu/2*log(1-2*ls))
#cc = exp(-ls*mu - mu/2*log(1-2*ls)+2*ls^2*mu/(1-2*ls)) #lovende
cc = exp(- mu/2*log(1-2*ls)+2*ls*mu/(1-2*ls)) #lovende
#cc = exp((mu)*pmax(-ls - 1/2*log(1-2*ls),ls^2/(1-2*ls))) #lovende
plot(ls,aa,
     ylim = c(min(c(aa, bb,cc)), max(c(aa, bb,cc))),type="l")

lines(ls, bb,col=2)
lines(ls, cc,col=3)
sum(bb<aa)
sum(cc<aa)


lines(ls, exp(ls*(thet^2 - nu_a - mu) +mu*ls/(1-2*ls) )/sqrt(1-2*ls), col=2)



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
lines(thet/a, rowmins, type="l",col=2)
lines(thet/a, rowmeans, type="l",col=3)
9*a^2

plot(thet/a, rowmeans, col=2,type="l")
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


EE = 1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)) + exp(-(a-thet)^2/2)*(a+thet)/sqrt(2*pi)
EE = EE + 1/2 * (thet^2-nu_a + 1)*erfc((a+thet)/sqrt(2)) + exp(-(a+thet)^2/2)*(a-thet)/sqrt(2*pi)
EE[EE<0] = 0
lines(thet/a, rowmaxes, type="l")
plot(thet/a, rowmeans, col=2,type="l")
abline(h = 9*log(n), col=5)
lines(thet/a, p*EE + 9*log(n) + 2*sqrt(p*pmax(EE,sqrt(EE))*9*log(n)) ,col=3)
p5 = 1/2*(thet^2)*erfc((a-thet)/sqrt(2))
lines(thet/a, rowmins, col=7)
lines(thet/a, p*EE - 9*log(n) -2*sqrt(p*EE*9*log(n)),col=5)
p1 = pmax(1/2*(thet^2 -nu_a +1)*erfc((a-thet)/sqrt(2)),0)
lines(thet/a, p*EE- 9*log(n) - 2*sqrt(p*p5*9*log(n)), col=6)
sum(rowmaxes > p*EE + 9*log(n) + 2*sqrt(p*pmax(EE,sqrt(EE))*9*log(n)))
sum(rowmins < p*EE - 9*log(n) - 2*sqrt(p*pmax(EE,sqrt(EE))*9*log(n)))
lines(thet/a,  p*EE + 9*log(n) + 2*sqrt(p*thet^2*9*log(n)) )
#lines(thet/a, p*EE, type="l", col=4)
plot(thet[1:10]/a, rowmaxes[1:10])
abline(h = 9*log(n), col=2)
lines(thet[1:10]/a, p*EE[1:10] + 9*log(n) )





# test of second statistic in combination
# ie if significant found, ensure that each is in mean larger than blabla
# this is basically a check on estimation of s

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
n = 100000
p = 20
N=n
maxs = floor(sqrt(p*log(n^(1/4))))
sums  =  matrix(NA, ncol = maxs, nrow = N)
counts  =  matrix(NA, ncol = maxs, nrow = N)
counts2  =  matrix(NA, ncol = maxs, nrow = N)
as = sums[1,]
nu_as = sums[1,]
delta_as = sums[1,]
ss = 1:maxs
for (t in 1:N) {
  y = rnorm(p)
  for (i in 1:maxs) {
    x = y[]
    counter = 0
    counter2 = 0
    a = sqrt(4*log(exp(1)*p*log(n^(1/4))/i^2))
    as[i] = a
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
    nu_as[i] = nu_a
    delta_a = qnorm(1/2*(1+pnorm(a)))
    delta_as[i] = delta_a
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
    
    sums[t,i] = sum(x)
    counts[t,i] = counter
    counts2[t,i] = counter2
    
  }
}
# plot(1:maxs, colmaxes(counts) ,type="l", col=2)
# lines(1:maxs, colMeans(counts))
# lines(1:maxs, 1:maxs,col=3)
# lines(1:maxs, colmins(counts), type="l", col=4)

maxx = colmaxes(counts) - 2*p*pnorm(as,lower.tail=FALSE)
minn = colmins(counts) - 2*p*pnorm(as,lower.tail=FALSE)
maxx / ss
minn / ss

maxx2 = colmaxes(counts2) 
minn2 = colmins(counts2)
maxx2 / ss
minn2 / ss

plot(maxx / ss)
plot(minn / ss)
plot(1:maxs, maxx ,type="l", col=2,
     ylim=c(min(minn), max(maxx)))
lines(1:maxs, colmins(counts) - 2*p*pnorm(as,lower.tail=FALSE), type="l", col=4)
lines(1:maxs, colMeans(counts) - 2*p*pnorm(as,lower.tail=FALSE))
lines(1:maxs, -(1:maxs),col=5)

plot(1:maxs, 28*sqrt(pmax(maxx,0)) ,type="l", col=2,
     ylim=c(min(minn), max(28*maxx)))
lines(1:maxs,-sqrt( -colmins(counts) + 2*p*pnorm(as,lower.tail=FALSE)), type="l", col=4)
lines(1:maxs, colMeans(counts) - 2*p*pnorm(as,lower.tail=FALSE))
lines(1:maxs, -(1:maxs),col=5)
lines(1:maxs, (1:maxs),col=5)





ss = 1:maxs
plot(ss, ss*log(exp(1)*p*log(n^4)/ss^2),type="l")
abline(h = sqrt(p*log(n^4)))
lines(ss, ss^2*log(exp(1)*p*log(n^4)/ss^2) )

plot(1:maxs, colmaxes(sums))
plot(1:maxs, colmins(sums))



# easier binomial trial: 
n = 10000
p = 100
N=n
maxs = floor(sqrt(p*log(n)))
ss = 1:maxs
s_hats  =  matrix(NA, ncol = maxs, nrow = N)
as = s_hats[1,]

for (i in ss) {
  a = sqrt(2*log(exp(1)*p*log(n)/i^2) + log(n)/i)
  as[i] = a

  s_hats[,i]= rbinom(N, p, 2*pnorm(a,lower.tail=FALSE))
}
maxx = colmaxes(s_hats) - 2*p*pnorm(as,lower.tail=FALSE)
minn = colmins(s_hats) - 2*p*pnorm(as,lower.tail=FALSE)
maxx / ss
minn / ss
plot(1:maxs, maxx ,type="l", col=2,
     ylim=c(min(minn), max(maxx)))
lines(1:maxs, colmins(s_hats) - 2*p*pnorm(as,lower.tail=FALSE), type="l", col=4)
lines(1:maxs, colMeans(s_hats) - 2*p*pnorm(as,lower.tail=FALSE))
lines(1:maxs, -(1:maxs),col=5)
lines(1:maxs, (1:maxs),col=5)



# easier binomial trial with bound on prob
n = 50000
p = 10
N=n
maxs = floor(sqrt(p*log(n)))
ss = 1:maxs
s_hats  =  matrix(NA, ncol = maxs, nrow = N)
as = s_hats[1,]

for (i in ss) {
  a = sqrt(4*log(exp(1)*p*log(n)/i^2))
  as[i] = a
  
  s_hats[,i]= rbinom(N, p, 2*exp(-a^2/2))
}
maxx = colmaxes(s_hats) - 2*p*exp(-as^2/2)
minn = colmins(s_hats) - 2*p*exp(-as^2/2)

plot(1:maxs, maxx ,type="l", col=2,
     ylim=c(min(c(minn, -maxs)), max(c(maxx, maxs))))

lines(1:maxs, minn, type="l", col=4)
lines(1:maxs, colMeans(s_hats) - 2*p*exp(-as^2/2))
lines(1:maxs, -(1:maxs),col=5)
lines(1:maxs, (1:maxs),col=5)
lines(1:maxs, (1:maxs) + log(n)/ss,col=5)


# new simu tail bound
#N = 100
p = 10
n = 10000000
N = n
s = floor(sqrt(p*log(n)))
#s = 500
#aa = sqrt(4*log(exp(1)*p*log(n)/s^2))
#nu_a = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))

#thet = seq(0, 1.2, by=0.05)*a
sums = matrix(0, nrow = N, ncol =s)
maxtots = matrix(0, nrow = N, ncol =s)
as = sums[1,]
for (iter in 1:N) {
  y = rnorm(p) 
  #y[1:s] = y[1:s] + as[1]
  for (t in 1:s) {
    maxzz = rep(0,p)
    x = y[]
    a = sqrt(4*log(exp(1)*p*log(n)/t^2))
    as[t] = a
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
      for (j in 1:p) {
        if(abs(x[j])>a){
        
          x[j] = x[j]^2 - nu_a -2*a^2
        
        }
        else{
          x[j]=0
        }
      }
    sums[iter,t] = sum(x)
    #x = x[x>0]
    #if(length(x)>t){
    #  maxtots[iter, t] = sum(x[1:t])
    #}else{
    #  maxtots[iter, t] = sum(x)
    #}
      
   
    
  }
}

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








# new simu tail bound
#N = 100
p = 1000
n = 1000
N = n
s = floor(sqrt(p*log(n)))
ss = 1:s
#s = 500
#aa = sqrt(4*log(exp(1)*p*log(n)/s^2))
#nu_a = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))

#thet = seq(0, 1.2, by=0.05)*a
sums = matrix(0, nrow = N, ncol =s)
maxtots = matrix(0, nrow = N, ncol =s)
as = sums[1,]
for (iter in 1:N) {
  y = rnorm(p) 
  #y[1:p] = y[1:p] + 3*as[1]
  
  for (t in 1:s) {
    maxzz = rep(0,p)
    a = sqrt(4*log(exp(1)*p*log(n)/t^2))
    as[t] = a
    #x = y[] + a*1.9
    x = y[]
    nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
    for (j in 1:p) {
      if(abs(x[j])>a){
        
        x[j] = x[j]^2 - nu_a -20*a^2
        
      }
      else{
        x[j]=0
      }
    }
    sums[iter,t] = sum(x) + 20*a^2*p*2*pnorm(a,lower.tail=FALSE)
    #x = x[x>0]
    #if(length(x)>t){
    #  maxtots[iter, t] = sum(x[1:t])
    #}else{
    #  maxtots[iter, t] = sum(x)
    #}
    
    
    
  }
}

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
plot(1:s, rowmaxes,
     ylim = c(min(c(rowmins,9*(- (1:s) * as^2 - log(n)))) ,
              max(c(rowmaxes, 9*(1:s)*as^2 +9*log(n)))) , type="l")
lines(1:s, rowmeans , type="l",col=2)
lines(1:s, rowmins, type="l",col=3)
lines(1:s, 9*(1:s * as^2 + log(n)),col=4,type="l")
lines(1:s, 9*(-(1:s) * as^2-log(n)),col=4,type="l")

lines(1:s, sqrt(s*aa^2*1:s*as^2))
lines(thet/a, rowmins, type="l",col=2)
lines(thet/a, rowmeans, type="l",col=3)
9*a^2












# check if it is actually true that sqrt(slog(n^4))< \xi(s)
n = 10
p = 1000000000
s = floor(sqrt(p*log(n^4)))
ss = 1:s
xis = ss*log(exp(1)*p*log(n^4)/ss^2) + log(n^4)
slog = sqrt(ss*log(n^4))
plot(xis,type="l")
lines(slog)
sum(slog>xis)

