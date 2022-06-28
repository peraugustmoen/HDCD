# checking idea of non-equal magnitudes in changes

n = 3
p = 1000000

maxs = floor(sqrt(p*log(n^4)))

ss = 1:maxs
as = sqrt(log(exp(1)*p*log(n^4)/ss^2))

xis = pmax(log(n^4), ss*as^2)
#xis = log(n^4)+ ss*as^2
len = length(xis)
actdiff = xis[2:len] / xis[1:(len-1)] -1
plot(actdiff,type="l",col=2)
lines(1/ss)

plot(1/ss[1:(len-1)] /(actdiff))

start = min(which(ss*as^2 > log(n^4)))
stop = max(which(xis<sqrt(p*log(n^4))))

plot(1/ss[start:stop] /(actdiff[start:stop]))

max(1/ss[start:stop] /(actdiff[start:stop]))


thet2s = rep(0,stop)
leadingthet = rep(0,stop)
secondleadingthet = rep(0,stop)
fracs = rep(0,stop)
incr = 1e-4
thet2stmp = rep(0,stop)
for (s in 2:stop) {
  thet2s = rep(1,s)/s*xis[s]
  thet2stmp = rep(1,s)/s*xis[s]
  right = TRUE
  while(right){
    thet2s[] = thet2stmp[]
    thet2stmp[1] = (1+incr)*thet2s[1]
    thet2stmp[2:s] = (xis[s] - thet2stmp[1])/(s-1)
    
    right = (which.max(cumsum(thet2stmp)/xis[1:s])==s)
  }
  fracs[s] = thet2s[1]/xis[s]
  leadingthet[s] = thet2s[1]
  
  if(s>2){
    thet2stmp[] = thet2s[]
    
    right = TRUE
    while(right){
      thet2s[] = thet2stmp[]
      thet2stmp[2] = (1+incr)*thet2s[2]
      thet2stmp[3:s] = (xis[s] - thet2s[1] -thet2stmp[2])/(s-2)
      
      right = (which.max(cumsum(thet2stmp)/xis[1:s])==s)
    }
    #fracs[s] = thet2s[1]/xis[s]
    secondleadingthet[s] = thet2s[2]
    
  }
  
  
}
plot(fracs,type="l",col=2)
lines(1/(1:maxs))

plot(fracs*(1:stop))

plot(leadingthet)
plot(secondleadingthet)

#conjecture: maximum uneven \thet = (\thet_1, ..., \thet_s) such that
#(\thet_1^2 + ... + \thet_t^2) / xi(t) 
# is maximized by t=s is
# thet_1^2 = xi(1)
# thet_2^2 = xi(2) - xi(1),
# . 
# . 
# . 

plot(c(0,cumsum(diff(xis[1:stop]))) + xis[1])

plot((diff(xis[1:stop])))
lines(as[2:stop]^2)
