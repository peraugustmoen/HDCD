p = 1000
n = 500
N = 1000

ss = 1:(floor(sqrt(p*log(n^4))))
res = matrix(NA, nrow = max(ss)+1, ncol=N)
as = rep(NA, length(ss)+1)
nu_as = rep(NA, length(ss)+1)
last = max(ss)+1
for (j in 1:N) {
  xo = rnorm(p)
  
  for (s in 1:(floor(sqrt(p*log(n^4)))+1)) {
    x = xo[]
    a = 0
    nu_a = 1
    if(s !=last){
      a = 2*sqrt(log(exp(1)*p*log(n^4)/s^2))
      nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
    }
    
    as[s] = a
    nu_as[s] = nu_a
    thresholded[,] = 0
    
    
    x[abs(x)<=a] = 0
    res[s,j] = sum(x[abs(x)>0]^2 - nu_a)
    
  }  
}
hist(res[1,])
maxes = apply(res, 1, max)
plot(maxes)
#lines(ss*log(exp(1)*p*log(n^4)/ss^2))
lines(sqrt(p*exp(-as^2/2)*log(n^4)) + log(n^4))
plot(sqrt(p*exp(-as^2/2)*log(n^4)) + log(n^4), type="l")

diffs = res[2:(max(ss)+1), ] - res[1:max(ss),]
maxdiffs = apply(diffs, 1, max)
mindiffs = apply(diffs, 1, min)
meandiffs = apply(diffs, 1, mean)
plot(maxdiffs)
plot(maxdiffs[1:(length(maxdiffs)-1)])
plot(mindiffs[1:(length(mindiffs)-1)])
plot(meandiffs)
## simulation of differences

