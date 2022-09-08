p = 1000
n = 500
N = 1000

ss = 1:(floor(sqrt(p*log(n^4))))
res = matrix(NA, nrow = max(ss)+1, ncol=N)
as = rep(NA, length(ss)+1)
nu_as = rep(NA, length(ss)+1)
last = max(ss)+1



as = sqrt(2*log(exp(1)*p*log(n^4)/ss^2))
nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))


thet = 2
expected = (as + thet)/(exp((as - thet)^2/2)* sqrt(2*pi)) + (1 - nu_as + thet^2)*pnorm((as - thet)/sqrt(2),lower.tail = FALSE)/2

expected = expected + (as + 3*thet)/(exp((as + thet)^2/2)*sqrt(2*pi)) + (1 - nu_as + thet^2)* pnorm((as + thet)/sqrt(2), lower.tail=FALSE)/2
#(a + 3 t)/(E^((a + t)^2/2) Sqrt[2 Pi]) + ((1 - n + t^2) Erfc[(a + t)/Sqrt[2]])/2

plot(as)
abline(h=thet)
plot(expected)



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