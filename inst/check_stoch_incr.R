

a = 2
#as = c(2)
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))

thets = seq(0,3,by=0.1)

N = 1000000

cdfs = list()
lambdas = seq(0, 0.5, by=0.001)
mgfs = matrix(data=NA, nrow = N, ncol = length(lambdas))
for (j in 1:length(thets)) {
  set.seed(100)
  z = rnorm(N) + thets[j]
  
  z[abs(z)<a] = 0
  z[z!=0] = z[z!=0]^2 - nu_a
  cdfs[[j]] = ecdf(z)
  for (h in 1:length(lambdas)) {
    mgfs[j,h] = mean( exp(-lambdas[h] * z))
  }
  
}

xs = seq(a^2-nu_a, 15,by=0.001)

plot(xs,cdfs[[1]](xs),type="l")
lines(xs,cdfs[[22]](xs),col=2)

plot(lambdas,mgfs[1,],type="l")
lines(lambdas,mgfs[2,],col=2)





a = 2
#as = c(2)
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))

thets = seq(0,3,by=0.1)

N = 100000

cdfs = list()
lambdas = seq(0, 0.5, by=0.001)
mgfs = matrix(data=NA, nrow = N, ncol = length(lambdas))
for (j in 1:length(thets)) {
  set.seed(100)
  z = rnorm(N) + thets[j]
  y = z^2-1
  z[abs(z)<a] = 0
  z[z!=0] = z[z!=0]^2 - nu_a
  cdfs[[j]] = ecdf(z-y)
  
}

xs = seq(-nu_a, 15,by=0.001)

t = 1
sum(cdfs[[t]](xs) > cdfs[[t+1]](xs) )
which(cdfs[[t]](xs) > cdfs[[t+1]](xs))
plot(xs,cdfs[[1]](xs),type="l")
lines(xs,cdfs[[30]](xs),col=2)

plot(lambdas,mgfs[1,],type="l")
lines(lambdas,mgfs[2,],col=2)


a = 2
#as = c(2)
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))

b = 0.1
#as = c(2)
nu_b = 1 + b*exp(dnorm(b, log=TRUE)-pnorm(b, lower.tail = FALSE, log=TRUE))


xx = seq(0,10,length.out=1000)

ff = function(x){
  aa = rep(0, length(x))
  bb = rep(0, length(x))
  
  aa[abs(x)>=a] = (x[abs(x)>=a]^2 - nu_a)
  bb[abs(x)>=b] = (x[abs(x)>=b]^2 - nu_b)
  
  return(aa - bb)
}
plot(xx, -ff(xx))

