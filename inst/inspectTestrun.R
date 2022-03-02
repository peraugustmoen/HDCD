#### check inspect

#multiple change-points:
p = 100
n =600
mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)

etas = c(20,60)
mus[,21:60] = mus[,21:60] + matrix(rep(c(rnorm(10)*10, rep(0,p-10)), 40) ,nrow=p)
mus[,41:60] = mus[,41:60] + matrix(rep(rnorm(p), 20) ,nrow=p)
plot(mus[1,])
X = mus+noise

# the inspect package cheats slightly, as it lowers lambda if necessary.. 
# lambda = 4*sqrt(log(p*n)) is the theoretically justified value of lambda
# while lambda = sqrt(log(p*log(n))/2) is what wang and samworth recommend in practice. 

xi = 4*sqrt(log(p*n))
lambda = sqrt(log(p*log(n))/2)
#lambda = 4*sqrt(log(p*n))

res = Inspect(X, xi=xi, alpha = 1.2, K = 10,eps=1e-10,
                lambda = lambda, maxiter=10000,debug=FALSE)
res$changepoints

library(InspectChangepoint)
res2= InspectChangepoint::inspect(X,lambda=lambda,threshold=xi,M=200)
res2 
