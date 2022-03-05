#### check inspect
library(HDCD)

#multiple change-points:
p = 1000
n =1000
mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
etas = c(round(n/4),round(2*n/4),round(3*n/4))
randomdense =rnorm(p)
randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
randomdense = randomdense/p^(1/4)/sqrt(n)*12
mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + matrix(rep(c(rnorm(10)*2/sqrt(n)*10, rep(0,p-10)), n-round(etas[1])) ,nrow=p)
mus[,round(etas[2]+1):n] = mus[,round(etas[2]+1):n] + matrix(rep(randomdense, n-round(etas[2])) ,nrow=p)
mus[,round(etas[3]+1):n] = mus[,round(etas[3]+1):n] + matrix(rep(rep(1,p)/p^(1/4)/sqrt(n)*12, n-round(etas[3])) ,nrow=p)
plot(mus[1,])
X = mus+noise

# the inspect package cheats slightly, as it lowers lambda if necessary..
# lambda = 4*sqrt(log(p*n)) is the theoretically justified value of lambda
# while lambda = sqrt(log(p*log(n))/2) is what wang and samworth recommend in practice.

xi = 4*sqrt(log(p*n))
lambda = sqrt(log(p*log(n)))
#lambda = 4*sqrt(log(p*n))

res = Inspect(X, xi=xi, alpha = 1.2, K = 10,eps=1e-10,
              lambda = lambda, maxiter=10000,debug=FALSE)
res$changepoints

res2 = HDCD (X, 4, 1, alpha = 1+1/6, K = 7, debug =FALSE)
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res2$changepoints
res2$coordinate
res2$CUSUMval


p = 100
n =500
mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)

etas = c(20,60)
#mus[,2:3] = mus[,2:3] + matrix(rep(c(rnorm(10)*10, rep(0,p-10)), 2) ,nrow=p)
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

res2 = HDCD (X, 8, 36, alpha = 1+1/6, K = 7, debug =TRUE)
res2$changepoints
s=-1
cusum = .Call(CUSUM_R, X , as.integer(s), as.integer(e), as.integer(P), as.integer(n))
