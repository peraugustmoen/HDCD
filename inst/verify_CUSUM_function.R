library(HDCD)

p = 5000
n =100
mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
etas = c(round(n/4),round(2*n/4),round(3*n/4))
X = noise

start = -1
stop = etas[2]-1
CUSUM_X = CUSUM(X, start, stop)
CUSUM_X = CUSUM_X[1:(p*(stop-start-1))]
CUSUM_X = matrix(data=CUSUM_X, nrow=p,ncol = stop-start-1)
#CUSUM_X = CUSUM_X[, 1:(stop-start-1)]
plot(colSums(CUSUM_X^2))

cusum_tt = InspectChangepoint::cusum.transform(X[,1:etas[2]])
plot(colSums(cusum_tt^2))
points(colSums(CUSUM_X^2), col=2)
