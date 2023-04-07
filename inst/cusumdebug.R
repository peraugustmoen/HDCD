
n = 2
p = 2000
X = matrix(rnorm(n*p), ncol = n)

c1 = CUSUM(X, -1, n-1)
c1 = matrix(c1,nrow = p, ncol = n-1)

k = 1
c2 = single_CUSUM(X,-1,n-1, k-1)
c2 = matrix(c2, nrow = p)


sum(abs(c1[,k] - c2))
 
c3 = InspectChangepoint::cusum.transform(X)
c3 = matrix(c3, nrow=p)

min(sum(abs(c3 - c1)), sum(abs(c3 +c1)))
min(sum(abs(c3[,k] - c2)), sum(abs(c3[,k] +c2)))
