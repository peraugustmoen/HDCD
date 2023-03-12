library(HDCD)

n = 100
p = 1000
X = matrix(rnorm(n*p), ncol = n, nrow = p)
Xold = X[,]
true = InspectChangepoint::rescale.variance(X)

Xold - X
myimp = rescale.variance(X)

sum(abs(myimp - true))

truescales = rescale.variance.scales(X)


myimp2 = rescale_variance(X)

sum(abs(true - myimp2$X))
sum(abs(myimp2$scales - truescales))

plot(myimp2$scales -truescales)

     