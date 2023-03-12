library(HDCD)

#multiple change-points:
p = 1000
n =200
mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
eta = round(0.25*n)
#eta = round(sample(1:floor(n/2),1))

k = p
sparse_const = 2.5
dense_const = 2.5
diff = 0
if(k<sqrt(p*log(n))){
  rootnorm = sparse_const/sqrt(eta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
  diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
}else{
  rootnorm = dense_const/sqrt(eta)*(p*log(n))^(1/4)
  diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
}

mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)
X = mus+noise

rescale_variance = TRUE


res1 = HDCD_test(X,droppartialsum = TRUE,fast=TRUE,rescale_variance = rescale_variance)
res1

res2 = HDCD_test(X,droppartialsum = TRUE,fast=FALSE,rescale_variance = rescale_variance)
res2

N = 1000
tol = 1/N
cc1 = HDCD_test_calibrate(n, p, N, tol, fast = TRUE, rescale_variance = rescale_variance )
cc2 = HDCD_test_calibrate(n, p, N, tol, fast = FALSE, rescale_variance = rescale_variance )


res3 = HDCD_test(X,droppartialsum = FALSE,fast=TRUE,rescale_variance = rescale_variance,
                 thresholds = cc1[[2]])
res3
res4 = HDCD_test(X,droppartialsum = FALSE,fast=FALSE,rescale_variance = rescale_variance,
                 thresholds = cc2[[2]])
res4



#checking pilliat

res5 = changepoint_test_Pilliat(X, rescale_variance = TRUE,debug =FALSE)
res5

cc3 = changepoint_test_Pilliat_calibrate(n,p, N=N, tol=tol,
                                                    rescale_variance = TRUE,debug=FALSE)

res6 = changepoint_test_Pilliat(X, rescale_variance = TRUE, empirical=TRUE, thresholds_partial = cc3[[1]], 
                                threshold_dense = cc3[[2]], thresholds_bj = cc3[[3]])
res6


#checking inspect


res7 = as.integer(single_Inspect(X)$cusumval>sqrt(log(p*log(n))/2))
res7
set.seed(122)
cc4 = Inspect_test_calibrate(n, p, N=100, tol,lambda = sqrt(log(p*log(n))/2), 
                                        maxiter=10000,rescale_variance = TRUE,debug =TRUE)

cc4
set.seed(122)
mm = -10000
for (i in 1:100) {
  # X = matrix(NA, nrow = p, ncol = n)

  # for(ii in 1:p){
  #   for(jj in 1:n){
  #     X[ii,jj] = .Call(stats:::C_rnorm, 1,0.0,1.0)
  #   }
  # }
  
  X = matrix(rnorm(n*p), nrow = p, ncol=n,byrow = FALSE)
  #rr = single_Inspect(rescale.variance(X))
  rr = single_Inspect(rescale.variance(X),debug=TRUE)
  #print(rr$cusumval)
  if(rr$cusumval>mm){
    mm = rr$cusumval
  }
}
mm
res8 = as.integer(single_Inspect(X)$cusumval>cc4[[1]])
res8
