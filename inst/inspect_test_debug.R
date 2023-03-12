library(HDCD)

#multiple change-points:
p = 50
n =100
mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
eta = round(0.25*n)
#eta = round(sample(1:floor(n/2),1))

k = p
sparse_const = 0
dense_const = 0
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

tt1 = single_Inspect(rescale.variance(X))
tt1



N = 1000
tol = 1/N*10


set.seed(300*32+1)
cc4 = Inspect_test_calibrate(n, p, N=N, tol=tol,lambda = sqrt(log(p*log(n))/2), 
                             maxiter=10000,rescale_variance = TRUE,debug =FALSE)

cc4
set.seed(300*32+1)
mm = -10000
for (i in 1:N) {
   # X = matrix(NA, nrow = p, ncol = n)
   # 
   # for(jj in 1:n){
   #   for(ii in 1:p){
   #     X[ii,jj] = .Call(stats:::C_rnorm, 1,0.0,1.0)
   #   }
   # }
  
  X = matrix(rnorm(n*p), nrow = p, ncol=n,byrow = FALSE)
  #print(X[p,n])
  #print(rescale.variance(X)[p,n])
  #rr = single_Inspect(rescale.variance(X))
  rr = single_Inspect(matrix(rescale_variance(X)$X, byrow= FALSE, ncol = n, nrow=p),debug=FALSE)
  #print(rr$cusumval)
  if(rr$cusumval>mm){
    mm = rr$cusumval
  }
}
mm
res8 = as.integer(single_Inspect(X)$cusumval>cc4[[1]])
res8
