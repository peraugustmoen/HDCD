library(HDCD)
library(doSNOW)
library(foreach)

n = 100
p = 10000

N = 5000
Nsim = 1000
tol = 1/2000

gridlen = 50
grid = seq(0,3, length.out = gridlen)

spars = c(1,10,30,100,500,1000)

tmp1 = rep(NA, Nsim) # without partial
tmp2 = rep(NA, Nsim) # with partial
tmp3 = rep(NA, Nsim) # pilliat
tmp4 = rep(NA, Nsim) # without partial, log(n) instead of log(n^4)
tmp5 = rep(NA, Nsim) # without partial, log(n) instead of log(n^4)

res = array(NA, dim= c(length(spars),5, gridlen))

thr = changepoint_test_HDCD_calibrate (n,p, as=NULL, nu_as=NULL, ts=NULL, twologn=NULL, N, tol,debug=TRUE)
#thr

thr2 = changepoint_test_Pilliat_calibrate(n,p,N = N, tol = tol ,debug=FALSE)
#thr2

thr3 = changepoint_test_HDCD_modified_calibrate(n,p, as=NULL, nu_as=NULL, ts=NULL, twologn=NULL, N, tol,debug=TRUE)

num_cores = 6

cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
result = foreach(i = 1:length(spars),.options.snow = opts) %dopar% {
  #for (i in 1:length(spars)) {
  library(HDCD)
  res = array(NA, dim = c(5,gridlen))
  t = spars[i]
  mus = matrix(0, nrow=p, ncol=n)
  #noise = matrix(rnorm(n*p), nrow=p, ncol=n)
  etas = c(round(n/4))
  #randomdense =rnorm(p)
  #randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
  #randomdense = randomdense/p^(1/4)/sqrt(n)*12
  k = t
  
  if(k<2*sqrt(p*log(n))){
    mus[,round(etas[1]+1):n] =  sqrt(max(log(exp(1)*p*log(n)/k^2), log(n)))*matrix(rep(c(rep(1,k)*1/sqrt(etas), rep(0,p-k)), n-round(etas[1])) ,nrow=p)
  }else{
    #randomdense =rnorm(p)
    #randomdense = randomdense/norm(randomdense,type="2")*p^(1/4)/sqrt(etas)
    mus[,round(etas[1]+1):n] =  matrix(rep(rep(1,p)/p^(1/4)/sqrt(etas), n-round(etas[1])) ,nrow=p)
  }
  
  for (j in 1:length(grid)) {
    gridconst = grid[j]
    mumu = mus * gridconst
    
    for (z in 1:Nsim) {
      X = matrix(rnorm(n*p), nrow=p, ncol=n) + mumu
      
      tmp1[z] = changepoint_test_HDCD(X, empirical =TRUE,thresholds=thr$without_partial,droppartialsum = TRUE,
                                      debug=FALSE)
      tmp2[z] = changepoint_test_HDCD(X, empirical =TRUE,thresholds=thr$with_partial,debug=FALSE)
      
      tmp3[z] = changepoint_test_Pilliat(X, empirical =TRUE,  thresholds_partial=thr2$thresholds_partial,
                                         threshold_dense = thr2$threshold_dense, thresholds_bj = thr2$thresholds_bj,
                                         debug =FALSE)
      
      tmp4[z] = changepoint_test_HDCD_modified(X, empirical =TRUE,thresholds=thr3$without_partial,droppartialsum = TRUE,
                                      debug=FALSE)
      tmp5[z] = changepoint_test_HDCD_modified(X, empirical =TRUE,thresholds=thr3$with_partial,debug=FALSE)
    }
    res[1,j] = mean(tmp1)
    res[2,j] = mean(tmp2)
    res[3,j] = mean(tmp3)
    res[4,j] = mean(tmp4)
    res[5,j] = mean(tmp5)
    
  }
  res
  #X = mus+noise
  
}
close(pb)
stopCluster(cl) 
for (i in 1:length(spars)) {
  res[i,,] = result[[i]]
}

plot(grid, res[4,1,],type="l") # without partial, logn4
lines(grid,res[4,2,],type="l",col=2) # with partial, logn4
lines(grid,res[4,3,],type="l", col=3) # pilliat
lines(grid,res[4,4,],type="l", col=4) # without partial, logn
lines(grid,res[4,5,],type="l", col=5) # with partial, logn





