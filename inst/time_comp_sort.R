NN = 20
library(HDCD)
zortz = 1:NN
nonzortz = 1:NN
for (nn in 1:NN) {
  
  
  #multiple change-points:
  p = 10000
  n =100
  mus = matrix(0, nrow=p, ncol=n)
  noise = matrix(rnorm(n*p), nrow=p, ncol=n)
  etas = c(round(n/4),round(2*n/4),round(3*n/4))
  randomdense =rnorm(p)
  randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
  randomdense = randomdense/p^(1/4)/sqrt(n)*12
  k = 1
  mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + sqrt(2*max(log(exp(1)*p*log(n)/k^2), log(n)))*matrix(rep(c(rep(1,k)*1/sqrt(n)*10, rep(0,p-k)), n-round(etas[1])) ,nrow=p)
  mus[,round(etas[2]+1):n] = mus[,round(etas[2]+1):n] + matrix(rep(randomdense, n-round(etas[2])) ,nrow=p)*1.2
  mus[,round(etas[3]+1):n] = mus[,round(etas[3]+1):n] + matrix(rep(rep(1,p)/p^(1/4)/sqrt(n)*12, n-round(etas[3])) ,nrow=p)*1.2
  
  #plot(mus[1,])
  X = mus+noise
  #plot(X[1,])
  
  #X[2,] = 10*X[2,]
  
  #sds = apply(X, MARGIN=1, FUN = function(x) median(abs(diff(x)))*1.05)
  #XX =t(apply(X, MARGIN=1, FUN = function(x) x/ (median(abs(diff(x)))*1.05)))
  
  # the inspect package cheats slightly, as it lowers lambda if necessary..
  # lambda = 4*sqrt(log(p*n)) is the theoretically justified value of lambda
  # while lambda = sqrt(log(p*log(n))/2) is what wang and samworth recommend in practice.
  

  
  a = system.time({res2 = HDCD (X[,], 2,2, alpha = 1+1/6, K = 7, threshold_d_test = 2, 
                            threshold_s_test = 2,debug= FALSE)})
  zortz[nn] = a[3]
  #NOTE: 8, 2 corresponds to EQUAL PENALTIES
  #res2$changepoints
  #res2$s
  #res2$coordinate
  
  a = system.time({res3 = HDCD (X[,], 2,2, alpha = 1+1/6, K = 7, threshold_d_test = 2, 
                            threshold_s_test = 1, droppartialsum = TRUE,debug= FALSE)})
  nonzortz[nn] = a[3]
  #NOTE: 8, 2 corresponds to EQUAL PENALTIES
  #res3$changepoints
  #res3$s
  #res2$coordinate
  
}






lenz = c(1000,10000,100000)
kz = c(2,5,10,100,1000)

insertsort = matrix(NA,nrow = length(lenz), ncol = length(kz))
qsort = matrix(NA,nrow = length(lenz), ncol = length(kz))
for (ll in 1:length(lenz)) {
  for (kk in 1:length(kz)) {
    tmp1  = 1:NN
    tmp2  = 1:NN
    for (nn in 1:NN) {
      x = rnorm(lenz[ll])
      a = system.time({sort_k_largest(x, kz[kk],0, lenz[ll])})
      tmp1[nn] = a[3]
      a = system.time({partial_quicksort(x, kz[kk],lenz[ll])})
      tmp2[nn] = a[3]
    }
    insertsort[ll, kk] = mean(tmp1)
    qsort[ll, kk] = mean(tmp2)
    
  }
}

plot( qsort[3,])
lines( insertsort[3,])

qsort
insertsort
