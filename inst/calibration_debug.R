library(HDCD)

#multiple change-points:
p = 50
n =50
mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
etas = c(round(n/4),round(2*n/4),round(3*n/4))
randomdense =rnorm(p)
randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
randomdense = randomdense/p^(1/4)/sqrt(n)*12
k = 8
#mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + 0.5*sqrt(2*max(log(exp(1)*p*log(n)/k^2), log(n)))*matrix(rep(c(rep(1,k)*1/sqrt(n)*10, rep(0,p-k)), n-round(etas[1])) ,nrow=p)
#mus[,round(etas[2]+1):n] = mus[,round(etas[2]+1):n] + matrix(rep(randomdense, n-round(etas[2])) ,nrow=p)*1.2
#mus[,round(etas[3]+1):n] = mus[,round(etas[3]+1):n] + matrix(rep(rep(1,p)/p^(1/4)/sqrt(n)*12, n-round(etas[3])) ,nrow=p)*1.2

plot(mus[1,])
X = mus+noise
plot(X[1,])

#X[2,] = 10*X[2,]

#sds = apply(X, MARGIN=1, FUN = function(x) median(abs(diff(x)))*1.05)
#XX =t(apply(X, MARGIN=1, FUN = function(x) x/ (median(abs(diff(x)))*1.05)))

# the inspect package cheats slightly, as it lowers lambda if necessary..
# lambda = 4*sqrt(log(p*n)) is the theoretically justified value of lambda
# while lambda = sqrt(log(p*log(n))/2) is what wang and samworth recommend in practice.


cc3 = changepoint_test_Pilliat_calibrate(n,p, N=10000, tol=1/100, debug=FALSE)
cc3

res12 = changepoint_test_Pilliat(X)

changepoint_test_Pilliat(X, empirical =TRUE,  thresholds_partial=cc3$thresholds_partial,
                         threshold_dense = cc3$threshold_dense, thresholds_bj = cc3$thresholds_bj,
                         debug =TRUE)









system.time({res11= Pilliat(X, 
                            K = 2, alpha = 1+1/6, empirical = TRUE, threshold_dense = cc3$threshold_dense, 
                            thresholds_partial = cc3$thresholds_partial, thresholds_bj = cc3$thresholds_bj,debug =FALSE)})


res11






# checking the modified (logn instead of logn4)


system.time({res2 = HDCD_modified (X[,], 4,1, alpha = 1+1/6, K = 7, threshold_d_test = 4, 
                                   threshold_s_test = 2, debug= FALSE)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res2$changepoints
res2$s
#res2$coordinate
ccc = HDCD_modified_calibrate(n,p, N=1000, tol=0.001,debug=FALSE)

system.time({res3 = HDCD_modified (X[,], 3,2, alpha = 1+1/6, K = 7, threshold_d_test = 3, 
                                   threshold_s_test = 1, droppartialsum = TRUE,debug= FALSE)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res3$changepoints
res3$s
#res2$coordinate

system.time({res4 = HDCD_modified (X[,], 3,2, alpha = 1+1/6, K = 2, threshold_d_test = 3, 
                                   threshold_s_test = 1, droppartialsum = FALSE, fast =TRUE,debug= FALSE)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res4$changepoints
res4$s
res4$startpoints
res4$endpoints
#res2$coordinate


# system.time({res3 = Scan (X[,],alpha = 1+1/6, K = 7,debug=FALSE)})
# #NOTE: 8, 2 corresponds to EQUAL PENALTIES
# res3$changepoints
# #res3$coordinate





rez = single_HDCD(X[,1:etas[2]], 4, 4)
rez





# 
# 
# res3 = HDCD (X, 4,sqrt(2), alpha = 1+1/6, K = 7, debug =FALSE)
# #NOTE: 8, 2 corresponds to EQUAL PENALTIES
# res3$changepoints
# 
# res2$s
# cbind(res2$startpoints, res2$endpoints)
# res2$thresholds
# 
# res2$coordinate
# res2$CUSUMval
# 
# 
# InspectChangepoint::inspect(X, lambda = lambda, threshold=xi, M = 1000)
# 
# p = 100
# n =500
# mus = matrix(0, nrow=p, ncol=n)
# noise = matrix(rnorm(n*p), nrow=p, ncol=n)
# 
# etas = c(20,60)
# #mus[,2:3] = mus[,2:3] + matrix(rep(c(rnorm(10)*10, rep(0,p-10)), 2) ,nrow=p)
# plot(mus[1,])
# X = mus+noise
# 
# # the inspect package cheats slightly, as it lowers lambda if necessary.. 
# # lambda = 4*sqrt(log(p*n)) is the theoretically justified value of lambda
# # while lambda = sqrt(log(p*log(n))/2) is what wang and samworth recommend in practice. 
# 
# xi = 4*sqrt(log(p*n))
# lambda = sqrt(log(p*log(n))/2)
# #lambda = 4*sqrt(log(p*n))
# 
# res = Inspect(X, xi=xi, alpha = 1.2, K = 10,eps=1e-10,
#               lambda = lambda, maxiter=10000,debug=FALSE)
# res$changepoints
# 
# res2 = HDCD (X, 8, 36, alpha = 1+1/6, K = 7, debug =TRUE)
# res2$changepoints
# s=-1
# cusum = .Call(CUSUM_R, X , as.integer(s), as.integer(e), as.integer(P), as.integer(n))
