library(HDCD)

#multiple change-points:
p = 10000
n =50

N = 1000
tol = 1/500
thr = changepoint_test_HDCD_calibrate (n,p, as=NULL, nu_as=NULL, ts=NULL, twologn=NULL, N, tol,debug=TRUE)
thr

thr2 = changepoint_test_Pilliat_calibrate(n,p,N = N, tol = tol ,debug=FALSE)
thr2

mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
etas = c(round(n/4),round(2*n/4),round(3*n/4))
randomdense =rnorm(p)
randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
randomdense = randomdense/p^(1/4)/sqrt(n)*12
k = 1
if(k<2*sqrt(p*log(n))){
  mus[,round(etas[1]+1):n] =  sqrt(max(log(exp(1)*p*log(n)/k^2), log(n)))*matrix(rep(c(rep(1,k)*1/sqrt(etas[1]), rep(0,p-k)), n-round(etas[1])) ,nrow=p)
}else{
  #randomdense =rnorm(p)
  #randomdense = randomdense/norm(randomdense,type="2")*p^(1/4)/sqrt(etas)
  mus[,round(etas[1]+1):n] =  matrix(rep(rep(1,p)/p^(1/4)/sqrt(etas[1]), n-round(etas[1])) ,nrow=p)
}
mus = mus*3
#mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + 0.3*sqrt(2*max(log(exp(1)*p*log(n)/k^2), log(n)))*matrix(rep(c(rep(1,k)*1/sqrt(n)*10, rep(0,p-k)), n-round(etas[1])) ,nrow=p)
#mus[,round(etas[2]+1):n] = mus[,round(etas[2]+1):n] + matrix(rep(randomdense, n-round(etas[2])) ,nrow=p)*1.2
#mus[,round(etas[3]+1):n] = mus[,round(etas[3]+1):n] + matrix(rep(rep(1,p)/p^(1/4)/sqrt(n)*12, n-round(etas[3])) ,nrow=p)*1.2

#plot(mus[1,])
X = mus+noise
#plot(X[1,])


changepoint_test_HDCD(X, empirical =TRUE,thresholds=thr$with_partial,debug=TRUE)

changepoint_test_HDCD(X, empirical =TRUE,thresholds=thr$without_partial,droppartialsum = TRUE,
                      debug=TRUE)


changepoint_test_Pilliat(X, empirical =TRUE,  thresholds_partial=thr2$thresholds_partial,
                         threshold_dense = thr2$threshold_dense, thresholds_bj = thr2$thresholds_bj,
                                           debug =TRUE)




cc = InspectChangepoint::cusum.transform(X)[,18]
sum(cc^2) -p
cc_s = rev(sort(cc^2))
cc_s[1]







# 
# 
# 
# 
# cc = InspectChangepoint::cusum.transform(X)[,18]
# sum(cc^2) -p
# cc_s = rev(sort(cc^2))
# cc_s[1]
# 
# 
# # Pilliat
# 
# 
# 
# #multiple change-points:
# # p = 1000
# # n =100
# # 
# # N = 2000
# # tol = 1/n
# 
# thr = changepoint_test_Pilliat_calibrate(n,p,N = N, tol = tol ,debug=FALSE)
# thr
# 
# 
# 
# mus = matrix(0, nrow=p, ncol=n)
# noise = matrix(rnorm(n*p), nrow=p, ncol=n)
# etas = c(round(n/4),round(2*n/4),round(3*n/4))
# randomdense =rnorm(p)
# randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
# randomdense = randomdense/p^(1/4)/sqrt(n)*12
# k = 1
# mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + 0.5*sqrt(2*max(log(exp(1)*p*log(n)/k^2), log(n)))*matrix(rep(c(rep(1,k)*1/sqrt(n)*10, rep(0,p-k)), n-round(etas[1])) ,nrow=p)
# #mus[,round(etas[2]+1):n] = mus[,round(etas[2]+1):n] + matrix(rep(randomdense, n-round(etas[2])) ,nrow=p)*1.2
# #mus[,round(etas[3]+1):n] = mus[,round(etas[3]+1):n] + matrix(rep(rep(1,p)/p^(1/4)/sqrt(n)*12, n-round(etas[3])) ,nrow=p)*1.2
# 
# #plot(mus[1,])
# X = mus+noise
# #plot(X[1,])
# 
# 
# changepoint_test_HDCD(X, empirical =TRUE,thresholds=thr$with_partial,debug=TRUE)
# 
# 
# cc = InspectChangepoint::cusum.transform(X)[,18]
# sum(cc^2) -p
# cc_s = rev(sort(cc^2))
# cc_s[1]
# 
#   
# 
# 

