library(HDCD)

#multiple change-points:
p = 100
n =40
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
mus = mus/1.3
plot(mus[1,])
X = mus+noise
plot(X[1,])

rescale_variance = FALSE

res2 = ESAC (X[,], 1.5,1, alpha = 1+1/6, K = 7, threshold_d_test = 2, fast =TRUE,droppartialsum = TRUE,
             threshold_s_test = 1,rescale_variance = rescale_variance, debug= TRUE)
res2$changepoints


# res2 = ESAC (X[,], 1.5,1, alpha = 1+1/6, K = 7, threshold_d_test = 2, fast =TRUE,droppartialsum = TRUE,
#              threshold_s_test = 1,rescale_variance = rescale_variance, NOT = FALSE,debug= FALSE)
# res2$changepoints

res3 = ESAC (X[,], 1.5,1, alpha = 1.2, K = 4, threshold_d_test = 2, fast =FALSE,droppartialsum =TRUE,
             threshold_s_test = 1,rescale_variance = rescale_variance, debug= TRUE)
res3$changepoints

res3 = ESAC (X[,], 1.5,1, alpha = 1.2, K = 4, threshold_d_test = 2, fast =FALSE,droppartialsum =TRUE,
             threshold_s_test = 1,rescale_variance = rescale_variance,trim=FALSE, midpoint=TRUE,debug= FALSE)
res3$changepoints
#X[2,] = 10*X[2,]

#sds = apply(X, MARGIN=1, FUN = function(x) median(abs(diff(x)))*1.05)
#XX =t(apply(X, MARGIN=1, FUN = function(x) x/ (median(abs(diff(x)))*1.05)))

# the inspect package cheats slightly, as it lowers lambda if necessary..
# lambda = 4*sqrt(log(p*n)) is the theoretically justified value of lambda
# while lambda = sqrt(log(p*log(n))/2) is what wang and samworth recommend in practice.

xi = 4*sqrt(log(p*n))
lambda = 2*sqrt(log(p*log(n)))
#lambda = xi
#lambda = 4*sqrt(log(p*n))

system.time({res = Inspect(X[,], xi=xi, alpha = 1+1/6, K = 7,eps=1e-10,
              lambda = lambda, maxiter=10000,rescale_variance = rescale_variance, debug=FALSE)})
res$changepoints


system.time({res2 = ESAC (X[,], 2,1, alpha = 1+1/6, K = 7, threshold_d_test = 2, 
                          threshold_s_test = 1,rescale_variance = rescale_variance,debug= FALSE)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res2$changepoints
res2$s
#res2$coordinate

system.time({res3 = ESAC (X[,], 2,1, alpha = 1+1/6, K = 7, threshold_d_test = 2, 
                          threshold_s_test = 1, droppartialsum = TRUE,debug= FALSE)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res3$changepoints
res3$s
#res2$coordinate

system.time({res4 = ESAC (X[,], 1.5,1, alpha = 2, K = 4, threshold_d_test = 2, 
                          threshold_s_test = 1, droppartialsum = FALSE, fast =TRUE,debug= FALSE)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res4$changepoints
res4$s
res4$startpoints
res4$endpoints
#res2$coordinate
# 
# N = 500
# tol = 1/n
# thr = changepoint_test_ESAC_calibrate (n,p, as=NULL, nu_as=NULL, ts=NULL, twologn=NULL, N, tol,debug=TRUE)
# thr
# system.time({res5 = ESAC (X[,], 2,2, empirical=TRUE,alpha = 1+1/6, K = 2,thresholds_test = 2*( thr$with_partial), threshold_d = 2, 
#                           threshold_s = 2, droppartialsum = FALSE, fast =TRUE,debug= FALSE)})
# #NOTE: 8, 2 corresponds to EQUAL PENALTIES
# res5$changepoints
# res5$s
# res5$startpoints
# res5$endpoints
# #res2$coordinate


system.time({res6= Pilliat(X, threshold_d_const=1.5, threshold_bj_const=6, threshold_partial_const=2,
                            K = 2, alpha = 1.5, test_all=TRUE,debug =TRUE)})


res6$changepoints

# checking empirical: 
cc = ESAC_calibrate(n,p, N=1000, tol=0.001,alpha = 2, K = 4, fast=TRUE, bonferroni = TRUE,
                    debug=FALSE, rescale_variance = rescale_variance)

system.time({res7 = ESAC (X[,], 2,1, empirical=TRUE,alpha = 2, K = 4, thresholds_test = cc[[2]], droppartialsum = FALSE, fast =TRUE,
                          rescale_variance = rescale_variance, debug= TRUE)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res7$changepoints
res7$s

system.time({res8 = ESAC (X[,], 2,1, empirical=TRUE,alpha = 2, K = 4, thresholds_test = cc[[1]], droppartialsum = TRUE, fast =TRUE,
                          rescale_variance = rescale_variance, debug= FALSE)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res8$changepoints
res8$s

cc2 = ESAC_calibrate(n,p, N=1000, tol=0.001,K=4,alpha = 1+1/5,fast=FALSE,
                     rescale_variance = rescale_variance, debug=FALSE)

system.time({res9 = ESAC (X[,], 2,1, empirical=TRUE,alpha = 1+1/5, K = 4, 
                          rescale_variance = rescale_variance, thresholds_test = cc2[[2]], droppartialsum = FALSE, fast =FALSE,debug= FALSE)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res9$changepoints
res9$s

system.time({res10 = ESAC (X[,], 1.5,1, empirical=TRUE,alpha = 1+1/6, K = 2, thresholds_test = cc2[[1]], droppartialsum = TRUE, fast =FALSE,
                           rescale_variance = rescale_variance, debug= FALSE)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res10$changepoints
res10$s

system.time({res101 = ESAC (X[,], 1.5,1, empirical=TRUE,alpha = 1+1/6, K = 2, thresholds_test = cc2[[1]], droppartialsum = TRUE, fast =FALSE,
                           rescale_variance = rescale_variance, debug= FALSE, cutoff = 0.2)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res101$changepoints
res101$s

set.seed(101)
cc3 = Pilliat_calibrate(n,p, N=1000, tol=0.001,K = 2, alpha = 2, bonferroni =TRUE,lens = NULL, rescale_variance =rescale_variance,
                        test_all = FALSE,debug=FALSE)
set.seed(101)
cc5 = Pilliat_calibrate(n,p, N=1000, tol=0.001,K = 2, alpha = 2, bonferroni = TRUE,lens = NULL, rescale_variance =rescale_variance,
                        test_all = TRUE,debug=FALSE)

# #multiple change-points:
# p = 100
# n =100
# mus = matrix(0, nrow=p, ncol=n)
# noise = matrix(rnorm(n*p), nrow=p, ncol=n)
# etas = c(round(n/4),round(2*n/4),round(3*n/4))
# randomdense =rnorm(p)
# randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
# randomdense = randomdense/p^(1/4)/sqrt(n)*12
# k = 1
# mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + sqrt(2*max(log(exp(1)*p*log(n)/k^2), log(n)))*matrix(rep(c(rep(1,k)*1/sqrt(n)*10, rep(0,p-k)), n-round(etas[1])) ,nrow=p)
# mus[,round(etas[2]+1):n] = mus[,round(etas[2]+1):n] + matrix(rep(randomdense, n-round(etas[2])) ,nrow=p)*1.2
# mus[,round(etas[3]+1):n] = mus[,round(etas[3]+1):n] + matrix(rep(rep(1,p)/p^(1/4)/sqrt(n)*12, n-round(etas[3])) ,nrow=p)*1.2
# mus = mus/1.2
# plot(mus[1,])
# X = mus+noise
# plot(X[1,])

system.time({res11= Pilliat(X, rescale_variance = rescale_variance,
                           K = 2, alpha =2, empirical = TRUE, threshold_dense = cc3$threshold_dense, 
                           thresholds_partial = cc3$thresholds_partial, thresholds_bj = cc3$thresholds_bj,debug =TRUE)})


res11$changepoints

system.time({res111= Pilliat(X, rescale_variance = rescale_variance,
                            K = 2, alpha = 2, empirical = TRUE, test_all = TRUE,threshold_dense = cc5$threshold_dense, 
                            thresholds_partial = cc5$thresholds_partial, thresholds_bj = cc5$thresholds_bj,debug =TRUE)})


res111$changepoints




cc4 = Inspect_calibrate(n, p, N=1000, tol=0.001,lambda = sqrt(log(p*log(n))/2) , alpha = 2, K = 4,eps=1e-10,
                                   maxiter=10000,debug =FALSE)


system.time({res12 = Inspect(X[,], xi=cc4[[1]], alpha = 2, K = 4,eps=1e-10,
                           lambda = sqrt(log(p*log(n))/2), maxiter=10000,debug=FALSE)})
res12$changepoints





source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")


Ncal = 1000
empirical_penalties = rep(NA, Ncal)
nind = 1
pind = 1
ps = c(p)
ns = c(n)
if(rescale_variance){
  for (i in 1:Ncal) {
    mynulldata <- matrix(rnorm(ns[nind]*ps[pind],0,1),nrow=ps[pind],ncol=ns[nind],byrow=FALSE) # 5 variates with 1000 time points, no change
    empirical_penalties[i] = wbs_penaltyfinder(InspectChangepoint::rescale.variance(mynulldata), SUBSET.normal_penalty, 100)
  }
}else{
  for (i in 1:Ncal) {
    mynulldata <- matrix(rnorm(ns[nind]*ps[pind],0,1),nrow=ps[pind],ncol=ns[nind],byrow=FALSE) # 5 variates with 1000 time points, no change
    empirical_penalties[i] = wbs_penaltyfinder(mynulldata, SUBSET.normal_penalty, 100)
  }
}

empirical_penalties = sort(empirical_penalties,decreasing = TRUE)
pen = empirical_penalties[1]

if(rescale_variance){
  res20 = change_main(InspectChangepoint::rescale.variance(X), SUBSET.normal, 100,penalties= pen)
}else{
  res20 = change_main(X, SUBSET.normal,100, penalties=pen)
}

res20[[2]]




















# checking the modified (logn instead of logn4)


system.time({res2 = ESAC_modified (X[,], 4,1, alpha = 1+1/6, K = 7, threshold_d_test = 4, 
                          threshold_s_test = 2, debug= FALSE)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res2$changepoints
res2$s
#res2$coordinate
ccc = ESAC_modified_calibrate(n,p, N=1000, tol=0.001,debug=FALSE)

system.time({res3 = ESAC_modified (X[,], 3,2, alpha = 1+1/6, K = 7, threshold_d_test = 3, 
                          threshold_s_test = 1, droppartialsum = TRUE,debug= FALSE)})
#NOTE: 8, 2 corresponds to EQUAL PENALTIES
res3$changepoints
res3$s
#res2$coordinate

system.time({res4 = ESAC_modified (X[,], 3,2, alpha = 1+1/6, K = 2, threshold_d_test = 3, 
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





rez = single_ESAC(X[,1:etas[2]], 4, 4)
rez





# 
# 
# res3 = ESAC (X, 4,sqrt(2), alpha = 1+1/6, K = 7, debug =FALSE)
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
# res2 = ESAC (X, 8, 36, alpha = 1+1/6, K = 7, debug =TRUE)
# res2$changepoints
# s=-1
# cusum = .Call(CUSUM_R, X , as.integer(s), as.integer(e), as.integer(P), as.integer(n))










# checking single SBS:
p = 500
n =100
mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
#etas = c(round(n/4),round(2*n/4),round(3*n/4))
etas = c(round(n/2))
randomdense =rnorm(p)
randomdense = randomdense/norm(randomdense,type="2")*sqrt(p)
randomdense = randomdense/p^(1/4)/sqrt(n)*12
k = 2
mus[,round(etas[1]+1):n] = mus[,round(etas[1]+1):n] + 0.5*sqrt(2*max(log(exp(1)*p*log(n)/k^2), log(n)))*matrix(rep(c(rep(1,k)*1/sqrt(n)*10, rep(0,p-k)), n-round(etas[1])) ,nrow=p)
#mus[,round(etas[2]+1):n] = mus[,round(etas[2]+1):n] + matrix(rep(randomdense, n-round(etas[2])) ,nrow=p)*1.2
#mus[,round(etas[3]+1):n] = mus[,round(etas[3]+1):n] + matrix(rep(rep(1,p)/p^(1/4)/sqrt(n)*12, n-round(etas[3])) ,nrow=p)*1.2

plot(mus[1,])
X = mus+noise
plot(X[1,])



max = single_SBS_calibrate(n,p,1000,100,debug=FALSE)


res = single_SBS(X, max, debug=FALSE)


# checking double CUSUM and SBS
library(hdbinseg)


res = dcbs.alg(X,height=1, thr = 0,cp.type=1,phi=-1,temporal=FALSE )

res = dcbs.alg(X)


# subset by tickle
source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")

empirical_penalty = 0
for (i in 1:10) {
  mynulldata <- InspectChangepoint::rescale.variance(matrix(rnorm(n*p,0,1),nrow=p,ncol=n,byrow=FALSE)) # 5 variates with 1000 time points, no change
  empirical_penalty = max(empirical_penalty, wbs_penaltyfinder(mynulldata, SUBSET.normal_penalty, 100))
}
 # finds an empirical penalty for the SUBSET method (with a normal likelihood cost function)

result <- change_main(InspectChangepoint::rescale.variance(X), SUBSET.normal, 100, empirical_penalty) # running the change detection wrapper with the SUBSET method, using the empirical penalty calculated above; result is an S3 list object with 3 elements. The first element gives the value of the test statistic (the difference in the cost) at the detected changepoints, the second element gives the location of the detected changepoints, and the final element is a matrix of binaries indicating which variates are affected by the change.


# checking double CUSUM and SBS
library(hdbinseg)
# YY = X[,]
# YY[,] = 0
# dcbs.thr(X, interval = c(1, dim(X)[2]), phi = 0.5, cp.type = 1,
#          do.clean.cp = FALSE, temporal = TRUE, scales = NULL, diag = FALSE,
#          sgn = NULL, B = 1000, q = 0.001, do.parallel = 1)
# 
# res = dcbs.alg(X,height=1, thr = 0,cp.type=1,phi=-1,temporal=FALSE )

res14 = dcbs.alg(X, cp.type = 1, phi = 0.5, thr = NULL, trim = NULL,
         height = NULL, temporal = TRUE, scales = NULL, diag = FALSE,
         B = 1000, q = 0.001, do.parallel = 4)

res15 = sbs.alg(X, cp.type = 1, thr = NULL, trim = NULL, height = NULL,
        temporal = TRUE, scales = NULL, diag = FALSE, B = 1000, q = 0.001,
        do.parallel = 4)



overlaps = function(a,b){
  return((a[1]+1) <= (b[2]-1) && (b[1]+1) <= (a[2]-1))
}

