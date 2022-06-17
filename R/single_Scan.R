#' @useDynLib HDCD cScan_single
#' @export
single_Scan = function(X,debug=FALSE){
  p = dim(X)[1]
  n = dim(X)[2]
  sequ = 1:p
  Ts = sqrt(2/sequ) *log(choose(p,sequ)*n^4) + sqrt(2*log(choose(p,sequ)*n^4))
  Ts[p] = sqrt(2*log(n^4)) + sqrt(2/p)*log(n^4)
  
  res = .Call(cScan_single, X,as.integer(n), as.integer(p), as.integer(debug),Ts)
  
  res$pos = res$pos+1
  return(res)
}


# single_Scan = function(X){
#   p = dim(X)[1]
#   n = dim(X)[2]
#   
#   start = -1
#   stop = n-1
#   CUSUM_X = CUSUM(X, start, stop)
#   CUSUM_X = CUSUM_X[1:(p*(stop-start-1))]
#   CUSUM_X = matrix(data=CUSUM_X, nrow=p,ncol = stop-start-1)
#   
#   maxes = rep(NA, n-1)
#   maxes_p = rep(NA, n-1)
#   sequ = 1:p
#   
#   linear = (colSums(CUSUM_X^2) - p)/sqrt(2*p)
#   H = sqrt(2*log(n)) + sqrt(2/p)*log(n)
#   linear = linear -H
#   
#   Ts = sqrt(2/sequ) *log(choose(p,sequ)*n) + sqrt(2*log(choose(p,sequ)*n))
#   for (t in 1:(n-1)) {
#     z_squared = CUSUM_X[,t]^2
#     z_squared_sorted = sort(z_squared,decreasing=TRUE)
#     cumsumm = cumsum(z_squared_sorted)
#     stats = 1/sqrt(2*sequ)*(cumsumm - sequ)-Ts
#     maxes_p[t] = which.max(stats)
#     maxes[t] = stats[maxes_p[t]]
#     #maxes[t] = max(stats[maxes_p[t]], linear)
#     if(maxes[t] < linear[t]){
#       maxes[t] = linear[t]
#       maxes_p[t] = p
#     }
#   }
#   
#   
#   
#   maxind = which.max(maxes)
#   ret = c()
#   ret$pos = maxind
#   ret$s = maxes_p[maxind]
#   return(ret)
# }