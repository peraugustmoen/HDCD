#' @useDynLib HDCD cScan
#' @export
Scan = function(X, threshold_const = 1, alpha = 1+1/6, K = 7, debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]
  
  lens = c(2)
  last = 2
  tmp = last
  while(alpha*last<n){
    #while(2*last<n){
    tmp = last
    last = floor(alpha*last)
    if(last==tmp){
      last = last+1
    }
    #last = last+1
    lens= c(lens, last)
  }
  
  
  sequ = 1:p
  Ts = sqrt(2/sequ) *log(choose(p,sequ)*n^4) + sqrt(2*log(choose(p,sequ)*n^4))
  Ts[p] = sqrt(2*log(n^4)) + sqrt(2/p)*log(n^4)
  Ts  = threshold_const*Ts
  
  res = .Call(cScan, X,as.integer(n), as.integer(p), Ts,
              as.integer(lens),as.integer(length(lens)), as.integer(K),
              as.integer(debug))
  
  
  if(res$changepointnum==0){
    return(NULL)
  }
  
  
  res$changepoints = as.integer(res$changepoints[1:res$changepointnum]+1)
  res$CUSUMval = res$CUSUMval[1:res$changepointnum]
  res$depth = as.integer(res$depth[1:res$changepointnum])
  res$coordinate = matrix(res$coordinate,nrow = p, ncol=n)
  #res$coordinate = res$coordinate[, 1:num_nonzero]
  srt_indices = as.integer(sort(res$changepoints, decreasing =FALSE, index.return=TRUE)$ix)
  res$changepoints = as.integer(res$changepoints[srt_indices])
  res$CUSUMval = res$CUSUMval[srt_indices]
  res$depth = as.integer(res$depth[srt_indices])
  res$coordinate = (res$coordinate[,srt_indices])
  res$startpoints = as.integer(res$startpoints[srt_indices])
  res$endpoints = as.integer(res$endpoints[srt_indices])
  res$s = as.integer(res$maxss[srt_indices])
  res$Ts = Ts
  
  
  return(res)
}