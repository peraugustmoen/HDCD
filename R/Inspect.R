#' @useDynLib Krafthack2022 myInspect
#' @export
myInspectR = function(X, threshold=2, adaptive_threshold=TRUE, alpha = 1.2, K = 20,eps=1e-10,
                     lambda = 1, maxiter=10000,debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]

  #maxjump = ceiling(2*sqrt(log(n)))
  #maxjump = ceiling(2*sqrt(log(n)))



  #last = as.integer(round(sqrt(log(n))))
  #last = as.integer(3)
  lens = c(1)
  last = 1
  tmp = last
  while(alpha*last<n){
    #while(2*last<n){
    tmp = last
    last = round(alpha*last)
    if(last==tmp){
      last = last+1
    }
    #last = last+1
    lens= c(lens, last)
  }
  #lens[7]=round(log(n)^2)
  #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
  # SEXP epsI, SEXP lambdaI, SEXP maxiterI)
  res = .Call(myInspect, X,as.integer(n), as.integer(p), threshold, as.integer(adaptive_threshold),
              as.integer(lens),as.integer(length(lens)), as.integer(K), eps, lambda, as.integer(maxiter),
              as.integer(debug))
  num_nonzero = sum(res$changepoints!=0)
  res$changepoints = res$changepoints[1:num_nonzero]
  res$CUSUMval = res$CUSUMval[1:num_nonzero]
  res$depth = res$depth[1:num_nonzero]
  srt_indices = sort(res$changepoints, decreasing =FALSE, index.return=TRUE)$ix
  res$changepoints = res$changepoints[srt_indices]
  res$CUSUMval = res$CUSUMval[srt_indices]
  res$depth = res$depth[srt_indices]
  if(num_nonzero==0){
    res$coordinate = NULL
  }
  else{
    res$coordinate = matrix(res$coordinate,nrow = p, ncol=n)
    res$coordinate = res$coordinate[,1:num_nonzero]
  }

  return(res)
}

##  Hausdorff distance
hausdorff = function(v1, v2){
  max = 0
  for (v  in v1) {
    tmp = min(abs(v-v2))
    if(tmp>max){
      max=tmp
    }
  }
  for (v  in v2) {
    tmp = min(abs(v-v1))
    if(tmp>max){
      max=tmp
    }
  }
  return(max)
}

hausdorff_avg = function(v1, v2){
  dists = rep(NA, length(v1))
  for (i in 1:length(v1)) {
    v = v1[i]
    dists[i] = min(abs(v-v2))

  }

  return(mean(dists))
}
