#' @useDynLib HDCD cInspect
#' @export
Inspect = function(X, lambda = sqrt(log(p*log(n))/2) , xi = sqrt(log(p*log(n))/2) , alpha = 1+1/6, K = 7,eps=1e-10,
                     maxiter=10000,rescale_variance = TRUE, debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]

  # FIX AUTOMATIC MONTE CARLO HERE



  #last = as.integer(round(sqrt(log(n))))
  #last = as.integer(3)
  lens = c(1)
  last = 1
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
  #lens[7]=round(log(n)^2)
  #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
  # SEXP epsI, SEXP lambdaI, SEXP maxiterI)
  res = .Call(cInspect, X[,],as.integer(n), as.integer(p), xi,
              as.integer(lens),as.integer(length(lens)), as.integer(K), eps, lambda, as.integer(maxiter),
              as.integer(rescale_variance),as.integer(debug))
  if(res$changepointnumber==0){
    return(list(scales = res$scales))
  }
  
  res$changepoints = as.integer(res$changepoints[1:res$changepointnumber]+1)
  res$CUSUMval = res$CUSUMval[1:res$changepointnumber]
  res$depth = as.integer(res$depth[1:res$changepointnumber])
  res$coordinate = matrix(res$coordinate,nrow = p, ncol=n)
  #res$coordinate = res$coordinate[, 1:num_nonzero]
  srt_indices = as.integer(sort(res$changepoints, decreasing =FALSE, index.return=TRUE)$ix)
  res$changepoints = as.integer(res$changepoints[srt_indices])
  res$CUSUMval = res$CUSUMval[srt_indices]
  res$depth = as.integer(res$depth[srt_indices])
  res$coordinate = res$coordinate[,srt_indices]

  return(res)
}

# #' @useDynLib HDCD sparse_svd
# #' @export
# check_svd = function(X, r1, c1, lambda, eps, maxiter){
#   return(.Call(sparse_svd, X, as.integer(r1), as.integer(c1), lambda, eps, as.integer(maxiter)))
# }


#' @useDynLib HDCD cInspect_calibrate
#' @export
Inspect_calibrate = function(n, p, N, tol,lambda = sqrt(log(p*log(n))/2) , alpha = 1+1/6, K = 7,eps=1e-10,
                   maxiter=10000,rescale_variance = TRUE,debug =FALSE){

  
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
    last = floor(alpha*last)
    if(last==tmp){
      last = last+1
    }
    #last = last+1
    lens= c(lens, last)
  }
  #lens[7]=round(log(n)^2)
  #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
  # SEXP epsI, SEXP lambdaI, SEXP maxiterI)
  toln = max(round(N*tol),1)

  res = .Call(cInspect_calibrate, as.integer(n), as.integer(p), as.integer(N),
              as.integer(toln),as.integer(lens), as.integer(length(lens)), as.integer(K),
              as.numeric(lambda), as.numeric(eps), as.integer(maxiter),
              as.integer(rescale_variance),as.integer(debug))
  
  return(res)
}
