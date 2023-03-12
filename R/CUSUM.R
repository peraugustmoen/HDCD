#' @useDynLib HDCD CUSUM_R
#' @export
CUSUM = function(X, start, stop){
  p = dim(X)[1]
  n = dim(X)[2]
  
  
  #lens[7]=round(log(n)^2)
  #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
  # SEXP epsI, SEXP lambdaI, SEXP maxiterI)
  
  #CUSUM_R(SEXP XI, SEXP sI, SEXP eI, SEXP pI, SEXP nI)
  res = .Call(CUSUM_R, X, as.integer(start), as.integer(stop), as.integer(p), as.integer(n))

  
  return(res)
}

#' @useDynLib HDCD single_CUSUM_R
#' @export
single_CUSUM = function(X, start, stop, pos){
  p = dim(X)[1]
  n = dim(X)[2]
  
  
  #lens[7]=round(log(n)^2)
  #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
  # SEXP epsI, SEXP lambdaI, SEXP maxiterI)
  
  #CUSUM_R(SEXP XI, SEXP sI, SEXP eI, SEXP pI, SEXP nI)
  res = .Call(single_CUSUM_R, X, as.integer(start), as.integer(stop), as.integer(p), as.integer(pos),as.integer(n))
  
  
  return(res)
}