#' @useDynLib HDCD cInspect_single
#' @export
single_Inspect = function(X, lambda, eps=1e-10,
                  maxiter=10000,debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]

  res = .Call(cInspect_single, X,as.integer(n), as.integer(p), 0,
              eps, lambda, as.integer(maxiter), as.integer(debug))
  
  
  
  res$pos = res$pos+1
  
  return(res)
}