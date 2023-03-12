#' @useDynLib HDCD rescale_variance_R
#' @export
rescale_variance = function(X,debug=FALSE){
  p = dim(X)[1]
  n = dim(X)[2]

  res = .Call(rescale_variance_R, as.numeric(X[,]),as.integer(n), as.integer(p), 
              as.integer(debug))
  
  
  return(res)
}
