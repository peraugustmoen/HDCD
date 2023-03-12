#' @useDynLib HDCD cSBS_single
#' @export
single_SBS = function(X, threshold, rescale_variance = TRUE,debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]
  
  if(n<2){
    return(NULL)
  }
  
  
  res = .Call(cSBS_single, X,as.integer(n), as.integer(p), as.numeric(threshold),
              as.integer(rescale_variance),as.integer(debug))
  
  res$pos = res$pos+1

  
  return(res)
}

#' @useDynLib HDCD cSBS_single_calibrate
#' @export
single_SBS_calibrate = function(n,p,N, tol, rescale_variance = TRUE,debug =FALSE){

  
  if(n<2){
    return(NULL)
  }
  toln = max(round(N*tol),1)
  
  res = .Call(cSBS_single_calibrate, as.integer(n), as.integer(p), as.integer(N),as.integer(toln),
          as.integer(rescale_variance),as.integer(debug))
  
  
  
  return(res)
}