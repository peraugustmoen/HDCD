#' @useDynLib HDCD cSBS_single
#' @export
single_SBS = function(X, threshold, debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]
  
  if(n<2){
    return(NULL)
  }
  
  
  res = .Call(cSBS_single, X,as.integer(n), as.integer(p), as.numeric(threshold),
              as.integer(debug))
  
  res$pos = res$pos+1

  
  return(res)
}

#' @useDynLib HDCD cSBS_single_calibrate
#' @export
single_SBS_calibrate = function(n,p,N, toln, debug =FALSE){

  
  if(n<2){
    return(NULL)
  }
  
  
  res = .Call(cSBS_single_calibrate, as.integer(n), as.integer(p), as.integer(N),as.integer(toln),
          as.integer(debug))
  
  
  
  return(res)
}