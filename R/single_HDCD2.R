#' @useDynLib HDCD cHDCD_single
#' @export
single_HDCD_two = function(X, threshold_d, threshold_s, debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]
  
  if(n<2){
    return(NULL)
  }
  max_s = sqrt(p*log(n))
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*4*p*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
  thresholds = nu_as[]
  
  thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n^4)/ss[2:length(ss)]^2), log(n^4))
  #thresholds[2:length(ss)] = threshold_s*(ss[2:length(ss)]*log(exp(1)-1 + sqrt(p*log(n^4))/ss[2:length(ss)])+ log(n^4))
  
  #thresholds[2:length(ss)] = threshold_s*(sqrt(p*exp(-as[2:length(ss)]^2/2)*log(n^4))+ log(n^4))
  
  
  #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n*log(sqrt(p*log(n)),2))/ss[2:length(ss)]^2), log(n*log(sqrt(p*log(n)),2)))
  thresholds[1] = threshold_d * max(sqrt(p*log(n^4)), log(n^4))
  
  #lens[7]=round(log(n)^2)
  #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
  # SEXP epsI, SEXP lambdaI, SEXP maxiterI)
  
  res = .Call(cHDCD_single, X,as.integer(n), as.integer(p), thresholds,
              as,nu_as, as.integer(length(as)), as.integer(debug))
  
  res$pos = res$pos+1
  res$s = ss[res$apos+1]
  
  
  return(res)
}