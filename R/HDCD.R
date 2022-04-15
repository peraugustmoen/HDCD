#' @useDynLib HDCD cHDCD
#' @export
HDCD = function(X, threshold_d, threshold_s, alpha = 1+1/6, K = 7, debug =FALSE){
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
  
  max_s = sqrt(4*p*log(n))
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*4*p*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
  thresholds = nu_as[]
  
  #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n^4)/ss[2:length(ss)]^2), log(n^4))
  thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)-1 + sqrt(p*log(n^4))/ss[2:length(ss)]), log(n^4))
  
  #thresholds[2:length(ss)] = threshold_s*(sqrt(p*exp(-as[2:length(ss)]^2/2)*log(n^4))+ log(n^4))
  
  
  #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n*log(sqrt(p*log(n)),2))/ss[2:length(ss)]^2), log(n*log(sqrt(p*log(n)),2)))
  thresholds[1] = threshold_d * sqrt(p*log(n^4))
  
  #lens[7]=round(log(n)^2)
  #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
  # SEXP epsI, SEXP lambdaI, SEXP maxiterI)

  res = .Call(cHDCD, X,as.integer(n), as.integer(p), thresholds,
              as.integer(lens),as.integer(length(lens)), as.integer(K), as,nu_as, 
              as.integer(length(as)), as.integer(debug))
  
  
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
  res$coordinate = as.integer(res$coordinate[,srt_indices])
  res$startpoints = as.integer(res$startpoints[srt_indices])
  res$endpoints = as.integer(res$endpoints[srt_indices])
  res$maxaposes = as.integer(res$maxaposes[srt_indices])
  res$s = as.integer(ss[res$maxaposes+1])
  res$thresholds = thresholds
  

  return(res)
}


# #' @useDynLib HDCD cHDCD
# #' @export
# HDCD = function(X, threshold_d, threshold_s, alpha = 1+1/6, K = 7, debug =FALSE){
#   p = dim(X)[1]
#   n = dim(X)[2]
#   
#   lens = c(2)
#   last = 2
#   tmp = last
#   while(alpha*last<n){
#     #while(2*last<n){
#     tmp = last
#     last = floor(alpha*last)
#     if(last==tmp){
#       last = last+1
#     }
#     #last = last+1
#     lens= c(lens, last)
#   }
#   
#   max_s = sqrt(4*p*log(n))
#   log2ss = 0:floor(log(max_s, base=2))
#   ss = 2^(log2ss)
#   ss = c(p, rev(ss))
#   as = ss[]
#   as[2:length(ss)] = 2*sqrt(log(4*p*log(n)/ss[2:length(ss)]^2))
#   as[1] = 0
#   nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
#   thresholds = nu_as[]
#   
#   #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n^4)/ss[2:length(ss)]^2), log(n^4))
#   thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]^2/exp(1)/sqrt(p*log(n^4)), log(n^4))
#   
#   
#   #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n*log(sqrt(p*log(n)),2))/ss[2:length(ss)]^2), log(n*log(sqrt(p*log(n)),2)))
#   thresholds[1] = threshold_d * sqrt(p*log(n))
#   
#   #lens[7]=round(log(n)^2)
#   #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
#   # SEXP epsI, SEXP lambdaI, SEXP maxiterI)
#   
#   res = .Call(cHDCD, X,as.integer(n), as.integer(p), thresholds,
#               as.integer(lens),as.integer(length(lens)), as.integer(K), as,nu_as, 
#               as.integer(length(as)), as.integer(debug))
#   
#   num_nonzero = sum(res$changepoints!=-1)
#   if(num_nonzero==0){
#     res$changepoints = c()
#     res$CUSUMval = c()
#     res$depth = c()
#     res$coordinate = c()
#     return(res)
#   }
#   
#   res$changepoints = res$changepoints[1:num_nonzero]+1
#   res$CUSUMval = res$CUSUMval[1:num_nonzero]
#   res$depth = res$depth[1:num_nonzero]
#   res$coordinate = matrix(res$coordinate,nrow = p, ncol=n)
#   #res$coordinate = res$coordinate[, 1:num_nonzero]
#   srt_indices = sort(res$changepoints, decreasing =FALSE, index.return=TRUE)$ix
#   res$changepoints = res$changepoints[srt_indices]
#   res$CUSUMval = res$CUSUMval[srt_indices]
#   res$depth = res$depth[srt_indices]
#   res$coordinate = res$coordinate[,srt_indices]
#   
#   return(res)
# }