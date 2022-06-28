#' @useDynLib HDCD cPilliat
#' @export
Pilliat = function(X, threshold_d_const=4, threshold_bj_const=6, threshold_partial_const=4,
                  K = 2, alpha = 1+1/6, debug =FALSE){
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
  
  # dense
  threshold_dense = threshold_d_const * (sqrt(p*log(2*n^2)) + log(2*n^2))
  
  #partial norm
  thresholds_partial = c()
  
  t = 1
  while(TRUE){
    thresholds_partial = c(thresholds_partial, t*log(2*exp(1)*p/t) + log(2*n^2))
    t = 2*t
    if(t>=p){
      break
    }
  }
  thresholds_partial = thresholds_partial*threshold_partial_const
  
  # berk jones
  
  thresholds_bj = c()
  xx = 1
  maxx = 1
  while(TRUE){
    #p0 = threshold_bj_const/(pi^2*xx^2*n^2)
    logp0 = log(threshold_bj_const) - 2*log(pi) - 2*log(xx) - 2*log(n)
    #qq = qbinom(1/n, p, p0,lower.tail=FALSE)
    qq = qbinom(p = logp0, size = p, prob = exp(log(2) + pnorm(xx,lower.tail = FALSE,log = TRUE)),
                lower.tail=FALSE,log.p = TRUE)
    if(debug){
      print(xx)
      print(p0)
      print(exp(log(2) + pnorm(xx,lower.tail = FALSE,log = TRUE)))
      print(qq)
    }
    if(qq==0){
      maxx = xx-1
      break
    }
    thresholds_bj = c(thresholds_bj,qq)
    xx = xx+1
  }
  
  if(debug){
    print(threshold_dense)
    print(thresholds_partial)
    print(thresholds_bj)
  }
  
  
  #res = .Call(cPilliat, X,as.integer(n), as.integer(p), thresholds,thresholds_test,
   #           as.integer(lens),as.integer(length(lens)), as.integer(K), as,nu_as, 
    #          as.integer(length(as)), as.integer(twologn), as.integer(ts),as.integer(debug))
  res = .Call(cPilliat, X, as.integer(n), as.integer(p), as.numeric(thresholds_partial), 
        as.numeric(threshold_dense), as.integer( thresholds_bj),
        as.integer(lens), as.integer(length(lens)), as.integer(K), as.integer(maxx),as.integer(debug))
  
  if(res$number_of_changepoints==0){
    return(NULL)
  }
  
  
  res$changepoints = as.integer(res$changepoints[1:res$number_of_changepoints]+1)

  #res$coordinate = res$coordinate[, 1:num_nonzero]
  srt_indices = as.integer(sort(res$changepoints, decreasing =FALSE, index.return=TRUE)$ix)
  res$changepoints = as.integer(res$changepoints[srt_indices])
  res$startpoints = as.integer(res$startpoints[srt_indices])
  res$endpoints = as.integer(res$endpoints[srt_indices])
  res$test_stat = as.integer(res$test_stat[srt_indices])
  
  
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