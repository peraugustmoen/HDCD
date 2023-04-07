#' #' @useDynLib HDCD cESAC_test
#' #' @useDynLib HDCD cESAC_test_calibrate
#' #' @export
#' changepoint_test_ESAC = function(X, empirical =FALSE, N = 100, tol = 0.05,thresholds=NULL, threshold_d, threshold_s,debug =FALSE,
#'                  droppartialsum = FALSE ){
#'   p = dim(X)[1]
#'   n = dim(X)[2]
#'   
#'   
#'   max_s = min(sqrt(4*p*log(n)), p)
#'   log2ss = 0:floor(log(max_s, base=2))
#'   ss = 2^(log2ss)
#'   ss = c(p, rev(ss))
#'   as = ss[]
#'   as[2:length(ss)] = sqrt(2*log(exp(1)*4*p*log(n)/ss[2:length(ss)]^2))
#'   as[1] = 0
#'   nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
#'   twologn = min(p, floor(2*log(n)))
#'   ts = ss[]
#'   
#'   if(empirical){
#'     if(is.null(thresholds)){
#'       res = changepoint_test_ESAC_calibrate(n,p, as, nu_as, ts, twologn,N,tol,debug)
#'       if(droppartialsum){
#'         thresholds = res[[1]]
#'       }
#'       else{
#'         thresholds = res[[2]]
#'       }
#'     }
#'   }else{
#'     thresholds = nu_as[]
#' 
#'     thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n^4)/ss[2:length(ss)]^2), log(n^4))
#'     #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)-1 + sqrt(p*log(n^4))/ss[2:length(ss)]), log(n^4))
#'     #thresholds[2:length(ss)] = threshold_s*(ss[2:length(ss)]*log(exp(1)-1 + sqrt(p*log(n^4))/ss[2:length(ss)])+ log(n^4))
#'     
#'     #thresholds[2:length(ss)] = threshold_s*(sqrt(p*exp(-as[2:length(ss)]^2/2)*log(n^4))+ log(n^4))
#'     
#'     
#'     #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n*log(sqrt(p*log(n)),2))/ss[2:length(ss)]^2), log(n*log(sqrt(p*log(n)),2)))
#'     thresholds[1] = threshold_d * (sqrt(p*log(n^4)) + log(n^4))
#'   }
#' 
#'   ind = thresholds > thresholds[1]
#'   thresholds = thresholds[!ind]
#'   as = as[!ind]
#'   nu_as = nu_as[!ind]
#'   #thresholds_test = threshold_test_const* thresholds
#'   ts = ss[!ind]
#'   twologn = min(p, floor(2*log(n)))
#'   if(droppartialsum){
#'     twologn=0
#'   }
#'   #print(twologn)
#'   #lens[7]=round(log(n)^2)
#'   #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
#'   # SEXP epsI, SEXP lambdaI, SEXP maxiterI)
#' 
#'   res = .Call(cESAC_test, as.numeric(X), as.integer(n), as.integer(p), as.numeric(thresholds), 
#'               as.numeric(as), as.numeric(nu_as), as.integer(length(as)), as.integer(droppartialsum), 
#'                                          as.integer(twologn), as.integer(ts), as.integer(debug))
#'   
#'   return(res)
#' }
#' 
#' #' @export
#' changepoint_test_ESAC_calibrate = function(n,p, as=NULL, nu_as=NULL, ts=NULL, twologn=NULL, N, tol,debug=FALSE){
#'   if(is.null(as) || is.null(nu_as) || is.null(ts) || is.null(twologn)){
#'     max_s = min(sqrt(4*p*log(n)), p)
#'     log2ss = 0:floor(log(max_s, base=2))
#'     ss = 2^(log2ss)
#'     ss = c(p, rev(ss))
#'     as = ss[]
#'     as[2:length(ss)] = sqrt(2*log(exp(1)*4*p*log(n)/ss[2:length(ss)]^2))
#'     as[1] = 0
#'     nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
#'     twologn = min(p, floor(2*log(n)))
#'     ts = ss[]
#'   }
#'   
#'   toln = max(round(N*tol),1)
#'   
#'   
#'   
#'   res= .Call(cESAC_test_calibrate, as.integer(n), as.integer(p), as.integer(N), 
#'                as.integer(toln), as.numeric(as), as.numeric(nu_as), as.integer(length(as)), 
#'                as.integer(twologn), as.integer(ts), as.integer(debug))
#'   #res[["as"]] = as
#'   #res[["nu_as"]] = nu_as
#'   return(res)
#' }
#' 
#' 
#' #' @export
#' changepoint_test_ESAC_modified = function(X, empirical =FALSE, N = 100, tol = 0.05,thresholds=NULL, threshold_d, threshold_s,debug =FALSE,
#'                                  droppartialsum = FALSE ){
#'   p = dim(X)[1]
#'   n = dim(X)[2]
#'   
#'   
#'   max_s = min(sqrt(p*log(n)), p)
#'   log2ss = 0:floor(log(max_s, base=2))
#'   ss = 2^(log2ss)
#'   ss = c(p, rev(ss))
#'   as = ss[]
#'   as[2:length(ss)] = sqrt(2*log(exp(1)*p*log(n)/ss[2:length(ss)]^2))
#'   as[1] = 0
#'   nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
#'   twologn = min(p, floor(log(n)/2))
#'   ts = ss[]
#'   
#'   if(empirical){
#'     if(is.null(thresholds)){
#'       res = changepoint_test_ESAC_calibrate(n,p, as, nu_as, ts, twologn,N,tol,debug)
#'       if(droppartialsum){
#'         thresholds = res[[1]]
#'       }
#'       else{
#'         thresholds = res[[2]]
#'       }
#'     }
#'   }else{
#'     thresholds = nu_as[]
#'     
#'     thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n)/ss[2:length(ss)]^2), log(n))
#'     #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)-1 + sqrt(p*log(n^4))/ss[2:length(ss)]), log(n^4))
#'     #thresholds[2:length(ss)] = threshold_s*(ss[2:length(ss)]*log(exp(1)-1 + sqrt(p*log(n^4))/ss[2:length(ss)])+ log(n^4))
#'     
#'     #thresholds[2:length(ss)] = threshold_s*(sqrt(p*exp(-as[2:length(ss)]^2/2)*log(n^4))+ log(n^4))
#'     
#'     
#'     #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n*log(sqrt(p*log(n)),2))/ss[2:length(ss)]^2), log(n*log(sqrt(p*log(n)),2)))
#'     thresholds[1] = threshold_d * (sqrt(p*log(n)) + log(n))
#'   }
#'   
#'   ind = thresholds > thresholds[1]
#'   thresholds = thresholds[!ind]
#'   as = as[!ind]
#'   nu_as = nu_as[!ind]
#'   #thresholds_test = threshold_test_const* thresholds
#'   ts = ss[!ind]
#'   twologn = min(p, floor(log(n)/2))
#'   if(droppartialsum){
#'     twologn=0
#'   }
#'   #print(twologn)
#'   #lens[7]=round(log(n)^2)
#'   #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
#'   # SEXP epsI, SEXP lambdaI, SEXP maxiterI)
#'   
#'   res = .Call(cESAC_test, as.numeric(X), as.integer(n), as.integer(p), as.numeric(thresholds), 
#'               as.numeric(as), as.numeric(nu_as), as.integer(length(as)), as.integer(droppartialsum), 
#'               as.integer(twologn), as.integer(ts), as.integer(debug))
#'   
#'   return(res)
#' }
#' 
#' #' @export
#' changepoint_test_ESAC_modified_calibrate = function(n,p, as=NULL, nu_as=NULL, ts=NULL, twologn=NULL, N, tol,debug=FALSE){
#'   if(is.null(as) || is.null(nu_as) || is.null(ts) || is.null(twologn)){
#'     max_s = min(sqrt(p*log(n)), p)
#'     log2ss = 0:floor(log(max_s, base=2))
#'     ss = 2^(log2ss)
#'     ss = c(p, rev(ss))
#'     as = ss[]
#'     as[2:length(ss)] = sqrt(2*log(exp(1)*p*log(n)/ss[2:length(ss)]^2))
#'     as[1] = 0
#'     nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
#'     twologn = min(p, floor(2*log(n)))
#'     ts = ss[]
#'   }
#'   
#'   toln = max(round(N*tol),1)
#'   
#'   
#'   
#'   res= .Call(cESAC_test_calibrate, as.integer(n), as.integer(p), as.integer(N), 
#'              as.integer(toln), as.numeric(as), as.numeric(nu_as), as.integer(length(as)), 
#'              as.integer(twologn), as.integer(ts), as.integer(debug))
#'   #res[["as"]] = as
#'   #res[["nu_as"]] = nu_as
#'   return(res)
#' }
#' 
#' 
#' 
