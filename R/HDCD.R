#' @useDynLib HDCD cHDCD
#' @useDynLib HDCD cHDCD_calibrate
#' @useDynLib HDCD cHDCD_test
#' @useDynLib HDCD cHDCD_test_calibrate
#' @export
HDCD = function(X, threshold_d=1.5, threshold_s=1.0, alpha = 1.5, K = 5, debug =FALSE,
                empirical = FALSE, thresholds = NULL, thresholds_test = NULL,
                threshold_d_test = threshold_d, threshold_s_test = threshold_s, droppartialsum = TRUE, fast = FALSE,
                rescale_variance = TRUE,trim = FALSE,NOT=TRUE, midpoint = FALSE){
  
  ## FIX THE AUTOMATIC MONTE CARLO HERE. tol and N are missing
  
  
  p = dim(X)[1]
  n = dim(X)[2]

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

  
  max_s = min(sqrt(p*log(n)), p)
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*p*4*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
  if(debug){
    print(ss)
  }
  
  ts = ss
  twologn = min(p, floor(log(n)/2))
  if(droppartialsum){
    twologn=0
  }

  if(is.null(thresholds_test)){
    if(empirical){
      ttt = HDCD_calibrate(n,p, as, nu_as, ts, twologn, lens, K,
                           alpha, K, N, tol,fast, rescale_variance,debug)
      if(droppartialsum){
        thresholds_test = ttt[[1]]
      }else{
        thresholds_test = ttt[[2]]
      }
    }
    else{
      thresholds_test = nu_as[] 
      thresholds_test[2:length(ss)] = threshold_s_test*(ss[2:length(ss)]*log(exp(1)*p*log(n^4)/ss[2:length(ss)]^2)+ log(n^4))
      thresholds_test[1] = threshold_d_test* (sqrt(p*log(n^4)) + log(n^4))
    }
    
    # if(debug){
    #   print(thresholds_test)
    # }
  }
  
  
  if(is.null(thresholds)){
    
    
    thresholds = nu_as[]
    
    thresholds[2:length(ss)] = threshold_s*(ss[2:length(ss)]*log(exp(1)*4*p*log(n)/ss[2:length(ss)]^2)+ log(n^4))
    #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)-1 + sqrt(p*log(n^4))/ss[2:length(ss)]), log(n^4))
    #thresholds[2:length(ss)] = threshold_s*(ss[2:length(ss)]*log(exp(1)-1 + sqrt(p*log(n^4))/ss[2:length(ss)])+ log(n^4))
    
    #thresholds[2:length(ss)] = threshold_s*(sqrt(p*exp(-as[2:length(ss)]^2/2)*log(n^4))+ log(n^4))
    
    
    #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n*log(sqrt(p*log(n)),2))/ss[2:length(ss)]^2), log(n*log(sqrt(p*log(n)),2)))
    thresholds[1] = threshold_d * (sqrt(p*log(n^4)) + log(n^4))
  }

  #ind = thresholds > thresholds[1]
  #thresholds = thresholds[!ind]
  #thresholds_test = thresholds_test[!ind]
  #as = as[!ind]
  #nu_as = nu_as[!ind]
  #thresholds_test = threshold_test_const* thresholds
  #ts = ss[!ind]

  if(debug){
    print("thresholds:")
    print(thresholds)
    print("thresholds_test:")
    print(thresholds_test)
    print("as:")
    print(as)
  }
  if(!empirical && twologn != 0){
    tt = 1
    count = 0
    h = length(thresholds_test)
    while(tt <=twologn){
      thresholds_test[h - count] = thresholds_test[h-count] + ts[h-count] * as[h-count]*as[h-count]  
      count = count+1
      tt = 2*tt
    }
    
  }
  
  #print(twologn)
  #lens[7]=round(log(n)^2)
  #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
  # SEXP epsI, SEXP lambdaI, SEXP maxiterI)

  res = .Call(cHDCD, X[,],as.integer(n), as.integer(p), thresholds,thresholds_test,
              as.integer(lens),as.integer(length(lens)), as.integer(K), as,nu_as, 
              as.integer(length(as)), as.integer(twologn), as.integer(ts), as.integer(fast),as.integer(rescale_variance),
              as.integer(trim),as.integer(NOT),as.integer(midpoint),as.integer(debug))
  
  
  if(res$changepointnum==0){
    return(list(scales = res$scales))
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
  res$coordinate = (res$coordinate[,srt_indices])
  res$startpoints = as.integer(res$startpoints[srt_indices])
  res$endpoints = as.integer(res$endpoints[srt_indices])
  res$maxaposes = as.integer(res$maxaposes[srt_indices])
  res$s = as.integer(ss[res$maxaposes+1])
  res$thresholds = thresholds
  res$thresholds_test = thresholds_test
  res$ts = ts
  

  return(res)
}

#' @export
HDCD_test = function(X, threshold_d=1.5, threshold_s=1.0,debug =FALSE,
                empirical = FALSE, thresholds = NULL, droppartialsum = TRUE, fast = FALSE,
                rescale_variance = TRUE){
  
  ## FIX THE AUTOMATIC MONTE CARLO HERE. tol and N are missing
  
  
  p = dim(X)[1]
  n = dim(X)[2]
  

  
  max_s = min(sqrt(p*log(n)), p)
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*p*4*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
  if(debug){
    print(ss)
  }
  
  ts = ss
  twologn = min(p, floor(log(n)/2))
  if(droppartialsum){
    twologn=0
  }
  
  if(is.null(thresholds)){
    if(empirical){
      ttt = HDCD_test_calibrate(n,p, as, nu_as, ts, twologn, lens, K,
                           alpha, K, N, tol,fast, rescale_variance,debug)
      if(droppartialsum){
        thresholds = ttt[[1]]
      }else{
        thresholds= ttt[[2]]
      }
    }
    else{
      thresholds = nu_as[] 
      thresholds[2:length(ss)] = threshold_s*(ss[2:length(ss)]*log(exp(1)*p*log(n^4)/ss[2:length(ss)]^2)+ log(n^4))
      thresholds[1] = threshold_d* (sqrt(p*log(n^4)) + log(n^4))
      if(!droppartialsum){
        tt = 1
        count = 0
        h = length(thresholds)
        while(tt <=twologn){
          thresholds[h - count] = thresholds[h-count] + ts[h-count] * as[h-count]*as[h-count]  
          count = count+1
          tt = 2*tt
        }
      }
    }
    
    # if(debug){
    #   print(thresholds_test)
    # }
  }
  
  
  
  
  #ind = thresholds > thresholds[1]
  #thresholds = thresholds[!ind]
  #thresholds_test = thresholds_test[!ind]
  #as = as[!ind]
  #nu_as = nu_as[!ind]
  #thresholds_test = threshold_test_const* thresholds
  #ts = ss[!ind]
  
  if(debug){
    print("thresholds:")
    print(thresholds)
    print("as:")
    print(as)
  }

  
  #print(twologn)
  #lens[7]=round(log(n)^2)
  #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
  # SEXP epsI, SEXP lambdaI, SEXP maxiterI)
  
  res = .Call(cHDCD_test, X[,],as.integer(n), as.integer(p), thresholds,
              as,nu_as, as.integer(length(as)),
              as.integer(twologn), as.integer(ts), as.integer(rescale_variance), 
              as.integer(fast), as.integer(debug))
              
            
  
  return(res)
}

#' @export
HDCD_test_calibrate = function(n, p, N, tol, fast = TRUE, rescale_variance = FALSE, droppartialsum = TRUE,
                               debug =FALSE){
  
  ## FIX THE AUTOMATIC MONTE CARLO HERE. tol and N are missing
  

  
  max_s = min(sqrt(p*log(n)), p)
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*p*4*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
  if(debug){
    print(ss)
  }
  
  ts = ss
  twologn = min(p, floor(log(n)/2))
  if(droppartialsum){
    twologn=0
  }
  
  
  toln = max(round(N*tol),1)
  
  
  res = .Call(cHDCD_test_calibrate, as.integer(n), as.integer(p), as.integer(N), as.integer(toln),
              as, nu_as, as.integer(length(as)), as.integer(twologn), as.integer(ts), as.integer(fast), 
              as.integer(rescale_variance), as.integer(debug))
  
  
  return(res)
}


#' @export
HDCD_calibrate = function(n,p, alpha = 1+1/6, K = 7, N, tol, bonferroni = TRUE, fast = FALSE,rescale_variance = TRUE,
                          droppartialsum = TRUE, debug=FALSE){
  
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
  
  max_s = min(sqrt(p*log(n)), p)
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*p*4*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
  if(debug){
    print(ss)
  }
  
  ts = ss
  twologn = min(p, floor(log(n)/2))
  if(droppartialsum){
    twologn=0
  }
  
  if(debug){
    print("as:")
    print(as)
  }
  
  #print(twologn)
  #lens[7]=round(log(n)^2)
  #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
  # SEXP epsI, SEXP lambdaI, SEXP maxiterI)
  
  
  if(bonferroni){
    hdcdtreshinfo =  HDCD_thresholds(n,p)
    h = hdcdtreshinfo[[2]]
    v = length(hdcdtreshinfo[[1]])
    if(h <= v && v >1){
      corr = 3
    }else if(v>1){
      corr = 2
    }else{
      corr = 1
    }
    tol = tol/corr
    N = corr*N
  }
  toln = max(round(N*tol),1)
  
  res= .Call(cHDCD_calibrate, as.integer(n), as.integer(p), as.integer(N),
             as.integer(toln), as.integer(lens), as.integer(length(lens)),
             as.integer(K), as.numeric(as), as.numeric(nu_as), as.integer(length(as)),
             as.integer(twologn), as.integer(ts), as.integer(fast),
             as.integer(rescale_variance),as.integer(debug))
  
  if(bonferroni){
    hdcdtreshinfo =  HDCD_thresholds(n,p)
    ttt = hdcdtreshinfo[[1]]
    h = hdcdtreshinfo[[2]]
    v = length(ttt)
    
    if(h <=length(hdcdtreshinfo[[1]]) && v>1){
      
      res[[1]] = c(res[[1]][1], max(res[[1]][2:(h-1)] / ttt[2:(h-1)])*ttt[2:(h-1)], max(res[[1]][h:v] / ttt[h:v])*ttt[h:v] )
      res[[2]] = c(res[[2]][1], max(res[[2]][2:(h-1)] / ttt[2:(h-1)])*ttt[2:(h-1)], max(res[[2]][h:v] / ttt[h:v])*ttt[h:v] )
    }else if(v>1){
      res[[1]] = c(res[[1]][1], max(res[[1]][2:v] / ttt[2:v])*ttt[2:v])
      res[[2]] = c(res[[2]][1], max(res[[2]][2:v] / ttt[2:v])*ttt[2:v])
    }else{
      res[[1]] = max(res[[1]] / ttt)*ttt
      res[[2]] = max(res[[2]] / ttt)*ttt
    }
  }
  return(res)
  #return(NULL)
}

# 
# #'
# HDCD_calibrate = function(n,p, as=NULL, nu_as=NULL, ts=NULL, twologn=NULL, lens = NULL,
#                           alpha = 1+1/6, K = 7, N, tol, fast = FALSE,rescale_variance = TRUE,debug=FALSE){
#   if(is.null(as) || is.null(nu_as) || is.null(ts) || is.null(twologn) ||
#      is.null(lens)){
#     lens = c(1)
#     last = 1
#     tmp = last
#     while(alpha*last<n){
#       #while(2*last<n){
#       tmp = last
#       last = floor(alpha*last)
#       if(last==tmp){
#         last = last+1
#       }
#       #last = last+1
#       lens= c(lens, last)
#     }
#     max_s = min(sqrt(p*log(n)), p)
#     log2ss = 0:floor(log(max_s, base=2))
#     ss = 2^(log2ss)
#     ss = c(p, rev(ss))
#     as = ss[]
#     as[2:length(ss)] = sqrt(2*log(exp(1)*4*p*log(n)/ss[2:length(ss)]^2))
#     as[1] = 0
#     nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
#     twologn = min(p, floor(log(n)/2))
#     ts = ss[]
#   }
#   if(debug){
#     print(n)
#     print(p)
#     print(N)
#     print(tol)
#     print("as:")
#     print(as)
#   }
#   
#   toln = max(round(N*tol),1)
#   
#  
#   res= .Call(cHDCD_calibrate, as.integer(n), as.integer(p), as.integer(N),
#              as.integer(toln), as.integer(lens), as.integer(length(lens)),
#              as.integer(K), as.numeric(as), as.numeric(nu_as), as.integer(length(as)),
#              as.integer(twologn), as.integer(ts), as.integer(fast),
#              as.integer(rescale_variance),as.integer(debug))
# 
#   return(res)
#   #return(NULL)
# }


#' @export
HDCD_modified = function(X, threshold_d=2, threshold_s=2, alpha = 1+1/6, K = 7, debug =FALSE,
                empirical = FALSE, thresholds_test = NULL,
                threshold_d_test = threshold_d, threshold_s_test = threshold_s, droppartialsum =TRUE, fast = FALSE){
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
  
  max_s = min(sqrt(p*log(n)), p)
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = 2*sqrt(log(exp(1)*p*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
  if(debug){
    print(ss)
  }
  
  if(is.null(thresholds_test)){
    thresholds_test = nu_as[] 
    thresholds_test[2:length(ss)] = threshold_s_test*(ss[2:length(ss)]*log(exp(1)*p*log(n)/ss[2:length(ss)]^2)+ 2*log(n))
    thresholds_test[1] = threshold_d_test* (sqrt(p*log(n)) + log(n))
    
    if(debug){
      print(thresholds_test)
    }
  }
  
  
  thresholds = nu_as[]
  
  thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n)/ss[2:length(ss)]^2)+ 2*log(n))
  #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)-1 + sqrt(p*log(n^4))/ss[2:length(ss)]), log(n^4))
  #thresholds[2:length(ss)] = threshold_s*(ss[2:length(ss)]*log(exp(1)-1 + sqrt(p*log(n^4))/ss[2:length(ss)])+ log(n^4))
  
  #thresholds[2:length(ss)] = threshold_s*(sqrt(p*exp(-as[2:length(ss)]^2/2)*log(n^4))+ log(n^4))
  
  
  #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n*log(sqrt(p*log(n)),2))/ss[2:length(ss)]^2), log(n*log(sqrt(p*log(n)),2)))
  thresholds[1] = threshold_d * (sqrt(p*log(n)) + log(n))
  
  #ind = thresholds > thresholds[1]
  #thresholds = thresholds[!ind]
  #thresholds_test = thresholds_test[!ind]
  #as = as[!ind]
  #nu_as = nu_as[!ind]
  #thresholds_test = threshold_test_const* thresholds
  #ts = ss[!ind]
  ts = ss
  twologn = min(p, floor(log(n)/2))
  if(droppartialsum){
    twologn=0
  }else if(!empirical){
    tt = 1
    count = 0
    h = length(thresholds_test)
    while(tt <=twologn){
      thresholds_test[h - count] = thresholds_test[h-count] + ts[h-count] * as[h-count]*as[h-count]  
      count = count+1
      tt = 2*tt
    }
    
  }
  if(debug){
    print(thresholds)
    print(thresholds_test)  
  }
  
  #print(twologn)
  #lens[7]=round(log(n)^2)
  #myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
  # SEXP epsI, SEXP lambdaI, SEXP maxiterI)
  
  res = .Call(cHDCD, X,as.integer(n), as.integer(p), thresholds,thresholds_test,
              as.integer(lens),as.integer(length(lens)), as.integer(K), as,nu_as, 
              as.integer(length(as)), as.integer(twologn), as.integer(ts), as.integer(fast),as.integer(debug))
  
  
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
  res$coordinate = (res$coordinate[,srt_indices])
  res$startpoints = as.integer(res$startpoints[srt_indices])
  res$endpoints = as.integer(res$endpoints[srt_indices])
  res$maxaposes = as.integer(res$maxaposes[srt_indices])
  res$s = as.integer(ss[res$maxaposes+1])
  res$thresholds = thresholds
  res$ts = ts
  
  
  return(res)
}

#' @export
HDCD_modified_calibrate = function(n,p, as=NULL, nu_as=NULL, ts=NULL, twologn=NULL, lens = NULL,
                          alpha = 1+1/6, K = 7, N, tol, fast = FALSE,debug=FALSE){
  if(is.null(as) || is.null(nu_as) || is.null(ts) || is.null(twologn) ||
     is.null(lens)){
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
    max_s = min(sqrt(p*log(n)), p)
    log2ss = 0:floor(log(max_s, base=2))
    ss = 2^(log2ss)
    ss = c(p, rev(ss))
    as = ss[]
    as[2:length(ss)] = sqrt(2*log(exp(1)*p*log(n)/ss[2:length(ss)]^2))
    as[1] = 0
    nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))
    twologn = min(p, floor(2*log(n)))
    ts = ss[]
  }
  if(debug){
    print(n)
    print(p)
    print(N)
    print(tol)
  }
  
  toln = max(round(N*tol),1)
  
  
  res= .Call(cHDCD_calibrate, as.integer(n), as.integer(p), as.integer(N),
             as.integer(toln), as.integer(lens), as.integer(length(lens)),
             as.integer(K), as.numeric(as), as.numeric(nu_as), as.integer(length(as)),
             as.integer(twologn), as.integer(ts), as.integer(fast), as.integer(debug))
  
  return(res)
  #return(NULL)
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

#' @export
HDCD_thresholds = function(n,p){
  
  ## FIX THE AUTOMATIC MONTE CARLO HERE. tol and N are missing
  
  
  max_s = min(sqrt(p*log(n)), p)
  log2ss = 0:floor(log(max_s, base=2))
  ss = 2^(log2ss)
  ss = c(p, rev(ss))
  as = ss[]
  as[2:length(ss)] = sqrt(2*log(exp(1)*p*4*log(n)/ss[2:length(ss)]^2))
  as[1] = 0
  nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))

  
  
  thresholds = nu_as[]
  
  thresholds[2:length(ss)] = (ss[2:length(ss)]*log(exp(1)*4*p*log(n)/ss[2:length(ss)]^2)+ log(n^4))
  #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)-1 + sqrt(p*log(n^4))/ss[2:length(ss)]), log(n^4))
  #thresholds[2:length(ss)] = threshold_s*(ss[2:length(ss)]*log(exp(1)-1 + sqrt(p*log(n^4))/ss[2:length(ss)])+ log(n^4))
  
  #thresholds[2:length(ss)] = threshold_s*(sqrt(p*exp(-as[2:length(ss)]^2/2)*log(n^4))+ log(n^4))
  
  
  #thresholds[2:length(ss)] = threshold_s*pmax(ss[2:length(ss)]*log(exp(1)*p*log(n*log(sqrt(p*log(n)),2))/ss[2:length(ss)]^2), log(n*log(sqrt(p*log(n)),2)))
  thresholds[1] =  (sqrt(p*log(n^4)) + log(n^4))
  
  twologn = min(p, floor(log(n)/2))
  tt = 1
  count = 0
  h = length(thresholds) +1
  while(tt <=twologn){
    #thresholds_test[h - count] = thresholds_test[h-count] + ts[h-count] * as[h-count]*as[h-count]  
    h = h-1
    tt = 2*tt
  }
  

  
  return(list(thresholds, h))
}


#' @export
rescale.variance <- function(x){
  p <- dim(x)[1]
  n <- dim(x)[2]
  for (j in 1:p){
    scale <- mad(diff(x[j,]))/sqrt(2)
    x[j,] <- x[j,] / scale
  }
  return(x)
}

#' @export
rescale.variance.scales <- function(x){
  p <- dim(x)[1]
  n <- dim(x)[2]
  scales = rep(NA, p)
  for (j in 1:p){
    scales[j] <- mad(diff(x[j,]))/sqrt(2)
    #x[j,] <- x[j,] / scale
  }
  return(scales)
}