#' @useDynLib HDCD cPilliat_test
#' @useDynLib HDCD cPilliat_test_calibrate
#' @useDynLib HDCD sort_test
#' @useDynLib HDCD partialsum_test
#' @export
changepoint_test_Pilliat = function(X, empirical =FALSE, N = 100, tol = 0.05,
                                    thresholds_partial=NULL, threshold_dense = NULL,
                                    thresholds_bj = NULL,
                                    threshold_d_const=4, threshold_bj_const=6, threshold_partial_const=4,
                                    rescale_variance = TRUE,fast = FALSE,
                                    debug =FALSE){
  p = dim(X)[1]
  n = dim(X)[2]
  
  log2p = floor(log(p,base=2))+1
  
  xx = 1
  maxx = 1
  while(TRUE){
    #p0 = threshold_bj_const/(pi^2*xx^2*n^2)
    logp0 = log(threshold_bj_const) - 2*log(pi) - 2*log(xx) - 2*log(n)
    #qq = qbinom(1/n, p, p0,lower.tail=FALSE)
    qq = qbinom(p = logp0, size = p, prob = exp(log(2) + pnorm(xx,lower.tail = FALSE,log = TRUE)),
                lower.tail=FALSE,log.p = TRUE)

    if(qq==0){
      maxx = xx-1
      break
    }
    xx = xx+1
  }
  
  if(empirical){
    if(is.null(thresholds_partial) || is.null(threshold_dense) || is.null(threshold_dense)){
      res = changepoint_test_Pilliat_calibrate(n,p, N=N, tol=tol,threshold_bj_const=threshold_bj_const,debug=debug)
      thresholds_partial = res$thresholds_partial
      thresholds_bj = res$thresholds_bj
      threshold_dense = res$threshold_dense
    }
  }else{
    
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
  }
  
  if(debug){
    print(threshold_dense)
    print(thresholds_partial)
    print(thresholds_bj)
  }
  
  
  #res = .Call(cPilliat, X,as.integer(n), as.integer(p), thresholds,thresholds_test,
  #           as.integer(lens),as.integer(length(lens)), as.integer(K), as,nu_as, 
  #          as.integer(length(as)), as.integer(twologn), as.integer(ts),as.integer(debug))
  #res = .Call(cPilliat, X, as.integer(n), as.integer(p), as.numeric(thresholds_partial), 
  #            as.numeric(threshold_dense), as.integer( thresholds_bj),
  #            as.integer(lens), as.integer(length(lens)), as.integer(K), as.integer(maxx),as.integer(debug))
  
  res = .Call(cPilliat_test, as.numeric(X),as.integer(n),as.integer(p),
              as.numeric(thresholds_partial), as.numeric(threshold_dense), 
              as.integer(thresholds_bj),as.integer(maxx), as.integer(rescale_variance),
              as.integer(fast),as.integer(debug))

  return(res)
}



#' @export
changepoint_test_Pilliat_calibrate = function(n,p, N=100, tol=0.05,threshold_bj_const=6,
                                              rescale_variance = TRUE,fast = FALSE,debug=FALSE){

  log2p = floor(log(p,base=2))+1
  
  xx = 1
  maxx = 1
  while(TRUE){
    #p0 = threshold_bj_const/(pi^2*xx^2*n^2)
    logp0 = log(threshold_bj_const) - 2*log(pi) - 2*log(xx) - 2*log(n)
    #qq = qbinom(1/n, p, p0,lower.tail=FALSE)
    qq = qbinom(p = logp0, size = p, prob = exp(log(2) + pnorm(xx,lower.tail = FALSE,log = TRUE)),
                lower.tail=FALSE,log.p = TRUE)
    
    if(qq==0){
      maxx = xx-1
      break
    }
    xx = xx+1
  }
  
  toln = max(round(N*tol),1)

  
  res= .Call(cPilliat_test_calibrate, as.integer(n), as.integer(p), as.integer(N), 
             as.integer(toln), as.integer(maxx), as.integer(log2p),as.integer(rescale_variance),
             as.integer(fast),as.integer(debug))
  #res[["as"]] = as
  #res[["nu_as"]] = nu_as
  return(res)
}





