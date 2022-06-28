
#' @export
hausdorff = function(v1, v2,n){
  if(is.null(v1)){
    if(!is.null(v2)){
      return(max(c(max(v2), max(n-v2))))
    }
    else{
      return(0)
    }
  }
  if(is.null(v2)){
    if(!is.null(v1)){
      return(max(c(max(v1), max(n-v1))))
    }
    else{
      return(0)
    }
  }
  # if(is.null(v1) | is.null(v2)){
  #   return(NULL)
  # }
  max = 0
  for (v  in v1) {
    tmp = min(abs(v-v2))
    if(tmp>max){
      max=tmp
    }
  }
  for (v  in v2) {
    tmp = min(abs(v-v1))
    if(tmp>max){
      max=tmp
    }
  }
  return(max)
}
#' @import mclust
#' @export
ARI = function(etas, eta_hats, n){
  
  # create label for each data point
  seg_etas = rep(1, n)
  counter = 1
  if(!is.null(etas)){
    for (i in 1:length(etas)) {
      counter = counter+1
      seg_etas[(etas[i]+1):n] = counter
    }
    
  }
  
  seg_eta_hats = rep(1, n)
  counter = 1
  if(!is.null(eta_hats)){
    for (i in 1:length(eta_hats)) {
      counter = counter+1
      seg_eta_hats[(eta_hats[i]+1):n] = counter
    }
    
  }
  return(mclust::adjustedRandIndex(seg_etas, seg_eta_hats))
  
}

hausdorff_avg = function(v1, v2){
  dists = rep(NA, length(v1))
  for (i in 1:length(v1)) {
    v = v1[i]
    dists[i] = min(abs(v-v2))
    
  }
  
  return(mean(dists))
}