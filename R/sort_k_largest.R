#' @useDynLib HDCD sort_k_largest_R
#' @export
sort_k_largest = function(vec, k,start, stop){
  return(.Call(sort_k_largest_R, as.numeric(vec), as.integer(k), as.integer(start), as.integer(stop)))
}


#' @useDynLib HDCD partial_quicksort_R
#' @export
partial_quicksort = function(vec, k,len){
  return(.Call(partial_quicksort_R, as.numeric(vec), as.integer(k), as.integer(len)))
}
