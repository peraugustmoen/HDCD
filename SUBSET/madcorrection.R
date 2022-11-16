madcorrection <- function(vector){#vector is the univariate time series in question
  return(vector/(sd(vector)))
}
