inspect.amoc_penalty = function(Y){#Y is the data matrix - each row corresponds to the set of observations made for a variate, each column to the set of observations at a given point in time.
  #thresh <- compute.threshold(n,p,nrep=100) #computing threshold

  n <- dim(Y)[2]
  d <- dim(Y)[1]
  
  values <- locate.change(Y)#, lambda=log(log(n)*d)) #changed from default lambda for sanity check
  
  return(list("Empty Slot","Another Empty Slot",values[[2]]))
  }
