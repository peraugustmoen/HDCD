library('InspectChangepoint')

############################################
##############inspect.AMOC##################
############################################

#inspect.AMOC carries out the base procedure for inspect
#assuming that there is at most one changepoint within 
#the system. As for the usual inspect, the best projection
#direction is found, a candidate changepoint location 
#computed, jointly with a test statistic, and we then 
#decide whether to accept this as a change location based
#on a threshold, which is calculated using simulations
#from the null.

inspect.amoc = function(Y,thresh){#Y is the data matrix with each column corresponding to observations at a time point, and each row corresponding to the observations of a variate. thresh is the penalty for admitted a changepoint to the model.
  n=dim(Y)[2]
  p=dim(Y)[1]
  
  #thresh <- compute.threshold(n,p,nrep=100) #computing threshold
  
  values <- locate.change(Y)#, lambda = log(log(n)*p)) #changed default lambda for sanity check
  
  if(values[[2]]>thresh){
    return(list(threshold=thresh,changepoint=as.numeric(values[1]),cost=10000,affected=rep(1,p)))
  }
  else{
    return(list(threshold=thresh,changepoint=NULL,cost=10000,affected=NULL))
  }
}

