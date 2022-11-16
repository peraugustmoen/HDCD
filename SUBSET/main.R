source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/BinWeight.R')
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/BinWeight_penalty.R')
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/inspect_AMOC.R')
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/inspect_AMOC_penalty.R')
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/madcorrection.R')
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/Max.R')
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/Max_penalty.R')
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/Mean.R')
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/Mean_penalty.R')
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/SUBSET_negbin.R')
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/SUBSET_negbin_penalty.R')
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/SUBSET_normal.R')
source('/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/SUBSET_normal_penalty.R')


#encases a multivariate test within the wild binary segmentation procedure for detecting multiple changepoints

change_main <- function(data,method,int_reps,penalties){#data - matrix corresponding to the data inputted, such that (i,j)th entry corresponds to the jth entry in the ith variate; method corresponds to the name of the test statistic to be applied (e.g. MV.single.CUSUM.max); int_reps is the number of random intervals to be simulated; penalties is a vector of length two - the first entry corresponds to beta and the second entry corresponds to alpha. 

#Get intervals of interest

new.number <- int_reps+int_reps

n <- dim(data)[2]
interval.points <- sample(1:n,new.number,replace=TRUE)
interval.mat <- matrix(interval.points,nrow=int_reps)
interval.mat <- t(apply(interval.mat,1,sort))

results <- wbs_act(data,method,interval.mat,penalties)

return(results)

}

wbs_act <- function(data,method,interval.mat,penalties){#data - matrix corresponding to the data inputted, such that (i,j)th entry corresponds to the jth entry in the ith variate; method corresponds to the name of the test statistic to be applied (e.g. MV.single.CUSUM.max); interval.mat is the set of random intervals entered as a matrix - the first column corresponds to the beginnings of the intervals, the second column corresponds to the ends; penalties is a vector of length two - the first entry corresponds to beta and the second entry corresponds to alpha. 

active.intervals <- matrix(c(1,dim(data)[2]),ncol=2) #start off with the interval c(1,n)

#the.number <- dim(active.intervals)[1]

current.result <- list(0,NULL,NULL)

while(dim(active.intervals)[1]>0.5){ #while there is still an active interval
if(active.intervals[dim(active.intervals)[1],2]-active.intervals[dim(active.intervals)[1],1]<0.5){
 if(dim(active.intervals)[1]>2.5){
	active.intervals <- active.intervals[-dim(active.intervals)[1],]
 }
 else if(dim(active.intervals)[1]>1.5){
	active.intervals <- t(as.matrix(active.intervals[-dim(active.intervals)[1],],nrow=1))
 }
 else{
	return(current.result)
 }
}
active.result <- do_wbs(data,method,interval.mat,penalties,active.intervals[dim(active.intervals)[1],]) #runs method on interval listed last in the list of active intervals
if(length(active.result[[2]])>0.5){#changepoint found
	value <- active.result[[2]]
  active.intervals <- rbind(c(active.intervals[dim(active.intervals)[1],1],value),c(value+1,active.intervals[dim(active.intervals)[1],2]),active.intervals) #splits up previous active interval into two about change, and places these at the top of the list of active intervals
	active.intervals <- active.intervals[-dim(active.intervals)[1],] #removes the previous interval from consideration
        if(current.result[[1]][length(current.result[[1]])]<0.0000000001){
		current.result[[1]] <- active.result[[1]]
		current.result[[2]] <- active.result[[2]]
		current.result[[3]] <- as.matrix(active.result[[3]],ncol=1)
	}
	else{
		current.result[[1]] <- c(current.result[[1]],active.result[[1]])
		critical.number <- sum(current.result[[2]]<active.result[[2]])
		current.result[[2]] <- c(current.result[[2]],active.result[[2]])
		if(critical.number<0.5){
			current.result[[3]] <- cbind(active.result[[3]],current.result[[3]])
                }
		else if(critical.number>((dim((current.result)[[3]])[2])-0.5)){
			current.result[[3]] <- cbind(current.result[[3]],active.result[[3]])
		}
		else{
			current.result[[3]] <- cbind(current.result[[3]][,1:critical.number],active.result[[3]],current.result[[3]][,(critical.number+1):(dim((current.result)[[3]])[2])])
		}
	}
}
else if(dim(active.intervals)[1]<1.5){
	return(current.result)
}
else if(dim(active.intervals)[1]<2.5){
	active.intervals <- t(as.matrix(active.intervals[-dim(active.intervals)[1],],nrow=1))
}
else{
	active.intervals <- active.intervals[-dim(active.intervals)[1],]
}
}
return(current.result)

}

do_wbs <- function(data,method,interval.mat,penalties,current.interval){#data - matrix corresponding to the data inputted, such that (i,j)th entry corresponds to the jth entry in the ith variate; method corresponds to the name of the test statistic to be applied (e.g. MV.single.CUSUM.max); interval.mat is the set of random intervals entered as a matrix - the first column corresponds to the beginnings of the intervals, the second column corresponds to the ends; penalties is a vector of length two - the first entry corresponds to beta and the second entry corresponds to alpha. current.interval is a vector of length 2 indicating the start and end points of the current analysis of interest.

relevant.intervals <- rbind(current.interval,interval.mat)

A <- dim(relevant.intervals)[1]

for(j in A:2){
	if(current.interval[1]>relevant.intervals[j,1]+0.5){
		relevant.intervals <- relevant.intervals[-j,]
	}
	else if(current.interval[2]<relevant.intervals[j,2]-0.5){
		relevant.intervals <- relevant.intervals[-j,]
	}
	else{}
}

B <- dim(relevant.intervals)[1]

if(length(B)<0.5){
	relevant.intervals<-matrix(relevant.intervals,ncol=2)
}

current.max <- 0
changepoint <- NULL
affected    <- NULL

if(length(B)<0.5){
	return(list(current.max,changepoint,affected))
}

for(i in 1:B){
  if(relevant.intervals[i,2]-relevant.intervals[i,1]<0.5){
  }
  else{
  new.result <- method(data[,relevant.intervals[i,1]:relevant.intervals[i,2]],penalties)
  value <- new.result[[3]]
  if(length(new.result[[2]])){
  	if(value > current.max){
		current.max <- value
        	changepoint <- current.interval[1]+new.result[[2]]-1
		affected <- new.result[[4]]
  	}
        else{}
  }
  else{}
  }
}

the.final.list <- list(current.max,changepoint,affected)

return(the.final.list)

}

wbs_penaltyfinder <- function(data,method,int_reps){#data - matrix corresponding to the data inputted, such that (i,j)th entry corresponds to the jth entry in the ith variate; method corresponds to the name of the test statistic to be applied (e.g. MV.single.CUSUM.max); penalties is a vector of length two - the first entry corresponds to beta and the second entry corresponds to alpha. 

n <- dim(data)[2]
interval.points <- sample(1:n,2*int_reps,replace=TRUE)
interval.mat <- matrix(interval.points,nrow=int_reps)
interval.mat <- t(apply(interval.mat,1,sort))

relevant.intervals <- interval.mat

B <- dim(relevant.intervals)[1]

current.max <- 0
changepoint <- NULL
affected    <- NULL

	for(i in 1:B){
	   if(relevant.intervals[i,2]-relevant.intervals[i,1]>0.5){
	    new.result <- method(data[,relevant.intervals[i,1]:relevant.intervals[i,2]])
  		value <- new.result[[3]]
  		if(value > current.max){
			current.max <- value
  		}
        	else{}
           }
           else{}
	}
return(current.max)
}
