library(HDCD)

#multiple change-points:
p = 1000
n =200

#calibrates = 
dense_const = 3
sparse_const = 3

mus = matrix(0, nrow=p, ncol=n)
etas = sort(sample(1:(n-1), 5))
sparsities = sample(c(1,2), length(etas), replace=TRUE)
for (j in 1:length(etas)) {
  Delta = 1
  if(j==1){
    Delta = min(etas[1], etas[2]-etas[1])
  }else if(j<length(etas)){
    Delta = min(etas[j] - etas[j-1], etas[j+1] - etas[j])
  }else{
    Delta = min(etas[j]-etas[j-1], n -etas[j])
  }
  if(sparsities[j]==1){
    # sparse
    sparsities[j] = sample(1:floor(sqrt(p*log(n))),1)
    k = sparsities[j]
    phi = sparse_const/sqrt(Delta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
  }else{
    #dense
    sparsities[j] = sample(ceiling(sqrt(p*log(n))):p,1)
    k = sparsities[j]
    phi = dense_const/sqrt(Delta)*(p*log(n))^(1/4)
  }
  coords = sample(1:p, sparsities[j])
  #diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
  diff = sample(c(-1,1), sparsities[j], replace=TRUE)
  diff = diff/norm(diff, type="2")
  mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
}

#set.seed(100)
X= matrix(rnorm(n*p), ncol = n) + mus


# theoretical lim
res = HDCD(X, alpha = 2, K=4,fast = TRUE, rescale_variance = FALSE)
etas
res$changepoints

res2 = HDCD(X, alpha = 1.5, K=5,fast = FALSE,rescale_variance=FALSE)
etas
res$changepoints

# empirical 
calibrated = TRUE
if(!calibrated){
  rr = HDCD_calibrate (n,p,alpha = 2, K = 4, N=1000, tol=1/1000, fast = TRUE,rescale_variance = FALSE,debug=TRUE)
  rr2 = HDCD_calibrate (n,p,alpha = 1.5, K = 5, N=1000, tol=1/1000, fast = FALSE,rescale_variance = FALSE,debug=TRUE)
  
  # source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")
  # 
  # Ncal = 1000
  # empirical_penalties = rep(NA, Ncal)
  # 
  # for (i in 1:Ncal) {
  #   mynulldata <- matrix(rnorm(n*p,0,1),nrow=p,ncol=n,byrow=FALSE) # 5 variates with 1000 time points, no change
  #   empirical_penalties[i] = wbs_penaltyfinder(mynulldata, SUBSET.normal_penalty, 100)
  #   print(i)
  # }
  # 
  # 
  # empirical_penalties = sort(empirical_penalties,decreasing = TRUE)
}

sink(file = "/Users/peraugust/Library/CloudStorage/OneDrive-UniversitetetiOslo/project1/simulations/Simulations_HDCD/sinks/sink.txt")
res3 = HDCD(X, alpha = 2, K=4,fast = TRUE, rescale_variance = FALSE, 
            thresholds_test = rr[[2]], debug=TRUE)

sink()

etas
res3$changepoints

res4 = HDCD(X, alpha = 2, K=4,fast = TRUE,rescale_variance=FALSE,droppartialsum = TRUE,
            thresholds_test = rr[[1]])
etas
res4$changepoints
sparsities

res5 = HDCD(X, alpha = 1.2, K=5,fast = FALSE,rescale_variance=FALSE,droppartialsum = TRUE,
            thresholds_test = rr2[[1]])
etas
res5$changepoints
sparsities

source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")
res6 = change_main(X, SUBSET.normal,100, penalties=(calibrates[[2]])[[9]])
etas
sort(res6[[2]])



# short simu study: 
N = 300
rez1 = 0
rez2 = 0
rez3 = 0
rez4 = 0
rez5 = 0
rez6 = 0

kerr1 = 0
kerr2 = 0
kerr3 = 0
kerr4 = 0
kerr5 = 0
kerr6 = 0

for (i in 1:N) {
  print(i)
  mus = matrix(0, nrow=p, ncol=n)
  etas = sort(sample(1:(n-1), 5))
  sparsities = sample(c(1,2), length(etas), replace=TRUE)
  for (j in 1:length(etas)) {
    Delta = 1
    if(j==1){
      Delta = min(etas[1], etas[2]-etas[1])
    }else if(j<length(etas)){
      Delta = min(etas[j] - etas[j-1], etas[j+1] - etas[j])
    }else{
      Delta = min(etas[j]-etas[j-1], n -etas[j])
    }
    if(sparsities[j]==1){
      # sparse
      sparsities[j] = sample(1:floor(sqrt(p*log(n))),1)
      k = sparsities[j]
      phi = sparse_const/sqrt(Delta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
    }else{
      #dense
      sparsities[j] = sample(ceiling(sqrt(p*log(n))):p,1)
      k = sparsities[j]
      phi = dense_const/sqrt(Delta)*(p*log(n))^(1/4)
    }
    coords = sample(1:p, sparsities[j])
    #diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
    diff = sample(c(-1,1), sparsities[j], replace=TRUE)
    diff = diff/norm(diff, type="2")
    mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
  }
  
  #set.seed(100)
  X= matrix(rnorm(n*p), ncol = n) + mus
  
  
  # theoretical lim
  res = HDCD(X, alpha = 2, K=4,fast = TRUE, rescale_variance = FALSE)
  etas
  #res$changepoints
  
  rez1 = rez1 + hausdorff(res$changepoints, etas,n) / N
  kerr1 = kerr1 +abs(length(etas)-length(res$changepoints))/ N
  
  res2 = HDCD(X, alpha = 1.5, K=5,fast = FALSE,rescale_variance=FALSE)
  etas
  #res$changepoints
  
  rez2 = rez2 + hausdorff(res2$changepoints, etas,n) / N
  kerr2 = kerr2 +abs(length(etas)-length(res2$changepoints))/ N
  # empirical 
  calibrated = TRUE
  if(!calibrated){
    rr = HDCD_calibrate (n,p,alpha = 2, K = 4, N=1000, tol=1/1000, fast = TRUE,rescale_variance = FALSE,debug=TRUE)
    rr2 = HDCD_calibrate (n,p,alpha = 1.5, K = 5, N=1000, tol=1/1000, fast = FALSE,rescale_variance = FALSE,debug=TRUE)
    
    # source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")
    # 
    # Ncal = 1000
    # empirical_penalties = rep(NA, Ncal)
    # 
    # for (i in 1:Ncal) {
    #   mynulldata <- matrix(rnorm(n*p,0,1),nrow=p,ncol=n,byrow=FALSE) # 5 variates with 1000 time points, no change
    #   empirical_penalties[i] = wbs_penaltyfinder(mynulldata, SUBSET.normal_penalty, 100)
    #   print(i)
    # }
    # 
    # 
    # empirical_penalties = sort(empirical_penalties,decreasing = TRUE)
  }
  
  #sink(file = "/Users/peraugust/Library/CloudStorage/OneDrive-UniversitetetiOslo/project1/simulations/Simulations_HDCD/sinks/sink.txt")
  res3 = HDCD(X, alpha = 2, K=4,fast = TRUE, rescale_variance = FALSE, 
              thresholds_test = rr[[2]], debug=FALSE)
  
  #sink()
  
  etas
  res3$changepoints
  
  rez3 = rez3 + hausdorff(res3$changepoints, etas,n) / N
  kerr3 = kerr3 +abs(length(etas)-length(res3$changepoints))/ N
  
  res4 = HDCD(X, alpha = 2, K=4,fast = TRUE,rescale_variance=FALSE,droppartialsum = TRUE,
              thresholds_test = rr[[1]])
  etas
  res4$changepoints
  rez4 = rez4 + hausdorff(res4$changepoints, etas,n) / N
  kerr4 = kerr4 +abs(length(etas)-length(res4$changepoints))/ N
  sparsities
  
  res5 = HDCD(X, alpha = 1.2, K=5,fast = FALSE,rescale_variance=FALSE,droppartialsum = TRUE,
              thresholds_test = rr2[[1]])
  etas
  res5$changepoints
  sparsities
  
  rez5 = rez5 + hausdorff(res5$changepoints, etas,n) / N
  kerr5 = kerr5 +abs(length(etas)-length(res5$changepoints))/ N
  
  source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")
  res6 = change_main(X, SUBSET.normal,100, penalties=(calibrates[[2]])[[9]])
  etas
  sort(res6[[2]])
  rez6 = rez6 + hausdorff(sort(res6[[2]]), etas,n) / N
  kerr6 = kerr6 +abs(length(etas)-length(sort(res6[[2]])))/ N
}
rez1
rez2
rez3
rez4
rez5
rez6

kerr1
kerr2
kerr3
kerr4
kerr5
kerr6
