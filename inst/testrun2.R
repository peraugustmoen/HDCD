library(HDCD)

#multiple change-points:
p = 1000
n =500

#calibrates = 
dense_const = 8
sparse_const = 8
 
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

calibrates = readRDS("/Users/peraugust/Library/CloudStorage/OneDrive-UniversitetetiOslo/project1/simulations/Simulations_HDCD/2022-11-16--11.56.35/multi/calibrates.RDA")

res3 = HDCD(X, alpha = 2, K=4,fast = TRUE, rescale_variance = FALSE, 
            thresholds_test = (calibrates[[5]])[[2]])
etas
res3$changepoints

res4 = HDCD(X, alpha = 1.5, K=5,fast = FALSE,rescale_variance=FALSE,droppartialsum = TRUE,
            thresholds_test = (calibrates[[5]])[[10]])
etas
res4$changepoints
 
source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")
res5 = change_main(X, SUBSET.normal,100, penalties=(calibrates[[5]])[[9]])
etas
sort(res5[[2]])
