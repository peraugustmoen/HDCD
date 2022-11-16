#### Simulation for multiple change-points

library(doSNOW)
library(HDCD)
library(foreach)
## same magnitude in all coords

maindir = "/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/Simulations_HDCD"
dateandtime = gsub(" ", "--",as.character(Sys.time()))
dateandtime = gsub(":", ".", dateandtime)
savedir = file.path(maindir, dateandtime)


save = TRUE

if(save){
  dir.create(savedir, showWarnings = FALSE)
  savedir = file.path(maindir, sprintf("%s/multi",dateandtime))
  dir.create(savedir, showWarnings = FALSE)
  
}



N = 12
num_cores = 6
sparse_const = 6
dense_const = 6
set.seed(1996)
rescale_variance = TRUE
Ncal = 12
tol = 1/Ncal
#tol = 0.001

#ns = c(500,1000,2000)
#ns = c(100,500,1000)
#ns = c(500,1000,2000)
#ns = c(100,500,1000)
#ps = c(1000,2000,10000)
#ns = c(100,500)
#ps = ns[]
#ps = c(500,1000,2000,5000)
#kvals=4
#totruns = length(ns)*length(ps)*kvals
#ns = c(100,200,500,2000)
ns = c(200,500)
#ns = c(500)
#ps = c(1000)
#ps = c(50,100,500,1000,4000)
ps = c(50,100,1000)
#ps = c(50,100,500,5000)
#ps = ns[]


# rez = config(7, n, p, mus)
# etas = rez[[1]]
# sparsities = rez[[2]]
# mus = rez[[3]]



# calibrate

cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = ((length(ns)*length(ps))), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
calibrates = foreach(z = 0:((length(ns)*length(ps))-1),.options.snow = opts) %dopar% {
  #  for(z in 0:((length(ns)*length(ps))-1)){
  library(HDCD)
  set.seed(300*z)
  rez = list()
  nind = floor(z/length(ps))+1
  pind = z%%length(ps)+1
  cc = HDCD_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,K=4, alpha = 2, fast = TRUE,
                      rescale_variance = rescale_variance, debug=FALSE)
  # cc2 = HDCD_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,debug=FALSE)
  lambda = sqrt(log(ps[pind]*log(ns[nind]))/2)
  cc2 = Inspect_calibrate(n = ns[nind], p = ps[pind], N=Ncal, tol=tol,lambda = lambda , alpha = 2, K = 4,eps=1e-10,
                          maxiter=10000,rescale_variance = rescale_variance, debug =FALSE)
  
  cc3 = Pilliat_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,K = 4, 
                          rescale_variance = rescale_variance, debug=FALSE)
  
  source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")
  
  
  empirical_penalties = rep(NA, Ncal)
  if(rescale_variance){
    for (i in 1:Ncal) {
      mynulldata <- matrix(rnorm(ns[nind]*ps[pind],0,1),nrow=ps[pind],ncol=ns[nind],byrow=FALSE) # 5 variates with 1000 time points, no change
      empirical_penalties[i] = wbs_penaltyfinder(InspectChangepoint::rescale.variance(mynulldata), SUBSET.normal_penalty, 100)
    }
  }else{
    for (i in 1:Ncal) {
      mynulldata <- matrix(rnorm(ns[nind]*ps[pind],0,1),nrow=ps[pind],ncol=ns[nind],byrow=FALSE) # 5 variates with 1000 time points, no change
      empirical_penalties[i] = wbs_penaltyfinder(mynulldata, SUBSET.normal_penalty, 100)
    }
  }
  
  empirical_penalties = sort(empirical_penalties,decreasing = TRUE)
  
  
  rez[[1]] = cc[[1]] #FAST threshoolds
  rez[[2]] = cc[[2]] #FAST thresholds
  rez[[3]] = cc2[[1]] #Inspect xi
  rez[[4]] = cc3[[1]] # pilliat
  rez[[5]] = cc3[[2]] # pilliat
  rez[[6]] = cc3[[3]] # pilliat
  rez[[7]] = ns[nind]
  rez[[8]] = ps[pind]
  rez[[9]] = empirical_penalties[max(1,round(tol*Ncal))] # subset penalty
  rez
  #calibrates[[z+1]] = rez
  
}
close(pb)
stopCluster(cl) 


numconfig = 7
config = function(i,n,p){
  #mus = matrix(0, nrow=p, ncol=n)
  mus = matrix(0, nrow=p, ncol=n)
  #k = 0
  etas = c()
  sparsities = c()
  sparsity = c()
  
  if(i==1){
    # no chgppts
  }
  else if(i==2){
    # two dense
    sparsity = "Dense"
    #etas = round(c(n/4, n*2/4, n*3/4))
    etas = round(c(n/3, n*2/3))
    #sparsities = c(p, round(p/2), round(p^(4/5)))
    sparsities = c(p, round(p^(4/5)))
    Delta = min(diff(c(0, etas, n)))
    
    for (j in 1:length(etas)) {
      phi = dense_const/sqrt(Delta)*(p*log(n))^(1/4)
      k = sparsities[j]
      coords = sample(1:p, sparsities[j]) 
      #diff = runif(k) *sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
    }
    
  }
  else if(i==3){
    # two sparse
    sparsity="Sparse"
    #etas = round(c(n/4, n*2/4, n*3/4))
    etas = round(c(n/3, n*2/3))
    #sparsities = round(c(1,p^(1/5), p^(2/5)))
    sparsities = round(c(1, p^(2/5)))
    Delta = min(diff(c(0, etas, n)))
    
    
    for (j in 1:length(etas)) {
      k = sparsities[j]
      phi = sparse_const/sqrt(Delta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
      coords = sample(1:p, sparsities[j])
      #diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
    }
    
  }
  else if(i==4){
    # two mixed
    sparsity = "Mixed"
    #etas = round(c(n/4, n*2/4, n*3/4))
    etas = round(c(n/3, n*2/3))
    Delta = min(diff(c(0, etas, n)))
    #sparsities = sample(c(1,2), 2, replace=TRUE)
    sparsities = c(1,2)
    phi = 0
    for (j in 1:length(etas)) {
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
    
  }
  else if(i==5){
    # 5 dense
    sparsity="Dense"
    etas = sort(sample(1:(n-1), 5))
    sparsities = sample(ceiling(sqrt(p*log(n))):p,length(etas), replace=TRUE)
    
    for (j in 1:length(etas)) {
      Delta = 1
      if(j==1){
        Delta = min(etas[1], etas[2]-etas[1])
      }else if(j<length(etas)){
        Delta = min(etas[j] - etas[j-1], etas[j+1] - etas[j])
      }else{
        Delta = min(etas[j]-etas[j-1], n -etas[j])
      }
      phi = dense_const/sqrt(Delta)*(p*log(n))^(1/4)
      k = sparsities[j]
      coords = sample(1:p, sparsities[j]) 
      #diff = runif(k) *sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] + matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
    }
  }
  else if(i==6){
    # 5 sparse
    sparsity = "Sparse"
    etas = sort(sample(1:(n-1), 5))
    sparsities = sample(1:floor(sqrt(p*log(n))),length(etas), replace=TRUE )
    
    
    for (j in 1:length(etas)) {
      Delta = 1
      if(j==1){
        Delta = min(etas[1], etas[2]-etas[1])
      }else if(j<length(etas)){
        Delta = min(etas[j] - etas[j-1], etas[j+1] - etas[j])
      }else{
        Delta = min(etas[j]-etas[j-1], n -etas[j])
      }
      k = sparsities[j]
      phi = sparse_const/sqrt(Delta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
      coords = sample(1:p, sparsities[j])
      #diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
    }
  }
  else if(i==7){
    # 5 mixed
    sparsity="Mixed"
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
  }
  return(list(etas, sparsities,mus,sparsity))
}


cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
#for(z in 1:N){
result = foreach(z = 1:N,.options.snow = opts) %dopar% {
  rez = list()
  set.seed(2*z)
  library(HDCD)
  library(hdbinseg)
  counter = 1
  source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")
  
  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      
      for (y in 1:numconfig) {
        noise = matrix(rnorm(n*p), nrow=p, ncol=n)
        conf = config(y, n, p)
        etas = conf[[1]]
        sparsities = conf[[2]]
        mus = conf[[3]]
        
        X = noise + mus
        
        
        rezi = list()
        rezi[["i"]] = i
        rezi[["j"]] = j
        rezi[["y"]] = y
        
        # hdcd_fast
        #xi = 4*sqrt(log(p*n))
        #lambda = sqrt(log(p*log(n)))
        lambda = sqrt(log(p*log(n))/2)
        
        
        a = proc.time()
        #res  = HDCD (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
        res = HDCD (X[,], 1.5,1, empirical=TRUE,alpha = 2, K = 4, thresholds_test = (calibrates[[j+(i-1)*length(ps)]])[[2]], droppartialsum = FALSE, fast =TRUE,
                    rescale_variance = rescale_variance, debug= FALSE)
        b=proc.time()
        rezi[["hdcd_fast_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["hdcd_fast_K"]]= res$changepointnumber
        rezi[["hdcd_fast_chgpts"]]= res$changepoints
        rezi[["hdcd_fast_hausd"]] = hausdorff(res$changepoints, etas,n)
        rezi[["hdcd_fast_K_error"]] = length(res$changepoints) - length(etas)
        rezi[["hdcd_fast_ari"]] = ARI(etas, res$changepoints, n)
        # hdcd
        
        # slow FAST:
        # a = proc.time()
        # res  = HDCD (X, 2,2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[4]] , droppartialsum = FALSE, fast =FALSE,debug= FALSE)
        # b=proc.time()
        # 
        # rezi[["hdcd_time"]] = (b-a)[1]+(b-a)[2]
        # rezi[["hdcd_K"]] = res$changepointnumber
        # rezi[["hdcd_chgpts"]]= res$changepoints
        # rezi[["hdcd_hausd"]] = hausdorff(res$changepoints, etas,n)
        # rezi[["hdcd_K_error"]] = length(res$changepoints) - length(etas)
        # rezi[["hdcd_ari"]] = ARI(etas, res$changepoints, n)
        # 
        # hdcd slow
        # if(dim(X)[2]>500){
        #   rezi[["pilliat_time"]] = NA
        #   rezi[["pilliat_K"]] = NA
        #   rezi[["pilliat_chgpts"]]= NA
        #   rezi[["pilliat_hausd"]] = NA
        #   rezi[["pilliat_K_error"]] = NA
        #   rezi[["pilliat_ari"]] = NA
        # }
        # else{
        #a = proc.time()
        #res = HDCD (X, 2,2, alpha = 1+1/6, K = 4, threshold_d_test = 1.5, 
        #            threshold_s_test = 1.5, droppartialsum = FALSE, fast =FALSE,debug= FALSE)
        #b=proc.time()
        
        a = proc.time()
        res  = Pilliat(X, K = 4, alpha = 2, empirical = TRUE, threshold_dense = (calibrates[[j+(i-1)*length(ps)]])[[5]], 
                       thresholds_partial = (calibrates[[j+(i-1)*length(ps)]])[[4]], thresholds_bj = (calibrates[[j+(i-1)*length(ps)]])[[6]],
                       rescale_variance = rescale_variance, debug =FALSE)
        b=proc.time()
        print(res)
        
        rezi[["pilliat_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["pilliat_K"]] = res$changepointnumber
        rezi[["pilliat_chgpts"]]= res$changepoints
        rezi[["pilliat_hausd"]] = hausdorff(res$changepoints, etas,n)
        rezi[["pilliat_K_error"]] = length(res$changepoints) - length(etas)
        rezi[["pilliat_ari"]] = ARI(etas, res$changepoints, n)  
        #}
        
        
        #inspect
        a = proc.time()
        #res  = HDCD (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
        res = Inspect(X[,], xi=(calibrates[[j+(i-1)*length(ps)]])[[3]], alpha = 2, K = 4,eps=1e-10,
                      lambda = lambda, maxiter=10000,
                      rescale_variance = rescale_variance, debug=FALSE)
        b=proc.time()
        rezi[["inspect_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["inspect_K"]]= res$changepointnumber
        rezi[["inspect_chgpts"]]= res$changepoints
        rezi[["inspect_hausd"]] = hausdorff(res$changepoints, etas,n)
        rezi[["inspect_K_error"]] = length(res$changepoints) - length(etas)
        rezi[["inspect_ari"]] = ARI(etas, res$changepoints, n)
        
        
        
        # subset
        a = proc.time()
        #res  = HDCD (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
        if(rescale_variance){
          res = change_main(InspectChangepoint::rescale.variance(X), SUBSET.normal, 100,penalties= (calibrates[[j+(i-1)*length(ps)]])[[9]])
        }else{
          res = change_main(X, SUBSET.normal,100, penalties=(calibrates[[j+(i-1)*length(ps)]])[[9]])
        }
        b=proc.time()
        rezi[["subset_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["subset_K"]]= length(res[[2]])
        rezi[["subset_chgpts"]]= sort(res[[2]])
        rezi[["subset_hausd"]] = hausdorff(sort(res[[2]]), etas,n)
        rezi[["subset_K_error"]] = length(res[[2]]) - length(etas)
        rezi[["subset_ari"]] = ARI(etas, sort(res[[2]]), n)
        
        
        
        # SBS
        a = proc.time()
        #res  = HDCD (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
        res = sbs.alg(X[,], cp.type = 1, thr = NULL, trim = NULL, height = NULL,
                      temporal = TRUE, scales = NULL, diag = FALSE, B = 100, q = 0.01,
                      do.parallel = 1)
        b=proc.time()
        if(length(res$ecp)==0){
          res$ecp = NULL
        }
        rezi[["sbs_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["sbs_K"]]= length(res$ecp)
        rezi[["sbs_chgpts"]]= sort(res$ecp)
        rezi[["sbs_hausd"]] = hausdorff(sort(res$ecp), etas,n)
        rezi[["sbs_K_error"]] = length(res$ecp) - length(etas)
        rezi[["sbs_ari"]] = ARI(etas, sort(res$ecp), n)
        
        # DC
        a = proc.time()
        res = dcbs.alg(X[,], cp.type = 1, phi = 0.5, thr = NULL, trim = NULL,
                       height = NULL, temporal = TRUE, scales = NULL, diag = FALSE,
                       B = 100, q = 0.01, do.parallel = 1)
        b=proc.time()
        if(length(res$ecp)==0){
          res$ecp = NULL
        }
        rezi[["dc_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["dc_K"]]= length(res$ecp)
        rezi[["dc_chgpts"]]= sort(res$ecp)
        rezi[["dc_hausd"]] = hausdorff(sort(res$ecp), etas,n)
        rezi[["dc_K_error"]] = length(res$ecp) - length(etas)
        rezi[["dc_ari"]] = ARI(etas, sort(res$ecp), n)
        
        rezi[["hehe"]] = res
        
        
        rezi[["true_K"]] = length(etas)
        rezi[["true_etas"]] = etas
        rezi[["true_sparsities"]] = sparsities
        
        rez[[counter]] = rezi
        counter = counter+1
        
        
        
      }
    }
  }
  rez
}
# bb = proc.time()
# print(bb-aa)
#parallel::stopCluster(cl)
close(pb)
stopCluster(cl) 

{
  hdcd_fast_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  hdcd_fast_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  hdcd_fast_time = array(0, dim= c(length(ns), length(ps), numconfig))
  hdcd_fast_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  inspect_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  inspect_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  inspect_time = array(0, dim= c(length(ns), length(ps), numconfig))
  inspect_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  pilliat_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  pilliat_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  pilliat_time = array(0, dim= c(length(ns), length(ps), numconfig))
  pilliat_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  subset_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  subset_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  subset_time = array(0, dim= c(length(ns), length(ps), numconfig))
  subset_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  sbs_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  sbs_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  sbs_time = array(0, dim= c(length(ns), length(ps), numconfig))
  sbs_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  dc_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  dc_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  dc_time = array(0, dim= c(length(ns), length(ps), numconfig))
  dc_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  pilliat_hausd[,,1] = NA
  hdcd_fast_hausd[,,1] = NA
  inspect_hausd[,,1] = NA
  subset_hausd[,,1] = NA
  sbs_hausd[,,1] = NA
  dc_hausd[,,1] = NA
  
  for (z in 1:N) {
    list = result[[z]]
    len = length(list)
    
    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      i = sublist[["i"]]
      j = sublist[["j"]]
      
      hdcd_fast_Kerr[i,j,y] = hdcd_fast_Kerr[i,j,y] + abs(sublist[["hdcd_fast_K_error"]])/N
      pilliat_Kerr[i,j,y] = pilliat_Kerr[i,j,y] + abs(sublist[["pilliat_K_error"]])/N
      inspect_Kerr[i,j,y] = inspect_Kerr[i,j,y] + abs(sublist[["inspect_K_error"]])/N
      subset_Kerr[i,j,y] = subset_Kerr[i,j,y] + abs(sublist[["subset_K_error"]])/N
      sbs_Kerr[i,j,y] = sbs_Kerr[i,j,y] + abs(sublist[["sbs_K_error"]])/N
      dc_Kerr[i,j,y] = dc_Kerr[i,j,y] + abs(sublist[["dc_K_error"]])/N
      
      hdcd_fast_time[i,j,y] = hdcd_fast_time[i,j,y] + sublist[["hdcd_fast_time"]]/N
      pilliat_time[i,j,y] = pilliat_time[i,j,y] + sublist[["pilliat_time"]]/N
      inspect_time[i,j,y] = inspect_time[i,j,y] + sublist[["inspect_time"]]/N
      subset_time[i,j,y] = subset_time[i,j,y] + sublist[["subset_time"]]/N
      sbs_time[i,j,y] = sbs_time[i,j,y] + sublist[["sbs_time"]]/N
      dc_time[i,j,y] = dc_time[i,j,y] + sublist[["dc_time"]]/N
      
      hdcd_fast_ari[i,j,y] = hdcd_fast_ari[i,j,y] + sublist[["hdcd_fast_ari"]]/N
      pilliat_ari[i,j,y] = pilliat_ari[i,j,y] + sublist[["pilliat_ari"]]/N
      inspect_ari[i,j,y] = inspect_ari[i,j,y] + sublist[["inspect_ari"]]/N
      subset_ari[i,j,y] = subset_ari[i,j,y] + sublist[["subset_ari"]]/N
      sbs_ari[i,j,y] = sbs_ari[i,j,y] + sublist[["sbs_ari"]]/N
      dc_ari[i,j,y] = dc_ari[i,j,y] + sublist[["dc_ari"]]/N
      
      # if(i==1 & j ==1 & y==1){
      #   print(sublist[["hdcd_ari"]])
      # }
      
      if(y!= 1){
        hdcd_fast_hausd[i,j,y] = hdcd_fast_hausd[i,j,y] + sublist[["hdcd_fast_hausd"]]/N
        pilliat_hausd[i,j,y] = pilliat_hausd[i,j,y] + sublist[["pilliat_hausd"]]/N
        inspect_hausd[i,j,y] = inspect_hausd[i,j,y] + sublist[["inspect_hausd"]]/N
        subset_hausd[i,j,y] = subset_hausd[i,j,y] + sublist[["subset_hausd"]]/N
        sbs_hausd[i,j,y] = sbs_hausd[i,j,y] + sublist[["sbs_hausd"]]/N
        dc_hausd[i,j,y] = dc_hausd[i,j,y] + sublist[["dc_hausd"]]/N
      }
      
      
      
    }
    
  }
  
  
}

if(save){
  saveRDS(result, file=sprintf("%s/result.RDA", savedir))
  saveRDS(hdcd_fast_hausd, file=sprintf("%s/hdcd_fast_haus.RDA", savedir))
  saveRDS(hdcd_fast_time, file=sprintf("%s/hdcd_fast_time.RDA", savedir))
  saveRDS(hdcd_fast_Kerr, file=sprintf("%s/hdcd_fast_Kerr.RDA", savedir))
  saveRDS(hdcd_fast_ari, file=sprintf("%s/hdcd_fast_ari.RDA", savedir))
  
  saveRDS(inspect_hausd, file=sprintf("%s/inspect_haus.RDA", savedir))
  saveRDS(inspect_time, file=sprintf("%s/inspect_time.RDA", savedir))
  saveRDS(inspect_Kerr, file=sprintf("%s/inspect_Kerr.RDA", savedir))
  saveRDS(inspect_ari, file=sprintf("%s/inspect_ari.RDA", savedir))
  
  saveRDS(pilliat_hausd, file=sprintf("%s/pilliat_haus.RDA", savedir))
  saveRDS(pilliat_time, file=sprintf("%s/pilliat_time.RDA", savedir))
  saveRDS(pilliat_Kerr, file=sprintf("%s/pilliat_Kerr.RDA", savedir))
  saveRDS(pilliat_ari, file=sprintf("%s/pilliat_ari.RDA", savedir))
  
  saveRDS(subset_hausd, file=sprintf("%s/subset_haus.RDA", savedir))
  saveRDS(subset_time, file=sprintf("%s/subset_time.RDA", savedir))
  saveRDS(subset_Kerr, file=sprintf("%s/subset_Kerr.RDA", savedir))
  saveRDS(subset_ari, file=sprintf("%s/subset_ari.RDA", savedir))
  
  saveRDS(sbs_hausd, file=sprintf("%s/sbs_haus.RDA", savedir))
  saveRDS(sbs_time, file=sprintf("%s/sbs_time.RDA", savedir))
  saveRDS(sbs_Kerr, file=sprintf("%s/sbs_Kerr.RDA", savedir))
  saveRDS(sbs_ari, file=sprintf("%s/sbs_ari.RDA", savedir))
  
  saveRDS(dc_hausd, file=sprintf("%s/dc_haus.RDA", savedir))
  saveRDS(dc_time, file=sprintf("%s/dc_time.RDA", savedir))
  saveRDS(dc_Kerr, file=sprintf("%s/dc_Kerr.RDA", savedir))
  saveRDS(dc_ari, file=sprintf("%s/dc_ari.RDA", savedir))
  
  
  infofile<-file(sprintf("%s/parameters.txt", savedir))
  writeLines(c(sprintf("N = %d", N),
               sprintf("n = %d", n),
               sprintf("p = %d", p),
               sprintf("Sparse constant = %f", sparse_const), 
               sprintf("Dense constant = %f", dense_const), 
               sprintf("Rescale variance = %d", as.integer(rescale_variance))),
             infofile)
  close(infofile)
}

# creating table:
if(save){
  # output latex table
  printlines = c("%%REMEMBER to use package \\usepackage{rotating}!!",
                 " \\begin{sidewaystable} \\centering",
                 "\\caption{Multiple changepoints}",
                 "\\label{}",
                 "\\small",
                 "\\begin{adjustbox}{scale=0.5,center}",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} cccc|ccccccc|ccccccc|ccccccc|ccccccc}",
                 "\\hline", 
                 "\\multicolumn{4}{c|}{Parameters} & \\multicolumn{7}{c|}{Hausdorff distance} &\\multicolumn{7}{c|}{$\\left | \\widehat{K}-K \\right |$} &\\multicolumn{7}{c|}{ARI} &\\multicolumn{7}{c}{Time in miliseconds} \\\\ \\hline ",
                 "$n$ & $p$ & Sparsity & K & \\text{HDCD} & \\text{HDCD'} & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} & \\text{HDCD} & \\text{HDCD'} & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} &\\text{HDCD} & \\text{HDCD'} & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} &\\text{HDCD}& \\text{HDCD'} & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} \\\\", 
                 "\\hline \\")
  
  
  for (i in 1:length(ns)) {
    
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:numconfig) {
        conf = config(y, n, p)
        etas = conf[[1]]
        #sparsities = conf[[2]]
        #mus = conf[[3]]
        sparsity = conf[[4]]
        if(is.null(sparsity)){
          sparsity="-"
        }
        
        #k = kfunc(y,n,p)
        #eta = round(0.4*n)
        #rootnorm = 0
        #if(k<sqrt(p*log(n^4))){
        #  rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
        #}else{
        #  rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
        #}
        string = sprintf("%d & %d & %s & %d ", n, p, sparsity, length(etas))
        
        
        if(y==1){
          for (t in 1:7) {
            string = sprintf("%s & - ", string)
            
          }
        }
        
        else{
          res = round(c(hdcd_fast_hausd[i,j,y],hdcd_long_hausd[i,j,y],  pilliat_hausd[i,j,y], inspect_hausd[i,j,y], sbs_hausd[i,j,y], subset_hausd[i,j,y],dc_hausd[i,j,y]),digits=3)
          minind = (res==min(na.omit(res)))
          res = c(hdcd_fast_hausd[i,j,y],hdcd_long_hausd[i,j,y], pilliat_hausd[i,j,y], inspect_hausd[i,j,y], sbs_hausd[i,j,y], subset_hausd[i,j,y],dc_hausd[i,j,y])
          
          for (t in 1:length(res)) {
            if(is.na(res[t])){
              string = sprintf("%s & - ", string)
            }
            else if(minind[t]){
              string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
            }else{
              string = sprintf("%s & %.3f", string, res[t])
            }
          }
        }
        
        res = round(c(hdcd_fast_Kerr[i,j,y],hdcd_long_Kerr[i,j,y], pilliat_Kerr[i,j,y], inspect_Kerr[i,j,y], sbs_Kerr[i,j,y], subset_Kerr[i,j,y],dc_Kerr[i,j,y]),digits=3)
        minind = (abs(res)==min(abs(res)))
        res = c(hdcd_fast_Kerr[i,j,y], hdcd_long_Kerr[i,j,y],pilliat_Kerr[i,j,y], inspect_Kerr[i,j,y], sbs_Kerr[i,j,y], subset_Kerr[i,j,y],dc_Kerr[i,j,y])
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        res = round(c(hdcd_fast_ari[i,j,y], hdcd_long_ari[i,j,y],pilliat_ari[i,j,y], inspect_ari[i,j,y], sbs_ari[i,j,y], subset_ari[i,j,y],dc_ari[i,j,y]),digits=3)
        minind = (abs(res)==max(abs(res)))
        res = c(hdcd_fast_ari[i,j,y], hdcd_long_ari[i,j,y],pilliat_ari[i,j,y], inspect_ari[i,j,y], sbs_ari[i,j,y], subset_ari[i,j,y],dc_ari[i,j,y])
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        
        res = round(1000*c(hdcd_fast_time[i,j,y],hdcd_long_time[i,j,y], pilliat_time[i,j,y], inspect_time[i,j,y], sbs_time[i,j,y], subset_time[i,j,y],dc_time[i,j,y]),digits=3)
        minind = (res==min(res))
        res = 1000*c(hdcd_fast_time[i,j,y], hdcd_long_time[i,j,y],pilliat_time[i,j,y], inspect_time[i,j,y], sbs_time[i,j,y], subset_time[i,j,y],dc_time[i,j,y])
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        string = sprintf("%s \\\\", string)
        printlines = c(printlines, string)
      }
    }
  }
  
  printlines = c(printlines, c("\\hline \\\\[-1.8ex]",
                               "\\end{tabular}",
                               "\\end{adjustbox}",
                               "\\end{sidewaystable}\\}"))
  texfile<-file(sprintf("%s/table.tex", savedir))
  writeLines(printlines, texfile)
  close(texfile)
  
}

