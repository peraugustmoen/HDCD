#### Simulation for multiple change-points
#install.packages("InspectChangepoint", repos = "http:cran.us.r-project.org")
#install.packages("hdbinseg",repos = "http:cran.us.r-project.org")
#install.packages("/mn/sarpanitu/ansatte-u2/pamoen/project_inspect/simulations_nov_22/HDCD", repos=NULL, type="source")

library(doSNOW)
library(HDCD)
library(foreach)
## same magnitude in all coords
#maindir = "/mn/sarpanitu/ansatte-u2/pamoen/project_inspect/simulations_jan_23/"
maindir = "/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/Simulations_HDCD"
dateandtime = gsub(" ", "--",as.character(Sys.time()))
dateandtime = gsub(":", ".", dateandtime)
savedir = file.path(maindir, dateandtime)
#calibrates_path = "/Users/peraugust/uioomrade/project_inspect/simulations_jan_23/2023-01-04--00.01.24/multi/calibrates.RDA"

save = TRUE

if(save){
  dir.create(savedir, showWarnings = FALSE)
  savedir = file.path(maindir, sprintf("%s/multi_timing",dateandtime))
  dir.create(savedir, showWarnings = FALSE)
  
}



N = 6
num_cores = 6
sparse_const = 3.5
dense_const = 3.5
set.seed(1996)
rescale_variance = FALSE
Ncal = 1000
tol = 1/Ncal
even_spread = TRUE
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
#ns = seq(100,2000,by=100)
ns = c(100)
#ns = c(200)
#ns = c(500)
#ps = c(1000)
ps = seq(100,1000,by=100)
#ps = c(100)
#ps = c(50,100)
#ps = c(50,100,500,5000)
#ps = ns[]


# rez = config(7, n, p, mus)
# etas = rez[[1]]
# sparsities = rez[[2]]
# mus = rez[[3]]





numconfig = 7
config = function(i,n,p,e){
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
    # 2 dense
    sparsity="Dense"
    etas = sort(sample(1:(n-1), 2))
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
      if(!e){
        diff = rnorm(sparsities[j])
      }
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] + matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
    }
    
  }
  else if(i==3){
    # 2 sparse
    sparsity = "Sparse"
    etas = sort(sample(1:(n-1), 2))
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
      if(!e){
        diff = rnorm(sparsities[j])
      }
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +matrix(rep(diff*phi, n-etas[j]), nrow = length(coords))
    }
    
  }
  else if(i==4){
    # 2 mixed
    sparsity="Mixed"
    etas = sort(sample(1:(n-1), 2))
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
      if(!e){
        diff = rnorm(sparsities[j])
      }
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
      if(!e){
        diff = rnorm(sparsities[j])
      }
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
      if(!e){
        diff = rnorm(sparsities[j])
      }
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
      if(!e){
        diff = rnorm(sparsities[j])
      }
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
  library(InspectChangepoint)
  counter = 1
  #source("/mn/sarpanitu/ansatte-u2/pamoen/project_inspect/simulations_nov_22/HDCD/SUBSET/main.R")
  source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")
  
  n = 1 
  p = 1
  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      
      #for (y in 1:numconfig) {
      y = 4
      noise = matrix(rnorm(n*p), nrow=p, ncol=n)
      conf = config(y, n, p,even_spread)
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
      res = HDCD (X[,], 1.5,1, empirical=FALSE,alpha = 1.5, K = 4,droppartialsum = TRUE, fast =FALSE,
                  rescale_variance = rescale_variance, debug= FALSE)
      b=proc.time()
      rezi[["hdcd_fast_time2"]] = (b-a)[3]
      rezi[["hdcd_fast_K"]]= res$changepointnumber
      rezi[["hdcd_fast_chgpts"]]= res$changepoints
      rezi[["hdcd_fast_hausd"]] = hausdorff(res$changepoints, etas,n)
      rezi[["hdcd_fast_K_error"]] = length(res$changepoints) - length(etas)
      rezi[["hdcd_fast_ari"]] = ARI(etas, res$changepoints, n)
      
      # a = proc.time()
      # #res  = HDCD (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
      # res = HDCD (X[,], 1.5,1, empirical=TRUE,alpha = 1.5, K = 5, thresholds_test = (calibrates[[j+(i-1)*length(ps)]])[[10]], droppartialsum = TRUE, fast =FALSE,
      #             rescale_variance = rescale_variance, debug= FALSE)
      # b=proc.time()
      # rezi[["hdcd_long_time2"]] = (b-a)[3]
      # rezi[["hdcd_long_K"]]= res$changepointnumber
      # rezi[["hdcd_long_chgpts"]]= res$changepoints
      # rezi[["hdcd_long_hausd"]] = hausdorff(res$changepoints, etas,n)
      # rezi[["hdcd_long_K_error"]] = length(res$changepoints) - length(etas)
      # rezi[["hdcd_long_ari"]] = ARI(etas, res$changepoints, n)
      # hdcd
      
      # slow FAST:
      # a = proc.time()
      # res  = HDCD (X, 2,2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[4]] , droppartialsum = FALSE, fast =FALSE,debug= FALSE)
      # b=proc.time()
      #
      # rezi[["hdcd_time2"]] = (b-a)[3]
      # rezi[["hdcd_K"]] = res$changepointnumber
      # rezi[["hdcd_chgpts"]]= res$changepoints
      # rezi[["hdcd_hausd"]] = hausdorff(res$changepoints, etas,n)
      # rezi[["hdcd_K_error"]] = length(res$changepoints) - length(etas)
      # rezi[["hdcd_ari"]] = ARI(etas, res$changepoints, n)
      #
      # hdcd slow
      # if(dim(X)[2]>500){
      #   rezi[["pilliat_time2"]] = NA
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
      res  = Pilliat(X, K = 2, alpha = 1.5, 
                     rescale_variance = rescale_variance, test_all = TRUE,debug =FALSE)
      b=proc.time()
      #print(res)
      
      rezi[["pilliat_time2"]] = (b-a)[3]
      rezi[["pilliat_K"]] = res$changepointnumber
      rezi[["pilliat_chgpts"]]= res$changepoints
      rezi[["pilliat_hausd"]] = hausdorff(res$changepoints, etas,n)
      rezi[["pilliat_K_error"]] = length(res$changepoints) - length(etas)
      rezi[["pilliat_ari"]] = ARI(etas, res$changepoints, n)
      #}
      
      
      #inspect
      a = proc.time()
      #res  = HDCD (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
      res = Inspect(X[,], alpha = 1.5, K = 4,eps=1e-10, lambda = 4*sqrt(log(n*p)), xi = 4*sqrt(log(n*p)),
                    maxiter=10000,
                    rescale_variance = rescale_variance, debug=FALSE)
      b=proc.time()
      rezi[["inspect_time2"]] = (b-a)[3]
      rezi[["inspect_K"]]= res$changepointnumber
      rezi[["inspect_chgpts"]]= res$changepoints
      rezi[["inspect_hausd"]] = hausdorff(res$changepoints, etas,n)
      rezi[["inspect_K_error"]] = length(res$changepoints) - length(etas)
      rezi[["inspect_ari"]] = ARI(etas, res$changepoints, n)
      
      
      
      # subset
      a = proc.time()
      #res  = HDCD (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
      if(rescale_variance){
        res = change_main(InspectChangepoint::rescale.variance(X), SUBSET.normal, 100,penalties= c(4*log(n),2*log(p)))
      }else{
        res = change_main(X, SUBSET.normal,100, penalties= c(4*log(n)))
      }
      b=proc.time()
      rezi[["subset_time2"]] = (b-a)[3]
      rezi[["subset_K"]]= length(res[[2]])
      rezi[["subset_chgpts"]]= sort(res[[2]])
      rezi[["subset_hausd"]] = hausdorff(sort(res[[2]]), etas,n)
      rezi[["subset_K_error"]] = length(res[[2]]) - length(etas)
      rezi[["subset_ari"]] = ARI(etas, sort(res[[2]]), n)
      
      
      
      # SBS
      a = proc.time()
      #res  = HDCD (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
      res = sbs.alg(X[,], cp.type = 1, thr = NULL, trim = NULL, height = NULL,
                    temporal = FALSE, scales = NULL, diag = FALSE, B = 100, q = 0.01,
                    do.parallel = 1)
      b=proc.time()
      if(length(res$ecp)==0){
        res$ecp = NULL
      }
      rezi[["sbs_time2"]] = (b-a)[3]
      rezi[["sbs_K"]]= length(res$ecp)
      rezi[["sbs_chgpts"]]= sort(res$ecp)
      rezi[["sbs_hausd"]] = hausdorff(sort(res$ecp), etas,n)
      rezi[["sbs_K_error"]] = length(res$ecp) - length(etas)
      rezi[["sbs_ari"]] = ARI(etas, sort(res$ecp), n)
      
      # DC
      a = proc.time()
      res = dcbs.alg(X[,], cp.type = 1, phi = -1, thr = NULL, trim = NULL,
                     height = NULL, temporal = TRUE, scales = NULL, diag = FALSE,
                     B = 100, q = 0.01, do.parallel = 1)
      b=proc.time()
      if(length(res$ecp)==0){
        res$ecp = NULL
      }
      rezi[["dc_time2"]] = (b-a)[3]
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
      
      
      
      #}
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
  hdcd_fast_time2 = array(0, dim= c(length(ns), length(ps), numconfig))
  hdcd_fast_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  # hdcd_long_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  # hdcd_long_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  # hdcd_long_time2 = array(0, dim= c(length(ns), length(ps), numconfig))
  # hdcd_long_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  inspect_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  inspect_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  inspect_time2 = array(0, dim= c(length(ns), length(ps), numconfig))
  inspect_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  pilliat_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  pilliat_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  pilliat_time2 = array(0, dim= c(length(ns), length(ps), numconfig))
  pilliat_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  subset_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  subset_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  subset_time2 = array(0, dim= c(length(ns), length(ps), numconfig))
  subset_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  sbs_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  sbs_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  sbs_time2 = array(0, dim= c(length(ns), length(ps), numconfig))
  sbs_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  dc_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  dc_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  dc_time2 = array(0, dim= c(length(ns), length(ps), numconfig))
  dc_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  pilliat_hausd[,,1] = NA
  #hdcd_long_hausd[,,1] = NA
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
      #hdcd_long_Kerr[i,j,y] = hdcd_long_Kerr[i,j,y] + abs(sublist[["hdcd_long_K_error"]])/N
      pilliat_Kerr[i,j,y] = pilliat_Kerr[i,j,y] + abs(sublist[["pilliat_K_error"]])/N
      inspect_Kerr[i,j,y] = inspect_Kerr[i,j,y] + abs(sublist[["inspect_K_error"]])/N
      subset_Kerr[i,j,y] = subset_Kerr[i,j,y] + abs(sublist[["subset_K_error"]])/N
      sbs_Kerr[i,j,y] = sbs_Kerr[i,j,y] + abs(sublist[["sbs_K_error"]])/N
      dc_Kerr[i,j,y] = dc_Kerr[i,j,y] + abs(sublist[["dc_K_error"]])/N
      
      hdcd_fast_time2[i,j,y] = hdcd_fast_time2[i,j,y] + sublist[["hdcd_fast_time2"]]/N
      #hdcd_long_time2[i,j,y] = hdcd_long_time2[i,j,y] + sublist[["hdcd_long_time2"]]/N
      pilliat_time2[i,j,y] = pilliat_time2[i,j,y] + sublist[["pilliat_time2"]]/N
      inspect_time2[i,j,y] = inspect_time2[i,j,y] + sublist[["inspect_time2"]]/N
      subset_time2[i,j,y] = subset_time2[i,j,y] + sublist[["subset_time2"]]/N
      sbs_time2[i,j,y] = sbs_time2[i,j,y] + sublist[["sbs_time2"]]/N
      dc_time2[i,j,y] = dc_time2[i,j,y] + sublist[["dc_time2"]]/N
      
      hdcd_fast_ari[i,j,y] = hdcd_fast_ari[i,j,y] + sublist[["hdcd_fast_ari"]]/N
      #hdcd_long_ari[i,j,y] = hdcd_long_ari[i,j,y] + sublist[["hdcd_long_ari"]]/N
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
        #hdcd_long_hausd[i,j,y] = hdcd_long_hausd[i,j,y] + sublist[["hdcd_long_hausd"]]/N
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
  saveRDS(hdcd_fast_time2, file=sprintf("%s/hdcd_fast_time2.RDA", savedir))
  saveRDS(hdcd_fast_Kerr, file=sprintf("%s/hdcd_fast_Kerr.RDA", savedir))
  saveRDS(hdcd_fast_ari, file=sprintf("%s/hdcd_fast_ari.RDA", savedir))
  
  # saveRDS(hdcd_long_hausd, file=sprintf("%s/hdcd_long_haus.RDA", savedir))
  # saveRDS(hdcd_long_time2, file=sprintf("%s/hdcd_long_time2.RDA", savedir))
  # saveRDS(hdcd_long_Kerr, file=sprintf("%s/hdcd_long_Kerr.RDA", savedir))
  # saveRDS(hdcd_long_ari, file=sprintf("%s/hdcd_long_ari.RDA", savedir))
  
  saveRDS(inspect_hausd, file=sprintf("%s/inspect_haus.RDA", savedir))
  saveRDS(inspect_time2, file=sprintf("%s/inspect_time2.RDA", savedir))
  saveRDS(inspect_Kerr, file=sprintf("%s/inspect_Kerr.RDA", savedir))
  saveRDS(inspect_ari, file=sprintf("%s/inspect_ari.RDA", savedir))
  
  saveRDS(pilliat_hausd, file=sprintf("%s/pilliat_haus.RDA", savedir))
  saveRDS(pilliat_time2, file=sprintf("%s/pilliat_time2.RDA", savedir))
  saveRDS(pilliat_Kerr, file=sprintf("%s/pilliat_Kerr.RDA", savedir))
  saveRDS(pilliat_ari, file=sprintf("%s/pilliat_ari.RDA", savedir))
  
  saveRDS(subset_hausd, file=sprintf("%s/subset_haus.RDA", savedir))
  saveRDS(subset_time2, file=sprintf("%s/subset_time2.RDA", savedir))
  saveRDS(subset_Kerr, file=sprintf("%s/subset_Kerr.RDA", savedir))
  saveRDS(subset_ari, file=sprintf("%s/subset_ari.RDA", savedir))
  
  saveRDS(sbs_hausd, file=sprintf("%s/sbs_haus.RDA", savedir))
  saveRDS(sbs_time2, file=sprintf("%s/sbs_time2.RDA", savedir))
  saveRDS(sbs_Kerr, file=sprintf("%s/sbs_Kerr.RDA", savedir))
  saveRDS(sbs_ari, file=sprintf("%s/sbs_ari.RDA", savedir))
  
  saveRDS(dc_hausd, file=sprintf("%s/dc_haus.RDA", savedir))
  saveRDS(dc_time2, file=sprintf("%s/dc_time2.RDA", savedir))
  saveRDS(dc_Kerr, file=sprintf("%s/dc_Kerr.RDA", savedir))
  saveRDS(dc_ari, file=sprintf("%s/dc_ari.RDA", savedir))
  
  
  infofile<-file(sprintf("%s/parameters.txt", savedir))
  writeLines(c(sprintf("N = %d", N),
               sprintf("n = %s", paste(ns, collapse=",")),
               sprintf("p = %s", paste(ps, collapse=",")),
               sprintf("Sparse constant = %f", sparse_const),
               sprintf("Dense constant = %f", dense_const),
               sprintf("Rescale variance = %d", as.integer(rescale_variance)),
               sprintf("Even spread= %d", even_spread)),
             infofile)
  close(infofile)
}





times = cbind(ps, ESAC= 1000*hdcd_fast_time2[1, ,4], Pilliat = 1000*pilliat_time2[1,,4], Inspect = 1000*inspect_time2[1,,4],
              SBS = 1000*sbs_time2[1,,4], SUBSET = 1000*subset_time2[1,,4], DC = 1000*dc_time2 [1,,4])

times = data.frame(times)
df <- melt(times ,  id.vars = 'ps', variable.name = 'Method')

ggplot(df, aes(ps,value))+ geom_line(mapping=aes(colour = Method))

#relative times:

times = cbind(ps, ESAC= hdcd_fast_time2[1, ,4] / hdcd_fast_time2[1,1 ,4], Pilliat = pilliat_time2[1,,4]/hdcd_fast_time2[1, 1,4], Inspect =inspect_time2[1,,4]/inspect_time2[1,1,4],
              SBS = sbs_time2[1,,4]/sbs_time2[1,1,4], SUBSET = subset_time2[1,,4]/subset_time2[1,1,4], DC = dc_time2 [1,,4]/dc_time2 [1,1,4])

times = data.frame(times)
df <- melt(times ,  id.vars = 'ps', variable.name = 'Method')

ggplot(df, aes(ps,value))+ geom_line(mapping=aes(colour = Method))

# in log:


times = cbind(ps, ESAC= log(1000*hdcd_fast_time2[1, ,4]), Pilliat = log(1000*pilliat_time2[1,,4]), Inspect = log(1000*inspect_time2[1,,4]),
              SBS = log(1000*sbs_time2[1,,4]), SUBSET = log(1000*subset_time2[1,,4]), DC = log(1000*dc_time2 [1,,4]))

times = data.frame(times)
df <- melt(times ,  id.vars = 'ps', variable.name = 'Method')

ggplot(df, aes(ps,value))+ geom_line(mapping=aes(colour = Method))







# creating table:
if(save){
  # output latex table
  printlines = c("%%REMEMBER to use package \\usepackage{rotating}!!",
                 " \\begin{table} \\centering",
                 "\\caption{Multiple changepoints}",
                 "\\label{}",
                 "\\small",
                 "\\begin{adjustbox}{scale=0.35,center}",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} cccc|cccccc|cccccc|cccccc}",
                 "\\hline",
                 "\\multicolumn{4}{c|}{Parameters} & \\multicolumn{6}{c|}{Hausdorff distance} &\\multicolumn{6}{c|}{$\\left | \\widehat{K}-K \\right |$}  &\\multicolumn{6}{c}{Time in miliseconds} \\\\ \\hline ",
                 "$n$ & $p$ & Sparsity & K & \\text{FAST}  & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} & \\text{FAST}  & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} &\\text{FAST} & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC}  \\\\",
                 "\\hline \\")
  
  
  for (i in 1:length(ns)) {
    
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:numconfig) {
        conf = config(y, n, p,even_spread)
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
          for (t in 1:6) {
            string = sprintf("%s & - ", string)
            
          }
        }
        
        else{
          res = round(c(hdcd_fast_hausd[i,j,y],  pilliat_hausd[i,j,y], inspect_hausd[i,j,y], sbs_hausd[i,j,y], subset_hausd[i,j,y],dc_hausd[i,j,y]),digits=3)
          minind = (res==min(na.omit(res)))
          #res = c(hdcd_fast_hausd[i,j,y],hdcd_long_hausd[i,j,y], pilliat_hausd[i,j,y], inspect_hausd[i,j,y], sbs_hausd[i,j,y], subset_hausd[i,j,y],dc_hausd[i,j,y])
          
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
        
        res = round(c(hdcd_fast_Kerr[i,j,y], pilliat_Kerr[i,j,y], inspect_Kerr[i,j,y], sbs_Kerr[i,j,y], subset_Kerr[i,j,y],dc_Kerr[i,j,y]),digits=3)
        minind = (abs(res)==min(abs(res)))
        #res = c(hdcd_fast_Kerr[i,j,y], hdcd_long_Kerr[i,j,y],pilliat_Kerr[i,j,y], inspect_Kerr[i,j,y], sbs_Kerr[i,j,y], subset_Kerr[i,j,y],dc_Kerr[i,j,y])
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        # res = round(c(hdcd_fast_ari[i,j,y], hdcd_long_ari[i,j,y],pilliat_ari[i,j,y], inspect_ari[i,j,y], sbs_ari[i,j,y], subset_ari[i,j,y],dc_ari[i,j,y]),digits=3)
        # minind = (abs(res)==max(abs(res)))
        # res = c(hdcd_fast_ari[i,j,y], hdcd_long_ari[i,j,y],pilliat_ari[i,j,y], inspect_ari[i,j,y], sbs_ari[i,j,y], subset_ari[i,j,y],dc_ari[i,j,y])
        # for (t in 1:length(res)) {
        #   if(minind[t]){
        #     string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
        #   }else{
        #     string = sprintf("%s & %.3f", string, res[t])
        #   }
        # }
        
        
        res = round(1000*c(hdcd_fast_time2[i,j,y], pilliat_time2[i,j,y], inspect_time2[i,j,y], sbs_time2[i,j,y], subset_time2[i,j,y],dc_time2[i,j,y]),digits=3)
        minind = (res==min(res))
        #res = 1000*c(hdcd_fast_time2[i,j,y], hdcd_long_time2[i,j,y],pilliat_time2[i,j,y], inspect_time2[i,j,y], sbs_time2[i,j,y], subset_time2[i,j,y],dc_time2[i,j,y])
        
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
                               "\\end{table}\\}"))
  texfile<-file(sprintf("%s/table.tex", savedir))
  writeLines(printlines, texfile)
  close(texfile)
  
}
