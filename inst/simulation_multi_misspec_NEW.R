#### Simulation for multiple change-points
#install.packages("InspectChangepoint", repos = "http:cran.us.r-project.org")
#install.packages("hdbinseg",repos = "http:cran.us.r-project.org")
#install.packages("/mn/sarpanitu/ansatte-u2/pamoen/project_inspect/simulations_nov_22/ESAC", repos=NULL, type="source")

library(doSNOW)
library(HDCD)
library(foreach)
## same magnitude in all coords

#maindir = "/mn/sarpanitu/ansatte-u2/pamoen/project_inspect/simulations_jan_23/"
maindir = "/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/Simulations_HDCD"
dateandtime = gsub(" ", "--",as.character(Sys.time()))
dateandtime = gsub(":", ".", dateandtime)
savedir = file.path(maindir, dateandtime)
calibrates_path = "/mn/sarpanitu/ansatte-u2/pamoen/project_inspect/simulations_jan_23/2023-01-04--00.01.24/multi/calibrates.RDA"

save = TRUE

if(save){
  dir.create(savedir, showWarnings = FALSE)
  savedir = file.path(maindir, sprintf("%s/multi_misspec",dateandtime))
  dir.create(savedir, showWarnings = FALSE)
  
}



N = 12
num_cores = 6
sparse_const = 3.5
dense_const = 3.5
set.seed(1996)
rescale_variance = TRUE
Ncal = 100
tol = 1/Ncal
even_spread = TRUE
calibrated = FALSE
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
ns = c(200)
#ns = c(200)
#ns = c(500)
#ps = c(1000)
ps = c(200)
#ps = c(100)
#ps = c(50,100)
#ps = c(50,100,500,5000)
#ps = ns[]


# rez = config(7, n, p, mus)
# etas = rez[[1]]
# sparsities = rez[[2]]
# mus = rez[[3]]



calibrates = NULL
# calibrate
if(calibrated){
  calibrates = readRDS(calibrates_path)
}else{
  #cl <- makeCluster(num_cores,type="SOCK")
  #registerDoSNOW(cl)
  #pb <- txtProgressBar(max = ((length(ns)*length(ps))), style = 3)
  #progress <- function(n) setTxtProgressBar(pb, n)
  #opts <- list(progress = progress)
  #calibrates = foreach(z = 0:((length(ns)*length(ps))-1),.options.snow = opts) %dopar% {
  for(z in 0:((length(ns)*length(ps))-1)){
    library(HDCD)
    #if(z%%10 == 0){
    #  print(z)
    #}
    set.seed(300*z+1)
    rez = list()
    nind = floor(z/length(ps))+1
    pind = z%%length(ps)+1
    cc = ESAC_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,K=4, alpha = 1.5, fast = FALSE,
                        rescale_variance = rescale_variance, 
                        bonferroni = TRUE, debug=FALSE)
    
    lambda = sqrt(log(ps[pind]*log(ns[nind]))/2)
    set.seed(300*z+1)
    cc2 = Inspect_calibrate(n = ns[nind], p = ps[pind], N=Ncal, tol=tol,lambda = lambda , alpha = 1.5, K = 4,eps=1e-10,
                            maxiter=10000,rescale_variance = rescale_variance, debug =FALSE)
    set.seed(300*z+1)
    cc3 = Pilliat_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,K = 2, alpha = 1.5, 
                            rescale_variance = rescale_variance, bonferroni=TRUE, test_all=TRUE, debug=FALSE)
    
    #source("/mn/sarpanitu/ansatte-u2/pamoen/project_inspect/simulations_nov_22/HDCD/SUBSET/main.R")
    source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")
    
    set.seed(300*z+1)
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
    calibrates = list(rez)
  }
  #close(pb)
  #stopCluster(cl)
  
  saveRDS(calibrates, file=sprintf("%s/calibrates.RDA", savedir))
}




numconfig = 12



## #' @import MASS
#' @import InspectChangepoint
config = function(i,h,n,p){
  #mus = matrix(0, nrow=p, ncol=n)
  noise = NULL
  mus = matrix(0, nrow=p, ncol=n)
  #X = NULL
  #k = 0
  etas = c()
  sparsities = c()
  sparsity = c()
  
  if(h==2){
    # generate changepoints
    # 3 mixed
    sparsity="Mixed"
    etas = round(c(n/4,2*n/4, 3*n/4))
    sparsities = sample(c(1,2), 3, replace=TRUE)
    for (j in 1:3){
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
        sparsities[j] = sample(1:floor(sqrt(p*log(n^4))),1)
        k = sparsities[j]
        phi = sparse_const/sqrt(Delta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
      }else{
        #dense
        sparsities[j] = sample(ceiling(sqrt(p*log(n^4))):p,1)
        k = sparsities[j]
        phi = dense_const/sqrt(Delta)*(p*log(n^4))^(1/4)
      }
      coords = sample(1:p, sparsities[j])
      diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +diff*phi
    }
  }
  
  if(i==1){
    #everything good
    noise = matrix(rnorm(n*p), nrow=p, ncol=n)
    
    
  }
  
  if(i==2){
    #uniform
    noise = matrix(runif(n*p,min=-sqrt(3), max=sqrt(3)), nrow=p, ncol=n)
    
    
  }
  
  if(i==3){
    #student t
    noise = matrix(rt(n*p,df=3), nrow=p, ncol=n)/sqrt(3)
    
    
  }
  
  if(i==4){
    #heavy tailed
    #noise = matrix(rcauchy(n*p),nrow=p,ncol=n)
    noise = matrix(rt(n*p,df=10), nrow=p, ncol=n)*sqrt(8/10)
    
  }
  
  
  if(i==5){
    #cs, loc 1
    rho = 0.1
    noise = t(mvrnorm(n = n, rep(0,p), Sigma = outer(1:p,1:p,function(a,b){rho^(abs(a-b))}), tol = 1e-6))
    
    
  }
  
  if(i==6){
    #cs, loc 2
    rho = 0.4
    noise = t(mvrnorm(n = n, rep(0,p), Sigma = outer(1:p,1:p,function(a,b){rho^(abs(a-b))}), tol = 1e-6))
    
    
  }
  
  if(i==7){
    #cs 1
    rho = 0.1
    noise = t(mvrnorm(n, mu = rep(0,p), Sigma = matrix(rho,p,p)+diag(p)*(1-rho)))
    
  }
  
  if(i==8){
    #cs 2
    rho = 0.4
    noise = t(mvrnorm(n, mu = rep(0,p), Sigma = matrix(rho,p,p)+diag(p)*(1-rho)))
    
  }
  
  if(i==9){
    #temp 1
    rho = 0.1
    noise = matrix(0,p,n)
    noise[,1] = rnorm(p)
    for (j in 1:p){
      for (t in 2:n) noise[j,t] = noise[j,t-1]*sqrt(rho) + rnorm(1)*sqrt(1-rho)
    }
    
    
  }
  
  if(i==10){
    #temp 2
    rho = 0.4
    noise = matrix(0,p,n)
    noise[,1] = rnorm(p)
    for (j in 1:p){
      for (t in 2:n) noise[j,t] = noise[j,t-1]*sqrt(rho) + rnorm(1)*sqrt(1-rho)
    }
    
    
  }
  
  if(i==11){
    #async 
    noise = matrix(rnorm(n*p), nrow=p, ncol=n)
    
    #sparsity="Mixed"
    #etas = sort(sample(1:(n-1), round(n/3)))
    #sparsities = sample(c(1,2), round(n/3), replace=TRUE)
    sparsity="Mixed"
    etas = round(c(n/4,2*n/4, 3*n/4))
    sparsities = sample(c(1,2), 3, replace=TRUE)
    mus = matrix(0, nrow=p, ncol=n)
    k = 0
    for (j in 1:3) {
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
        sparsities[j] = sample(1:floor(sqrt(p*log(n^4))),1)
        k = sparsities[j]
        phi = sparse_const/sqrt(Delta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
      }else{
        #dense
        sparsities[j] = sample(ceiling(sqrt(p*log(n^4))):p,1)
        k = sparsities[j]
        phi = dense_const/sqrt(Delta)*(p*log(n^4))^(1/4)
      }
      coords = sample(1:p, sparsities[j])
      diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      L = round(Delta/2)
      
      for (zz in 1:length(coords)) {
        start = sample(seq(etas[j]-L,etas[j]+L ),1)
        coord = coords[zz]
        mus[coord, start:n] = mus[coord, start:n] +diff[zz]*phi
        
      }
      
    }
    
    
    
    
    
  }
  
  if(i==12){
    #gradual
    noise = matrix(rnorm(n*p), nrow=p, ncol=n)
    #chgpt
    #sparsity="Mixed"
    #etas = sort(sample(1:(n-1), round(n/3)))
    #sparsities = sample(c(1,2), round(n/3), replace=TRUE)
    sparsity="Mixed"
    etas = round(c(n/4,2*n/4, 3*n/4))
    sparsities = sample(c(1,2), 3, replace=TRUE)
    mus = matrix(0, nrow=p, ncol=n)
    k=0
    for (j in 1:3) {
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
        sparsities[j] = sample(1:floor(sqrt(p*log(n^4))),1)
        k = sparsities[j]
        phi = sparse_const/sqrt(Delta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
      }else{
        #dense
        sparsities[j] = sample(ceiling(sqrt(p*log(n^4))):p,1)
        k = sparsities[j]
        phi = dense_const/sqrt(Delta)*(p*log(n^4))^(1/4)
      }
      coords = sample(1:p, sparsities[j])
      diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      
      for (zz in 1:length(coords)) {
        start = etas[j]
        stop = n
        if(j<3){
          stop = etas[j+1]
        }
        coord = coords[zz]
        mus[coord, start:stop] = mus[coord, start:stop] +diff[zz]*phi*(1:(stop-start+1))/(stop-start+1)
        if(j<3){
          mus[coord, (stop+1):n]= mus[coord, (stop+1):n]+diff[zz]*phi
        }
      }
      
    }
    
  }
  
  X = mus + noise
  
  #sds = apply(X, MARGIN=1, FUN = function(x) median(abs(diff(x)))*1.05)
  #X =t(apply(X, MARGIN=1, FUN = function(x) x/ (median(abs(diff(x)))*1.05)))
  #X = rescale.variance(X)
  
  
  return(list(etas, sparsities,mus, X,sparsity,noise))
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
  library(MASS)
  library(InspectChangepoint)
  counter = 1
  #source("/mn/sarpanitu/ansatte-u2/pamoen/project_inspect/simulations_nov_22/HDCD/SUBSET/main.R")
  source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")
  n = 1 
  p = 1
  i = 1
  j=1
  n = ns[i]
    
  p = ps[j]
  
  for(h in 1:2){
      for (y in 1:numconfig) {
        if(y %in% c(11,12) & h==1){
          next
        }
        #noise = matrix(rnorm(n*p), nrow=p, ncol=n)
        conf = config(y, h,n, p)
        etas = conf[[1]]
        sparsities = conf[[2]]
        #mus = conf[[3]]
        X = conf[[4]]
        sparsity = conf[[5]]
        #noise = conf[[6]]
        
        #X = noise + mus
        
        rezi = list()
        rezi[["i"]] = i
        rezi[["j"]] = j
        rezi[["y"]] = y
        rezi[["h"]] = h
        
        # ESAC_fast
        #xi = 4*sqrt(log(p*n))
        #lambda = sqrt(log(p*log(n)))
        lambda = sqrt(log(p*log(n))/2)
        
        
        a = proc.time()
        #res  = ESAC (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
        res = ESAC (X[,], 1.5,1, empirical=TRUE,alpha = 1.5, K = 4, thresholds_test = (calibrates[[j+(i-1)*length(ps)]])[[1]], droppartialsum = TRUE, fast =FALSE,
                    rescale_variance = rescale_variance, debug= FALSE)
        b=proc.time()
        rezi[["ESAC_fast_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["ESAC_fast_K"]]= res$changepointnumber
        rezi[["ESAC_fast_chgpts"]]= res$changepoints
        rezi[["ESAC_fast_hausd"]] = hausdorff(res$changepoints, etas,n)
        rezi[["ESAC_fast_K_error"]] = length(res$changepoints) - length(etas)
        rezi[["ESAC_fast_ari"]] = ARI(etas, res$changepoints, n)
        
        # a = proc.time()
        # #res  = ESAC (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
        # res = ESAC (X[,], 1.5,1, empirical=TRUE,alpha = 1.5, K = 5, thresholds_test = (calibrates[[j+(i-1)*length(ps)]])[[10]], droppartialsum = TRUE, fast =FALSE,
        #             rescale_variance = rescale_variance, debug= FALSE)
        # b=proc.time()
        # rezi[["ESAC_long_time"]] = (b-a)[1]+(b-a)[2]
        # rezi[["ESAC_long_K"]]= res$changepointnumber
        # rezi[["ESAC_long_chgpts"]]= res$changepoints
        # rezi[["ESAC_long_hausd"]] = hausdorff(res$changepoints, etas,n)
        # rezi[["ESAC_long_K_error"]] = length(res$changepoints) - length(etas)
        # rezi[["ESAC_long_ari"]] = ARI(etas, res$changepoints, n)
        # ESAC
        
        # slow FAST:
        # a = proc.time()
        # res  = ESAC (X, 2,2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[4]] , droppartialsum = FALSE, fast =FALSE,debug= FALSE)
        # b=proc.time()
        #
        # rezi[["ESAC_time"]] = (b-a)[1]+(b-a)[2]
        # rezi[["ESAC_K"]] = res$changepointnumber
        # rezi[["ESAC_chgpts"]]= res$changepoints
        # rezi[["ESAC_hausd"]] = hausdorff(res$changepoints, etas,n)
        # rezi[["ESAC_K_error"]] = length(res$changepoints) - length(etas)
        # rezi[["ESAC_ari"]] = ARI(etas, res$changepoints, n)
        #
        # ESAC slow
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
        #res = ESAC (X, 2,2, alpha = 1+1/6, K = 4, threshold_d_test = 1.5,
        #            threshold_s_test = 1.5, droppartialsum = FALSE, fast =FALSE,debug= FALSE)
        #b=proc.time()
        
        a = proc.time()
        res  = Pilliat(X, K = 2, alpha = 1.5, empirical = TRUE, threshold_dense = (calibrates[[j+(i-1)*length(ps)]])[[5]],
                       thresholds_partial = (calibrates[[j+(i-1)*length(ps)]])[[4]], thresholds_bj = (calibrates[[j+(i-1)*length(ps)]])[[6]],
                       rescale_variance = rescale_variance, test_all = TRUE,debug =FALSE)
        b=proc.time()
        #print(res)
        
        rezi[["pilliat_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["pilliat_K"]] = res$changepointnumber
        rezi[["pilliat_chgpts"]]= res$changepoints
        rezi[["pilliat_hausd"]] = hausdorff(res$changepoints, etas,n)
        rezi[["pilliat_K_error"]] = length(res$changepoints) - length(etas)
        rezi[["pilliat_ari"]] = ARI(etas, res$changepoints, n)
        #}
        
        
        #inspect
        a = proc.time()
        #res  = ESAC (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
        res = Inspect(X[,], xi=(calibrates[[j+(i-1)*length(ps)]])[[3]], alpha = 1.5, K = 4,eps=1e-10,
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
        #res  = ESAC (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
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
        #res  = ESAC (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
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
        res = dcbs.alg(X[,], cp.type = 1, phi = -1, thr = NULL, trim = NULL,
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
  rez
}
# bb = proc.time()
# print(bb-aa)
#parallel::stopCluster(cl)
close(pb)
stopCluster(cl)

{
  ESAC_fast_hausd = array(0, dim = c(2, numconfig) )
  ESAC_fast_Kerr = array(0, dim = c(2, numconfig) )
  ESAC_fast_time = array(0, dim= c(2, numconfig))
  ESAC_fast_ari = array(0, dim= c(2, numconfig))
  
  # ESAC_long_hausd = array(0, dim = c(2, numconfig) )
  # ESAC_long_Kerr = array(0, dim = c(2, numconfig) )
  # ESAC_long_time = array(0, dim= c(2, numconfig))
  # ESAC_long_ari = array(0, dim= c(2, numconfig))
  
  inspect_hausd = array(0, dim = c(2, numconfig) )
  inspect_Kerr = array(0, dim = c(2, numconfig) )
  inspect_time = array(0, dim= c(2, numconfig))
  inspect_ari = array(0, dim= c(2, numconfig))
  
  pilliat_hausd = array(0, dim = c(2, numconfig) )
  pilliat_Kerr = array(0, dim = c(2, numconfig) )
  pilliat_time = array(0, dim= c(2, numconfig))
  pilliat_ari = array(0, dim= c(2, numconfig))
  
  subset_hausd = array(0, dim = c(2, numconfig) )
  subset_Kerr = array(0, dim = c(2, numconfig) )
  subset_time = array(0, dim= c(2, numconfig))
  subset_ari = array(0, dim= c(2, numconfig))
  
  sbs_hausd = array(0, dim = c(2, numconfig) )
  sbs_Kerr = array(0, dim = c(2, numconfig) )
  sbs_time = array(0, dim= c(2, numconfig))
  sbs_ari = array(0, dim= c(2, numconfig))
  
  dc_hausd = array(0, dim = c(2, numconfig) )
  dc_Kerr = array(0, dim = c(2, numconfig) )
  dc_time = array(0, dim= c(2, numconfig))
  dc_ari = array(0, dim= c(2, numconfig))
  
  pilliat_hausd[1,] = NA
  #ESAC_long_hausd[,,1] = NA
  ESAC_fast_hausd[1,] = NA
  inspect_hausd[1,] = NA
  subset_hausd[1,] = NA
  sbs_hausd[1,] = NA
  dc_hausd[1,] = NA
  
  for (z in 1:N) {
    list = result[[z]]
    len = length(list)
    
    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      h = sublist[["h"]]
      
      ESAC_fast_Kerr[h,y] = ESAC_fast_Kerr[h,y] + abs(sublist[["ESAC_fast_K_error"]])/N
      #ESAC_long_Kerr[h,y] = ESAC_long_Kerr[h,y] + abs(sublist[["ESAC_long_K_error"]])/N
      pilliat_Kerr[h,y] = pilliat_Kerr[h,y] + abs(sublist[["pilliat_K_error"]])/N
      inspect_Kerr[h,y] = inspect_Kerr[h,y] + abs(sublist[["inspect_K_error"]])/N
      subset_Kerr[h,y] = subset_Kerr[h,y] + abs(sublist[["subset_K_error"]])/N
      sbs_Kerr[h,y] = sbs_Kerr[h,y] + abs(sublist[["sbs_K_error"]])/N
      dc_Kerr[h,y] = dc_Kerr[h,y] + abs(sublist[["dc_K_error"]])/N
      
      ESAC_fast_time[h,y] = ESAC_fast_time[h,y] + sublist[["ESAC_fast_time"]]/N
      #ESAC_long_time[h,y] = ESAC_long_time[h,y] + sublist[["ESAC_long_time"]]/N
      pilliat_time[h,y] = pilliat_time[h,y] + sublist[["pilliat_time"]]/N
      inspect_time[h,y] = inspect_time[h,y] + sublist[["inspect_time"]]/N
      subset_time[h,y] = subset_time[h,y] + sublist[["subset_time"]]/N
      sbs_time[h,y] = sbs_time[h,y] + sublist[["sbs_time"]]/N
      dc_time[h,y] = dc_time[h,y] + sublist[["dc_time"]]/N
      
      ESAC_fast_ari[h,y] = ESAC_fast_ari[h,y] + sublist[["ESAC_fast_ari"]]/N
      #ESAC_long_ari[h,y] = ESAC_long_ari[h,y] + sublist[["ESAC_long_ari"]]/N
      pilliat_ari[h,y] = pilliat_ari[h,y] + sublist[["pilliat_ari"]]/N
      inspect_ari[h,y] = inspect_ari[h,y] + sublist[["inspect_ari"]]/N
      subset_ari[h,y] = subset_ari[h,y] + sublist[["subset_ari"]]/N
      sbs_ari[h,y] = sbs_ari[h,y] + sublist[["sbs_ari"]]/N
      dc_ari[h,y] = dc_ari[h,y] + sublist[["dc_ari"]]/N
      
      # if(i==1 & j ==1 & y==1){
      #   print(sublist[["ESAC_ari"]])
      # }
      
      if(h!= 1){
        ESAC_fast_hausd[h,y] = ESAC_fast_hausd[h,y] + sublist[["ESAC_fast_hausd"]]/N
        #ESAC_long_hausd[h,y] = ESAC_long_hausd[h,y] + sublist[["ESAC_long_hausd"]]/N
        pilliat_hausd[h,y] = pilliat_hausd[h,y] + sublist[["pilliat_hausd"]]/N
        inspect_hausd[h,y] = inspect_hausd[h,y] + sublist[["inspect_hausd"]]/N
        subset_hausd[h,y] = subset_hausd[h,y] + sublist[["subset_hausd"]]/N
        sbs_hausd[h,y] = sbs_hausd[h,y] + sublist[["sbs_hausd"]]/N
        dc_hausd[h,y] = dc_hausd[h,y] + sublist[["dc_hausd"]]/N
      }
      
      
      
    }
    
  }
  
  
}

if(save){
  saveRDS(result, file=sprintf("%s/result.RDA", savedir))
  saveRDS(ESAC_fast_hausd, file=sprintf("%s/ESAC_fast_haus.RDA", savedir))
  saveRDS(ESAC_fast_time, file=sprintf("%s/ESAC_fast_time.RDA", savedir))
  saveRDS(ESAC_fast_Kerr, file=sprintf("%s/ESAC_fast_Kerr.RDA", savedir))
  saveRDS(ESAC_fast_ari, file=sprintf("%s/ESAC_fast_ari.RDA", savedir))
  
  # saveRDS(ESAC_long_hausd, file=sprintf("%s/ESAC_long_haus.RDA", savedir))
  # saveRDS(ESAC_long_time, file=sprintf("%s/ESAC_long_time.RDA", savedir))
  # saveRDS(ESAC_long_Kerr, file=sprintf("%s/ESAC_long_Kerr.RDA", savedir))
  # saveRDS(ESAC_long_ari, file=sprintf("%s/ESAC_long_ari.RDA", savedir))
  
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
               sprintf("n = %s", paste(ns, collapse=",")),
               sprintf("p = %s", paste(ps, collapse=",")),
               sprintf("Sparse constant = %f", sparse_const),
               sprintf("Dense constant = %f", dense_const),
               sprintf("Rescale variance = %d", as.integer(rescale_variance)),
               sprintf("Even spread= %d", even_spread)),
             infofile)
  close(infofile)
}


models = c("$M$", 
           "$M_{\\text{Unif}}$", 
           "$M_{t_{3}}$",
           "$M_{t_{10}}$",
           "$M_{\\text{cs, loc}}(\\rho = 0.1)$",
           "$M_{\\text{cs, loc}}(\\rho = 0.4)$",
           "$M_{\\text{cs}}(\\rho = 0.1)$",
           "$M_{\\text{cs}}(\\rho = 0.4)$",
           "$M_{\\text{AR}}(\\rho = 0.1)$",
           "$M_{\\text{AR}}(\\rho = 0.4)$",
           "$M_{\\text{async}}$",
           "$M_{\\text{grad}}$")

# creating table:
if(save){
  # output latex table
  printlines = c("%%REMEMBER to use package \\usepackage{rotating}!!",
                 " \\begin{table} \\centering",
                 "\\caption{Multiple changepoints under misspecified models}",
                 "\\label{tablemulti_misspec}",
                 "\\small",
                 "\\begin{adjustbox}{width=\\columnwidth}",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} cc|cccccc|cccccc}",
                 "\\hline",
                 "\\multicolumn{2}{c|}{Parameters} & \\multicolumn{6}{c|}{Hausdorff distance} &\\multicolumn{6}{c|}{$\\left | \\widehat{K}-K \\right |$}   \\\\ \\hline ",
                 "Model  & K & \\text{FAST}  & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} & \\text{FAST}  & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET}  & \\text{DC} \\\\",
                 "\\hline \\")
  
  
  n = ns
  p = ps
    
    for (y in 1:numconfig) {
      for(h in c(1,2)){
        
        if(h==1 & y %in% c(11,12)){
          next
        }
        conf = config(y, h,n, p)
        etas = conf[[1]]
        sparsities = conf[[2]]
        #mus = conf[[3]]
        X = conf[[4]]
        sparsity = conf[[5]]
        # if(is.null(sparsity)){
        #   sparsity="-"
        # }
        
        
        #k = kfunc(y,n,p)
        #eta = round(0.4*n)
        #rootnorm = 0
        #if(k<sqrt(p*log(n^4))){
        #  rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
        #}else{
        #  rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
        #}
        K = 0
        if(h>1){
          K=3
        }
        string = sprintf("%s & %d ", models[y], K)
        
        
        if(h==1){
          for (t in 1:6) {
            string = sprintf("%s & - ", string)
            
          }
        }else{
          res = round(c(ESAC_fast_hausd[h,y],  pilliat_hausd[h,y], inspect_hausd[h,y], sbs_hausd[h,y], subset_hausd[h,y],dc_hausd[h,y]),digits=3)
          minind = (res==min(na.omit(res)))
          #res = c(ESAC_fast_hausd[h,y],ESAC_long_hausd[h,y], pilliat_hausd[h,y], inspect_hausd[h,y], sbs_hausd[h,y], subset_hausd[h,y],dc_hausd[h,y])
          
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
        
        res = round(c(ESAC_fast_Kerr[h,y], pilliat_Kerr[h,y], inspect_Kerr[h,y], sbs_Kerr[h,y], subset_Kerr[h,y],dc_Kerr[h,y]),digits=3)
        minind = (abs(res)==min(abs(res)))
        #res = c(ESAC_fast_Kerr[i,j,y], ESAC_long_Kerr[i,j,y],pilliat_Kerr[i,j,y], inspect_Kerr[i,j,y], sbs_Kerr[i,j,y], subset_Kerr[i,j,y],dc_Kerr[i,j,y])
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        # res = round(c(ESAC_fast_ari[i,j,y], ESAC_long_ari[i,j,y],pilliat_ari[i,j,y], inspect_ari[i,j,y], sbs_ari[i,j,y], subset_ari[i,j,y],dc_ari[i,j,y]),digits=3)
        # minind = (abs(res)==max(abs(res)))
        # res = c(ESAC_fast_ari[i,j,y], ESAC_long_ari[i,j,y],pilliat_ari[i,j,y], inspect_ari[i,j,y], sbs_ari[i,j,y], subset_ari[i,j,y],dc_ari[i,j,y])
        # for (t in 1:length(res)) {
        #   if(minind[t]){
        #     string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
        #   }else{
        #     string = sprintf("%s & %.3f", string, res[t])
        #   }
        # }
        
        
        string = sprintf("%s \\\\", string)
        printlines = c(printlines, string)
      }
    }
  
  
  printlines = c(printlines, c("\\hline \\\\[-1.8ex]",
                               "\\end{tabular}",
                               "\\end{adjustbox}",
                               "\\end{table}"))
  texfile<-file(sprintf("%s/table.tex", savedir))
  writeLines(printlines, texfile)
  close(texfile)
  
}
