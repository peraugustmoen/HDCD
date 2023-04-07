#### Simulation for a single change-point
library(doSNOW)
library(HDCD)
library(foreach)
## same magnitude in all coords

maindir = "/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/Simulations_HDCD"
dateandtime = gsub(" ", "--",as.character(Sys.time()))
dateandtime = gsub(":", ".", dateandtime)
savedir = file.path(maindir, dateandtime)

source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/SUBSET_normal.R")
save = TRUE

if(save){
  dir.create(savedir, showWarnings = FALSE)
  savedir = file.path(maindir, sprintf("%s/single",dateandtime))
  dir.create(savedir, showWarnings =FALSE )
  plotdir = file.path(maindir, sprintf("%s/single/plots",dateandtime))
  dir.create(plotdir, showWarnings =FALSE )
  
}
#ns = c(200,500)
ns = c(200,500)
ps = c(100,1000,5000)
#ps = c(200)
#ps = ns[]
#ps = c(500,1000,2000,5000)
kvals=9
totruns = length(ns)*length(ps)*kvals

# simsbsN = 1000
# simsbstoln = 1

#for sbs:
# pis = matrix(nrow = length(ns), ncol =length(ps))
# for (i in 1:length(ns)) {
#   for (j in 1:length(ps)) {
#     pis[i,j] =  single_SBS_calibrate(ns[i],ps[j],simsbsN,simsbstoln,debug=FALSE)
# 
#   }
# }


N = 1000
num_cores = 6
sparse_const = 1.8
dense_const = 1.8
set.seed(1996)
rescale_variance = TRUE
#ns = c(500,1000,2000)
#ns = c(100,500,1000)
#ns = c(500,1000,2000)
#ns = c(100,500,1000)
#ps = c(1000,2000,10000)
#ns = c(200,500)



kfunc = function(y, n, p){
  k=1 
  eta = 1
  change = 1
  etadiv = 2
  if(y==1){
    k = 1
    eta = round(n/2)
    change = 0 
  }else if(y==2){
    k = 1
    eta = round(n/2)
  }
  else if(y==3){
    k = 1
    eta = round(n/5)
    etadiv = 5
  }else if(y==4){
    k = ceiling(p^(1/3))
    eta = round(n/2)
  }else if(y==5){
    k = ceiling(p^(1/3))
    eta = round(n/5)
    etadiv = 5
  }else if(y==6){
    k = ceiling(sqrt(p*log(n)))
    eta = round(n/2)
  }else if(y==7){
    k = ceiling(sqrt(p*log(n)))
    eta = round(n/5)
    etadiv = 5
  }else if(y==8){
    k = p
    eta = round(n/2)
  }else{
    k = p
    eta = round(n/5)
    etadiv = 5
  }
  
  return(c(k,eta,change,etadiv))
}

# calibrate
calibrated = FALSE
if(!calibrated){
  Ncal = 1000
  tol = 1/Ncal*10
  cl <- makeCluster(num_cores,type="SOCK")
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = ((length(ns)*length(ps))), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  calibrates = foreach(z = 0:((length(ns)*length(ps))-1),.options.snow = opts) %dopar% {
    #for(z in 0:((length(ns)*length(ps))-1)){
    library(HDCD)
    library(InspectChangepoint)
    library(hdbinseg)
    
    
    rez = list()
    nind = floor(z/length(ps))+1
    pind = z%%length(ps)+1
    set.seed(z+1)
    cc1 = ESAC_test_calibrate(ns[nind], ps[pind], N=3*Ncal, tol=tol/3, fast = TRUE, rescale_variance = rescale_variance )
    set.seed(z+1)
    cc4 = ESAC_test_calibrate(ns[nind], ps[pind], N=3*Ncal, tol=tol/3, fast = FALSE, rescale_variance = rescale_variance )
    
    
   
    # cc2 = ESAC_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,debug=FALSE)
    lambda = sqrt(log(ps[pind]*log(ns[nind]))/2)
    set.seed(z+1)
    cc2 = Inspect_test_calibrate(n=ns[nind], p = ps[pind], N=Ncal, tol=tol,lambda = lambda, 
                           maxiter=10000,rescale_variance = rescale_variance,debug =FALSE)
    set.seed(z+1)
    cc3 = Pilliat_test_calibrate(n = ns[nind],p=ps[pind], N=Ncal*3, tol=tol/3,
                                       rescale_variance = rescale_variance,debug=FALSE)
    
    ESACtreshinfo =  ESAC_thresholds(ns[nind],ps[pind])
    corr = length(ESACtreshinfo[[1]])
    set.seed(z+1)
    cc_bon= ESAC_test_calibrate(ns[nind], ps[pind], N=1, tol=tol/corr, fast = FALSE, rescale_variance = rescale_variance )
    
    pilliatthreshinfo = Pilliat_thresholds(ns[nind],ps[pind],tol)
    corr = length(pilliatthreshinfo[[1]]) + length(pilliatthreshinfo[[2]])
    set.seed(z+1)
    cc_bon2 = Pilliat_test_calibrate(n = ns[nind],p=ps[pind], N=1, tol=tol/corr,
                                             rescale_variance = rescale_variance,debug=FALSE)
    
    
    
    source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/main.R")
    
    set.seed(z+1)
    empirical_penalties = rep(NA, Ncal)
    empirical_penalties2 = rep(NA, Ncal)
    empirical_penalties3 = rep(NA, Ncal)
    
    
    if(rescale_variance){
      for (i in 1:(Ncal)) {
        mynulldata <- matrix(rnorm(ns[nind]*ps[pind],0,1),nrow=ps[pind],ncol=ns[nind],byrow=FALSE) # 5 variates with 1000 time points, no change
        #empirical_penalties[i] = wbs_penaltyfinder(InspectChangepoint::rescale.variance(mynulldata), SUBSET.normal_penalty, 100)
        empirical_penalties[i] = SUBSET.normal_penalty(InspectChangepoint::rescale.variance(mynulldata))[[3]]
        empirical_penalties2[i] = dcbs.alg(InspectChangepoint::rescale.variance(mynulldata),height=1, thr = 0,cp.type=1,phi=-1,temporal=FALSE )$mat[6]
        #empirical_penalties3[i] = single_Inspect(InspectChangepoint::rescale.variance(mynulldata),debug=FALSE)$cusumval
        
      }
    }else{
      for (i in 1:(Ncal)) {
        mynulldata <- matrix(rnorm(ns[nind]*ps[pind],0,1),nrow=ps[pind],ncol=ns[nind],byrow=FALSE) # 5 variates with 1000 time points, no change
        empirical_penalties[i] = SUBSET.normal_penalty((mynulldata))[[3]]
        empirical_penalties2[i] = dcbs.alg((mynulldata),height=1, thr = 0,cp.type=1,phi =-1,temporal=FALSE )$mat[6]
        #empirical_penalties3[i] = single_Inspect(matrix(mynulldata, byrow= TRUE, ncol = n, nrow=p),debug=FALSE)$cusumval
        
      }
    }
    
    empirical_penalties = sort(empirical_penalties,decreasing = TRUE)
    empirical_penalties2 = sort(empirical_penalties2,decreasing = TRUE)
    empirical_penalties3 = sort(empirical_penalties3,decreasing = TRUE)
    
    
    rez[[1]] = cc1[[1]] #FAST threshoolds
    rez[[2]] = cc1[[2]] #FAST thresholds
    #rez[[3]] = cc2[[1]] #Inspect xi
    rez[[3]] = cc2 #Inspect xi
    rez[[4]] = cc_bon2[[1]] # pilliat
    rez[[5]] = cc_bon2[[2]] # pilliat
    rez[[6]] = pilliatthreshinfo$thresholds_bj # pilliat
    rez[[7]] = ns[nind]
    rez[[8]] = ps[pind]
    rez[[9]] = empirical_penalties[max(1,round(tol*Ncal))] # subset penalty
    rez[[10]] = cc_bon[[1]] #FAST threshoolds
    rez[[11]] = cc_bon[[2]] #FAST thresholds
    rez[[12]] = empirical_penalties2[max(1,round(tol*Ncal))] 
    rez[[13]] = single_SBS_calibrate(ns[nind],ps[pind],Ncal,tol,rescale_variance= rescale_variance,
                                     debug=FALSE)
    
    ESACtreshinfo =  ESAC_thresholds(ns[nind],ps[pind])
    ttt = ESACtreshinfo[[1]]
    h = ESACtreshinfo[[2]]
    end = length(cc1[[2]])
    v = length(cc1[[2]])
    

    ## THIS IS OUTDATED. THIS IS DONE IN ESAC_TEST_CALIBRATE NOW
    if(h <=length(ESACtreshinfo[[1]])){

      rez[[14]] = c(cc1[[1]][1], max(cc1[[1]][2:(h-1)] / ttt[2:(h-1)])*ttt[2:(h-1)], max(cc1[[2]][h:end] / ttt[h:end])*ttt[h:end] )
      rez[[15]] = c(cc4[[1]][1], max(cc4[[1]][2:(h-1)] / ttt[2:(h-1)])*ttt[2:(h-1)], max(cc4[[2]][h:end] / ttt[h:end])*ttt[h:end] )
      rez[[19]] = c(cc4[[1]][1], max(cc4[[1]][2:(h-1)] / ttt[2:(h-1)])*ttt[2:(h-1)], max(cc4[[1]][h:end] / ttt[h:end])*ttt[h:end] )
    }else if(v>1){
      rez[[14]] = c(cc1[[1]][1], max(cc1[[1]][2:v] / ttt[2:v])*ttt[2:v])
      rez[[15]] = c(cc4[[1]][1], max(cc4[[1]][2:v] / ttt[2:v])*ttt[2:v])
      rez[[19]] = c(cc4[[1]][1], max(cc4[[1]][2:v] / ttt[2:v])*ttt[2:v])
    }else{
      rez[[14]] = max(cc1[[1]] / ttt)*ttt
      rez[[15]] = max(cc4[[1]] / ttt)*ttt
      rez[[19]] = max(cc4[[1]] / ttt)*ttt
    }
    
    
    pilliatthreshinfo = Pilliat_thresholds(ns[nind],ps[pind],tol/3)
    
    
    rez[[16]] = max(cc3[[1]] / pilliatthreshinfo$thresholds_partial) * pilliatthreshinfo$thresholds_partial
    rez[[17]] = cc3[[2]]
    rez[[18]] = pilliatthreshinfo$thresholds_bj
    
    #calibrates[[z+1]] = rez
    rez
  }
  
  close(pb)
  stopCluster(cl) 
  
  if(save){
    saveRDS(calibrates, file=sprintf("%s/calibrates.RDA", savedir))
  }
}else{
  calibrates = readRDS("/Users/peraugust/Library/CloudStorage/OneDrive-UniversitetetiOslo/project1/simulations/Simulations_HDCD/2022-12-04--13.47.30/single/calibrates.RDA")
}



# FIRST WITH EVEN CHANGES

#cl <- parallel::makeForkCluster(6)
#doParallel::registerDoParallel(cl)


cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
result = foreach(z = 1:N,.options.snow = opts) %dopar% {
  rez = list()
  counter = 1
  set.seed(z+1337)
  library(HDCD)
  library(hdbinseg)
  source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/SUBSET_normal.R")
  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:kvals) {
        rezi = list()
        # determine k!!
        tt = kfunc(y,n,p)
        k = tt[1]
        eta = tt[2]
        chg = tt[3]
        mus = matrix(0, nrow=p, ncol=n)
        noise = matrix(rnorm(n*p), nrow=p, ncol=n)
        #eta = round(0.2*n)
        #eta = round(sample(1:floor(n/2),1))
        
        if(chg!=0){
          diff = 0
          if(k<sqrt(p*log(n))){
            rootnorm = sparse_const/sqrt(eta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
            diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
          }else{
            rootnorm = dense_const/sqrt(eta)*(p*log(n))^(1/4)
            diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
          }
          
          mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)
        }
        X = mus+noise
        
        
        # DC
        a = proc.time()
        res_dc = dcbs.alg(X[,],height=1, thr = (calibrates[[j+(i-1)*length(ps)]])[[12]],cp.type=1,phi=-1,temporal=FALSE )
        b=proc.time()
        
        #res_dc$ecp
        rezi["dc_time"] = (b-a)[1]+(b-a)[2]
        if(!is.numeric(res_dc$ecp)){
          rezi["dc_res"] = 0
        }else{
          rezi["dc_res"]= 1
        }
        #then sbs
        a = proc.time()
        res_sbs = sbs.alg(X[,],height=1, thr = rep(sqrt((calibrates[[j+(i-1)*length(ps)]])[[13]]),p),cp.type=1,temporal=FALSE )
        #ress = single_SBS(X,(calibrates[[j+(i-1)*length(ps)]])[[13]]) )
        #ress$pos
        #res_sbs$ecp
        b=proc.time()
        
        if(!rescale_variance){
          res_sbs = single_SBS(X[,],(calibrates[[j+(i-1)*length(ps)]])[[13]]) 
        }else{
          res_sbs = single_SBS(rescale.variance(X),(calibrates[[j+(i-1)*length(ps)]])[[13]]) 
        }
        
        rezi["sbs_time"] = (b-a)[1]+(b-a)[2]
        if(res_sbs$maxval <= 0){
          #rezi["sbs_res"] = round(n/2)
          #ress = single_SBS(X,pis[i,j] )
          rezi["sbs_res"] = 0
          
        }else{
          rezi["sbs_res"] = 1
        }
        
        if(rescale_variance){
          X = rescale.variance(X)
        }
        
        # # now rescale variance
        # X = rescale.variance(X)
        # first inspect
        lambda = sqrt(log(p*log(n))/2) 
        
        a = proc.time()
        #if(rescale_variance){
        #  res_inspect = single_Inspect(InspectChangepoint::rescale.variance(X), lambda)
        #}else{
          res_inspect = single_Inspect(X[,], lambda)
        #}
        
        b=proc.time()
        rezi["inspect_time"] = (b-a)[1]+(b-a)[2]
        rezi["inspect_res"]= as.integer(res_inspect$cusumval>(calibrates[[j+(i-1)*length(ps)]])[[3]])
        
        # then ESAC
        a = proc.time()
        #res_ESAC  = single_ESAC(X[,], 1.5,1, debug =FALSE)
        res_ESAC = ESAC_test(X,droppartialsum = TRUE,fast=FALSE,rescale_variance = rescale_variance,
                   thresholds = (calibrates[[j+(i-1)*length(ps)]])[[19]])
        b=proc.time()
        rezi["ESAC_time"] = (b-a)[1]+(b-a)[2]
        rezi["ESAC_res"] = res_ESAC
        
        
        a = proc.time()
        #res_ESAC  = single_ESAC(X[,], 1.5,1, debug =FALSE)
        res_ESAC_slow2 = ESAC_test(X,droppartialsum = FALSE,fast=FALSE,rescale_variance = rescale_variance,
                             thresholds = (calibrates[[j+(i-1)*length(ps)]])[[11]])
        b=proc.time()
        
        rezi["ESAC_slow2_time"] = (b-a)[1]+(b-a)[2]
        rezi["ESAC_slow2_res"] = res_ESAC_slow2
        
        a = proc.time()
        #res_ESAC  = single_ESAC(X[,], 1.5,1, debug =FALSE)
        res_ESAC_slow = ESAC_test(X,droppartialsum = FALSE,fast=FALSE,rescale_variance = rescale_variance,
                                  thresholds = (calibrates[[j+(i-1)*length(ps)]])[[15]])
        b=proc.time()
        # a = proc.time()
        # #res_ESAC  = single_ESAC(X[,], 1.5,1, debug =FALSE)
        # res_ESAC_slow = single_ESAC(X,threshold_s = 1.5, threshold_d=1.0)
        # b=proc.time()
        
        rezi["ESAC_slow_time"] = (b-a)[1]+(b-a)[2]
        rezi["ESAC_slow_res"] = res_ESAC_slow
        #rezi["ESAC_s"] = res_ESAC$s
        
        # then scan
        # a = proc.time()
        # #res  = single_Scan (X[,], debug =FALSE)
        # b=proc.time()
        # scan_time = (b-a)[1]+(b-a)[2]
        # scan_res = res$pos
        # scan_s = res$s
        
        
        
        #rez = cbind(rez, c(z,i,j,y,k,inspect_res, inspect_time, ESAC_res, ESAC_time,ESAC_s, scan_res, scan_time, scan_s))
        
        # then subset
        a = proc.time()
        #if(rescale_variance){
        #  res_subset  = SUBSET.normal(rescale.variance(X),(calibrates[[j+(i-1)*length(ps)]])[[9]])
        #}else{
          res_subset  = SUBSET.normal(X[,],(calibrates[[j+(i-1)*length(ps)]])[[9]])
        #}
        
        b=proc.time()
        rezi["subset_time"] = (b-a)[1]+(b-a)[2]
        
        if(is.null(res_subset$cpt)){
          rezi["subset_res"] = 0
        }else{
          rezi["subset_res"] = 1
        }
        
        # then Pilliat
        a = proc.time()
        res_pilliat = Pilliat_test(X[,], rescale_variance = FALSE, empirical=TRUE, thresholds_partial = (calibrates[[j+(i-1)*length(ps)]])[[16]], 
                                 threshold_dense = (calibrates[[j+(i-1)*length(ps)]])[[17]], thresholds_bj = (calibrates[[j+(i-1)*length(ps)]])[[18]])
        b=proc.time()
        
        
        rezi["pilliat_time"] = (b-a)[1]+(b-a)[2]
        rezi["pilliat_res"] = res_pilliat
        
        
        rezi["z"] = z
        rezi["i"] = i
        rezi["j"] = j
        rezi["y"] = y
        rezi["k"] = k
        rezi["eta"] = eta
        rezi["chg"] = chg
        
        # rezi = list(z,i,j,y,k,inspect_res, inspect_time, ESAC_res, ESAC_time,ESAC_s, scan_res, scan_time, scan_s, sbs_res, sbs_time, 
        #                   subset_res, subset_time, dc_res, dc_time)
        
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
  inspect_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  inspect_time = array(0, dim= c(length(ns), length(ps), kvals))
  inspect_det = array(0, dim= c(length(ns), length(ps), kvals))
  ESAC_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  ESAC_time = array(0, dim= c(length(ns), length(ps), kvals))
  ESAC_det = array(0, dim= c(length(ns), length(ps), kvals))
  ESAC_slow_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  ESAC_slow_time = array(0, dim= c(length(ns), length(ps), kvals))
  ESAC_slow_det = array(0, dim= c(length(ns), length(ps), kvals))
  scan_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  scan_time = array(0, dim= c(length(ns), length(ps), kvals))
  scan_det = array(NA, dim= c(length(ns), length(ps), kvals))  
  sbs_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  sbs_time = array(0, dim= c(length(ns), length(ps), kvals))
  sbs_det = array(0, dim= c(length(ns), length(ps), kvals))
  subset_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  subset_time = array(0, dim= c(length(ns), length(ps), kvals))
  subset_det = array(0, dim= c(length(ns), length(ps), kvals))  
  dc_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  dc_time = array(0, dim= c(length(ns), length(ps), kvals))
  dc_det = array(0, dim= c(length(ns), length(ps), kvals))
  ESAC_slow2_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  ESAC_slow2_time = array(0, dim= c(length(ns), length(ps), kvals))
  ESAC_slow2_det = array(0, dim= c(length(ns), length(ps), kvals))
  
  pilliat_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  pilliat_time = array(0, dim= c(length(ns), length(ps), kvals))
  pilliat_det = array(0, dim= c(length(ns), length(ps), kvals))
  
  for (z in 1:N) {
    list = result[[z]]
    len = length(list)
    
    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      i = sublist[["i"]]
      j = sublist[["j"]]
      eta = sublist[["eta"]]
      
      inspect_res[i,j,y,z] = sublist["inspect_res"][[1]]
      inspect_time[i,j,y] = inspect_time[i,j,y] +  sublist["inspect_time"][[1]]/N
      inspect_det[i,j,y] = inspect_det[i,j,y] + (inspect_res[i,j,y,z])/N
      
      ESAC_res[i,j,y,z] = sublist["ESAC_res"][[1]]
      ESAC_time[i,j,y] = ESAC_time[i,j,y] + sublist["ESAC_time"][[1]]/N
      ESAC_det[i,j,y] = ESAC_det[i,j,y] + (ESAC_res[i,j,y,z] )/N
      
      ESAC_slow_res[i,j,y,z] = sublist["ESAC_slow_res"][[1]]
      ESAC_slow_time[i,j,y] = ESAC_slow_time[i,j,y] + sublist["ESAC_slow_time"][[1]]/N
      ESAC_slow_det[i,j,y] = ESAC_slow_det[i,j,y] + (ESAC_slow_res[i,j,y,z] )/N
      
      ESAC_slow2_res[i,j,y,z] = sublist["ESAC_slow2_res"][[1]]
      ESAC_slow2_time[i,j,y] = ESAC_slow2_time[i,j,y] + sublist["ESAC_slow2_time"][[1]]/N
      ESAC_slow2_det[i,j,y] = ESAC_slow2_det[i,j,y] + (ESAC_slow2_res[i,j,y,z] )/N
      
      sbs_res[i,j,y,z] = sublist["sbs_res"][[1]]
      sbs_time[i,j,y] = sbs_time[i,j,y] +  sublist["sbs_time"][[1]]/N
      sbs_det[i,j,y] = sbs_det[i,j,y] + (sbs_res[i,j,y,z] )/N
      
      subset_res[i,j,y,z] = sublist["subset_res"][[1]]
      subset_time[i,j,y] = subset_time[i,j,y] +  sublist["subset_time"][[1]]/N
      subset_det[i,j,y] = subset_det[i,j,y] + (subset_res[i,j,y,z])/N
      
      dc_res[i,j,y,z] = sublist["dc_res"][[1]]
      dc_time[i,j,y] = dc_time[i,j,y] + sublist["dc_time"][[1]]/N
      dc_det[i,j,y] = dc_det[i,j,y] + (dc_res[i,j,y,z])/N
      
      pilliat_res[i,j,y,z] = sublist["pilliat_res"][[1]]
      pilliat_time[i,j,y] = pilliat_time[i,j,y] + sublist["pilliat_time"][[1]]/N
      pilliat_det[i,j,y] = pilliat_det[i,j,y] + (pilliat_res[i,j,y,z])/N
      
      
    }
  }
  
  # for (i in 1:length(ns)) {
  #   n = ns[i]
  #   for(j in 1:length(ps)){
  #     #p = ps[j]
  #     for (y in 1:kvals) {
  #       subset = (result[2,] == i& result[3,] ==j & result[4,] == y )
  #       inspect_res[i,j,y,] = result[6, subset]
  #       inspect_time[i,j,y] = sum(result[7, subset])
  #       eta = round(n*0.2)
  #       inspect_det[i,j,y] = mean((inspect_res[i,j,y,] - eta)^2)
  #       ESAC_res[i,j,y,] = result[8, subset]
  #       ESAC_time[i,j,y] = sum(result[9, subset])
  #       ESAC_s[i,j,y,] = result[10, subset]
  #       ESAC_det[i,j,y] = mean((ESAC_res[i,j,y,] - eta)^2)
  #       
  #       scan_res[i,j,y,] = result[11, subset]
  #       scan_time[i,j,y] = sum(result[12, subset])
  #       scan_s[i,j,y,] = result[13, subset]
  #       scan_det[i,j,y] = mean((scan_res[i,j,y,] - eta)^2)
  #       
  #       sbs_res[i,j,y,] = result[14, subset]
  #       sbs_time[i,j,y] = sum(result[15, subset])
  #       sbs_det[i,j,y] = mean((sbs_res[i,j,y,] - eta)^2)
  #       
  #       subset_res[i,j,y,] = result[16, subset]
  #       subset_time[i,j,y] = sum(result[17, subset])
  #       subset_det[i,j,y] = mean((subset_res[i,j,y,] - eta)^2)
  #       
  #       dc_res[i,j,y,] = result[18, subset]
  #       dc_time[i,j,y] = sum(result[19, subset])
  #       dc_det[i,j,y] = mean((dc_res[i,j,y,] - eta)^2)
  #     }
  #   }
  # }
}







inspect_det - ESAC_det
inspect_det
ESAC_det

par(mfrow=c(1,2))
hist(inspect_res[1,2,4,])
hist(ESAC_res[1,2,4,])

#saving: 
if(save){
  saveRDS(inspect_det, file=sprintf("%s/inspect_det.RDA", savedir))
  saveRDS(inspect_time, file=sprintf("%s/inspect_time.RDA", savedir))
  saveRDS(inspect_res, file=sprintf("%s/inspect_res.RDA", savedir))
  saveRDS(ESAC_res, file=sprintf("%s/ESAC_res.RDA", savedir))
  saveRDS(ESAC_time, file=sprintf("%s/ESAC_time.RDA", savedir))
  saveRDS(ESAC_det, file=sprintf("%s/ESAC_det.RDA", savedir))
  saveRDS(ESAC_slow_res, file=sprintf("%s/ESAC_slow_res.RDA", savedir))
  saveRDS(ESAC_slow_time, file=sprintf("%s/ESAC_slow_time.RDA", savedir))
  saveRDS(ESAC_slow_det, file=sprintf("%s/ESAC_slow_det.RDA", savedir))
  saveRDS(ESAC_slow2_res, file=sprintf("%s/ESAC_slow_res.RDA", savedir))
  saveRDS(ESAC_slow2_time, file=sprintf("%s/ESAC_slow_time.RDA", savedir))
  saveRDS(ESAC_slow2_det, file=sprintf("%s/ESAC_slow_det.RDA", savedir))
  saveRDS(scan_res, file=sprintf("%s/scan_res.RDA", savedir))
  saveRDS(scan_time, file=sprintf("%s/scan_time.RDA", savedir))
  saveRDS(scan_det, file=sprintf("%s/scan_det.RDA", savedir))
  saveRDS(sbs_det, file=sprintf("%s/sbs_det.RDA", savedir))
  saveRDS(sbs_res, file=sprintf("%s/sbs_res.RDA", savedir))
  saveRDS(sbs_time, file=sprintf("%s/sbs_time.RDA", savedir))
  saveRDS(subset_det, file=sprintf("%s/subset_det.RDA", savedir))
  saveRDS(subset_res, file=sprintf("%s/subset_res.RDA", savedir))
  saveRDS(subset_time, file=sprintf("%s/subset_time.RDA", savedir))
  saveRDS(pilliat_det, file=sprintf("%s/pilliat_det.RDA", savedir))
  saveRDS(pilliat_res, file=sprintf("%s/pilliat_res.RDA", savedir))
  saveRDS(pilliat_time, file=sprintf("%s/pilliat_time.RDA", savedir))
  
  infofile<-file(sprintf("%s/parameters.txt", savedir))
  writeLines(c(sprintf("N = %d", N),
               sprintf("ns = %s", paste(ns, sep=" ", collapse=" ")),
               sprintf("ps = %s", paste(ps, sep=" ", collapse=" ")), 
               sprintf("Sparse constant = %f", sparse_const), 
               sprintf("Dense constant = %f", dense_const)),
             infofile)
  close(infofile)
}

# creating table:
if(save){
  # output latex table
  printlines = c("\\begin{table}[H] \\centering",
                 "\\caption{Single change-point detection}",
                 "\\label{tablesingletest}",
                 "\\small",
                 "\\begin{adjustbox}{scale=0.75,center}",
                 #"\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|ccccccc|ccccccc}",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|cccccc|cccccc}",
                 "\\hline", 
                 #"\\multicolumn{5}{c|}{Parameters} & \\multicolumn{7}{c|}{Detection rate} &\\multicolumn{7}{c|}{Time in milliseconds} \\\\ \\hline ",
                 "\\multicolumn{5}{c|}{Parameters} & \\multicolumn{6}{c|}{Detection rate} &\\multicolumn{6}{c|}{Time in milliseconds} \\\\ \\hline ",
                 #"$n$ & $p$ & $k$ &  $\\eta$ & $\\phi$ & \\text{FAST} & \\text{FAST'} & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET} & \\text{DC} & \\text{FAST} & \\text{FAST'} & \\text{Pilliat} & \\text{Inspect} &\\text{SBS} & \\text{SUBSET}& \\text{DC} \\\\", 
                 "$n$ & $p$ & $k$ &  $\\eta$ & $\\phi$ & \\text{FAST} & \\text{Pilliat} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET} & \\text{DC} & \\text{FAST}  & \\text{Pilliat} & \\text{Inspect} &\\text{SBS} & \\text{SUBSET}& \\text{DC} \\\\", 
                 "\\hline \\")
  
  for (i in 1:length(ns)) {
    
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:kvals) {
        tt = kfunc(y,n,p)
        k = tt[1]
        eta = tt[2]
        chg = tt[3]
        etadiv = tt[4]
        if(etadiv == 2 && y != 1){
          next
        }
        rootnorm = 0
        if(k<sqrt(p*log(n))){
          rootnorm = sparse_const/sqrt(eta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n))^(1/4)
        }
        if(chg==1){
          string = sprintf("%d & %d & %d & $n/%d$ & %.2f", n, p, k, as.integer(etadiv), rootnorm)
        }else{
          string = sprintf("%d & %d & - & - & -", n, p)
        }
        
        
        #res = round(c(ESAC_det[i,j,y],inspect_det[i,j,y], scan_det[i,j,y]),digits=3)
        #res = round(c(ESAC_slow_det[i,j,y],ESAC_det[i,j,y],pilliat_det[i,j,y],inspect_det[i,j,y], sbs_det[i,j,y], subset_det[i,j,y], dc_det[i,j,y]),digits=3)
        res = round(c(ESAC_slow_det[i,j,y],pilliat_det[i,j,y],inspect_det[i,j,y], sbs_det[i,j,y], subset_det[i,j,y], dc_det[i,j,y]),digits=3)
        minind = (res==max(res))
        if(chg==0){
          minind = (res==min(res))
        }
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        
        #res = round(c(ESAC_slow_time[i,j,y],ESAC_time[i,j,y],pilliat_time[i,j,y],inspect_time[i,j,y], sbs_time[i,j,y], subset_time[i,j,y], dc_time[i,j,y])*1000,digits=1)
        res = round(c(ESAC_slow_time[i,j,y],pilliat_time[i,j,y],inspect_time[i,j,y], sbs_time[i,j,y], subset_time[i,j,y], dc_time[i,j,y])*1000,digits=1)
        #res = round(c(ESAC_time[i,j,y],inspect_time[i,j,y], scan_time[i,j,y])/N*1000,digits=3)
        minind = (res==min(res))
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.1f} ", string, res[t])
          }else{
            string = sprintf("%s & %.1f", string, res[t])
          }
        }
        string = sprintf("%s \\\\", string)
        printlines = c(printlines, string)
        
        
        
        
        
      }
    }
  }
  printlines = c(printlines, "\\hline \\multicolumn{5}{c|}{Average detection rate}")
  #res = round(c(mean(ESAC_slow_det[,,2:kvals]),mean(ESAC_det[,,2:kvals]),mean(pilliat_det[,,2:kvals]),mean(inspect_det[,,2:kvals]), mean(sbs_det[,,2:kvals]), mean(subset_det[,,2:kvals]), mean(dc_det[,,2:kvals])),digits=3)
  res = round(c(mean(ESAC_slow_det[,,c(3,5,7,9)]),mean(pilliat_det[,,c(3,5,7,9)]),mean(inspect_det[,,c(3,5,7,9)]), mean(sbs_det[,,c(3,5,7,9)]), mean(subset_det[,,c(3,5,7,9)]), mean(dc_det[,,c(3,5,7,9)])),digits=3)
  minind = (res==max(res))
  string =""
  for (t in 1:length(res)) {
    if(minind[t]){
      string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
    }else{
      string = sprintf("%s & %.3f", string, res[t])
    }
  }
  string = sprintf("%s \\\\", string)
  printlines = c(printlines, string)
  
  printlines = c(printlines, c("\\hline \\\\[-1.8ex]",
                               "\\end{tabular}",
                               "\\end{adjustbox}",
                               "\\end{table}"))
  texfile<-file(sprintf("%s/table_even.tex", savedir))
  writeLines(printlines, texfile)
  close(texfile)
  
}
