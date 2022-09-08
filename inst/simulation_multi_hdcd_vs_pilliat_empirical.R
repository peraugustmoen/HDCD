#### Simulation for a single change-point
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



N = 1000
num_cores = 6
sparse_const = 3
dense_const = 3
set.seed(1996)

Ncal = 1000
tol = 1/Ncal

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
ns = c(100,200)
#ps = c(50,100,500,1000,4000)
ps = c(50,100,500,5000)
#ps = ns[]


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
    # three dense
    sparsity = "Dense"
    etas = round(c(n/4, n*2/4, n*3/4))
    sparsities = c(p, round(p/2), round(p^(4/5)))
    Delta = min(diff(c(0, etas, n)))
    
    for (j in 1:length(etas)) {
      phi = dense_const/sqrt(Delta)*(p*log(n^4))^(1/4)
      k = sparsities[j]
      coords = sample(1:p, sparsities[j]) 
      diff = runif(k) *sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +diff*phi
    }
    
  }
  else if(i==3){
    # three sparse
    sparsity="Sparse"
    etas = round(c(n/4, n*2/4, n*3/4))
    sparsities = round(c(1,p^(1/5), p^(2/5)))
    Delta = min(diff(c(0, etas, n)))
    
    
    for (j in 1:length(etas)) {
      k = sparsities[j]
      phi = sparse_const/sqrt(Delta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
      coords = sample(1:p, sparsities[j])
      diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +diff*phi
    }
    
  }
  else if(i==4){
    # three mixed
    sparsity = "Mixed"
    etas = round(c(n/4, n*2/4, n*3/4))
    Delta = min(diff(c(0, etas, n)))
    sparsities = sample(c(1,2), 3, replace=TRUE)
    phi = 0
    for (j in 1:length(etas)) {
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
  else if(i==5){
    # n/3 dense
    sparsity="Dense"
    etas = sort(sample(1:(n-1), round(n/3)))
    sparsities = sample(ceiling(sqrt(p*log(n^4))):p,round(n/3), replace=TRUE)
    
    for (j in 1:length(etas)) {
      Delta = 1
      if(j==1){
        Delta = min(etas[1], etas[2]-etas[1])
      }else if(j<length(etas)){
        Delta = min(etas[j] - etas[j-1], etas[j+1] - etas[j])
      }else{
        Delta = min(etas[j]-etas[j-1], n -etas[j])
      }
      phi = dense_const/sqrt(Delta)*(p*log(n^4))^(1/4)
      k = sparsities[j]
      coords = sample(1:p, sparsities[j]) 
      diff = runif(k) *sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +diff*phi
    }
  }
  else if(i==6){
    # n/3 sparse
    sparsity = "Sparse"
    etas = sort(sample(1:(n-1), round(n/3)))
    sparsities = sample(1:floor(sqrt(p*log(n^4))),round(n/3), replace=TRUE )
    
    
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
      phi = sparse_const/sqrt(Delta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
      coords = sample(1:p, sparsities[j])
      diff = runif(k)*sample(c(-1,1), sparsities[j], replace=TRUE)
      diff = diff/norm(diff, type="2")
      mus[coords, (etas[j]+1):n] = mus[coords, (etas[j]+1):n] +diff*phi
    }
  }
  else if(i==7){
    # n/3 mixed
    sparsity="Mixed"
    etas = sort(sample(1:(n-1), round(n/3)))
    sparsities = sample(c(1,2), round(n/3), replace=TRUE)
    for (j in 1:round(n/3)) {
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
  return(list(etas, sparsities,mus,sparsity))
}


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
  #for(z in 0:((length(ns)*length(ps))-1)){
  library(HDCD)
  rez = list()
  nind = floor(z/length(ps))+1
  pind = z%%length(ps)+1
  cc = HDCD_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,K=2,fast = TRUE,debug=FALSE)
  cc2 = HDCD_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,debug=FALSE)
  cc3 = Pilliat_calibrate(ns[nind],ps[pind], N=Ncal, tol=tol,K = 2, debug=FALSE)
  rez[[1]] = cc[[1]]
  rez[[2]] = cc[[2]]
  rez[[3]] = cc2[[1]]
  rez[[4]] = cc2[[2]]
  rez[[5]] = cc3[[1]]
  rez[[6]] = cc3[[2]]
  rez[[7]] = cc3[[3]]
  rez[[8]] = ns[nind]
  rez[[9]] = ps[pind]
  rez
  #calibrates[[z+1]] = rez
  
}
close(pb)
stopCluster(cl) 



cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
result = foreach(z = 1:N,.options.snow = opts) %dopar% {
  rez = list()
  set.seed(z)
  library(HDCD)
  counter = 1
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
        xi = 4*sqrt(log(p*n))
        lambda = sqrt(log(p*log(n)))
        
        a = proc.time()
        res  = HDCD (X, 2,2, K=2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[2]] , droppartialsum = FALSE, fast =TRUE,debug= FALSE)
        b=proc.time()
        rezi[["hdcd_fast_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["hdcd_fast_K"]]= res$changepointnumber
        rezi[["hdcd_fast_chgpts"]]= res$changepoints
        rezi[["hdcd_fast_hausd"]] = hausdorff(res$changepoints, etas,n)
        rezi[["hdcd_fast_K_error"]] = length(res$changepoints) - length(etas)
        rezi[["hdcd_fast_ari"]] = ARI(etas, res$changepoints, n)
        # hdcd
        
        a = proc.time()
        res  = HDCD (X, 2,2, empirical=TRUE, thresholds_test =(calibrates[[j+(i-1)*length(ps)]])[[4]] , droppartialsum = FALSE, fast =FALSE,debug= FALSE)
        b=proc.time()
        
        rezi[["hdcd_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["hdcd_K"]] = res$changepointnumber
        rezi[["hdcd_chgpts"]]= res$changepoints
        rezi[["hdcd_hausd"]] = hausdorff(res$changepoints, etas,n)
        rezi[["hdcd_K_error"]] = length(res$changepoints) - length(etas)
        rezi[["hdcd_ari"]] = ARI(etas, res$changepoints, n)
        
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
        res  = Pilliat(X, 
                       K = 2, alpha = 1+1/6, empirical = TRUE, threshold_dense = (calibrates[[j+(i-1)*length(ps)]])[[6]], 
                       thresholds_partial = (calibrates[[j+(i-1)*length(ps)]])[[5]], thresholds_bj = (calibrates[[j+(i-1)*length(ps)]])[[7]],debug =FALSE)
        b=proc.time()
        print(res)
        
        rezi[["pilliat_time"]] = (b-a)[1]+(b-a)[2]
        rezi[["pilliat_K"]] = res$changepointnumber
        rezi[["pilliat_chgpts"]]= res$changepoints
        rezi[["pilliat_hausd"]] = hausdorff(res$changepoints, etas,n)
        rezi[["pilliat_K_error"]] = length(res$changepoints) - length(etas)
        rezi[["pilliat_ari"]] = ARI(etas, res$changepoints, n)  
        #}
        
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
  
  hdcd_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  hdcd_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  hdcd_time = array(0, dim= c(length(ns), length(ps), numconfig))
  hdcd_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  pilliat_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  pilliat_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  pilliat_time = array(0, dim= c(length(ns), length(ps), numconfig))
  pilliat_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  pilliat_hausd[,,1] = NA
  hdcd_fast_hausd[,,1] = NA
  hdcd_hausd[,,1] = NA
  
  for (z in 1:N) {
    list = result[[z]]
    len = length(list)
    
    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      i = sublist[["i"]]
      j = sublist[["j"]]
      
      hdcd_fast_Kerr[i,j,y] = hdcd_fast_Kerr[i,j,y] + sublist[["hdcd_fast_K_error"]]/N
      pilliat_Kerr[i,j,y] = pilliat_Kerr[i,j,y] + sublist[["pilliat_K_error"]]/N
      hdcd_Kerr[i,j,y] = hdcd_Kerr[i,j,y] + sublist[["hdcd_K_error"]]/N
      
      hdcd_fast_time[i,j,y] = hdcd_fast_time[i,j,y] + sublist[["hdcd_fast_time"]]/N
      pilliat_time[i,j,y] = pilliat_time[i,j,y] + sublist[["pilliat_time"]]/N
      hdcd_time[i,j,y] = hdcd_time[i,j,y] + sublist[["hdcd_time"]]/N
      
      hdcd_fast_ari[i,j,y] = hdcd_fast_ari[i,j,y] + sublist[["hdcd_fast_ari"]]/N
      pilliat_ari[i,j,y] = pilliat_ari[i,j,y] + sublist[["pilliat_ari"]]/N
      hdcd_ari[i,j,y] = hdcd_ari[i,j,y] + sublist[["hdcd_ari"]]/N
      
      # if(i==1 & j ==1 & y==1){
      #   print(sublist[["hdcd_ari"]])
      # }
      
      if(y!= 1){
        hdcd_fast_hausd[i,j,y] = hdcd_fast_hausd[i,j,y] + sublist[["hdcd_fast_hausd"]]/N
        pilliat_hausd[i,j,y] = pilliat_hausd[i,j,y] + sublist[["pilliat_hausd"]]/N
        hdcd_hausd[i,j,y] = hdcd_hausd[i,j,y] + sublist[["hdcd_hausd"]]/N
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
  
  saveRDS(hdcd_hausd, file=sprintf("%s/hdcd_haus.RDA", savedir))
  saveRDS(hdcd_time, file=sprintf("%s/hdcd_time.RDA", savedir))
  saveRDS(hdcd_Kerr, file=sprintf("%s/hdcd_Kerr.RDA", savedir))
  saveRDS(hdcd_ari, file=sprintf("%s/hdcd_ari.RDA", savedir))
  
  saveRDS(pilliat_hausd, file=sprintf("%s/pilliat_haus.RDA", savedir))
  saveRDS(pilliat_time, file=sprintf("%s/pilliat_time.RDA", savedir))
  saveRDS(pilliat_Kerr, file=sprintf("%s/pilliat_Kerr.RDA", savedir))
  saveRDS(pilliat_ari, file=sprintf("%s/pilliat_ari.RDA", savedir))
  
  
  infofile<-file(sprintf("%s/parameters.txt", savedir))
  writeLines(c(sprintf("N = %d", N),
               sprintf("n = %d", n),
               sprintf("p = %d", p),
               sprintf("Sparse constant = %f", sparse_const), 
               sprintf("Dense constant = %f", dense_const)),
             infofile)
  close(infofile)
}

# creating table:
if(save){
  # output latex table
  printlines = c("%%REMEMBER to use package \\usepackage{rotating}!!",
                 "\\begin{table}[!htbp] \\centering",
                 "\\caption{Multiple changepoints}",
                 "\\label{}",
                 "\\small",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} cccc|ccc|ccc|ccc|ccc}",
                 "\\hline", 
                 "\\multicolumn{4}{c|}{Parameters} & \\multicolumn{3}{c|}{Hausdorff distance} &\\multicolumn{3}{c|}{$\\widehat{K}-K$} &\\multicolumn{3}{c|}{ARI} &\\multicolumn{3}{c}{Time in miliseconds} \\\\ \\hline ",
                 "$n$ & $p$ & Sparsity & K & \\text{HDCD} & \\text{HDCD fast} & \\text{Pilliat}& \\text{HDCD} & \\text{HDCD fast} & \\text{Pilliat} & \\text{HDCD} & \\text{HDCD fast} & \\text{Pilliat}  & \\text{HDCD} & \\text{HDCD fast} & \\text{Pilliat} \\\\", 
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
        
        
        res = round(c(hdcd_hausd[i,j,y],hdcd_fast_hausd[i,j,y], pilliat_hausd[i,j,y]),digits=3)
        minind = (res==min(na.omit(res)))
        res = c(hdcd_hausd[i,j,y],hdcd_fast_hausd[i,j,y], pilliat_hausd[i,j,y])
        
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
        
        res = round(c(hdcd_Kerr[i,j,y],hdcd_fast_Kerr[i,j,y], pilliat_Kerr[i,j,y]),digits=3)
        minind = (abs(res)==min(abs(res)))
        res = c(hdcd_Kerr[i,j,y],hdcd_fast_Kerr[i,j,y], pilliat_Kerr[i,j,y])
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        res = round(c(hdcd_ari[i,j,y],hdcd_fast_ari[i,j,y], pilliat_ari[i,j,y]),digits=3)
        minind = (abs(res)==max(abs(res)))
        res = c(hdcd_ari[i,j,y],hdcd_fast_ari[i,j,y], pilliat_ari[i,j,y])
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        
        res = round(1000*c(hdcd_time[i,j,y],hdcd_fast_time[i,j,y], pilliat_time[i,j,y]),digits=3)
        minind = (res==min(res))
        res = 1000*c(hdcd_time[i,j,y],hdcd_fast_time[i,j,y], pilliat_time[i,j,y])
        
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
                               "\\end{table}"))
  texfile<-file(sprintf("%s/table.tex", savedir))
  writeLines(printlines, texfile)
  close(texfile)
  
}

