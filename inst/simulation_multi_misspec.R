#### Simulation for a single change-point
library(doSNOW)
library(HDCD)
library(foreach)
## same magnitude in all coords

maindir = "/Users/peraugust/OneDrive - Universitetet i Oslo/project_inspect/simulations/Simulations_HDCD"
dateandtime = gsub(" ", "--",as.character(Sys.time()))
dateandtime = gsub(":", ".", dateandtime)
savedir = file.path(maindir, dateandtime)


save = TRUE

if(save){
  dir.create(savedir, showWarnings = FALSE)
  savedir = file.path(maindir, sprintf("%s/multi_misspec",dateandtime))
  dir.create(savedir, showWarnings = FALSE)
  
}


# TODO: FIX I KODEN -- er litt kopi av inspectkoden nå

N = 1000
num_cores = 6
sparse_const = 3
dense_const = 3
set.seed(1996)

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
ns = c(500)
#ps = c(50,100,500,1000,4000)
ps = c(500)
#ps = ns[]

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
        mus[coord, start:n] = mus[coord, start:n] +diff[coord]*phi
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
        mus[coord, start:stop] = mus[coord, start:stop] +diff[coord]*phi*(1:(stop-start+1))/(stop-start+1)
        if(j<3){
          mus[coord, (stop+1):n]= mus[coord, (stop+1):n]+diff[coord]*phi
        }
      }
      
    }
  
  }
  
  X = mus + noise
  
  #sds = apply(X, MARGIN=1, FUN = function(x) median(abs(diff(x)))*1.05)
  #X =t(apply(X, MARGIN=1, FUN = function(x) x/ (median(abs(diff(x)))*1.05)))
  X = rescale.variance(X)
  
  
  return(list(etas, sparsities,mus, X,sparsity,noise))
}


# rezz = config(12, h=2,n=100,p=100)
# etas = rezz[[1]]
# sparsities = rezz[[2]]
# mus = rezz[[3]]
# X = rezz[[4]]
# sparsity = rezz[[5]]
# noise = rezz[[6]]
# plot(X[2,])
# plot(mus[2,])
# plot(noise[2,])

# rez = config(7, n, p, mus)
# etas = rez[[1]]
# sparsities = rez[[2]]
# mus = rez[[3]]









cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
result = foreach(z = 1:N,.options.snow = opts) %dopar% {
  rez = list()
  set.seed(z)
  library(HDCD)
  library(MASS)
  counter = 1
  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
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
          
          # inspect
          xi = 4*sqrt(log(p*n))
          lambda = sqrt(log(p*log(n)))
          
          a = proc.time()
          res = Inspect(X, lambda, xi)
          b=proc.time()
          rezi[["inspect_time"]] = (b-a)[1]+(b-a)[2]
          rezi[["inspect_K"]]= res$changepointnumber
          rezi[["inspect_chgpts"]]= res$changepoints
          rezi[["inspect_hausd"]] = hausdorff(res$changepoints, etas,n)
          rezi[["inspect_K_error"]] = length(res$changepoints) - length(etas)
          rezi[["inspect_ari"]] = ARI(etas, res$changepoints, n)
          # hdcd
          a = proc.time()
          res = HDCD (X, 2,1)
          b=proc.time()
          rezi[["hdcd_time"]] = (b-a)[1]+(b-a)[2]
          rezi[["hdcd_K"]] = res$changepointnumber
          rezi[["hdcd_chgpts"]]= res$changepoints
          rezi[["hdcd_hausd"]] = hausdorff(res$changepoints, etas,n)
          rezi[["hdcd_K_error"]] = length(res$changepoints) - length(etas)
          rezi[["hdcd_ari"]] = ARI(etas, res$changepoints, n)
          # scan
          # hdcd
          a = proc.time()
          res = Scan(X)
          b=proc.time()
          rezi[["scan_time"]] = (b-a)[1]+(b-a)[2]
          rezi[["scan_K"]] = res$changepointnumber
          rezi[["scan_chgpts"]]= res$changepoints
          rezi[["scan_hausd"]] = hausdorff(res$changepoints, etas,n)
          rezi[["scan_K_error"]] = length(res$changepoints) - length(etas)
          rezi[["scan_ari"]] = ARI(etas, res$changepoints, n)
          
          rezi[["true_K"]] = length(etas)
          rezi[["true_etas"]] = etas
          rezi[["true_sparsities"]] = sparsities
          
          rez[[counter]] = rezi
          counter = counter+1
          
          
          
        }
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
  inspect_hausd = array(0, dim = c(length(ns), length(ps), numconfig,2) )
  inspect_Kerr = array(0, dim = c(length(ns), length(ps), numconfig,2) )
  inspect_time = array(0, dim= c(length(ns), length(ps), numconfig,2))
  inspect_ari = array(0, dim= c(length(ns), length(ps), numconfig,2))
  
  hdcd_hausd = array(0, dim = c(length(ns), length(ps), numconfig,2) )
  hdcd_Kerr = array(0, dim = c(length(ns), length(ps), numconfig,2) )
  hdcd_time = array(0, dim= c(length(ns), length(ps), numconfig,2))
  hdcd_ari = array(0, dim= c(length(ns), length(ps), numconfig,2))
  
  scan_hausd = array(0, dim = c(length(ns), length(ps), numconfig,2) )
  scan_Kerr = array(0, dim = c(length(ns), length(ps), numconfig,2) )
  scan_time = array(0, dim= c(length(ns), length(ps), numconfig,2))
  scan_ari = array(0, dim= c(length(ns), length(ps), numconfig,2))
  
  scan_hausd[,,,1] = NA
  inspect_hausd[,,,1] = NA
  hdcd_hausd[,,,1] = NA
  
  for (z in 1:N) {
    list = result[[z]]
    len = length(list)
    
    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      i = sublist[["i"]]
      j = sublist[["j"]]
      h = sublist[["h"]]
      
      inspect_Kerr[i,j,y,h] = inspect_Kerr[i,j,y,h] + sublist[["inspect_K_error"]]/N
      scan_Kerr[i,j,y,h] = scan_Kerr[i,j,y,h] + sublist[["scan_K_error"]]/N
      hdcd_Kerr[i,j,y,h] = hdcd_Kerr[i,j,y,h] + sublist[["hdcd_K_error"]]/N
      
      inspect_time[i,j,y,h] = inspect_time[i,j,y,h] + sublist[["inspect_time"]]/N
      scan_time[i,j,y,h] = scan_time[i,j,y,h] + sublist[["scan_time"]]/N
      hdcd_time[i,j,y,h] = hdcd_time[i,j,y,h] + sublist[["hdcd_time"]]/N
      
      inspect_ari[i,j,y,h] = inspect_ari[i,j,y,h] + sublist[["inspect_ari"]]/N
      scan_ari[i,j,y,h] = scan_ari[i,j,y,h] + sublist[["scan_ari"]]/N
      hdcd_ari[i,j,y,h] = hdcd_ari[i,j,y,h] + sublist[["hdcd_ari"]]/N
      
      # if(i==1 & j ==1 & y==1){
      #   print(sublist[["hdcd_ari"]])
      # }
      
      if(h!= 1){
        inspect_hausd[i,j,y,h] = inspect_hausd[i,j,y,h] + sublist[["inspect_hausd"]]/N
        scan_hausd[i,j,y,h] = scan_hausd[i,j,y,h] + sublist[["scan_hausd"]]/N
        hdcd_hausd[i,j,y,h] = hdcd_hausd[i,j,y,h] + sublist[["hdcd_hausd"]]/N
      }
      
      
      
    }
    
  }
  
  
}

if(save){
  saveRDS(result, file=sprintf("%s/result.RDA", savedir))
  saveRDS(inspect_hausd, file=sprintf("%s/inspect_haus.RDA", savedir))
  saveRDS(inspect_time, file=sprintf("%s/inspect_time.RDA", savedir))
  saveRDS(inspect_Kerr, file=sprintf("%s/inspect_Kerr.RDA", savedir))
  saveRDS(inspect_ari, file=sprintf("%s/inspect_ari.RDA", savedir))
  
  saveRDS(hdcd_hausd, file=sprintf("%s/hdcd_haus.RDA", savedir))
  saveRDS(hdcd_time, file=sprintf("%s/hdcd_time.RDA", savedir))
  saveRDS(hdcd_Kerr, file=sprintf("%s/hdcd_Kerr.RDA", savedir))
  saveRDS(hdcd_ari, file=sprintf("%s/hdcd_ari.RDA", savedir))
  
  saveRDS(scan_hausd, file=sprintf("%s/scan_haus.RDA", savedir))
  saveRDS(scan_time, file=sprintf("%s/scan_time.RDA", savedir))
  saveRDS(scan_Kerr, file=sprintf("%s/scan_Kerr.RDA", savedir))
  saveRDS(scan_ari, file=sprintf("%s/scan_ari.RDA", savedir))
  
  
  infofile<-file(sprintf("%s/parameters.txt", savedir))
  writeLines(c(sprintf("N = %d", N),
               sprintf("n = %d", n),
               sprintf("p = %d", p),
               sprintf("Sparse constant = %f", sparse_const), 
               sprintf("Dense constant = %f", dense_const)),
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
                 "\\begin{table}[!htbp] \\centering",
                 "\\caption{Multiple changepoints, misspecified model}",
                 "\\label{}",
                 "\\small",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|ccc|ccc|ccc}",
                 "\\hline", 
                 "\\multicolumn{5}{c|}{Parameters} & \\multicolumn{3}{c|}{Hausdorff distance} &\\multicolumn{3}{c|}{$\\widehat{K}-K$} &\\multicolumn{3}{c|}{ARI}  \\\\ \\hline ",
                 "Model & $n$ & $p$ & Sparsity & K & \\text{HDCD} & \\text{Inspect} & \\text{Scan}& \\text{HDCD} & \\text{Inspect} & \\text{Scan} & \\text{HDCD} & \\text{Inspect} & \\text{Scan}  \\\\", 
                 "\\hline \\")
  
  
  for (i in 1:length(ns)) {
    
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (h in 1:2) {
        for (y in 1:numconfig) {
          #conf = config(y, n, p)
          #etas = conf[[1]]
          #sparsities = conf[[2]]
          #mus = conf[[3]]
          #sparsity = conf[[4]]
          #if(is.null(sparsity)){
          #  sparsity="-"
          #}
          if(h==1 & y %in% c(11,12)){
            next
          }
          sparsity = "-"
          etas = c()
          if(h==2){
            sparsity="Mixed"
            etas = round(c(n/4, 2*n/4, 3*n/4))
          }
          
          #k = kfunc(y,n,p)
          #eta = round(0.4*n)
          #rootnorm = 0
          #if(k<sqrt(p*log(n^4))){
          #  rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
          #}else{
          #  rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
          #}
          string = sprintf("%s & %d & %d & %s & %d ", models[y], n, p, sparsity, length(etas))
          
          
          res = round(c(hdcd_hausd[i,j,y,h],inspect_hausd[i,j,y,h], scan_hausd[i,j,y,h]),digits=3)
          minind = (res==min(na.omit(res)))
          res = c(hdcd_hausd[i,j,y,h],inspect_hausd[i,j,y,h], scan_hausd[i,j,y,h])
          
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
          
          res = round(c(hdcd_Kerr[i,j,y,h],inspect_Kerr[i,j,y,h], scan_Kerr[i,j,y,h]),digits=3)
          minind = (abs(res)==min(abs(res)))
          res = c(hdcd_Kerr[i,j,y,h],inspect_Kerr[i,j,y,h], scan_Kerr[i,j,y,h])
          
          for (t in 1:length(res)) {
            if(minind[t]){
              string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
            }else{
              string = sprintf("%s & %.3f", string, res[t])
            }
          }
          
          res = round(c(hdcd_ari[i,j,y,h],inspect_ari[i,j,y,h], scan_ari[i,j,y,h]),digits=3)
          minind = (abs(res)==max(abs(res)))
          res = c(hdcd_ari[i,j,y,h],inspect_ari[i,j,y,h], scan_ari[i,j,y,h])
          for (t in 1:length(res)) {
            if(minind[t]){
              string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
            }else{
              string = sprintf("%s & %.3f", string, res[t])
            }
          }
          
          
          # res = round(1000*c(hdcd_time[i,j,y,h],inspect_time[i,j,y,h], scan_time[i,j,y,h]),digits=3)
          # minind = (res==min(res))
          # res = 1000*c(hdcd_time[i,j,y,h],inspect_time[i,j,y,h], scan_time[i,j,y,h])
          # 
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
    }
  }
  
  printlines = c(printlines, c("\\hline \\\\[-1.8ex]",
                               "\\end{tabular}",
                               "\\end{table}"))
  texfile<-file(sprintf("%s/table.tex", savedir))
  writeLines(printlines, texfile)
  close(texfile)
  
}

