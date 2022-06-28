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
  savedir = file.path(maindir, sprintf("%s/multi",dateandtime))
  dir.create(savedir, showWarnings = FALSE)
  
}



N = 100
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
ns = c(100,200)
#ps = c(50,100,500,1000,4000)
ps = c(50,100,500)
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
        res = HDCD (X, 2,2)
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
  rez
}
# bb = proc.time()
# print(bb-aa)
#parallel::stopCluster(cl)
close(pb)
stopCluster(cl) 

{
  inspect_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  inspect_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  inspect_time = array(0, dim= c(length(ns), length(ps), numconfig))
  inspect_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  hdcd_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  hdcd_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  hdcd_time = array(0, dim= c(length(ns), length(ps), numconfig))
  hdcd_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  scan_hausd = array(0, dim = c(length(ns), length(ps), numconfig) )
  scan_Kerr = array(0, dim = c(length(ns), length(ps), numconfig) )
  scan_time = array(0, dim= c(length(ns), length(ps), numconfig))
  scan_ari = array(0, dim= c(length(ns), length(ps), numconfig))
  
  scan_hausd[,,1] = NA
  inspect_hausd[,,1] = NA
  hdcd_hausd[,,1] = NA
  
  for (z in 1:N) {
    list = result[[z]]
    len = length(list)
    
    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      i = sublist[["i"]]
      j = sublist[["j"]]
      
      inspect_Kerr[i,j,y] = inspect_Kerr[i,j,y] + sublist[["inspect_K_error"]]/N
      scan_Kerr[i,j,y] = scan_Kerr[i,j,y] + sublist[["scan_K_error"]]/N
      hdcd_Kerr[i,j,y] = hdcd_Kerr[i,j,y] + sublist[["hdcd_K_error"]]/N
      
      inspect_time[i,j,y] = inspect_time[i,j,y] + sublist[["inspect_time"]]/N
      scan_time[i,j,y] = scan_time[i,j,y] + sublist[["scan_time"]]/N
      hdcd_time[i,j,y] = hdcd_time[i,j,y] + sublist[["hdcd_time"]]/N
      
      inspect_ari[i,j,y] = inspect_ari[i,j,y] + sublist[["inspect_ari"]]/N
      scan_ari[i,j,y] = scan_ari[i,j,y] + sublist[["scan_ari"]]/N
      hdcd_ari[i,j,y] = hdcd_ari[i,j,y] + sublist[["hdcd_ari"]]/N
      
      # if(i==1 & j ==1 & y==1){
      #   print(sublist[["hdcd_ari"]])
      # }
      
      if(y!= 1){
        inspect_hausd[i,j,y] = inspect_hausd[i,j,y] + sublist[["inspect_hausd"]]/N
        scan_hausd[i,j,y] = scan_hausd[i,j,y] + sublist[["scan_hausd"]]/N
        hdcd_hausd[i,j,y] = hdcd_hausd[i,j,y] + sublist[["hdcd_hausd"]]/N
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
                 "$n$ & $p$ & Sparsity & K & \\text{HDCD} & \\text{Inspect} & \\text{Scan}& \\text{HDCD} & \\text{Inspect} & \\text{Scan} & \\text{HDCD} & \\text{Inspect} & \\text{Scan}  & \\text{HDCD} & \\text{Inspect} & \\text{Scan} \\\\", 
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
        
        
        res = round(c(hdcd_hausd[i,j,y],inspect_hausd[i,j,y], scan_hausd[i,j,y]),digits=3)
        minind = (res==min(na.omit(res)))
        res = c(hdcd_hausd[i,j,y],inspect_hausd[i,j,y], scan_hausd[i,j,y])
        
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
        
        res = round(c(hdcd_Kerr[i,j,y],inspect_Kerr[i,j,y], scan_Kerr[i,j,y]),digits=3)
        minind = (abs(res)==min(abs(res)))
        res = c(hdcd_Kerr[i,j,y],inspect_Kerr[i,j,y], scan_Kerr[i,j,y])
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        res = round(c(hdcd_ari[i,j,y],inspect_ari[i,j,y], scan_ari[i,j,y]),digits=3)
        minind = (abs(res)==max(abs(res)))
        res = c(hdcd_ari[i,j,y],inspect_ari[i,j,y], scan_ari[i,j,y])
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        
        res = round(1000*c(hdcd_time[i,j,y],inspect_time[i,j,y], scan_time[i,j,y]),digits=3)
        minind = (res==min(res))
        res = 1000*c(hdcd_time[i,j,y],inspect_time[i,j,y], scan_time[i,j,y])
        
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

