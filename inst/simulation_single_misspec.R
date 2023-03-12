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
#calibrates_path = "/mn/sarpanitu/ansatte-u2/pamoen/project_inspect/simulations_jan_23/2023-01-04--00.01.24/multi/calibrates.RDA"

save = TRUE

if(save){
  dir.create(savedir, showWarnings = FALSE)
  savedir = file.path(maindir, sprintf("%s/single_misspec",dateandtime))
  dir.create(savedir, showWarnings = FALSE)
  
}



N = 1000
num_cores = 6
sparse_const = 3
dense_const = 3
set.seed(1996)
rescale_variance = TRUE

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
ns = c(200)
#ns = c(200)
#ns = c(500)
#ps = c(1000)
ps = c(200)
#ps = c(100)
#ps = c(50,100)
#ps = c(50,100,500,5000)
#ps = ns[]
n = ns
p = ps

# rez = config(7, n, p, mus)
# etas = rez[[1]]
# sparsities = rez[[2]]
# mus = rez[[3]]

etas = round(n/5)

simsbsN = 1000
simsbstoln = 1/simsbsN
#for sbs: 
set.seed(1996)
pis = matrix(nrow = length(ns), ncol =length(ps))
for (i in 1:length(ns)) {
  for (j in 1:length(ps)) {
    pis[i,j] =  single_SBS_calibrate(ns[i],ps[j],simsbsN,simsbstoln,rescale_variance = rescale_variance,debug=FALSE)
    
  }
}

numconfig = 12



## #' @import MASS
#' @import InspectChangepoint
config = function(i,h,n,p){
  # h = 1 sparse
  # h = 2 dense
  #mus = matrix(0, nrow=p, ncol=n)
  noise = NULL
  mus = matrix(0, nrow=p, ncol=n)
  #X = NULL
  #k = 0
  #etas = c()
  sparsities = c()
  sparsity = c()
  Delta = etas
  k = 1
  phi = 1
  
  if(h==1){
    sparsity = "Sparse"
    Delta = etas
    sparsities = sample(1:floor(sqrt(p*log(n))),1)
    k = sparsities[j]
    phi = sparse_const/sqrt(Delta)*sqrt(((k*log(exp(1)*p*log(n)/k^2)+ log(n))))
    
  }else{
    sparsity = "Dense"
    Delta = etas
    sparsities = sample((ceiling(sqrt(p*log(n)))):p,1)
    k = sparsities[j]
    phi = dense_const/sqrt(Delta)*(p*log(n))^(1/4)
  }
  
  coords = sample(1:p, sparsities)
  diff = sample(c(-1,1), sparsities, replace=TRUE)
  diff = diff/norm(diff, type="2")
  mus[coords, (etas+1):n] = mus[coords, (etas+1):n] +diff*phi
  
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
    L = floor(Delta/2)
    
    for (zz in 1:length(coords)) {
      start = sample(seq(etas-L,etas+L ),1)
      coord = coords[zz]
      hh = mus[coord,n]
      mus[coord,] = 0
      mus[coord, (start+1):n] = hh
      #mus[coord, start:n] = mus[coord, start:n] +diff[zz]*phi
    }
    
    
    
  }
  
  if(i==12){
    #gradual
    noise = matrix(rnorm(n*p), nrow=p, ncol=n)
    L = floor(Delta/2)
    for (zz in 1:length(coords)) {
      start = sample(seq(etas-L,etas+L ),1)
      coord = coords[zz]
      hh = mus[coord,n]
      mus[coord,] = 0
      mus[coord, (etas+L+1):n] = hh
      start = etas - L+1
      stop = etas+L+1
      mus[coord, start:stop] = hh*(1:(stop-start+1))/(stop-start+1)
      #mus[coord, start:n] = mus[coord, start:n] +diff[zz]*phi
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
      #if(y %in% c(11,12) & h==1){
      #  next
      #}
      #noise = matrix(rnorm(n*p), nrow=p, ncol=n)
      conf = config(y, h,n, p)
      eta = conf[[1]]
      sparsities = conf[[2]]
      mus = conf[[3]]
      X = conf[[4]]
      sparsity = conf[[5]]
      #noise = conf[[6]]
      
      #X = noise + mus
      
      rezi = list()
      #rezi[["i"]] = i
      #rezi[["j"]] = j
      rezi[["y"]] = y
      rezi[["h"]] = h
      
      
      # DC
      a = proc.time()
      res_dc = dcbs.alg(X[,],phi = -1,height=1, thr = 0,cp.type=1,temporal=FALSE )
      b=proc.time()
      
      #res_dc$ecp
      rezi["dc_time"] = (b-a)[1]+(b-a)[2]
      if(is.null(res_dc$ecp)){
        rezi["dc_res"] = round(n/2)
      }else{
        rezi["dc_res"]= res_dc$ecp
      }
      #then sbs
      a = proc.time()
      res_sbs = sbs.alg(X[,],height=1, thr = rep(sqrt(pis[i,j]),p),cp.type=1,temporal=FALSE )
      #ress = single_SBS(X,pis[i,j] )
      #ress$pos
      #res_sbs$ecp
      b=proc.time()
      
      rezi["sbs_time"] = (b-a)[1]+(b-a)[2]
      if(!is.numeric(res_sbs$ecp)){
        #rezi["sbs_res"] = round(n/2)
        if(!rescale_variance){
          ress = single_SBS(X[,],pis[i,j] )
        }else{
          ress = single_SBS(rescale.variance(X),pis[i,j] )
        }
        
        rezi["sbs_res"] = ress$pos
        
      }else{
        rezi["sbs_res"] = res_sbs$ecp
      }
      
      
      
      # now rescale variance
      #X = rescale.variance(X)
      X = matrix(rescale_variance(X)$X,byrow = FALSE, ncol = n, nrow = p)
      # first inspect
      lambda = sqrt(log(p*log(n))/2) 
      
      a = proc.time()
      res_inspect = single_Inspect(X[,], lambda)
      b=proc.time()
      rezi["inspect_time"] = (b-a)[1]+(b-a)[2]
      rezi["inspect_res"]= res_inspect$pos
      
      # then hdcd
      a = proc.time()
      res_hdcd  = single_HDCD(X[,], 1.5,1, debug =FALSE)
      b=proc.time()
      rezi["hdcd_time"] = (b-a)[1]+(b-a)[2]
      rezi["hdcd_res"] = res_hdcd$pos
      rezi["hdcd_s"] = res_hdcd$s
      
      # then scan
      # a = proc.time()
      # #res  = single_Scan (X[,], debug =FALSE)
      # b=proc.time()
      # scan_time = (b-a)[1]+(b-a)[2]
      # scan_res = res$pos
      # scan_s = res$s
      
      
      
      #rez = cbind(rez, c(z,i,j,y,k,inspect_res, inspect_time, hdcd_res, hdcd_time,hdcd_s, scan_res, scan_time, scan_s))
      
      # then subset
      a = proc.time()
      res_subset  = SUBSET.normal(X[,])
      b=proc.time()
      rezi["subset_time"] = (b-a)[1]+(b-a)[2]
      
      if(is.null(res_subset$cpt)){
        rezi["subset_res"] = round(n/2)
      }else{
        rezi["subset_res"] = res_subset$cpt
      }
      
      rezi["z"] = z
      rezi["h"] = h
      rezi["y"] = y
      rezi["k"] = k
      #rezi["X"] = X[,]
      #rezi["mu"] = mus
      rezi["eta"] = eta
      
      # rezi = list(z,i,j,y,k,inspect_res, inspect_time, hdcd_res, hdcd_time,hdcd_s, scan_res, scan_time, scan_s, sbs_res, sbs_time, 
      #                   subset_res, subset_time, dc_res, dc_time)
      
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
  inspect_res = array(NA, dim = c(2,numconfig,N) )
  inspect_time = array(0, dim= c(2,numconfig))
  inspect_mse = array(0, dim= c(2,numconfig))
  hdcd_res = array(NA, dim = c(2,numconfig,N) )
  hdcd_s = array(NA, dim = c(2,numconfig,N) )
  hdcd_time = array(0, dim= c(2,numconfig))
  hdcd_mse = array(0, dim= c(2,numconfig))
  scan_res = array(NA, dim = c(2,numconfig,N) )
  scan_s = array(NA, dim = c(2,numconfig,N) )
  scan_time = array(0, dim= c(2,numconfig))
  scan_mse = array(NA, dim= c(2,numconfig))  
  sbs_res = array(NA, dim = c(2,numconfig,N) )
  sbs_time = array(0, dim= c(2,numconfig))
  sbs_mse = array(0, dim= c(2,numconfig))
  subset_res = array(NA, dim = c(2,numconfig,N) )
  subset_time = array(0, dim= c(2,numconfig))
  subset_mse = array(0, dim= c(2,numconfig))  
  dc_res = array(NA, dim = c(2,numconfig,N) )
  dc_time = array(0, dim= c(2,numconfig))
  dc_mse = array(0, dim= c(2,numconfig))
  
  for (z in 1:N) {
    list = result[[z]]
    len = length(list)
    
    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      h = sublist[["h"]]
      eta = sublist[["eta"]]
      
      inspect_res[h,y,z] = sublist["inspect_res"][[1]]
      inspect_time[h,y] = inspect_time[h,y] +  sublist["inspect_time"][[1]]/N
      inspect_mse[h,y] = inspect_mse[h,y] + (inspect_res[h,y,z] - eta)^2/N
      
      hdcd_res[h,y,z] = sublist["hdcd_res"][[1]]
      hdcd_time[h,y] = hdcd_time[h,y] + sublist["hdcd_time"][[1]]/N
      hdcd_mse[h,y] = hdcd_mse[h,y] + (hdcd_res[h,y,z] - eta)^2/N
      
      sbs_res[h,y,z] = sublist["sbs_res"][[1]]
      sbs_time[h,y] = sbs_time[h,y] +  sublist["sbs_time"][[1]]/N
      sbs_mse[h,y] = sbs_mse[h,y] + (sbs_res[h,y,z] - eta)^2/N
      
      subset_res[h,y,z] = sublist["subset_res"][[1]]
      subset_time[h,y] = subset_time[h,y] +  sublist["subset_time"][[1]]/N
      subset_mse[h,y] = subset_mse[h,y] + (subset_res[h,y,z] - eta)^2/N
      
      dc_res[h,y,z] = sublist["dc_res"][[1]]
      dc_time[h,y] = dc_time[h,y] + sublist["dc_time"][[1]]/N
      dc_mse[h,y] = dc_mse[h,y] + (dc_res[h,y,z] - eta)^2/N
      
      
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
  #       inspect_mse[i,j,y] = mean((inspect_res[i,j,y,] - eta)^2)
  #       hdcd_res[i,j,y,] = result[8, subset]
  #       hdcd_time[i,j,y] = sum(result[9, subset])
  #       hdcd_s[i,j,y,] = result[10, subset]
  #       hdcd_mse[i,j,y] = mean((hdcd_res[i,j,y,] - eta)^2)
  #       
  #       scan_res[i,j,y,] = result[11, subset]
  #       scan_time[i,j,y] = sum(result[12, subset])
  #       scan_s[i,j,y,] = result[13, subset]
  #       scan_mse[i,j,y] = mean((scan_res[i,j,y,] - eta)^2)
  #       
  #       sbs_res[i,j,y,] = result[14, subset]
  #       sbs_time[i,j,y] = sum(result[15, subset])
  #       sbs_mse[i,j,y] = mean((sbs_res[i,j,y,] - eta)^2)
  #       
  #       subset_res[i,j,y,] = result[16, subset]
  #       subset_time[i,j,y] = sum(result[17, subset])
  #       subset_mse[i,j,y] = mean((subset_res[i,j,y,] - eta)^2)
  #       
  #       dc_res[i,j,y,] = result[18, subset]
  #       dc_time[i,j,y] = sum(result[19, subset])
  #       dc_mse[i,j,y] = mean((dc_res[i,j,y,] - eta)^2)
  #     }
  #   }
  # }
}





inspect_mse - hdcd_mse
inspect_mse
hdcd_mse

par(mfrow=c(1,2))
hist(inspect_res[1,2,4,])
hist(hdcd_res[1,2,4,])

#saving: 
if(save){
  saveRDS(inspect_mse, file=sprintf("%s/inspect_mse.RDA", savedir))
  saveRDS(inspect_time, file=sprintf("%s/inspect_time.RDA", savedir))
  saveRDS(inspect_res, file=sprintf("%s/inspect_res.RDA", savedir))
  saveRDS(hdcd_res, file=sprintf("%s/hdcd_res.RDA", savedir))
  saveRDS(hdcd_s, file=sprintf("%s/hdcd_s.RDA", savedir))
  saveRDS(hdcd_time, file=sprintf("%s/hdcd_time.RDA", savedir))
  saveRDS(hdcd_mse, file=sprintf("%s/hdcd_mse.RDA", savedir))
  saveRDS(scan_res, file=sprintf("%s/scan_res.RDA", savedir))
  saveRDS(scan_s, file=sprintf("%s/scan_s.RDA", savedir))
  saveRDS(scan_time, file=sprintf("%s/scan_time.RDA", savedir))
  saveRDS(scan_mse, file=sprintf("%s/scan_mse.RDA", savedir))
  saveRDS(sbs_mse, file=sprintf("%s/sbs_mse.RDA", savedir))
  saveRDS(sbs_res, file=sprintf("%s/sbs_res.RDA", savedir))
  saveRDS(sbs_time, file=sprintf("%s/sbs_time.RDA", savedir))
  saveRDS(subset_mse, file=sprintf("%s/subset_mse.RDA", savedir))
  saveRDS(subset_res, file=sprintf("%s/subset_res.RDA", savedir))
  saveRDS(subset_time, file=sprintf("%s/subset_time.RDA", savedir))
  saveRDS(dc_mse, file=sprintf("%s/dc_mse.RDA", savedir))
  saveRDS(dc_res, file=sprintf("%s/dc_res.RDA", savedir))
  saveRDS(dc_time, file=sprintf("%s/dc_time.RDA", savedir))
  
  infofile<-file(sprintf("%s/parameters.txt", savedir))
  writeLines(c(sprintf("N = %d", N),
               sprintf("ns = %s", paste(ns, sep=" ", collapse=" ")),
               sprintf("ps = %s", paste(ps, sep=" ", collapse=" ")), 
               sprintf("Sparse constant = %f", sparse_const), 
               sprintf("Dense constant = %f", dense_const),
               sprintf("Even spread = %d", even_spread)),
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
  printlines = c("\\begin{table}[H] \\centering",
                 "\\caption{Single change-point estimation under misspecified model}",
                 "\\label{tablesinglelocmisspec}",
                 "\\small",
                 "\\begin{adjustbox}{width=\\columnwidth}",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} cc|ccccc}",
                 "\\hline", 
                 "\\multicolumn{2}{c|}{Parameters} & \\multicolumn{5}{c|}{MSE} \\\\ \\hline ",
                 "Model & Sparsity & \\text{FAST} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET} & \\text{DC} \\\\", 
                 "\\hline \\")
  
  #for (i in 1:length(ns)) {
    
   # n = ns[i]
   # for(j in 1:length(ps)){
    #  p = ps[j]
  
    for (y in 1:numconfig) {
      for(h in 1:2){
        #k = kfunc(y,n,p)
        if(h==1){
          dens = "Sparse"
        }else{
          dens = "Dense"
        }
        eta = round(0.2*n)
        rootnorm = 0
        if(k<sqrt(p*log(n))){
          rootnorm = sparse_const/sqrt(eta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n))^(1/4)
        }
        string = sprintf("%s & %s", models[y], dens)
        
        #res = round(c(hdcd_mse[i,j,y],inspect_mse[i,j,y], scan_mse[i,j,y]),digits=3)
        res = round(c(hdcd_mse[h,y],inspect_mse[h,y], sbs_mse[h,y], subset_mse[h,y], dc_mse[h,y]),digits=1)
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
  
  #printlines = c(printlines, "\\hline \\multicolumn{5}{c|}{Average MSE}")
  #res = round(c(mean(hdcd_slow_det[,,2:kvals]),mean(hdcd_det[,,2:kvals]),mean(pilliat_det[,,2:kvals]),mean(inspect_det[,,2:kvals]), mean(sbs_det[,,2:kvals]), mean(subset_det[,,2:kvals]), mean(dc_det[,,2:kvals])),digits=3)
  #res = round(c(mean(hdcd_mse),mean(inspect_mse), mean(sbs_mse), mean(subset_mse), mean(dc_mse)),digits=3)
  #minind = (res==min(res))
  #string =""
  #for (t in 1:length(res)) {
  #  if(minind[t]){
  #    string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
  #  }else{
  #    string = sprintf("%s & %.3f", string, res[t])
  #  }
  #}
  #string = sprintf("%s", string)
  printlines = c(printlines, string)
  
  printlines = c(printlines, c("\\hline \\\\[-1.8ex]",
                               "\\end{tabular}",
                               "\\end{adjustbox}",
                               "\\end{table}"))
  texfile<-file(sprintf("%s/table_even.tex", savedir))
  writeLines(printlines, texfile)
  close(texfile)
  
}
