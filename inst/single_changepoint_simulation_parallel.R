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
  savedir = file.path(maindir, sprintf("%s/single",dateandtime))
  dir.create(savedir, showWarnings = FALSE)
  
}


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
ns = c(100)
ps = c(100,1000)
#ps = ns[]
#ps = c(500,1000,2000,5000)
kvals=4
totruns = length(ns)*length(ps)*kvals




kfunc = function(y, n, p){
  k=1
  if(y==1){
    k = 1
  }else if(y==2){
    k = ceiling(p^(1/3))
  }else if(y==3){
    k = ceiling(sqrt(p*log(n^4)))
  }else{
    k = p
  }
  return(k)
}


# FIRST WITH EVEN CHANGES

#cl <- parallel::makeForkCluster(6)
#doParallel::registerDoParallel(cl)


cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
result = foreach(z = 1:N, .combine="cbind",.options.snow = opts) %dopar% {
  rez = c()
  set.seed(z)
  library(HDCD)
  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:kvals) {
        
        # determine k!!
        k = kfunc(y,n,p)
        
        mus = matrix(0, nrow=p, ncol=n)
        noise = matrix(rnorm(n*p), nrow=p, ncol=n)
        eta = round(0.4*n)
        
        
        diff = 0
        if(k<sqrt(p*log(n^4))){
          rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
          diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
          diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
        }
        
        mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)
        
        X = mus+noise
        
        # first inspect
        lambda = sqrt(log(p*log(n))/2) 
        
        a = proc.time()
        res = single_Inspect(X[,], lambda)
        b=proc.time()
        inspect_time = (b-a)[1]+(b-a)[2]
        inspect_res= res$pos
        
        # then hdcd
        a = proc.time()
        res  = single_HDCD (X[,], 2,1, debug =FALSE)
        b=proc.time()
        hdcd_time = (b-a)[1]+(b-a)[2]
        hdcd_res = res$pos
        hdcd_s = res$s
        
        # then scan
        a = proc.time()
        res  = single_Scan (X[,], debug =FALSE)
        b=proc.time()
        scan_time = (b-a)[1]+(b-a)[2]
        scan_res = res$pos
        scan_s = res$s
        
        rez = cbind(rez, c(z,i,j,y,k,inspect_res, inspect_time, hdcd_res, hdcd_time,hdcd_s, scan_res, scan_time, scan_s))

        
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
  inspect_time = array(NA, dim= c(length(ns), length(ps), kvals))
  inspect_mse = array(NA, dim= c(length(ns), length(ps), kvals))
  hdcd_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  hdcd_s = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  hdcd_time = array(NA, dim= c(length(ns), length(ps), kvals))
  hdcd_mse = array(NA, dim= c(length(ns), length(ps), kvals))
  scan_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  scan_s = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  scan_time = array(NA, dim= c(length(ns), length(ps), kvals))
  scan_mse = array(NA, dim= c(length(ns), length(ps), kvals))
  
  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
      #p = ps[j]
      for (y in 1:kvals) {
        subset = (result[2,] == i& result[3,] ==j & result[4,] == y )
        inspect_res[i,j,y,] = result[6, subset]
        inspect_time[i,j,y] = sum(result[7, subset])
        eta = round(n*0.4)
        inspect_mse[i,j,y] = mean((inspect_res[i,j,y,] - eta)^2)
        hdcd_res[i,j,y,] = result[8, subset]
        hdcd_time[i,j,y] = sum(result[9, subset])
        hdcd_s[i,j,y,] = result[10, subset]
        hdcd_mse[i,j,y] = mean((hdcd_res[i,j,y,] - eta)^2)
        
        scan_res[i,j,y,] = result[11, subset]
        scan_time[i,j,y] = sum(result[12, subset])
        scan_s[i,j,y,] = result[13, subset]
        scan_mse[i,j,y] = mean((scan_res[i,j,y,] - eta)^2)
      }
    }
  }
}





inspect_mse - hdcd_mse
inspect_mse
hdcd_mse

#par(mfrow=c(1,2))
#hist(inspect_res[3,3,1,])
#hist(hdcd_res[3,3,1,])

#saving: 
if(save){
  saveRDS(inspect_mse, file=sprintf("%s/inspect_mse.RDA", savedir))
  saveRDS(inspect_time, file=sprintf("%s/inspect_time.RDA", savedir))
  saveRDS(inspect_mse, file=sprintf("%s/inspect_mse.RDA", savedir))
  saveRDS(hdcd_res, file=sprintf("%s/hdcd_res.RDA", savedir))
  saveRDS(hdcd_s, file=sprintf("%s/hdcd_s.RDA", savedir))
  saveRDS(hdcd_time, file=sprintf("%s/hdcd_time.RDA", savedir))
  saveRDS(hdcd_mse, file=sprintf("%s/hdcd_mse.RDA", savedir))
  saveRDS(scan_res, file=sprintf("%s/scan_res.RDA", savedir))
  saveRDS(scan_s, file=sprintf("%s/scan_s.RDA", savedir))
  saveRDS(scan_time, file=sprintf("%s/scan_time.RDA", savedir))
  saveRDS(scan_mse, file=sprintf("%s/scan_mse.RDA", savedir))
  
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
  printlines = c("\\begin{table}[!htbp] \\centering",
                 "\\caption{Change spread evenly across coordinates}",
                 "\\label{}",
                 "\\small",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|ccc|ccc}",
                 "\\hline", 
                 "\\multicolumn{5}{c|}{Parameters} & \\multicolumn{3}{c|}{MSE} &\\multicolumn{3}{c}{Time in miliseconds} \\\\ \\hline ",
                 "$n$ & $p$ & $k$ &  $\\eta$ & $\\phi$ & \\text{HDCD} & \\text{Inspect} & \\text{Scan} & \\text{HDCD} & \\text{Inspect} & \\text{Scan} \\\\", 
                 "\\hline \\")
  
  for (i in 1:length(ns)) {
    
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:kvals) {
        k = kfunc(y,n,p)
        eta = round(0.4*n)
        rootnorm = 0
        if(k<sqrt(p*log(n^4))){
          rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
        }
        string = sprintf("%d & %d & %d & %d & %.2f", n, p, k, eta, rootnorm)
        
        res = round(c(hdcd_mse[i,j,y],inspect_mse[i,j,y], scan_mse[i,j,y]),digits=3)
        minind = (res==min(res))
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        
        res = round(c(hdcd_time[i,j,y],inspect_time[i,j,y], scan_time[i,j,y])/N*1000,digits=3)
        minind = (res==min(res))
        
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
  texfile<-file(sprintf("%s/table_even.tex", savedir))
  writeLines(printlines, texfile)
  close(texfile)
  
}


# THEN UEVEN CHANGES



cl <- makeCluster(num_cores,type="SOCK")
registerDoSNOW(cl)
pb <- txtProgressBar(max = N, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

result_uneven = foreach(z = 1:N, .combine="cbind",.options.snow = opts) %dopar% {
  set.seed(z)
  rez = c()
  library(HDCD)
  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:kvals) {
        
        # determine k!!
        k = kfunc(y,n,p)
        
        mus = matrix(0, nrow=p, ncol=n)
        noise = matrix(rnorm(n*p), nrow=p, ncol=n)
        eta = round(0.4*n)
        
        
        diff = c(rep(1, k)/(sqrt(1:k)), rep(0, p-k))
        diff = diff / norm(diff, type="2")
        if(k<sqrt(p*log(n^4))){
          rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
          diff = rootnorm * diff
        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
          diff = rootnorm * diff
        }
        
        mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)
        
        X = mus+noise
        
        # first inspect
        lambda = sqrt(log(p*log(n))/2) 
        
        a = proc.time()
        res = single_Inspect(X[,], lambda)
        b=proc.time()
        inspect_time = (b-a)[1]+(b-a)[2]
        inspect_res= res$pos
        
        # then hdcd
        a = proc.time()
        res  = single_HDCD (X[,], 2,1, debug =FALSE)
        b=proc.time()
        hdcd_time = (b-a)[1]+(b-a)[2]
        hdcd_res = res$pos
        hdcd_s = res$s
        
        # then scan
        a = proc.time()
        res  = single_Scan (X[,], debug =FALSE)
        b=proc.time()
        scan_time = (b-a)[1]+(b-a)[2]
        scan_res = res$pos
        scan_s = res$s
        

        rez = cbind(rez, c(z,i,j,y,k,inspect_res, inspect_time, hdcd_res, hdcd_time,hdcd_s, scan_res, scan_time, scan_s))        
        
      }
    }
  }
  
  rez
}
# bb = proc.time()
# print(bb-aa)
close(pb)
stopCluster(cl) 


{
  inspect_res_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  inspect_time_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
  inspect_mse_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
  hdcd_res_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  hdcd_s_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  hdcd_time_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
  hdcd_mse_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
  scan_res_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  scan_s_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  scan_time_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
  scan_mse_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
  
  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
      #p = ps[j]
      for (y in 1:kvals) {
        subset = (result_uneven[2,] == i& result_uneven[3,] ==j & result_uneven[4,] == y )
        inspect_res_uneven[i,j,y,] = result_uneven[6, subset]
        inspect_time_uneven[i,j,y] = sum(result_uneven[7, subset])
        eta = round(n*0.4)
        inspect_mse_uneven[i,j,y] = mean((inspect_res_uneven[i,j,y,] - eta)^2)
        hdcd_res_uneven[i,j,y,] = result_uneven[8, subset]
        hdcd_time_uneven[i,j,y] = sum(result_uneven[9, subset])
        hdcd_s_uneven[i,j,y,] = result_uneven[10, subset]
        hdcd_mse_uneven[i,j,y] = mean((hdcd_res_uneven[i,j,y,] - eta)^2)
        
        scan_res_uneven[i,j,y,] = result_uneven[11, subset]
        scan_time_uneven[i,j,y] = sum(result_uneven[12, subset])
        scan_s_uneven[i,j,y,] = result_uneven[13, subset]
        scan_mse_uneven[i,j,y] = mean((scan_res_uneven[i,j,y,] - eta)^2)
      }
    }
  }
}





inspect_mse_uneven - hdcd_mse_uneven
inspect_mse_uneven
hdcd_mse_uneven

# par(mfrow=c(1,2))
# hist(inspect_res_uneven[3,3,1,])
# hist(hdcd_res_uneven[3,3,1,])

#saving: 
if(save){
  saveRDS(inspect_mse_uneven, file=sprintf("%s/inspect_mse_uneven.RDA", savedir))
  saveRDS(inspect_time_uneven, file=sprintf("%s/inspect_time_uneven.RDA", savedir))
  saveRDS(inspect_mse_uneven, file=sprintf("%s/inspect_mse_uneven.RDA", savedir))
  saveRDS(hdcd_res_uneven, file=sprintf("%s/hdcd_res_uneven.RDA", savedir))
  saveRDS(hdcd_s_uneven, file=sprintf("%s/hdcd_s_uneven.RDA", savedir))
  saveRDS(hdcd_time_uneven, file=sprintf("%s/hdcd_time_uneven.RDA", savedir))
  saveRDS(hdcd_mse_uneven, file=sprintf("%s/hdcd_mse_uneven.RDA", savedir))
  saveRDS(scan_res_uneven, file=sprintf("%s/scan_res_uneven.RDA", savedir))
  saveRDS(scan_s_uneven, file=sprintf("%s/scan_s_uneven.RDA", savedir))
  saveRDS(scan_time_uneven, file=sprintf("%s/scan_time_uneven.RDA", savedir))
  saveRDS(scan_mse_uneven, file=sprintf("%s/scan_mse_uneven.RDA", savedir))
  

}

# creating table:
if(save){
  # output latex table
  printlines = c("\\begin{table}[!htbp] \\centering",
                 "\\caption{Change spread unevenly across coordinates}",
                 "\\label{}",
                 "\\small",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|ccc|ccc}",
                 "\\hline", 
                 "\\multicolumn{5}{c|}{Parameters} & \\multicolumn{3}{c|}{MSE} &\\multicolumn{3}{c}{Time in miliseconds} \\\\ \\hline ",
                 "$n$ & $p$ & $k$ &  $\\eta$ & $\\phi$ & \\text{HDCD} & \\text{Inspect} & \\text{Scan} & \\text{HDCD} & \\text{Inspect} & \\text{Scan} \\\\", 
                 "\\hline \\")
  
  for (i in 1:length(ns)) {
    
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:kvals) {
        k = kfunc(y,n,p)
        eta = round(0.4*n)
        rootnorm = 0
        if(k<sqrt(p*log(n^4))){
          rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
        }
        string = sprintf("%d & %d & %d & %d & %.2f", n, p, k, eta, rootnorm)
        
        res = round(c(hdcd_mse_uneven[i,j,y],inspect_mse_uneven[i,j,y], scan_mse_uneven[i,j,y]),digits=3)
        minind = (res==min(res))
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }

        
        res = round(c(hdcd_time_uneven[i,j,y],inspect_time_uneven[i,j,y], scan_time_uneven[i,j,y])/N*1000,digits=3)
        minind = (res==min(res))
        
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
  texfile<-file(sprintf("%s/table_uneven.tex", savedir))
  writeLines(printlines, texfile)
  close(texfile)
  
}








# 
# 
# 
# # decreasing magnitude
# # 
# # N = 100
# # 
# # set.seed(100)
# # 
# # #ns = c(500,1000,2000)
# # ns = c(50,100,500)
# # #ps = c(1000,2000,10000)
# # ps = ns[]
# # kvals=4
# # totruns = length(ns)*length(ps)*kvals
# 
# {
#   inspect_res_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
#   inspect_time_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
#   inspect_mse_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
#   hdcd_res_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
#   hdcd_s_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
#   hdcd_time_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
#   hdcd_mse_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
# }
# for (i in 1:length(ns)) {
#   n = ns[i]
#   for(j in 1:length(ps)){
#     p = ps[j]
#     for (y in 1:kvals) {
#       
#       # determine k!!
#       k= 1
#       if(y==1){
#         k = 2
#       }else if(y==2){
#         k = ceiling(p^(1/3))
#       }else if(y==3){
#         k = ceiling(sqrt(p*log(n^4)))
#       }else{
#         k = p
#       }
#       
#       
#       
#       tottimeinspect = 0
#       tottimehdcd = 0
#       
#       for (z in 1:N) {
#         mus = matrix(0, nrow=p, ncol=n)
#         noise = matrix(rnorm(n*p), nrow=p, ncol=n)
#         eta = round(0.4*n)
#         
#         
#         diff = c(rep(1, k)/(sqrt(1:k)), rep(0, p-k))
#         diff = diff / norm(diff, type="2")
#         
#         if(k<sqrt(p*log(n^4))){
#           rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
#           diff = rootnorm *diff
#         }else{
#           rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
#           diff = rootnorm * diff
#         }
#         
#         mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)
#         
#         X = mus+noise
#         
#         # first inspect
#         lambda = sqrt(log(p*log(n))/2) 
#         
#         a = proc.time()
#         res = single_Inspect(X, lambda)
#         b=proc.time()
#         tottimeinspect = tottimeinspect + (b-a)[1]+(b-a)[2]
#         inspect_res_uneven[i,j,y,z] = res$pos
#         
#         # then hdcd
#         a = proc.time()
#         res  = single_HDCD (X, 2,1, debug =FALSE)
#         b=proc.time()
#         tottimehdcd = tottimehdcd + (b-a)[1]+(b-a)[2]
#         hdcd_res_uneven[i,j,y,z] = res$pos
#         hdcd_s_uneven[i,j,y,z] = res$s
#         
#       }
#       inspect_time_uneven[i,j,y] = tottimeinspect
#       inspect_mse_uneven[i,j,y] = mean((inspect_res_uneven[i,j,y,]-eta)^2)
#       
#       
#       hdcd_time_uneven[i,j,y] = tottimehdcd
#       hdcd_mse_uneven[i,j,y] = mean((hdcd_res_uneven[i,j,y,]-eta)^2)
#       
#       
#       
#     }
#   }
# }
# 
# #saving
# if(save){
#   saveRDS(inspect_mse_uneven, file=sprintf("%s/inspect_mse_uneven.RDA", savedir))
#   saveRDS(inspect_time_uneven, file=sprintf("%s/inspect_time_uneven.RDA", savedir))
#   saveRDS(inspect_mse_uneven, file=sprintf("%s/inspect_mse_uneven.RDA", savedir))
#   saveRDS(hdcd_res_uneven, file=sprintf("%s/hdcd_res_uneven.RDA", savedir))
#   saveRDS(hdcd_s_uneven, file=sprintf("%s/hdcd_s_uneven.RDA", savedir))
#   saveRDS(hdcd_time_uneven, file=sprintf("%s/hdcd_time_uneven.RDA", savedir))
#   saveRDS(hdcd_mse_uneven, file=sprintf("%s/hdcd_mse_uneven.RDA", savedir))
#   
# }
# # saving table:
# if(save){
#   # output latex table
#   printlines = c("\\begin{table}[!htbp] \\centering",
#                  "\\caption{Change spread unevenly across coordinates}",
#                  "\\label{}",
#                  "\\small",
#                  "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|cc|cc}",
#                  "\\\\[-1.8ex]\\hline", 
#                  "\\hline \\\\[-1.8ex]",
#                  "$n$ & $p$ & $k$ &  $\\eta$ & $\\phi$ & \\text{HDCD MSE} & \\text{Inspect MSE} & \\text{HDCD (ms)} & \\text{Inspect (ms)} \\\\", 
#                  "\\hline \\\\[-1.8ex]")
#   
#   for (i in 1:length(ns)) {
#     
#     n = ns[i]
#     for(j in 1:length(ps)){
#       p = ps[j]
#       for (y in 1:kvals) {
#         k = 1
#         if(y==1){
#           k = 2
#         }else if(y==2){
#           k = ceiling(p^(1/3))
#         }else if(y==3){
#           k = ceiling(sqrt(p*log(n^4)))
#         }else{
#           k = p
#         }
#         eta = round(0.4*n)
#         rootnorm = 0
#         if(k<sqrt(p*log(n^4))){
#           rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
#         }else{
#           rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
#         }
#         string = sprintf("%d & %d & %d & %d & %.2f", n, p, k, eta, rootnorm)
#         
#         res = c(hdcd_mse_uneven[i,j,y],inspect_mse_uneven[i,j,y])
#         minind = which.min(res)
#         if(minind>1){
#           for (t in 1:(minind-1)) {
#             string = sprintf("%s & %.4f", string, res[t])
#           }
#         }
#         string = sprintf("%s & \\textbf{%.4f} ", string, res[minind])
#         if(minind<length(res)){
#           for (t in (minind+1):length(res)) {
#             string = sprintf("%s & %.4f", string, res[t])
#           }
#         }
#         
#         res = c(hdcd_time_uneven[i,j,y],inspect_time_uneven[i,j,y])/N*1000
#         minind = which.min(res)
#         if(minind>1){
#           for (t in 1:(minind-1)) {
#             string = sprintf("%s & %.2f", string, res[t])
#           }
#         }
#         string = sprintf("%s & \\textbf{%.2f} ", string, res[minind])
#         if(minind<length(res)){
#           for (t in (minind+1):length(res)) {
#             string = sprintf("%s & %.2f", string, res[t])
#           }
#         }
#         string = sprintf("%s \\\\", string)
#         printlines = c(printlines, string)
#       }
#     }
#   }
#   
#   printlines = c(printlines, c("\\hline \\\\[-1.8ex]",
#                                "\\end{tabular}",
#                                "\\end{table}"))
#   texfile<-file(sprintf("%s/table_uneven.tex", savedir))
#   writeLines(printlines, texfile)
#   close(texfile)
#   
# }
# 
# 
# 
# inspect_mse_uneven - hdcd_mse_uneven
# 
# inspect_time_uneven / hdcd_time_uneven
# 
# 
# 
# 
# 
# 
# # p = 1000
# # n =500
# # mus = matrix(0, nrow=p, ncol=n)
# # noise = matrix(rnorm(n*p), nrow=p, ncol=n)
# # eta = round(0.4*n)
# # k = p
# # 
# # diff = 0
# # if(k<sqrt(p*log(n^4))){
# #   rootnorm = 2/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
# #   diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
# # }else{
# #   rootnorm = 2/sqrt(eta)*(p*log(n^4))^(1/4)
# #   diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
# # }
# # 
# # mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)
# # 
# # plot(mus[1,])
# # X = mus+noise
# # plot(X[1,])
# # 
# # xi = 4*sqrt(log(p*n))
# # lambda = sqrt(log(p*log(n))/2) # double check this is what they recc!!!!!!!
# # 
# # res = single_Inspect(X, lambda)
# # res$pos
# # 
# # res2 = single_HDCD (X, 2,2, debug =FALSE)
# # res2$pos
# # res2$s
# # 
# 
