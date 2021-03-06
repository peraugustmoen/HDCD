#### Simulation for a single change-point

## same magnitude in all coords

maindir = "/Users/peraugust/OneDrive - Universitetet i Oslo/project_inspect/simulations/Simulations_HDCD"
dateandtime = gsub(" ", "--",as.character(Sys.time()))
dateandtime = gsub(":", ".", dateandtime)
savedir = file.path(maindir, dateandtime)
save = TRUE

if(save){
  dir.create(savedir, showWarnings = FALSE)
  
}

library(HDCD)

N = 100
sparse_const = 3
dense_const = 3
set.seed(1996)

#ns = c(500,1000,2000)
ns = c(100,500,1000)
#ps = c(1000,2000,10000)
ps = ns[]
kvals=4
totruns = length(ns)*length(ps)*kvals

{
inspect_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
inspect_time = array(NA, dim= c(length(ns), length(ps), kvals))
inspect_mse = array(NA, dim= c(length(ns), length(ps), kvals))
hdcd_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
hdcd_s = array(NA, dim = c(length(ns), length(ps), kvals,N) )
hdcd_time = array(NA, dim= c(length(ns), length(ps), kvals))
hdcd_mse = array(NA, dim= c(length(ns), length(ps), kvals))
}

for (i in 1:length(ns)) {
  n = ns[i]
  for(j in 1:length(ps)){
    p = ps[j]
    for (y in 1:kvals) {
      
      # determine k!!
      k=1
      if(y==1){
        k = 2
      }else if(y==2){
        k = ceiling(p^(1/3))
      }else if(y==3){
        k = ceiling(sqrt(p*log(n^4)))
      }else{
        k = p
      }
      
      
      
      tottimeinspect = 0
      tottimehdcd = 0
      
      for (z in 1:N) {
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
        res = single_Inspect(X, lambda)
        b=proc.time()
        tottimeinspect = tottimeinspect + (b-a)[1]+(b-a)[2]
        inspect_res[i,j,y,z] = res$pos
        
        # then hdcd
        a = proc.time()
        res  = single_HDCD (X, 2,1, debug =FALSE)
        b=proc.time()
        tottimehdcd = tottimehdcd + (b-a)[1]+(b-a)[2]
        hdcd_res[i,j,y,z] = res$pos
        hdcd_s[i,j,y,z] = res$s
        
      }
      inspect_time[i,j,y] = tottimeinspect
      inspect_mse[i,j,y] = mean((inspect_res[i,j,y,]-eta)^2)
      
      
      hdcd_time[i,j,y] = tottimehdcd
      hdcd_mse[i,j,y] = mean((hdcd_res[i,j,y,]-eta)^2)
      
      
      
    }
  }
}

inspect_mse - hdcd_mse
inspect_mse
hdcd_mse

par(mfrow=c(1,2))
hist(inspect_res[3,3,1,])
hist(hdcd_res[3,3,1,])

#saving: 
if(save){
  saveRDS(inspect_mse, file=sprintf("%s/inspect_mse.RDA", savedir))
  saveRDS(inspect_time, file=sprintf("%s/inspect_time.RDA", savedir))
  saveRDS(inspect_mse, file=sprintf("%s/inspect_mse.RDA", savedir))
  saveRDS(hdcd_res, file=sprintf("%s/hdcd_res.RDA", savedir))
  saveRDS(hdcd_s, file=sprintf("%s/hdcd_s.RDA", savedir))
  saveRDS(hdcd_time, file=sprintf("%s/hdcd_time.RDA", savedir))
  saveRDS(hdcd_mse, file=sprintf("%s/hdcd_mse.RDA", savedir))
  
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
                 "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|cc|cc}",
                 "\\\\[-1.8ex]\\hline", 
                 "\\hline \\\\[-1.8ex]",
                 "$n$ & $p$ & $k$ &  $\\eta$ & $\\phi$ & \\text{HDCD MSE} & \\text{Inspect MSE} & \\text{HDCD Time (ms)} & \\text{Inspect Time (ms)} \\\\", 
                 "\\hline \\\\[-1.8ex]")
  
  for (i in 1:length(ns)) {
    
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:kvals) {
        k = 1
        if(y==1){
          k = 2
        }else if(y==2){
          k = ceiling(p^(1/3))
        }else if(y==3){
          k = ceiling(sqrt(p*log(n^4)))
        }else{
          k = p
        }
        eta = round(0.4*n)
        rootnorm = 0
        if(k<sqrt(p*log(n^4))){
          rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
        }
        string = sprintf("%d & %d & %d & %d & %.2f", n, p, k, eta, rootnorm)
        
        res = c(hdcd_mse[i,j,y],inspect_mse[i,j,y])
        minind = which.min(res)
        if(minind>1){
          for (t in 1:(minind-1)) {
            string = sprintf("%s & %.4f", string, res[t])
          }
        }
        string = sprintf("%s & \\textbf{%.4f} ", string, res[minind])
        if(minind<length(res)){
          for (t in (minind+1):length(res)) {
            string = sprintf("%s & %.4f", string, res[t])
          }
        }
        
        res = c(hdcd_time[i,j,y],inspect_time[i,j,y])/N*1000
        minind = which.min(res)
        if(minind>1){
          for (t in 1:(minind-1)) {
            string = sprintf("%s & %.2f", string, res[t])
          }
        }
        string = sprintf("%s & \\textbf{%.2f} ", string, res[minind])
        if(minind<length(res)){
          for (t in (minind+1):length(res)) {
            string = sprintf("%s & %.2f", string, res[t])
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


# decreasing magnitude
# 
# N = 100
# 
# set.seed(100)
# 
# #ns = c(500,1000,2000)
# ns = c(50,100,500)
# #ps = c(1000,2000,10000)
# ps = ns[]
# kvals=4
# totruns = length(ns)*length(ps)*kvals

{
inspect_res_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
inspect_time_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
inspect_mse_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
hdcd_res_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
hdcd_s_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
hdcd_time_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
hdcd_mse_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
}
for (i in 1:length(ns)) {
  n = ns[i]
  for(j in 1:length(ps)){
    p = ps[j]
    for (y in 1:kvals) {
      
      # determine k!!
      k= 1
      if(y==1){
        k = 2
      }else if(y==2){
        k = ceiling(p^(1/3))
      }else if(y==3){
        k = ceiling(sqrt(p*log(n^4)))
      }else{
        k = p
      }
      
      
      
      tottimeinspect = 0
      tottimehdcd = 0
      
      for (z in 1:N) {
        mus = matrix(0, nrow=p, ncol=n)
        noise = matrix(rnorm(n*p), nrow=p, ncol=n)
        eta = round(0.4*n)
        
        
        diff = c(rep(1, k)/(sqrt(1:k)), rep(0, p-k))
        diff = diff / norm(diff, type="2")
      
        if(k<sqrt(p*log(n^4))){
          rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
          diff = rootnorm *diff
        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
          diff = rootnorm * diff
        }
        
        mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)
        
        X = mus+noise
        
        # first inspect
        lambda = sqrt(log(p*log(n))/2) 
        
        a = proc.time()
        res = single_Inspect(X, lambda)
        b=proc.time()
        tottimeinspect = tottimeinspect + (b-a)[1]+(b-a)[2]
        inspect_res_uneven[i,j,y,z] = res$pos
        
        # then hdcd
        a = proc.time()
        res  = single_HDCD (X, 2,1, debug =FALSE)
        b=proc.time()
        tottimehdcd = tottimehdcd + (b-a)[1]+(b-a)[2]
        hdcd_res_uneven[i,j,y,z] = res$pos
        hdcd_s_uneven[i,j,y,z] = res$s
        
      }
      inspect_time_uneven[i,j,y] = tottimeinspect
      inspect_mse_uneven[i,j,y] = mean((inspect_res_uneven[i,j,y,]-eta)^2)
      
      
      hdcd_time_uneven[i,j,y] = tottimehdcd
      hdcd_mse_uneven[i,j,y] = mean((hdcd_res_uneven[i,j,y,]-eta)^2)
      
      
      
    }
  }
}

#saving
if(save){
  saveRDS(inspect_mse_uneven, file=sprintf("%s/inspect_mse_uneven.RDA", savedir))
  saveRDS(inspect_time_uneven, file=sprintf("%s/inspect_time_uneven.RDA", savedir))
  saveRDS(inspect_mse_uneven, file=sprintf("%s/inspect_mse_uneven.RDA", savedir))
  saveRDS(hdcd_res_uneven, file=sprintf("%s/hdcd_res_uneven.RDA", savedir))
  saveRDS(hdcd_s_uneven, file=sprintf("%s/hdcd_s_uneven.RDA", savedir))
  saveRDS(hdcd_time_uneven, file=sprintf("%s/hdcd_time_uneven.RDA", savedir))
  saveRDS(hdcd_mse_uneven, file=sprintf("%s/hdcd_mse_uneven.RDA", savedir))

}
# saving table:
if(save){
  # output latex table
  printlines = c("\\begin{table}[!htbp] \\centering",
                 "\\caption{Change spread unevenly across coordinates}",
                 "\\label{}",
                 "\\small",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|cc|cc}",
                 "\\\\[-1.8ex]\\hline", 
                 "\\hline \\\\[-1.8ex]",
                 "$n$ & $p$ & $k$ &  $\\eta$ & $\\phi$ & \\text{HDCD MSE} & \\text{Inspect MSE} & \\text{HDCD Time (ms)} & \\text{Inspect Time (ms)} \\\\", 
                 "\\hline \\\\[-1.8ex]")
  
  for (i in 1:length(ns)) {
    
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:kvals) {
        k = 1
        if(y==1){
          k = 2
        }else if(y==2){
          k = ceiling(p^(1/3))
        }else if(y==3){
          k = ceiling(sqrt(p*log(n^4)))
        }else{
          k = p
        }
        eta = round(0.4*n)
        rootnorm = 0
        if(k<sqrt(p*log(n^4))){
          rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
        }
        string = sprintf("%d & %d & %d & %d & %.2f", n, p, k, eta, rootnorm)
        
        res = c(hdcd_mse_uneven[i,j,y],inspect_mse_uneven[i,j,y])
        minind = which.min(res)
        if(minind>1){
          for (t in 1:(minind-1)) {
            string = sprintf("%s & %.4f", string, res[t])
          }
        }
        string = sprintf("%s & \\textbf{%.4f} ", string, res[minind])
        if(minind<length(res)){
          for (t in (minind+1):length(res)) {
            string = sprintf("%s & %.4f", string, res[t])
          }
        }
        
        res = c(hdcd_time_uneven[i,j,y],inspect_time_uneven[i,j,y])/N*1000
        minind = which.min(res)
        if(minind>1){
          for (t in 1:(minind-1)) {
            string = sprintf("%s & %.2f", string, res[t])
          }
        }
        string = sprintf("%s & \\textbf{%.2f} ", string, res[minind])
        if(minind<length(res)){
          for (t in (minind+1):length(res)) {
            string = sprintf("%s & %.2f", string, res[t])
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



inspect_mse_uneven - hdcd_mse_uneven

inspect_time_uneven / hdcd_time_uneven






p = 1000
n =500
mus = matrix(0, nrow=p, ncol=n)
noise = matrix(rnorm(n*p), nrow=p, ncol=n)
eta = round(0.4*n)
k = 30

diff = 0
if(k<sqrt(p*log(n^4))){
  rootnorm = 3/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
  diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
}else{
  rootnorm = 3/sqrt(eta)*(p*log(n^4))^(1/4)
  diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
}

mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)

plot(mus[1,])
X = mus+noise
plot(X[1,])

xi = 4*sqrt(log(p*n))
lambda = sqrt(log(p*log(n))/2) # double check this is what they recc!!!!!!!

system.time({res = single_Inspect(X, lambda)})
res$pos

system.time({res2 = single_HDCD (X, 2,1, debug =FALSE)})
res2$pos
res2$s

system.time({res3 = single_Scan(X)})
res3$pos
res3$s

