# hdcd calibration
# this is to find the best constants and penalty funcs to estimate a single change point
# seems like method 1 (with n^4 and + on the log) with 
# dense_const = 1.5 (const1)
# sparse_const =  1.5 (const2)
# works very well

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
ns = c(200,500)
ps = c(100,1000,5000)
#ps = ns[]
#ps = c(500,1000,2000,5000)
kvals=4
totruns = length(ns)*length(ps)*kvals

const1 = seq(1,2, by = 0.1)
const2 = seq(0.5,2.5, by = 0.1)

# const1 = c(1,2)
# const2 = c(1,2)


N = 36*10
num_cores = 6
sparse_const = 2.5
dense_const = 2.5
set.seed(1996)

#ns = c(500,1000,2000)
#ns = c(100,500,1000)
#ns = c(500,1000,2000)
#ns = c(100,500,1000)
#ps = c(1000,2000,10000)
#ns = c(200,500)



kfunc = function(y, n, p){
  k=1
  if(y==1){
    k = 1
  }else if(y==2){
    k = ceiling(p^(1/3))
  }else if(y==3){
    k = ceiling(sqrt(p*log(n)))
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
# result = foreach(z = 1:N, .combine="cbind",.options.snow = opts) %dopar% {
result = foreach(z = 1:N,.options.snow = opts) %dopar% {
  rez = list()
  count = 1
  set.seed(z)
  library(HDCD)
  source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/SUBSET_normal.R")
  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:kvals) {
        
        # determine k!!
        k = kfunc(y,n,p)
        
        mus = matrix(0, nrow=p, ncol=n)
        noise = matrix(rnorm(n*p), nrow=p, ncol=n)
        #eta = round(0.2*n)
        eta = max(1, round(0.5*runif(1)*n))
        
        
        diff = 0
        if(k<sqrt(p*log(n))){
          rootnorm = sparse_const/sqrt(eta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
          diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n))^(1/4)
          diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
        }
        
        mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)
        
        X = mus+noise
        
        X = rescale.variance(X)
        # first inspect
        
        # logn4 and plus
        mat = array(NA, dim = c(4,length(const1),length(const2)))
        # # logn4 and max
        # mat2 = matrix(NA, nrow = length(const1), ncol = length(const2))
        # # logn and plus
        # mat3 = matrix(NA, nrow = length(const1), ncol = length(const2))
        # # logn and max
        # mat4 = matrix(NA, nrow = length(const1), ncol = length(const2))
        
        # 5 & 6 -- other formula??
        
        for (cc1 in 1:length(const1)) {
          for (cc2 in 1:length(const2)) {
            res  = single_HDCD (X[,], const1[cc1],const2[cc2], debug =FALSE)
            mat[1,cc1, cc2] = res$pos 
            res  = single_HDCD_two(X[,], const1[cc1],const2[cc2], debug =FALSE)
            mat[2,cc1, cc2] = res$pos 
            # res  = single_HDCD_three(X[,], const1[cc1],const2[cc2], debug =FALSE)
            # mat[3,cc1, cc2] = res$pos 
            # res  = single_HDCD_four(X[,], const1[cc1],const2[cc2], debug =FALSE)
            # mat[4,cc1, cc2] = res$pos 
            
            
          }
        }
        
        
        #rez = cbind(rez, c(z,i,j,y,k,inspect_res, inspect_time, hdcd_res, hdcd_time,hdcd_s, scan_res, scan_time, scan_s, sbs_res, sbs_time, 
        #                   subset_res, subset_time, dc_res, dc_time))
        rez[[count]] = list()
        rez[[count]][[1]] = c(z,i,j,y,k, eta)
        rez[[count]][[2]] = mat
        count = count+1
        

        
        
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
  
  mse1 = array(0, dim=c(length(ns), length(ps), kvals, length(const1), length(const2)))
  mse2 = array(0, dim=c(length(ns), length(ps), kvals, length(const1), length(const2)))
  mse3 = array(0, dim=c(length(ns), length(ps), kvals, length(const1), length(const2)))
  mse4 = array(0, dim=c(length(ns), length(ps), kvals, length(const1), length(const2)))
  
  mse_avg1 = array(0, dim=c(length(const1), length(const2)))
  mse_avg2 = array(0, dim=c(length(const1), length(const2)))
  mse_avg3 = array(0, dim=c(length(const1), length(const2)))
  mse_avg4 = array(0, dim=c(length(const1), length(const2)))
  
  est1 = array(0, dim=c(length(ns), length(ps), kvals, length(const1), length(const2),N))
  est2 = array(0, dim=c(length(ns), length(ps), kvals, length(const1), length(const2),N))
  est3 = array(0, dim=c(length(ns), length(ps), kvals, length(const1), length(const2),N))
  est4 = array(0, dim=c(length(ns), length(ps), kvals, length(const1), length(const2),N))
  
  
  for(rr1 in result){
    for (run in rr1) {
      
      
      eta = run[[1]][6]
      ii = run[[1]][2]
      jj = run[[1]][3]
      yy = run[[1]][4]
      zz = run[[1]][1]
      eta = run[[1]][6]
      
      mat = run[[2]]
      
      for (cc1 in 1:length(const1)) {
        for (cc2 in 1:length(const2)) {
          mse1[ii, jj, yy, cc1,cc2] = mse1[ii, jj, yy, cc1,cc2] + (mat[1, cc1, cc2] - eta)^2
          mse2[ii, jj, yy, cc1,cc2] = mse2[ii, jj, yy, cc1,cc2] + (mat[2, cc1, cc2] - eta)^2
          mse3[ii, jj, yy, cc1,cc2] = mse3[ii, jj, yy, cc1,cc2] + (mat[3, cc1, cc2] - eta)^2
          mse4[ii, jj, yy, cc1,cc2] = mse4[ii, jj, yy, cc1,cc2] + (mat[4, cc1, cc2] - eta)^2
          
          mse_avg1[cc1,cc2] = mse_avg1[cc1,cc2] + (mat[1, cc1, cc2] - eta)^2/ns[ii]
          mse_avg2[cc1,cc2] = mse_avg2[cc1,cc2] + (mat[2, cc1, cc2] - eta)^2/ns[ii]
          mse_avg3[cc1,cc2] = mse_avg3[cc1,cc2] + (mat[3, cc1, cc2] - eta)^2/ns[ii]
          mse_avg4[cc1,cc2] = mse_avg4[cc1,cc2] + (mat[4, cc1, cc2] - eta)^2/ns[ii]
          
          # if(yy ==1 || yy==2 ){
          #   mse_avg1[cc1,cc2] = mse_avg1[cc1,cc2] + (mat[1, cc1, cc2] - eta)^2
          # }
          
          est1[ii,jj,yy,cc1,cc2,zz] = mat[1, cc1, cc2]
          est2[ii,jj,yy,cc1,cc2,zz] = mat[2, cc1, cc2]
          est3[ii,jj,yy,cc1,cc2,zz] = mat[3, cc1, cc2]
          est4[ii,jj,yy,cc1,cc2,zz] = mat[4, cc1, cc2]
        }
        
      }
    }
  }
  mse1 = mse1/N
  mse2 = mse2/N
  mse3 = mse3/N
  mse4 = mse4/N
  
  mse_avg1 = mse_avg1 / (N * length(ns)*length(ps)*kvals)
  mse_avg2 = mse_avg2 / (N * length(ns)*length(ps)*kvals)
  mse_avg3 = mse_avg3 / (N * length(ns)*length(ps)*kvals)
  mse_avg4 = mse_avg4 / (N * length(ns)*length(ps)*kvals)
}


saveRDS(mse_avg1, file=sprintf("%s/mse_avg1.RDA", savedir))
saveRDS(mse_avg2, file=sprintf("%s/mse_avg2.RDA", savedir))
saveRDS(mse_avg3, file=sprintf("%s/mse_avg3.RDA", savedir))
saveRDS(mse_avg4, file=sprintf("%s/mse_avg4.RDA", savedir))
saveRDS(ns, file=sprintf("%s/ns.RDA", savedir))
saveRDS(ps, file=sprintf("%s/ps.RDA", savedir))
saveRDS(const1, file=sprintf("%s/const1.RDA", savedir))
saveRDS(const2, file=sprintf("%s/const2.RDA", savedir))


# mse1 minimized at i = 11, j = 14
min(mse_avg1)
min(mse_avg2)
which(mse_avg1 == min(mse_avg1), arr.ind=TRUE)
const1[5]
const2[14]

which(mse_avg2 == min(mse_avg2), arr.ind=TRUE)
#const1[11]
#const2[14]






ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()




















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
  sbs_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  sbs_time = array(NA, dim= c(length(ns), length(ps), kvals))
  sbs_mse = array(NA, dim= c(length(ns), length(ps), kvals))
  subset_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  subset_time = array(NA, dim= c(length(ns), length(ps), kvals))
  subset_mse = array(NA, dim= c(length(ns), length(ps), kvals))  
  dc_res = array(NA, dim = c(length(ns), length(ps), kvals,N) )
  dc_time = array(NA, dim= c(length(ns), length(ps), kvals))
  dc_mse = array(NA, dim= c(length(ns), length(ps), kvals))
  
  for (i in 1:length(ns)) {
    n = ns[i]
    for(j in 1:length(ps)){
      #p = ps[j]
      for (y in 1:kvals) {
        subset = (result[2,] == i& result[3,] ==j & result[4,] == y )
        inspect_res[i,j,y,] = result[6, subset]
        inspect_time[i,j,y] = sum(result[7, subset])
        eta = round(n*0.2)
        inspect_mse[i,j,y] = mean((inspect_res[i,j,y,] - eta)^2)
        hdcd_res[i,j,y,] = result[8, subset]
        hdcd_time[i,j,y] = sum(result[9, subset])
        hdcd_s[i,j,y,] = result[10, subset]
        hdcd_mse[i,j,y] = mean((hdcd_res[i,j,y,] - eta)^2)
        
        scan_res[i,j,y,] = result[11, subset]
        scan_time[i,j,y] = sum(result[12, subset])
        scan_s[i,j,y,] = result[13, subset]
        scan_mse[i,j,y] = mean((scan_res[i,j,y,] - eta)^2)
        
        sbs_res[i,j,y,] = result[14, subset]
        sbs_time[i,j,y] = sum(result[15, subset])
        sbs_mse[i,j,y] = mean((sbs_res[i,j,y,] - eta)^2)
        
        subset_res[i,j,y,] = result[16, subset]
        subset_time[i,j,y] = sum(result[17, subset])
        subset_mse[i,j,y] = mean((subset_res[i,j,y,] - eta)^2)
        
        dc_res[i,j,y,] = result[18, subset]
        dc_time[i,j,y] = sum(result[19, subset])
        dc_mse[i,j,y] = mean((dc_res[i,j,y,] - eta)^2)
      }
    }
  }
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
                 "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|ccccc|ccccc}",
                 "\\hline", 
                 "\\multicolumn{5}{c|}{Parameters} & \\multicolumn{5}{c|}{MSE} &\\multicolumn{5}{c|}{Time in miliseconds} \\\\ \\hline ",
                 "$n$ & $p$ & $k$ &  $\\eta$ & $\\phi$ & \\text{FAST} & \\text{Inspect} & \\text{SBS} & \\text{SUBSET} & \\text{DC} & \\text{FAST} & \\text{Inspect} &\\text{SBS} & \\text{SUBSET}& \\text{DC} \\\\", 
                 "\\hline \\")
  
  for (i in 1:length(ns)) {
    
    n = ns[i]
    for(j in 1:length(ps)){
      p = ps[j]
      for (y in 1:kvals) {
        k = kfunc(y,n,p)
        eta = round(0.2*n)
        rootnorm = 0
        if(k<sqrt(p*log(n))){
          rootnorm = sparse_const/sqrt(eta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
        }else{
          rootnorm = dense_const/sqrt(eta)*(p*log(n))^(1/4)
        }
        string = sprintf("%d & %d & %d & %d & %.2f", n, p, k, eta, rootnorm)
        
        #res = round(c(hdcd_mse[i,j,y],inspect_mse[i,j,y], scan_mse[i,j,y]),digits=3)
        res = round(c(hdcd_mse[i,j,y],inspect_mse[i,j,y], sbs_mse[i,j,y], subset_mse[i,j,y], dc_mse[i,j,y]),digits=3)
        minind = (res==min(res))
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
          }else{
            string = sprintf("%s & %.3f", string, res[t])
          }
        }
        
        
        res = round(c(hdcd_time[i,j,y],inspect_time[i,j,y], sbs_time[i,j,y], subset_time[i,j,y], dc_time[i,j,y])/N*1000,digits=3)
        #res = round(c(hdcd_time[i,j,y],inspect_time[i,j,y], scan_time[i,j,y])/N*1000,digits=3)
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


library(logcondens)

rezz_inspect = logConDens(inspect_res[1,3,1,],smoothed=FALSE, print=FALSE)
rezz_hdcd = logConDens(hdcd_res[1,3,1,],smoothed=FALSE, print=FALSE)

xs = seq(1,n, by=1)
dens_hdcd = rep(NA, length(xs))
dens_inspect = rep(NA, length(xs))

for (i in 1:length(xs)) {
  xx = xs[i]
  dens_hdcd[i] = (evaluateLogConDens(xx, rezz_hdcd, which = 1:3)[, c("log-density", "density", "CDF")])[2]
  dens_inspect[i] = (evaluateLogConDens(xx, rezz_inspect, which = 1:3)[, c("log-density", "density", "CDF")])[2]
}

matplot(xs, cbind(dens_hdcd, dens_inspect), type="l",
        col = c(2,3), lty = c(1,1), lwd = c(2,2))

library(ggplot2)

par(mfrow = c(2,2))
for (i in 1:length(ns)) {
  for(j in 1:length(ps)){
    for (kk in 1:4) {
      n = ns[i]
      p = ps[i]
      print(sprintf("i = %d, n = %d, j = %d, p = %d, kk = %d, k = %d", 
                    i, ns[i], j, ps[j], kk, kfunc(kk, ns[i], ps[i])))
      rezz_inspect = logConDens(inspect_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      rezz_hdcd = logConDens(hdcd_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      rezz_sbs = NULL
      if(var(sbs_res[i,j,kk,])>0.0001){
        rezz_sbs = logConDens(sbs_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      }
      
      rezz_subset = logConDens(subset_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      rezz_scan = logConDens(scan_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      rezz_dc = logConDens(dc_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      
      
      xs = seq(1,n, by=1)
      dens_hdcd = rep(NA, length(xs))
      dens_inspect = rep(NA, length(xs))
      dens_sbs = rep(NA, length(xs))
      dens_subset = rep(NA, length(xs))
      dens_scan = rep(NA, length(xs))
      dens_dc = rep(NA, length(xs))
      
      for (ii in 1:length(xs)) {
        xx = xs[ii]
        dens_hdcd[ii] = (evaluateLogConDens(xx, rezz_hdcd, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        dens_inspect[ii] = (evaluateLogConDens(xx, rezz_inspect, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        if(var(sbs_res[i,j,kk,])>0.0001){
          dens_sbs[ii] = (evaluateLogConDens(xx, rezz_sbs, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        }else{dens_sbs[ii] = 1/length(xs)}
        dens_subset[ii] = (evaluateLogConDens(xx, rezz_subset, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        dens_scan[ii] = (evaluateLogConDens(xx, rezz_scan, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        dens_dc[ii] = (evaluateLogConDens(xx, rezz_dc, which = 1:3)[, c("log-density", "density", "CDF")])[2]
      }
      
      # hdcddf <- data.frame(
      #   xs,dens_hdcd
      # )
      # 
      # inspectdf <- data.frame(
      #   xs,dens_inspect
      # )
      # 
      # cols = c("x","y")
      # 
      # colnames(inspectdf) = cols
      # colnames(hdcddf) = cols
      # p = ggplot() +
      #   geom_line(data = hdcddf, aes(x = x, y = y), color = "blue") +
      #   geom_line(data = inspectdf, aes(x = x, y = y), color = "red") +
      #   xlab('estimated change-point location') +
      #   ylab('density') 
      
      dff = data.frame(x = xs, values = c(dens_hdcd, dens_inspect, dens_sbs, dens_subset,dens_dc),
                       fun = c(rep("FAST", length(xs)), rep("Inspect", length(xs)), 
                               rep("SBS", length(xs)),rep("SUBSET", length(xs)),rep("DC", length(xs)))
      )
      
      p = ggplot(dff,aes(x,values,col=fun))+geom_line(size=1.3, aes(linetype=fun)) +
        theme(legend.position="top",legend.title=element_blank(),text = element_text(size = 30)) +
        xlab('estimated change-point location') 
      #   ylab('density') 
      
      
      #print(p)
      ggsave(filename = sprintf("%s/i=%d_j=%d_kk=%d.pdf", plotdir,i, j, kk),
             plot = p,device = "pdf", dpi = 300)
      matplot(xs, cbind(dens_hdcd, dens_inspect, dens_sbs, dens_subset,dens_dc), type="l",
              col = c(2,3,4,5,6), lty =rep(1,5), lwd = rep(1,5))
      legend(x = "topright", legend = 
               c("FAST", "Inspect", "SBS", "SUBSET","DC"), lty=rep(1,5), col = c(2,3,4,5,6), 
             lwd = rep(2,6))
      
    }    
  }
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
        eta = round(n*0.2)
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
        eta = round(0.2*n)
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
