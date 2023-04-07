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
  savedir = file.path(maindir, sprintf("%s/single_rate",dateandtime))
  dir.create(savedir, showWarnings =FALSE )
  plotdir = file.path(maindir, sprintf("%s/single_rate/plots",dateandtime))
  dir.create(plotdir, showWarnings =FALSE )
  
}
#ns = c(200,500)
ns = c(200)
ps = c(1000)
#ps = c(100,1000,5000)
#ps = ns[]
#ps = c(500,1000,2000,5000)
kvals=4
totruns = length(ns)*length(ps)*kvals

simsbsN = 1000
simsbstoln = 1/simsbsN

#for sbs: 
set.seed(1996)
pis = matrix(nrow = length(ns), ncol =length(ps))
for (i in 1:length(ns)) {
  for (j in 1:length(ps)) {
    pis[i,j] =  single_SBS_calibrate(ns[i],ps[j],simsbsN,simsbstoln,rescale_variance = TRUE,debug=FALSE)
    
  }
}


N = 1000
num_cores = 6
sparse_const = 2.5
dense_const = 2.5
even_spread = TRUE
NN = 40
consts = seq(2.0, 5, length.out = NN)
#ns = c(500,1000,2000)
#ns = c(100,500,1000)
#ns = c(500,1000,2000)
#ns = c(100,500,1000)
#ps = c(1000,2000,10000)
#ns = c(200,500)

rescale_variance=TRUE


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
result = foreach(z = 1:N,.options.snow = opts) %dopar% {
  rez = list()
  counter = 1
  set.seed(z)
  library(HDCD)
  library(hdbinseg)
  source("/Users/peraugust/OneDrive - Universitetet i Oslo/project1/simulations/HDCD/SUBSET/SUBSET_normal.R")
  #for (i in 1:length(ns)) {
    n = ns[1]
    for(j in 1:NN){
      p = ps[1]
      for (y in 1:kvals) {
        rezi = list()
        # determine k!!
        k = kfunc(y,n,p)
        
        mus = matrix(0, nrow=p, ncol=n)
        noise = matrix(rnorm(n*p), nrow=p, ncol=n)
        eta = round(0.2*n)
        #eta = round(sample(1:floor(n/2),1))
        
        
        diff = 0
        if(k<sqrt(p*log(n))){
          rootnorm = consts[j]/sqrt(eta)*sqrt((c(k*log(exp(1)*p*log(n)/k^2)+ log(n))))
          diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
          
        }else{
          rootnorm = consts[j]/sqrt(eta)*(p*log(n))^(1/4)
          diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
        }
        if(!even_spread){
          diff = c(rnorm(k), rep(0,p-k))
          diff = diff/norm(diff, type="2")
          diff = diff*rootnorm
        }
        
        mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)
        X = mus+noise
        
        
        # DC
        a = proc.time()
        res_dc = dcbs.alg(X[,],phi = -1, height=1, thr = 0,cp.type=1,temporal=FALSE )
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
        res_sbs = sbs.alg(X[,],height=1, thr = rep(sqrt(pis[1,1]),p),cp.type=1,temporal=FALSE )
        #ress = single_SBS(X,pis[i,j] )
        #ress$pos
        #res_sbs$ecp
        b=proc.time()
        
        rezi["sbs_time"] = (b-a)[1]+(b-a)[2]
        if(!is.numeric(res_sbs$ecp)){
          #rezi["sbs_res"] = round(n/2)
          if(!rescale_variance){
            ress = single_SBS(X[,],pis[1,1] )
          }else{
            ress = single_SBS(rescale.variance(X),pis[1,1] )
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
        
        # then ESAC
        a = proc.time()
        res_ESAC  = single_ESAC(X[,], 1.5,1, debug =FALSE)
        b=proc.time()
        rezi["ESAC_time"] = (b-a)[1]+(b-a)[2]
        rezi["ESAC_res"] = res_ESAC$pos
        rezi["ESAC_s"] = res_ESAC$s
        
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
        res_subset  = SUBSET.normal(X[,])
        b=proc.time()
        rezi["subset_time"] = (b-a)[1]+(b-a)[2]
        
        if(is.null(res_subset$cpt)){
          rezi["subset_res"] = round(n/2)
        }else{
          rezi["subset_res"] = res_subset$cpt
        }
        
        rezi["z"] = z
        #rezi["i"] = i
        rezi["j"] = j
        rezi["y"] = y
        rezi["k"] = k
        rezi["eta"] = eta
        
        # rezi = list(z,i,j,y,k,inspect_res, inspect_time, ESAC_res, ESAC_time,ESAC_s, scan_res, scan_time, scan_s, sbs_res, sbs_time, 
        #                   subset_res, subset_time, dc_res, dc_time)
        
        rez[[counter]] = rezi
        counter = counter+1
        
      }
    }
  #}
  
  rez
}
# bb = proc.time()
# print(bb-aa)
#parallel::stopCluster(cl)
close(pb)
stopCluster(cl) 

{
  inspect_res = array(NA, dim = c( NN, kvals,N) )
  #inspect_time = array(0, dim= c( NN, kvals))
  inspect_mae = array(0, dim= c( NN, kvals))
  ESAC_res = array(NA, dim = c( NN, kvals,N) )
  #ESAC_s = array(NA, dim = c( NN, kvals,N) )
  #ESAC_time = array(0, dim= c( NN, kvals))
  ESAC_mae = array(0, dim= c( NN, kvals))
  #scan_res = array(NA, dim = c( NN, kvals,N) )
  #scan_s = array(NA, dim = c( NN, kvals,N) )
  #scan_time = array(0, dim= c( NN, kvals))
  #scan_mae = array(NA, dim= c( NN, kvals))  
  sbs_res = array(NA, dim = c( NN, kvals,N) )
  #sbs_time = array(0, dim= c( NN, kvals))
  sbs_mae = array(0, dim= c( NN, kvals))
  subset_res = array(NA, dim = c( NN, kvals,N) )
  #subset_time = array(0, dim= c( NN, kvals))
  subset_mae = array(0, dim= c( NN, kvals))  
  dc_res = array(NA, dim = c( NN, kvals,N) )
  #dc_time = array(0, dim= c( NN, kvals))
  dc_mae = array(0, dim= c( NN, kvals))
  
  for (z in 1:N) {
    list = result[[z]]
    len = length(list)
    
    for (t in 1:len) {
      sublist = list[[t]]
      y = sublist[["y"]]
      #i = sublist[["i"]]
      j = sublist[["j"]]
      eta = sublist[["eta"]]
      
      inspect_res[j,y,z] = sublist["inspect_res"][[1]]
      #inspect_time[j,y] = inspect_time[j,y] +  sublist["inspect_time"][[1]]/N
      inspect_mae[j,y] = inspect_mae[j,y] + abs((inspect_res[j,y,z] - eta)^2)/N
      
      ESAC_res[j,y,z] = sublist["ESAC_res"][[1]]
      #ESAC_time[j,y] = ESAC_time[j,y] + sublist["ESAC_time"][[1]]/N
      ESAC_mae[j,y] = ESAC_mae[j,y] + abs((ESAC_res[j,y,z] - eta)^2)/N
      
      sbs_res[j,y,z] = sublist["sbs_res"][[1]]
      #sbs_time[j,y] = sbs_time[j,y] +  sublist["sbs_time"][[1]]/N
      sbs_mae[j,y] = sbs_mae[j,y] + abs((sbs_res[j,y,z] - eta)^2)/N
      
      subset_res[j,y,z] = sublist["subset_res"][[1]]
      #subset_time[j,y] = subset_time[j,y] +  sublist["subset_time"][[1]]/N
      subset_mae[j,y] = subset_mae[j,y] + abs((subset_res[j,y,z] - eta)^2)/N
      
      dc_res[j,y,z] = sublist["dc_res"][[1]]
      #dc_time[j,y] = dc_time[j,y] + sublist["dc_time"][[1]]/N
      dc_mae[j,y] = dc_mae[j,y] + abs((dc_res[j,y,z] - eta)^2)/N
      
      
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
  #       inspect_mae[i,j,y] = mean((inspect_res[i,j,y,] - eta)^2)
  #       ESAC_res[i,j,y,] = result[8, subset]
  #       ESAC_time[i,j,y] = sum(result[9, subset])
  #       ESAC_s[i,j,y,] = result[10, subset]
  #       ESAC_mae[i,j,y] = mean((ESAC_res[i,j,y,] - eta)^2)
  #       
  #       scan_res[i,j,y,] = result[11, subset]
  #       scan_time[i,j,y] = sum(result[12, subset])
  #       scan_s[i,j,y,] = result[13, subset]
  #       scan_mae[i,j,y] = mean((scan_res[i,j,y,] - eta)^2)
  #       
  #       sbs_res[i,j,y,] = result[14, subset]
  #       sbs_time[i,j,y] = sum(result[15, subset])
  #       sbs_mae[i,j,y] = mean((sbs_res[i,j,y,] - eta)^2)
  #       
  #       subset_res[i,j,y,] = result[16, subset]
  #       subset_time[i,j,y] = sum(result[17, subset])
  #       subset_mae[i,j,y] = mean((subset_res[i,j,y,] - eta)^2)
  #       
  #       dc_res[i,j,y,] = result[18, subset]
  #       dc_time[i,j,y] = sum(result[19, subset])
  #       dc_mae[i,j,y] = mean((dc_res[i,j,y,] - eta)^2)
  #     }
  #   }
  # }
}

kk = 1
plot(consts^2, inspect_mae[1:NN, kk],type="l")
abline(h=0)
lines(consts^2, ESAC_mae[1:NN, kk], type="l",col=2)
lines(consts^2, subset_mae[1:NN, kk], type="l",col=3)
lines(consts^2, dc_mae[1:NN, kk], type="l",col=4)
lines(consts^2, sbs_mae[1:NN, kk], type="l",col=5)




inspect_mae - ESAC_mae
inspect_mae
ESAC_mae

par(mfrow=c(1,2))
hist(inspect_res[1,2,4,])
hist(ESAC_res[1,2,4,])

#saving: 
if(save){
  saveRDS(inspect_mae, file=sprintf("%s/inspect_mae.RDA", savedir))
  #saveRDS(inspect_time, file=sprintf("%s/inspect_time.RDA", savedir))
  saveRDS(inspect_res, file=sprintf("%s/inspect_res.RDA", savedir))
  saveRDS(ESAC_res, file=sprintf("%s/ESAC_res.RDA", savedir))
  saveRDS(ESAC_s, file=sprintf("%s/ESAC_s.RDA", savedir))
  #saveRDS(ESAC_time, file=sprintf("%s/ESAC_time.RDA", savedir))
  saveRDS(ESAC_mae, file=sprintf("%s/ESAC_mae.RDA", savedir))
  saveRDS(scan_res, file=sprintf("%s/scan_res.RDA", savedir))
  saveRDS(scan_s, file=sprintf("%s/scan_s.RDA", savedir))
  #saveRDS(scan_time, file=sprintf("%s/scan_time.RDA", savedir))
  saveRDS(scan_mae, file=sprintf("%s/scan_mae.RDA", savedir))
  saveRDS(sbs_mae, file=sprintf("%s/sbs_mae.RDA", savedir))
  saveRDS(sbs_res, file=sprintf("%s/sbs_res.RDA", savedir))
  #saveRDS(sbs_time, file=sprintf("%s/sbs_time.RDA", savedir))
  saveRDS(subset_mae, file=sprintf("%s/subset_mae.RDA", savedir))
  saveRDS(subset_res, file=sprintf("%s/subset_res.RDA", savedir))
  #saveRDS(subset_time, file=sprintf("%s/subset_time.RDA", savedir))
  saveRDS(dc_mae, file=sprintf("%s/dc_mae.RDA", savedir))
  saveRDS(dc_res, file=sprintf("%s/dc_res.RDA", savedir))
  #saveRDS(dc_time, file=sprintf("%s/dc_time.RDA", savedir))
  
  infofile<-file(sprintf("%s/parameters.txt", savedir))
  writeLines(c(sprintf("N = %d", N),
               sprintf("ns = %s", paste(ns, sep=" ", collapse=" ")),
               sprintf("ps = %s", paste(ps, sep=" ", collapse=" ")), 
               sprintf("Sparse constant = %f", sparse_const), 
               sprintf("Dense constant = %f", dense_const),
               sprintf("Even spread = %d", as.integer(even_spread))),
             infofile)
  close(infofile)
}




library(ggplot2)
library(latex2exp)
par(mfrow = c(2,2))
#for (i in 1:length(ns)) {
#  for(j in 1:length(ps)){
    for (kk in 1:4) {
      n = ns[1]
      p = ps[1]
      print(sprintf("i = %d, n = %d, j = %d, p = %d, kk = %d, k = %d", 
                    i, ns[i], j, ps[j], kk, kfunc(kk, ns[i], ps[i])))
      
      rezz_inspect = logConDens(inspect_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      rezz_ESAC = logConDens(ESAC_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      rezz_sbs = NULL
      sbs_res2 = sbs_res[i,j,kk,]
      sbs_res2[sbs_res2 == 1] = sample(1:(n-1), size=sum(sbs_res2==1),replace=TRUE)
      #if(var(sbs_res[i,j,kk,])>0.0001){
      rezz_sbs = logConDens(sbs_res2,smoothed=FALSE, print=FALSE)
      #}
      
      rezz_subset = logConDens(subset_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      #rezz_scan = logConDens(scan_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      rezz_dc = logConDens(dc_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      
      eta = round(0.2*n)
      xs = seq(eta-20,eta+20, by=1)
      dens_ESAC = rep(NA, length(xs))
      dens_inspect = rep(NA, length(xs))
      dens_sbs = rep(NA, length(xs))
      dens_subset = rep(NA, length(xs))
      dens_scan = rep(NA, length(xs))
      dens_dc = rep(NA, length(xs))
      
      for (ii in 1:length(xs)) {
        xx = xs[ii]
        dens_ESAC[ii] = (evaluateLogConDens(xx, rezz_ESAC, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        dens_inspect[ii] = (evaluateLogConDens(xx, rezz_inspect, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        if(var(sbs_res[i,j,kk,])>0.0001){
          dens_sbs[ii] = (evaluateLogConDens(xx, rezz_sbs, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        }else{dens_sbs[ii] = 1/length(xs)}
        dens_subset[ii] = (evaluateLogConDens(xx, rezz_subset, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        #dens_scan[ii] = (evaluateLogConDens(xx, rezz_scan, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        dens_dc[ii] = (evaluateLogConDens(xx, rezz_dc, which = 1:3)[, c("log-density", "density", "CDF")])[2]
      }
      
      # ESACdf <- data.frame(
      #   xs,dens_ESAC
      # )
      # 
      # inspectdf <- data.frame(
      #   xs,dens_inspect
      # )
      # 
      # cols = c("x","y")
      # 
      # colnames(inspectdf) = cols
      # colnames(ESACdf) = cols
      # p = ggplot() +
      #   geom_line(data = ESACdf, aes(x = x, y = y), color = "blue") +
      #   geom_line(data = inspectdf, aes(x = x, y = y), color = "red") +
      #   xlab('estimated change-point location') +
      #   ylab('density') 
      
      dff = data.frame(x = xs, values = c(dens_ESAC, dens_inspect, dens_sbs, dens_subset,dens_dc),
                       fun = factor(c(rep("FAST", length(xs)), rep("Inspect", length(xs)),
                                      rep("SBS", length(xs)),rep("SUBSET", length(xs)),rep("DC", length(xs))),levels=c("FAST", "Inspect", "SBS", "SUBSET","DC"))
      )
      # dff = data.frame(x = xs, values = c(dens_ESAC, dens_inspect, dens_sbs, dens_subset,dens_dc),
      #                 fun = as.factor(c(rep("FAST", length(xs)), rep("Inspect", length(xs)),
      #                         rep("SBS", length(xs)),rep("SUBSET", length(xs)),rep("DC", length(xs))))
      #                 )
      #levels(dff$fun) <- c("FAST", "Inspect", "SUBSET", "DC","SBS")
      
      kexp = c("", "=\\lceil p^{1/3}\\rceil", "=\\lceil (p\\log n)^{1/2} \\rceil","=p" )[kk]
      #p = ggplot(dff,aes(x,values,col=fun))+geom_line(size=1.3, aes(linetype=fun)) +
      p = ggplot(dff,aes(x,values,col=fun))+geom_line(size=1.3, linetype="solid") +
        theme(legend.position="top",legend.title=element_blank(),text = element_text(size = 40)) +
        xlab('estimated change-point location')  + scale_color_manual(values = c("red1", "steelblue","gold4","springgreen1","plum3")) +
        ylab('density')  + 
        theme(axis.text=element_text(size=20), axis.title=element_text(size=25),plot.title = element_text(size=40,hjust = 0.5)) +
        ggtitle(TeX(sprintf("Sparsity $k %s = %d$\n",kexp, kfunc(kk,n,p))))
      
      
      #print(p)
      ggsave(filename = sprintf("%s/i=%d_j=%d_kk=%d.pdf", plotdir,i, j, kk),
             plot = p,device = "pdf", dpi = 300,height=10,width=14)
      matplot(xs, cbind(dens_ESAC, dens_inspect, dens_sbs, dens_subset,dens_dc), type="l",
              col = c(2,3,4,5,6), lty =rep(1,5), lwd = rep(1,5))
      legend(x = "topright", legend = 
               c("FAST", "Inspect", "SBS", "SUBSET","DC"), lty=rep(1,5), col = c(2,3,4,5,6), 
             lwd = rep(2,6))
      
    }    
#  }
#}























# creating table:
if(save){
  # output latex table
  printlines = c("\\begin{table}[H] \\centering",
                 "\\caption{Single change-point estimation}",
                 "\\label{tablesingleloc}",
                 "\\small",
                 "\\begin{adjustbox}{scale=0.9,center}",
                 "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|ccccc|ccccc}",
                 "\\hline", 
                 "\\multicolumn{5}{c|}{Parameters} & \\multicolumn{5}{c|}{mae} &\\multicolumn{5}{c|}{Time in milliseconds} \\\\ \\hline ",
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
        
        #res = round(c(ESAC_mae[i,j,y],inspect_mae[i,j,y], scan_mae[i,j,y]),digits=3)
        res = round(c(ESAC_mae[i,j,y],inspect_mae[i,j,y], sbs_mae[i,j,y], subset_mae[i,j,y], dc_mae[i,j,y]),digits=1)
        minind = (res==min(res))
        
        for (t in 1:length(res)) {
          if(minind[t]){
            string = sprintf("%s & \\textbf{%.1f} ", string, res[t])
          }else{
            string = sprintf("%s & %.1f", string, res[t])
          }
        }
        
        
        res = round(c(ESAC_time[i,j,y],inspect_time[i,j,y], sbs_time[i,j,y], subset_time[i,j,y], dc_time[i,j,y])*1000,digits=1)
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
  printlines = c(printlines, "\\hline \\multicolumn{5}{c|}{Average mae}")
  #res = round(c(mean(ESAC_slow_det[,,2:kvals]),mean(ESAC_det[,,2:kvals]),mean(pilliat_det[,,2:kvals]),mean(inspect_det[,,2:kvals]), mean(sbs_det[,,2:kvals]), mean(subset_det[,,2:kvals]), mean(dc_det[,,2:kvals])),digits=3)
  res = round(c(mean(ESAC_mae),mean(inspect_mae), mean(sbs_mae), mean(subset_mae), mean(dc_mae)),digits=3)
  minind = (res==min(res))
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


library(logcondens)

rezz_inspect = logConDens(inspect_res[1,3,1,],smoothed=FALSE, print=FALSE)
rezz_ESAC = logConDens(ESAC_res[1,3,1,],smoothed=FALSE, print=FALSE)

xs = seq(1,n, by=1)
dens_ESAC = rep(NA, length(xs))
dens_inspect = rep(NA, length(xs))

for (i in 1:length(xs)) {
  xx = xs[i]
  dens_ESAC[i] = (evaluateLogConDens(xx, rezz_ESAC, which = 1:3)[, c("log-density", "density", "CDF")])[2]
  dens_inspect[i] = (evaluateLogConDens(xx, rezz_inspect, which = 1:3)[, c("log-density", "density", "CDF")])[2]
}

matplot(xs, cbind(dens_ESAC, dens_inspect), type="l",
        col = c(2,3), lty = c(1,1), lwd = c(2,2))

library(ggplot2)
library(latex2exp)
par(mfrow = c(2,2))
for (i in 1:length(ns)) {
  for(j in 1:length(ps)){
    for (kk in 1:4) {
      n = ns[i]
      p = ps[i]
      print(sprintf("i = %d, n = %d, j = %d, p = %d, kk = %d, k = %d", 
                    i, ns[i], j, ps[j], kk, kfunc(kk, ns[i], ps[i])))
      
      rezz_inspect = logConDens(inspect_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      rezz_ESAC = logConDens(ESAC_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      rezz_sbs = NULL
      sbs_res2 = sbs_res[i,j,kk,]
      sbs_res2[sbs_res2 == 1] = sample(1:(n-1), size=sum(sbs_res2==1),replace=TRUE)
      #if(var(sbs_res[i,j,kk,])>0.0001){
      rezz_sbs = logConDens(sbs_res2,smoothed=FALSE, print=FALSE)
      #}
      
      rezz_subset = logConDens(subset_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      #rezz_scan = logConDens(scan_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      rezz_dc = logConDens(dc_res[i,j,kk,],smoothed=FALSE, print=FALSE)
      
      eta = round(0.2*n)
      xs = seq(eta-20,eta+20, by=1)
      dens_ESAC = rep(NA, length(xs))
      dens_inspect = rep(NA, length(xs))
      dens_sbs = rep(NA, length(xs))
      dens_subset = rep(NA, length(xs))
      dens_scan = rep(NA, length(xs))
      dens_dc = rep(NA, length(xs))
      
      for (ii in 1:length(xs)) {
        xx = xs[ii]
        dens_ESAC[ii] = (evaluateLogConDens(xx, rezz_ESAC, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        dens_inspect[ii] = (evaluateLogConDens(xx, rezz_inspect, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        if(var(sbs_res[i,j,kk,])>0.0001){
          dens_sbs[ii] = (evaluateLogConDens(xx, rezz_sbs, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        }else{dens_sbs[ii] = 1/length(xs)}
        dens_subset[ii] = (evaluateLogConDens(xx, rezz_subset, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        #dens_scan[ii] = (evaluateLogConDens(xx, rezz_scan, which = 1:3)[, c("log-density", "density", "CDF")])[2]
        dens_dc[ii] = (evaluateLogConDens(xx, rezz_dc, which = 1:3)[, c("log-density", "density", "CDF")])[2]
      }
      
      # ESACdf <- data.frame(
      #   xs,dens_ESAC
      # )
      # 
      # inspectdf <- data.frame(
      #   xs,dens_inspect
      # )
      # 
      # cols = c("x","y")
      # 
      # colnames(inspectdf) = cols
      # colnames(ESACdf) = cols
      # p = ggplot() +
      #   geom_line(data = ESACdf, aes(x = x, y = y), color = "blue") +
      #   geom_line(data = inspectdf, aes(x = x, y = y), color = "red") +
      #   xlab('estimated change-point location') +
      #   ylab('density') 
      
      dff = data.frame(x = xs, values = c(dens_ESAC, dens_inspect, dens_sbs, dens_subset,dens_dc),
                       fun = factor(c(rep("FAST", length(xs)), rep("Inspect", length(xs)),
                                      rep("SBS", length(xs)),rep("SUBSET", length(xs)),rep("DC", length(xs))),levels=c("FAST", "Inspect", "SBS", "SUBSET","DC"))
      )
      # dff = data.frame(x = xs, values = c(dens_ESAC, dens_inspect, dens_sbs, dens_subset,dens_dc),
      #                 fun = as.factor(c(rep("FAST", length(xs)), rep("Inspect", length(xs)),
      #                         rep("SBS", length(xs)),rep("SUBSET", length(xs)),rep("DC", length(xs))))
      #                 )
      #levels(dff$fun) <- c("FAST", "Inspect", "SUBSET", "DC","SBS")
      
      kexp = c("", "=\\lceil p^{1/3}\\rceil", "=\\lceil (p\\log n)^{1/2} \\rceil","=p" )[kk]
      #p = ggplot(dff,aes(x,values,col=fun))+geom_line(size=1.3, aes(linetype=fun)) +
      p = ggplot(dff,aes(x,values,col=fun))+geom_line(size=1.3, linetype="solid") +
        theme(legend.position="top",legend.title=element_blank(),text = element_text(size = 40)) +
        xlab('estimated change-point location')  + scale_color_manual(values = c("red1", "steelblue","gold4","springgreen1","plum3")) +
        ylab('density')  + 
        theme(axis.text=element_text(size=20), axis.title=element_text(size=25),plot.title = element_text(size=40,hjust = 0.5)) +
        ggtitle(TeX(sprintf("Sparsity $k %s = %d$\n",kexp, kfunc(kk,n,p))))
      
      
      #print(p)
      ggsave(filename = sprintf("%s/i=%d_j=%d_kk=%d.pdf", plotdir,i, j, kk),
             plot = p,device = "pdf", dpi = 300,height=10,width=14)
      matplot(xs, cbind(dens_ESAC, dens_inspect, dens_sbs, dens_subset,dens_dc), type="l",
              col = c(2,3,4,5,6), lty =rep(1,5), lwd = rep(1,5))
      legend(x = "topright", legend = 
               c("FAST", "Inspect", "SBS", "SUBSET","DC"), lty=rep(1,5), col = c(2,3,4,5,6), 
             lwd = rep(2,6))
      
    }    
  }
}



# 
# 
# # THEN UEVEN CHANGES
# 
# 
# 
# cl <- makeCluster(num_cores,type="SOCK")
# registerDoSNOW(cl)
# pb <- txtProgressBar(max = N, style = 3)
# progress <- function(n) setTxtProgressBar(pb, n)
# opts <- list(progress = progress)
# 
# result_uneven = foreach(z = 1:N, .combine="cbind",.options.snow = opts) %dopar% {
#   set.seed(z)
#   rez = c()
#   library(HDCD)
#   for (i in 1:length(ns)) {
#     n = ns[i]
#     for(j in 1:length(ps)){
#       p = ps[j]
#       for (y in 1:kvals) {
#         
#         # determine k!!
#         k = kfunc(y,n,p)
#         
#         mus = matrix(0, nrow=p, ncol=n)
#         noise = matrix(rnorm(n*p), nrow=p, ncol=n)
#         eta = round(0.4*n)
#         
#         
#         diff = c(rep(1, k)/(sqrt(1:k)), rep(0, p-k))
#         diff = diff / norm(diff, type="2")
#         if(k<sqrt(p*log(n^4))){
#           rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
#           diff = rootnorm * diff
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
#         res = single_Inspect(X[,], lambda)
#         b=proc.time()
#         inspect_time = (b-a)[1]+(b-a)[2]
#         inspect_res= res$pos
#         
#         # then ESAC
#         a = proc.time()
#         res  = single_ESAC (X[,], 2,1, debug =FALSE)
#         b=proc.time()
#         ESAC_time = (b-a)[1]+(b-a)[2]
#         ESAC_res = res$pos
#         ESAC_s = res$s
#         
#         # then scan
#         a = proc.time()
#         res  = single_Scan (X[,], debug =FALSE)
#         b=proc.time()
#         scan_time = (b-a)[1]+(b-a)[2]
#         scan_res = res$pos
#         scan_s = res$s
#         
# 
#         rez = cbind(rez, c(z,i,j,y,k,inspect_res, inspect_time, ESAC_res, ESAC_time,ESAC_s, scan_res, scan_time, scan_s))        
#         
#       }
#     }
#   }
#   
#   rez
# }
# # bb = proc.time()
# # print(bb-aa)
# close(pb)
# stopCluster(cl) 
# 
# 
# {
#   inspect_res_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
#   inspect_time_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
#   inspect_mae_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
#   ESAC_res_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
#   ESAC_s_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
#   ESAC_time_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
#   ESAC_mae_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
#   scan_res_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
#   scan_s_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
#   scan_time_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
#   scan_mae_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
#   
#   for (i in 1:length(ns)) {
#     n = ns[i]
#     for(j in 1:length(ps)){
#       #p = ps[j]
#       for (y in 1:kvals) {
#         subset = (result_uneven[2,] == i& result_uneven[3,] ==j & result_uneven[4,] == y )
#         inspect_res_uneven[i,j,y,] = result_uneven[6, subset]
#         inspect_time_uneven[i,j,y] = sum(result_uneven[7, subset])
#         eta = round(n*0.2)
#         inspect_mae_uneven[i,j,y] = mean((inspect_res_uneven[i,j,y,] - eta)^2)
#         ESAC_res_uneven[i,j,y,] = result_uneven[8, subset]
#         ESAC_time_uneven[i,j,y] = sum(result_uneven[9, subset])
#         ESAC_s_uneven[i,j,y,] = result_uneven[10, subset]
#         ESAC_mae_uneven[i,j,y] = mean((ESAC_res_uneven[i,j,y,] - eta)^2)
#         
#         scan_res_uneven[i,j,y,] = result_uneven[11, subset]
#         scan_time_uneven[i,j,y] = sum(result_uneven[12, subset])
#         scan_s_uneven[i,j,y,] = result_uneven[13, subset]
#         scan_mae_uneven[i,j,y] = mean((scan_res_uneven[i,j,y,] - eta)^2)
#       }
#     }
#   }
# }
# 
# 
# 
# 
# 
# inspect_mae_uneven - ESAC_mae_uneven
# inspect_mae_uneven
# ESAC_mae_uneven
# 
# # par(mfrow=c(1,2))
# # hist(inspect_res_uneven[3,3,1,])
# # hist(ESAC_res_uneven[3,3,1,])
# 
# #saving: 
# if(save){
#   saveRDS(inspect_mae_uneven, file=sprintf("%s/inspect_mae_uneven.RDA", savedir))
#   saveRDS(inspect_time_uneven, file=sprintf("%s/inspect_time_uneven.RDA", savedir))
#   saveRDS(inspect_mae_uneven, file=sprintf("%s/inspect_mae_uneven.RDA", savedir))
#   saveRDS(ESAC_res_uneven, file=sprintf("%s/ESAC_res_uneven.RDA", savedir))
#   saveRDS(ESAC_s_uneven, file=sprintf("%s/ESAC_s_uneven.RDA", savedir))
#   saveRDS(ESAC_time_uneven, file=sprintf("%s/ESAC_time_uneven.RDA", savedir))
#   saveRDS(ESAC_mae_uneven, file=sprintf("%s/ESAC_mae_uneven.RDA", savedir))
#   saveRDS(scan_res_uneven, file=sprintf("%s/scan_res_uneven.RDA", savedir))
#   saveRDS(scan_s_uneven, file=sprintf("%s/scan_s_uneven.RDA", savedir))
#   saveRDS(scan_time_uneven, file=sprintf("%s/scan_time_uneven.RDA", savedir))
#   saveRDS(scan_mae_uneven, file=sprintf("%s/scan_mae_uneven.RDA", savedir))
#   
# 
# }
# 
# # creating table:
# if(save){
#   # output latex table
#   printlines = c("\\begin{table}[!htbp] \\centering",
#                  "\\caption{Change spread unevenly across coordinates}",
#                  "\\label{}",
#                  "\\small",
#                  "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|ccc|ccc}",
#                  "\\hline", 
#                  "\\multicolumn{5}{c|}{Parameters} & \\multicolumn{3}{c|}{mae} &\\multicolumn{3}{c}{Time in miliseconds} \\\\ \\hline ",
#                  "$n$ & $p$ & $k$ &  $\\eta$ & $\\phi$ & \\text{ESAC} & \\text{Inspect} & \\text{Scan} & \\text{ESAC} & \\text{Inspect} & \\text{Scan} \\\\", 
#                  "\\hline \\")
#   
#   for (i in 1:length(ns)) {
#     
#     n = ns[i]
#     for(j in 1:length(ps)){
#       p = ps[j]
#       for (y in 1:kvals) {
#         k = kfunc(y,n,p)
#         eta = round(0.2*n)
#         rootnorm = 0
#         if(k<sqrt(p*log(n^4))){
#           rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
#         }else{
#           rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
#         }
#         string = sprintf("%d & %d & %d & %d & %.2f", n, p, k, eta, rootnorm)
#         
#         res = round(c(ESAC_mae_uneven[i,j,y],inspect_mae_uneven[i,j,y], scan_mae_uneven[i,j,y]),digits=3)
#         minind = (res==min(res))
#         
#         for (t in 1:length(res)) {
#           if(minind[t]){
#             string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
#           }else{
#             string = sprintf("%s & %.3f", string, res[t])
#           }
#         }
# 
#         
#         res = round(c(ESAC_time_uneven[i,j,y],inspect_time_uneven[i,j,y], scan_time_uneven[i,j,y])/N*1000,digits=3)
#         minind = (res==min(res))
#         
#         for (t in 1:length(res)) {
#           if(minind[t]){
#             string = sprintf("%s & \\textbf{%.3f} ", string, res[t])
#           }else{
#             string = sprintf("%s & %.3f", string, res[t])
#           }
#         }
#         
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
# 
# 
# 
# 
# 
# # 
# # 
# # 
# # # decreasing magnitude
# # # 
# # # N = 100
# # # 
# # # set.seed(100)
# # # 
# # # #ns = c(500,1000,2000)
# # # ns = c(50,100,500)
# # # #ps = c(1000,2000,10000)
# # # ps = ns[]
# # # kvals=4
# # # totruns = length(ns)*length(ps)*kvals
# # 
# # {
# #   inspect_res_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
# #   inspect_time_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
# #   inspect_mae_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
# #   ESAC_res_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
# #   ESAC_s_uneven = array(NA, dim = c(length(ns), length(ps), kvals,N) )
# #   ESAC_time_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
# #   ESAC_mae_uneven = array(NA, dim= c(length(ns), length(ps), kvals))
# # }
# # for (i in 1:length(ns)) {
# #   n = ns[i]
# #   for(j in 1:length(ps)){
# #     p = ps[j]
# #     for (y in 1:kvals) {
# #       
# #       # determine k!!
# #       k= 1
# #       if(y==1){
# #         k = 2
# #       }else if(y==2){
# #         k = ceiling(p^(1/3))
# #       }else if(y==3){
# #         k = ceiling(sqrt(p*log(n^4)))
# #       }else{
# #         k = p
# #       }
# #       
# #       
# #       
# #       tottimeinspect = 0
# #       tottimeESAC = 0
# #       
# #       for (z in 1:N) {
# #         mus = matrix(0, nrow=p, ncol=n)
# #         noise = matrix(rnorm(n*p), nrow=p, ncol=n)
# #         eta = round(0.4*n)
# #         
# #         
# #         diff = c(rep(1, k)/(sqrt(1:k)), rep(0, p-k))
# #         diff = diff / norm(diff, type="2")
# #         
# #         if(k<sqrt(p*log(n^4))){
# #           rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
# #           diff = rootnorm *diff
# #         }else{
# #           rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
# #           diff = rootnorm * diff
# #         }
# #         
# #         mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)
# #         
# #         X = mus+noise
# #         
# #         # first inspect
# #         lambda = sqrt(log(p*log(n))/2) 
# #         
# #         a = proc.time()
# #         res = single_Inspect(X, lambda)
# #         b=proc.time()
# #         tottimeinspect = tottimeinspect + (b-a)[1]+(b-a)[2]
# #         inspect_res_uneven[i,j,y,z] = res$pos
# #         
# #         # then ESAC
# #         a = proc.time()
# #         res  = single_ESAC (X, 2,1, debug =FALSE)
# #         b=proc.time()
# #         tottimeESAC = tottimeESAC + (b-a)[1]+(b-a)[2]
# #         ESAC_res_uneven[i,j,y,z] = res$pos
# #         ESAC_s_uneven[i,j,y,z] = res$s
# #         
# #       }
# #       inspect_time_uneven[i,j,y] = tottimeinspect
# #       inspect_mae_uneven[i,j,y] = mean((inspect_res_uneven[i,j,y,]-eta)^2)
# #       
# #       
# #       ESAC_time_uneven[i,j,y] = tottimeESAC
# #       ESAC_mae_uneven[i,j,y] = mean((ESAC_res_uneven[i,j,y,]-eta)^2)
# #       
# #       
# #       
# #     }
# #   }
# # }
# # 
# # #saving
# # if(save){
# #   saveRDS(inspect_mae_uneven, file=sprintf("%s/inspect_mae_uneven.RDA", savedir))
# #   saveRDS(inspect_time_uneven, file=sprintf("%s/inspect_time_uneven.RDA", savedir))
# #   saveRDS(inspect_mae_uneven, file=sprintf("%s/inspect_mae_uneven.RDA", savedir))
# #   saveRDS(ESAC_res_uneven, file=sprintf("%s/ESAC_res_uneven.RDA", savedir))
# #   saveRDS(ESAC_s_uneven, file=sprintf("%s/ESAC_s_uneven.RDA", savedir))
# #   saveRDS(ESAC_time_uneven, file=sprintf("%s/ESAC_time_uneven.RDA", savedir))
# #   saveRDS(ESAC_mae_uneven, file=sprintf("%s/ESAC_mae_uneven.RDA", savedir))
# #   
# # }
# # # saving table:
# # if(save){
# #   # output latex table
# #   printlines = c("\\begin{table}[!htbp] \\centering",
# #                  "\\caption{Change spread unevenly across coordinates}",
# #                  "\\label{}",
# #                  "\\small",
# #                  "\\begin{tabular}{@{\\extracolsep{1pt}} ccccc|cc|cc}",
# #                  "\\\\[-1.8ex]\\hline", 
# #                  "\\hline \\\\[-1.8ex]",
# #                  "$n$ & $p$ & $k$ &  $\\eta$ & $\\phi$ & \\text{ESAC mae} & \\text{Inspect mae} & \\text{ESAC (ms)} & \\text{Inspect (ms)} \\\\", 
# #                  "\\hline \\\\[-1.8ex]")
# #   
# #   for (i in 1:length(ns)) {
# #     
# #     n = ns[i]
# #     for(j in 1:length(ps)){
# #       p = ps[j]
# #       for (y in 1:kvals) {
# #         k = 1
# #         if(y==1){
# #           k = 2
# #         }else if(y==2){
# #           k = ceiling(p^(1/3))
# #         }else if(y==3){
# #           k = ceiling(sqrt(p*log(n^4)))
# #         }else{
# #           k = p
# #         }
# #         eta = round(0.4*n)
# #         rootnorm = 0
# #         if(k<sqrt(p*log(n^4))){
# #           rootnorm = sparse_const/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
# #         }else{
# #           rootnorm = dense_const/sqrt(eta)*(p*log(n^4))^(1/4)
# #         }
# #         string = sprintf("%d & %d & %d & %d & %.2f", n, p, k, eta, rootnorm)
# #         
# #         res = c(ESAC_mae_uneven[i,j,y],inspect_mae_uneven[i,j,y])
# #         minind = which.min(res)
# #         if(minind>1){
# #           for (t in 1:(minind-1)) {
# #             string = sprintf("%s & %.4f", string, res[t])
# #           }
# #         }
# #         string = sprintf("%s & \\textbf{%.4f} ", string, res[minind])
# #         if(minind<length(res)){
# #           for (t in (minind+1):length(res)) {
# #             string = sprintf("%s & %.4f", string, res[t])
# #           }
# #         }
# #         
# #         res = c(ESAC_time_uneven[i,j,y],inspect_time_uneven[i,j,y])/N*1000
# #         minind = which.min(res)
# #         if(minind>1){
# #           for (t in 1:(minind-1)) {
# #             string = sprintf("%s & %.2f", string, res[t])
# #           }
# #         }
# #         string = sprintf("%s & \\textbf{%.2f} ", string, res[minind])
# #         if(minind<length(res)){
# #           for (t in (minind+1):length(res)) {
# #             string = sprintf("%s & %.2f", string, res[t])
# #           }
# #         }
# #         string = sprintf("%s \\\\", string)
# #         printlines = c(printlines, string)
# #       }
# #     }
# #   }
# #   
# #   printlines = c(printlines, c("\\hline \\\\[-1.8ex]",
# #                                "\\end{tabular}",
# #                                "\\end{table}"))
# #   texfile<-file(sprintf("%s/table_uneven.tex", savedir))
# #   writeLines(printlines, texfile)
# #   close(texfile)
# #   
# # }
# # 
# # 
# # 
# # inspect_mae_uneven - ESAC_mae_uneven
# # 
# # inspect_time_uneven / ESAC_time_uneven
# # 
# # 
# # 
# # 
# # 
# # 
# # # p = 1000
# # # n =500
# # # mus = matrix(0, nrow=p, ncol=n)
# # # noise = matrix(rnorm(n*p), nrow=p, ncol=n)
# # # eta = round(0.4*n)
# # # k = p
# # # 
# # # diff = 0
# # # if(k<sqrt(p*log(n^4))){
# # #   rootnorm = 2/sqrt(eta)*sqrt(max(c(k*log(exp(1)*p*log(n^4)/k^2), log(n^4))))
# # #   diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
# # # }else{
# # #   rootnorm = 2/sqrt(eta)*(p*log(n^4))^(1/4)
# # #   diff = rootnorm * c(sample(c(-1,1), replace=TRUE,k), rep(0,p-k))/sqrt(k)
# # # }
# # # 
# # # mus[,round(eta+1):n] = mus[,round(eta+1):n] + matrix(rep(diff, n-round(eta)) ,nrow=p)
# # # 
# # # plot(mus[1,])
# # # X = mus+noise
# # # plot(X[1,])
# # # 
# # # xi = 4*sqrt(log(p*n))
# # # lambda = sqrt(log(p*log(n))/2) # double check this is what they recc!!!!!!!
# # # 
# # # res = single_Inspect(X, lambda)
# # # res$pos
# # # 
# # # res2 = single_ESAC (X, 2,2, debug =FALSE)
# # # res2$pos
# # # res2$s
# # # 
# # 
