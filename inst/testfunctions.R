#ESAC

library(HDCD)
n = 50
p = 50
set.seed(100)
# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# Adding a single sparse change-point:
X[1:5, 26:n] = X[1:5, 26:n] +2

# Vanilla ESAC:
res = ESAC(X)
res$changepoints

# Manually setting leading constants for \lambda(t) and \gamma(t)
res = ESAC(X, 
           threshold_d = 2, threshold_s = 2, #leading constants for \lambda(t)
           threshold_d_test = 2, threshold_s_test = 2 #leading constants for \gamma(t)
           )
res$changepoints #estimated change-point locations

# Empirical choice of thresholds:
res = ESAC(X, empirical = TRUE, N = 100, tol = 1/100)
res$changepoints

# Manual empirical choice of thresholds (equivalent to the above)
thresholds_emp = ESAC_calibrate(n,p, N=100, tol=1/100)
res = ESAC(X, thresholds_test = thresholds_emp[[1]])
res$changepoints







# ESAC_calibrate
library(HDCD)
n = 50
p = 50

set.seed(100)
thresholds_emp = ESAC_calibrate(n,p, N=100, tol=1/100)
set.seed(100)
thresholds_emp_without_bonferroni = ESAC_calibrate(n,p, N=100, tol=1/100,bonferroni=FALSE)
thresholds_emp[[1]] # vector of \gamma(t) for t = p,...,1
thresholds_emp_without_bonferroni[[1]] # vector of \gamma(t) for t = p,...,1

# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# Adding a single sparse change-point:
X[1:5, 26:n] = X[1:5, 26:n] +2

res = ESAC(X, thresholds_test = thresholds_emp[[1]])
res$changepoints

#ESAC single change-point test
library(HDCD)
n = 50
p = 50

# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
Y = matrix(rnorm(n*p), ncol = n, nrow=p)

# Adding a single sparse change-point to X (and not Y):
X[1:5, 26:n] = X[1:5, 26:n] +2

# Vanilla ESAC:
resX = ESAC_test(X)
resX
resY = ESAC_test(Y)
resY

# Manually setting leading constants for \lambda(t) and \gamma(t)
resX = ESAC_test(X, 
           threshold_d = 2, threshold_s = 2, #leading constants for \gamma(t)
           )
resX 
resY = ESAC_test(Y, 
                 threshold_d = 2, threshold_s = 2, #leading constants for \gamma(t)
)
resY

# Empirical choice of thresholds:
resX = ESAC_test(X, empirical = TRUE, N = 100, tol = 1/100)
resX
resY = ESAC_test(Y, empirical = TRUE, N = 100, tol = 1/100)
resY

# Manual empirical choice of thresholds (equivalent to the above)
thresholds_test_emp = ESAC_test_calibrate(n,p, N=100, tol=1/100,bonferroni=TRUE)
resX = ESAC_test(X, thresholds = thresholds_test_emp[[1]])
resX
resY = ESAC_test(Y, thresholds = thresholds_test_emp[[1]])
resY


# ESAC_test_calibrate
library(HDCD)
n = 50
p = 50

set.seed(100)
thresholds_emp = ESAC_test_calibrate(n,p, bonferroni=TRUE,N=100, tol=1/100)
set.seed(100)
thresholds_emp_without_bonferroni = ESAC_test_calibrate(n,p, bonferroni=FALSE,N=100, tol=1/100)
thresholds_emp[[1]] # vector of \gamma(t) for t = p,...,1
thresholds_emp_without_bonferroni[[1]] # vector of \gamma(t) for t = p,...,1

# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
Y = matrix(rnorm(n*p), ncol = n, nrow=p)

# Adding a single sparse change-point to X (and not Y):
X[1:5, 26:n] = X[1:5, 26:n] +2
resX = ESAC_test(X, thresholds = thresholds_emp[[1]])
resX
resY = ESAC_test(Y,  thresholds = thresholds_emp[[1]])
resY











# Pilliat

library(HDCD)
n = 50
p = 50
set.seed(100)
# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# Adding a single sparse change-point:
X[1:5, 26:n] = X[1:5, 26:n] +2

# Vanilla Pilliat:
res = Pilliat(X)
res$changepoints

# Manually setting leading constants for detection thresholds
res = Pilliat(X, threshold_d_const = 4, threshold_bj_const = 6, threshold_partial_const=4)
res$changepoints #estimated change-point locations

# Empirical choice of thresholds:
res = Pilliat(X, empirical = TRUE, N = 100, tol = 1/100)
res$changepoints

# Manual empirical choice of thresholds (equivalent to the above)
thresholds_emp = Pilliat_calibrate(n,p, N=100, tol=1/100)
res = Pilliat(X, threshold_dense =thresholds_emp$threshold_dense, 
              thresholds_bj = thresholds_emp$thresholds_bj,
              thresholds_partial =thresholds_emp$thresholds_partial )
res$changepoints



# pilliat calibrate
library(HDCD)
n = 50
p = 50

set.seed(100)
thresholds_emp = Pilliat_calibrate(n,p, N=100, tol=1/100)
thresholds_emp$thresholds_partial # thresholds for partial sum statistic
thresholds_emp$thresholds_bj # thresholds for Berk-Jones statistic
thresholds_emp$threshold_dense # thresholds for Berk-Jones statistic
set.seed(100)
thresholds_emp_without_bonferroni = Pilliat_calibrate(n,p, N=100, tol=1/100,bonferroni = FALSE)
thresholds_emp_without_bonferroni$thresholds_partial # thresholds for partial sum statistic
thresholds_emp_without_bonferroni$thresholds_bj # thresholds for Berk-Jones statistic
thresholds_emp_without_bonferroni$threshold_dense # thresholds for Berk-Jones statistic

# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# Adding a single sparse change-point:
X[1:5, 26:n] = X[1:5, 26:n] +2

res = Pilliat(X, threshold_dense =thresholds_emp$threshold_dense, 
              thresholds_bj = thresholds_emp$thresholds_bj,
              thresholds_partial =thresholds_emp$thresholds_partial )
res$changepoints






#Pilliat single change-point test
library(HDCD)
n = 50
p = 50

# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
Y = matrix(rnorm(n*p), ncol = n, nrow=p)

# Adding a single sparse change-point to X (and not Y):
X[1:5, 26:n] = X[1:5, 26:n] +2

# Vanilla Pilliat test:
resX = Pilliat_test(X)
resX
resY = Pilliat_test(Y)
resY

# Manually setting leading constants for the theoretical thresholds for the three test statistics used
resX = Pilliat_test(X, 
                 threshold_d_const=4, 
                 threshold_bj_const=6, 
                 threshold_partial_const=4
)
resX 
resY = Pilliat_test(Y, 
                    threshold_d_const=4, 
                    threshold_bj_const=6, 
                    threshold_partial_const=4
)
resY

# Empirical choice of thresholds:
resX = Pilliat_test(X, empirical = TRUE, N = 100, tol = 1/100)
resX
resY = Pilliat_test(Y, empirical = TRUE, N = 100, tol = 1/100)
resY

# Manual empirical choice of thresholds (equivalent to the above)
thresholds_test_emp = Pilliat_test_calibrate(n,p, N=100, tol=1/100,bonferroni=TRUE)
resX = Pilliat_test(X, 
                    threshold_dense=thresholds_test_emp$threshold_dense, 
                    thresholds_bj = thresholds_test_emp$thresholds_bj, 
                    thresholds_partial = thresholds_test_emp$thresholds_partial
)
resX
resY = Pilliat_test(Y, 
                    threshold_dense=thresholds_test_emp$threshold_dense, 
                    thresholds_bj = thresholds_test_emp$thresholds_bj, 
                    thresholds_partial = thresholds_test_emp$thresholds_partial
)
resY


# ESAC_test_calibrate
library(HDCD)
n = 50
p = 50

set.seed(100)
thresholds_test_emp = Pilliat_test_calibrate(n,p, bonferroni=TRUE,N=100, tol=1/100)
set.seed(100)
thresholds_test_emp_without_bonferroni = Pilliat_test_calibrate(n,p, bonferroni=FALSE,N=100, tol=1/100)
thresholds_test_emp # thresholds with bonferroni correction
thresholds_test_emp_without_bonferroni # thresholds without bonferroni correction

# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
Y = matrix(rnorm(n*p), ncol = n, nrow=p)

# Adding a single sparse change-point to X (and not Y):
X[1:5, 26:n] = X[1:5, 26:n] +2
resX = Pilliat_test(X, 
                    threshold_dense=thresholds_test_emp$threshold_dense, 
                    thresholds_bj = thresholds_test_emp$thresholds_bj, 
                    thresholds_partial = thresholds_test_emp$thresholds_partial
)
resX
resY = Pilliat_test(Y, 
                    threshold_dense=thresholds_test_emp$threshold_dense, 
                    thresholds_bj = thresholds_test_emp$thresholds_bj, 
                    thresholds_partial = thresholds_test_emp$thresholds_partial
)
resY

















#Inspect
library(HDCD)
n = 50
p = 50
set.seed(100)
# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# Adding a single sparse change-point:
X[1:5, 26:n] = X[1:5, 26:n] +2

# Vanilla Inspect:
res = Inspect(X)
res$changepoints

# Manually setting \lambda and \xi
res = Inspect(X, 
           lambda = sqrt(log(p*log(n))/2),
           xi = 4*sqrt(log(n*p))
)
res$changepoints #estimated change-point locations

# Empirical choice of thresholds:
res = Inspect(X, empirical=TRUE, N = 100, tol = 1/100)
res$changepoints

# Manual empirical choice of thresholds (equivalent to the above)
thresholds_emp = Inspect_calibrate(n,p, N=100, tol=1/100)
res = Inspect(X, xi = thresholds_emp$max_value)
res$changepoints







# Inspect_calibrate
library(HDCD)
n = 50
p = 50

set.seed(100)
thresholds_emp = Inspect_calibrate(n,p, N=100, tol=1/100)
thresholds_emp$max_value # xi

# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# Adding a single sparse change-point:
X[1:5, 26:n] = X[1:5, 26:n] +2

res = Inspect(X, xi = thresholds_emp$max_value)
res$changepoints


#Inspect single change-point test
library(HDCD)
n = 50
p = 50

# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
Y = matrix(rnorm(n*p), ncol = n, nrow=p)

# Adding a single sparse change-point to X (and not Y):
X[1:5, 26:n] = X[1:5, 26:n] +2

# Vanilla Inspect:
resX = Inspect_test(X)
resX
resY = Inspect_test(Y)
resY

# Manually setting \lambda and \xi:
resX = Inspect_test(X, 
                 lambda = sqrt(log(p*log(n))/2),
                 xi = 4*sqrt(log(n*p))
)
resX 
resY = Inspect_test(Y, 
                 lambda = sqrt(log(p*log(n))/2),
                 xi = 4*sqrt(log(n*p))
)
resY

# Empirical choice of thresholds:
resX = Inspect_test(X, empirical = TRUE, N = 100, tol = 1/100)
resX
resY = Inspect_test(Y, empirical = TRUE, N = 100, tol = 1/100)
resY

# Manual empirical choice of thresholds (equivalent to the above)
thresholds_test_emp = Inspect_test_calibrate(n,p, N=100, tol=1/100)
resX = Inspect_test(X, xi = thresholds_test_emp$max_value)
resX
resY = Inspect_test(Y, xi = thresholds_test_emp$max_value)
resY


# Inspect_test_calibrate
library(HDCD)
n = 50
p = 50

set.seed(100)
thresholds_emp = Inspect_test_calibrate(n,p,N=100, tol=1/100)
thresholds_emp


# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
Y = matrix(rnorm(n*p), ncol = n, nrow=p)

# Adding a single sparse change-point to X (and not Y):
X[1:5, 26:n] = X[1:5, 26:n] +2
resX = Inspect_test(X, xi = thresholds_emp$max_value)
resX
resY = Inspect_test(Y,  xi = thresholds_emp$max_value)
resY





















# Single ESAC
library(HDCD)
n = 500
p = 500
set.seed(101)
# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# Adding a single sparse change-point:
X[1:5, 201:500] = X[1:5, 201:500] +2

res = single_ESAC(X,rescale_variance=TRUE)
res$pos

# Manually setting the leading constants for \lambda(t):
res = single_ESAC(X, threshold_d = 2, threshold_s = 2)
res$pos


# Single Inspect
library(HDCD)
n = 500
p = 500
set.seed(101)
# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# Adding a single sparse change-point:
X[1:5, 201:500] = X[1:5, 201:500] +2

res = single_Inspect(X,rescale_variance=TRUE)
res$pos

# Manually setting the leading constants for \lambda(t):
res = single_Inspect(X, lambda = 2*sqrt(log(p*log(n))/2))
res$pos




# Single SBS
library(HDCD)
n = 500
p = 500
set.seed(101)
# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# Adding a single sparse change-point:
X[1:5, 201:500] = X[1:5, 201:500] +2

res = single_SBS(X,threshold=7,rescale_variance=TRUE)
res$pos

# Choose threhsold by Monte Carlo:
res = single_SBS(X,empirical=TRUE,rescale_variance=TRUE)
res$pos


# Single SBS threshold simu
library(HDCD)
n = 50
p = 50
set.seed(101)

# Simulate threshold
pi_T_squared = single_SBS_calibrate(n=n,p=p,N=100, tol=1/100, rescale_variance = TRUE)
pi_T_squared


# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# Adding a single sparse change-point:
X[1:5, 26:n] = X[1:5, 26:n] +2

# Run SBS
res = single_SBS(X,threshold=sqrt(pi_T_squared),rescale_variance=TRUE)
res$pos




# Rescale variance
library(HDCD)
n = 50
p = 500
set.seed(101)
# Generating data
X = matrix(rnorm(n*p), ncol = n, nrow=p)

ret = rescale_variance(X)
res$X #rescaled matrix
ret$scales #estimated noise level for each time series (each row)

# Note that the rescaled matrix is in (column wise) vector form. To transform it back to a matrix,
# do the following:
rescaled_X = matrix(ret$X, nrow = p, ncol=n)



# Hausdorff distance
n = 400
true_changepoints = c(50,100)
est_changepoints = c(51,110)
hausdorff(true_changepoints, est_changepoints,n)
hausdorff(true_changepoints, NULL,n)
hausdorff(NULL, est_changepoints,n)
hausdorff(NULL,NULL)

# ARI
library(HDCD)
n = 400
true_changepoints = c(50,100)
est_changepoints = c(51,110)
ARI(true_changepoints, est_changepoints,n)


# CUSUM
n = 10
p = 10
set.seed(101)
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# CUSUM over the full data (s,e] = (0,n]
X_cusum = CUSUM(X)
X_cusum

# CUSUM over (s,e] = (3,9]:
s = 3
e = 9
X_cusum = CUSUM(X, start = s-1, stop = e-1)


# single CUSUM
n = 10
p = 10
set.seed(101)
X = matrix(rnorm(n*p), ncol = n, nrow=p)
# CUSUM over the full data (s,e] = (0,n] evaluated at position v=4
position = 4
X_cusum_single = single_CUSUM(X,pos = position-1)
X_cusum_single

# verifying that this corresponds to the 4-th row of output of CUSUM():
X_cusum = CUSUM(X)
X_cusum[,4]
