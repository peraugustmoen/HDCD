# SUBSET
Sparse and Dense Changepoint Detection for High-Dimensional Time Series

SUBSET (Sparse and Ubiquitous Binary Segmentation in Efficient Time) is a changepoint detection method designed for the offline, multivariate setting. For full details of the method, and an application to the detection of change in terror levels in various global regions, see "A computationally efficient, high-dimensional multiple changepoint procedure with application to global terrorism incidence", by Tickle, Eckley and Fearnhead, which can be found here: https://rss.onlinelibrary.wiley.com/doi/10.1111/rssa.12695. The mechanism by which multiple changepoints are detected is called Wild Binary Segmentation, more details of which can be found in "Wild Binary Segmentation for multiple change-point detection" by Fryzlewicz, which can be found here: https://projecteuclid.org/journals/annals-of-statistics/volume-42/issue-6/Wild-binary-segmentation-for-multiple-change-point-detection/10.1214/14-AOS1245.full.

The central functions of interest here can be found in main.R, in particular change_main and wbs_penaltyfinder. 

change_main performs the sparse and dense change detection using a test statistic of interest. It expects four entries. The first of these is the data itself, given as a d by n matrix in which d represents the number of variates and n represents the number of time points. The second is the name of the method, given as a string; currently supported are "BinWeight" (corresponding to the Sparsified Binary Segmentation, or SBS, approach of Cho and Fryzlewicz, which can be found here: https://www.jstor.org/stable/24774746); "inspect.amoc" (corresponding to the Inspect method of Wang and Samworth, which can be found here: https://rss.onlinelibrary.wiley.com/doi/pdf/10.1111/rssb.12243); "Max"; "Mean" (both corresponding to cross-variate aggregations of the CUSUM statistics given in, for instance: https://onlinelibrary.wiley.com/doi/epdf/10.1002/jae.1272); "SUBSET.negbin" and "SUBSET.normal" (both of which are utilisations of the SUBSET method, the first on data with abrupt changes in the success probability parameter (with a possibly changing over-dispersion) and the second on data with abrupt changes in mean under Gaussian noise). The third input to change_main is the number of intervals which are simulated at the beginning of the procedure to allow for Wild Binary Segmentation to detect multiple changepoints. The final input is the vector of (if appropriate - some of the above methods may only require a single value) penalties. These/this can be some value from theory (see, for instance, the penalties derived for the SUBSET method in the change-in-mean setting in the article above), but the penalties can also be computed using the second function in main.R, namely wbs_penaltyfinder.

wbs_penaltyfinder computes an empirical value of the threshold/penalty needed to admit a single change; it takes three arguments - the data as a d by n matrix, the name of the method (any of the methods above appended with "_penalty", for example "Max_penalty") and the number of intervals simulated for the Wild Binary Segmentation multiple change detector.

For instance, if you were interested in detecting changes in mean under standard Gaussian noise in a 100-variate setting, you could run something like the following:

source('main.R')

mydata <- matrix(c(rnorm(2500,0,1),rnorm(2500,1,1)),nrow=5,ncol=1000,byrow=FALSE) # 5 variates with 1000 time points, change in all at time 500

mynulldata <- matrix(rnorm(5000,0,1),nrow=5,ncol=1000,byrow=FALSE) # 5 variates with 1000 time points, no change

empirical_penalty <- wbs_penaltyfinder(mynulldata, SUBSET.normal_penalty, 500) # finds an empirical penalty for the SUBSET method (with a normal likelihood cost function)

result <- change_main(mydata, SUBSET.normal, 500, empirical_penalty) # running the change detection wrapper with the SUBSET method, using the empirical penalty calculated above; result is an S3 list object with 3 elements. The first element gives the value of the test statistic (the difference in the cost) at the detected changepoints, the second element gives the location of the detected changepoints, and the final element is a matrix of binaries indicating which variates are affected by the change.
