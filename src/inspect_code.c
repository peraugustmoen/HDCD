/*#define A(r, c) (A[(r) + DEPTH*(c)])
#define AA(r, c) (AA[(r) + DEPTH*(c)])
#define B(r, c) (B[(r) + DEPTH*(c)])*/

#define cord(r,c) ((r) + (DEPTH)*(c))
#define cord_spec(r,c, D) ((r) + (D)*(c))

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
//#include <Rinternals.h> // RK addition
#include <R_ext/RS.h>	// RK addition
#include <R_ext/Lapack.h> // RK addition
#include <R_ext/BLAS.h> // RK addition
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

// function to multiply two matrices
void internal_matmult(double * first,
                      double *second,
                      double * result,
                      int r1, int c1, int r2, int c2) {

   // Initializing elements of matrix mult to 0.
   /*for (int i = 0; i < r1; ++i) {
      for (int j = 0; j < c2; ++j) {
         result[i][j] = 0;
      }
   }*/
   memset(result, 0, r1*c2*sizeof(double));

   // Multiplying first and second matrices and storing it in result
   for (int i = 0; i < r1; ++i) {
      for (int j = 0; j < c2; ++j) {
         for (int k = 0; k < c1; ++k) {
            result[cord_spec(i,j,r1)] += first[cord_spec(i,k,r1)] * second[cord_spec(k,j,r2)];
         }
      }
   }
}


SEXP matmult(SEXP AI, SEXP BI, SEXP r1I, SEXP c1I, SEXP r2I, SEXP c2I){
    PROTECT(AI);
    PROTECT(BI);
    PROTECT(c1I);
    PROTECT(r1I);
    PROTECT(r2I);
    PROTECT(c2I);

    double * A = REAL(AI);
    double * B = REAL(BI);
    int r1 = *INTEGER(r1I);
    int c1 = *INTEGER(c1I);
    int c2 = *INTEGER(c2I);
    int r2 = *INTEGER(r2I);

    UNPROTECT(4);
    SEXP out = (allocVector(NILSXP,1));

    if(c1 != r2){
        printf("matrix dims do not match");
        return out;
    }

    out = PROTECT(allocVector(REALSXP, r1*c2));
    double * res = REAL(out);

    internal_matmult(A,B, res, r1,c1,r2,c2);

    UNPROTECT(3);
    return(out);
}

// computes A %*% t(A)
void internal_matmultrightT(double * first,
                      double * result,
                      int r1, int c1) {


   memset(result, 0, r1*r1*sizeof(double));

   // Multiplying first and second matrices and storing it in result
   for (int i = 0; i < r1; ++i) {
      for (int j = 0; j < r1; ++j) {
         for (int k = 0; k < c1; ++k) {
            result[cord_spec(i,j,r1)] += first[cord_spec(i,k,r1)] * first[cord_spec(j,k,r1)];
         }
      }
   }
}


SEXP matmultrightT(SEXP AI, SEXP r1I, SEXP c1I){
    PROTECT(AI);
    PROTECT(c1I);
    PROTECT(r1I);


    double * A = REAL(AI);
    int r1 = *INTEGER(r1I);
    int c1 = *INTEGER(c1I);

    UNPROTECT(2);
    SEXP out = (allocVector(NILSXP,1));


    out = PROTECT(allocVector(REALSXP, r1*r1));
    double * res = REAL(out);

    internal_matmultrightT(A,res, r1,c1);

    UNPROTECT(2);
    return(out);
}

// computes t(A)%*% A
void internal_matmultleftT(double * first,
                      double * result,
                      int r1, int c1) {


   memset(result, 0, c1*c1*sizeof(double));

   // Multiplying first and second matrices and storing it in result
   for (int i = 0; i < c1; ++i) {
      for (int j = 0; j < c1; ++j) {
         for (int k = 0; k < r1; ++k) {
            result[cord_spec(i,j,c1)] += first[cord_spec(k,i,r1)] * first[cord_spec(k,j,r1)];
         }
      }
   }


}


SEXP matmultleftT(SEXP AI, SEXP r1I, SEXP c1I){
    PROTECT(AI);
    PROTECT(c1I);
    PROTECT(r1I);


    double * A = REAL(AI);
    int r1 = *INTEGER(r1I);
    int c1 = *INTEGER(c1I);

    UNPROTECT(2);
    SEXP out = (allocVector(NILSXP,1));


    out = PROTECT(allocVector(REALSXP, c1*c1));
    double * res = REAL(out);

    internal_matmultleftT(A,res, r1,c1);

    UNPROTECT(2);
    return(out);
}

void internal_soft_thresh(double * x, int len, double lambda){
    double tmp;
    for (int i = 0; i < len; ++i)
    {
        tmp = fabs(x[i]) - lambda;
        if(tmp<0){
            tmp=0.0;
        }
        else{
            if(x[i]<0){
                tmp = -tmp;
            }
        }
        x[i] = tmp;
    }
}

SEXP soft_thresh(SEXP xI, SEXP lenI, SEXP lambdaI){
    PROTECT(xI);
    PROTECT(lenI);
    PROTECT(lambdaI);


    int len = *INTEGER(lenI);
    double lambda = *REAL(lambdaI);
    double * x = REAL(xI);

    UNPROTECT(2);


    //SEXP out = PROTECT(allocVector(REALSXP, 1));

    internal_soft_thresh(x, len, lambda);

    UNPROTECT(1);
    return(xI);
}

double * internal_power_method(double * A, int n, double eps, int maxiter, double * vec1, double * vec2,int debug ){
    if(maxiter == 0){
        maxiter = 10000;
    }

    // need vec1, vec2 to be len n
/*    SEXP out = PROTECT(allocVector(REALSXP, n));
    SEXP out2 = PROTECT(allocVector(REALSXP, n));
    SEXP tmpsexp;
    double * v = REAL(out);
    double * v2 = REAL(out2);
*/
    double * v = vec1;
    double * v2 = vec2;
    double * tmp;

    //printf("n = %d\n", n);
    //printf("eps = %f\n", eps);
    //printf("maxiter = %d\n", maxiter);
    // init v
    double zum =0.0;
    GetRNGstate();

    for (int i = 0; i < n; ++i)
    {
        v[i] = norm_rand();
    }
    PutRNGstate();
    for (int i = 0; i < n; ++i)
    {
        zum += v[i]*v[i];
    }
    double sumsq = sqrt(zum);
    for (int i = 0; i < n; ++i)
    {
        v[i] = v[i] /sumsq;
    }
    int i;
    double diff = 0.0;
    zum=0.0;
    for (i = 0; i < maxiter; ++i)
    {
        //void internal_matmult(double * first, double *second, double * result,int r1, int c1, int r2, int c2)
        internal_matmult(A, v, v2, n,n,n,1);
        zum = 0.0;
        for (int j = 0; j < n; ++j)
        {
            zum+= v2[j]*v2[j];
        }
        sumsq = sqrt(zum);
        //printf("ZUM = %f\n", zum);
        if (fabs(zum)<1e-15)
        {
            if(debug){
                printf("ERROR IN POWERMETHOD: REACHED 0 VECTOR\n");
            }

            return NULL;
        }
        diff = 0.0;
        for (int j = 0; j < n; ++j)
        {
            v2[j] = v2[j]/sumsq;
            diff+= (v2[j] - v[j])*(v2[j] - v[j]);
        }

        // switching
        tmp = v;
        //tmpsexp = out;
        v = v2;
        //out = out2;
        v2 = tmp;
        //out2 = tmpsexp;
        //printf("diff=%f\n", diff);
        if(diff<eps){
            break;
        }

    }
    if(i==(maxiter-1)){
        printf("WARNING: power method did not converge");
    }
    if (debug) {
        printf("num iter: %d\n", i);
    }
    //UNPROTECT(2);
    return v;



}

//double* internal_power_method(double * A, int n, double eps, int maxiter){
SEXP power_method(SEXP AI, SEXP nI, SEXP epsI, SEXP maxiterI){
    PROTECT(AI);
    PROTECT(nI);
    PROTECT(epsI);
    PROTECT(maxiterI);

    double * A = REAL(AI);
    int n = *(INTEGER(nI));
    double eps = *(REAL(epsI));
    int maxiter = *(INTEGER(maxiterI));
    UNPROTECT(3);
    SEXP vSEXP = PROTECT(allocVector(REALSXP, n));
    SEXP v2SEXP= PROTECT(allocVector(REALSXP, n));
    double * v = REAL(vSEXP);
    double * v2 = REAL(v2SEXP);

    double * res = internal_power_method(A, n, eps, maxiter, v, v2,0);

    SEXP resSEXP = v2SEXP;
    if (res == v)
    {
        resSEXP = vSEXP;
    }
    UNPROTECT(3);
    return(resSEXP);


}

double * internal_sparse_svd(double * Z, int r1, int c1, double lambda, double eps, int maxiter,
                        double * mhat, double * mhatprod, double * v, double * v2,int debug){
    // only schatten=2 for now
    //
    // call soft thresholding

    // need mhat to be r1*c1
    // and mhatprod to be (min(r1, c1))^2
    // result needs to be p

    //SEXP mhatSEXP = PROTECT(allocVector(REALSXP, r1*c1));
    //double *mhat = REAL(mhatSEXP);
    //SEXP mhatprodSEXP;
    memcpy(mhat, Z, r1*c1*sizeof(double));
    internal_soft_thresh(mhat, r1*c1, lambda);
    // call power method

    //SEXP ret;
    //SEXP tmpret;
    double * projection;
    double * tmproj;
    double summ=0;
    double sumsq=0;
    if(r1<c1){
        //internal_power_method(double * A, int n, double eps, int maxiter)
        //mhatprodSEXP = PROTECT(allocVector(REALSXP, r1*r1));
        //double *mhatprod = REAL(mhatprodSEXP);
        internal_matmultrightT(mhat, mhatprod,r1, c1);
        //internal_power_method(mhatprod, r1, eps, maxiter, double * vec1, double * vec2 )
        projection = internal_power_method(mhatprod, r1, eps, maxiter, v, v2,debug);
        if(projection==NULL){
          return NULL;
        }
    }
    else{

        //mhatprodSEXP = PROTECT(allocVector(REALSXP, c1*c1));
        //double *mhatprod = REAL(mhatprodSEXP);
        internal_matmultleftT(mhat, mhatprod,r1, c1);

        tmproj = internal_power_method(mhatprod, c1, eps, maxiter,v,v2,debug);
        if(tmproj==NULL){
          return NULL;
        }
        //PROTECT(tmpret);
        //tmpretarr = REAL(tmpret);
        //ret = PROTECT(allocVector(REALSXP, c1));
        //tmparr = REAL(ret);
        projection = v;
        if(tmproj ==v){
            projection = v2;
        }
        internal_matmult(mhat, tmproj, projection, r1, c1, r1, 1);

        for (int i = 0; i < r1; ++i)
        {
            summ+=projection[i]*projection[i];
        }
        sumsq=sqrt(summ);
        for (int i = 0; i < r1; ++i)
        {
            projection[i] = projection[i]/sumsq;
        }


    }

    /*UNPROTECT(3);
    if (r1>=c1)
    {
        UNPROTECT(1);
    }*/
    return projection;
}

SEXP sparse_svd(SEXP ZI, SEXP r1I, SEXP c1I, SEXP lambdaI, SEXP epsI, SEXP maxiterI){
    PROTECT(ZI);
    PROTECT(r1I);
    PROTECT(c1I);
    PROTECT(lambdaI);
    PROTECT(epsI);
    PROTECT(maxiterI);

    double *Z= REAL(ZI);
    int r1 = *(INTEGER(r1I));
    int c1 = *(INTEGER(c1I));
    double lambda= *(REAL(lambdaI));
    double eps = *(REAL(epsI));
    int maxiter= *(INTEGER(maxiterI));
    UNPROTECT(5);

    int maxlen = r1;
    int minlen = c1;
    if(c1>r1){
        maxlen = c1;
        minlen = r1;
    }
    SEXP vec1SEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP vec2SEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP mhatSEXP = PROTECT(allocVector(REALSXP, r1*c1));
    SEXP mhatprodSEXP = PROTECT(allocVector(REALSXP, minlen*minlen));

    double * vec1 = REAL(vec1SEXP);
    double * vec2 = REAL(vec2SEXP);
    double * mhat = REAL(mhatSEXP);
    double * mhatprod = REAL(mhatprodSEXP);

    //internal_sparse_svd(double * Z, int r1, int c1, double lambda, double eps, int maxiter,
                        //double * mhat, double * mhatprod, double * v, double * v2)

    double* retarr = internal_sparse_svd(Z, r1, c1, lambda, eps, maxiter,
                        mhat, mhatprod, vec1, vec2,0);
    SEXP ret = vec2SEXP;
    if(retarr==vec1){
        ret = vec1SEXP;
    }
    UNPROTECT(5);
    return ret;

}





void inspectCUSUM(double * cumsums, double * cusum, int s, int e, int p){
    if(e-s<2){
        return;
    }
    int n = e-s;
    int t;
    //printf("%f", cumsum[5]);
    for (int i = 0; i < p; ++i)
    {
        for (int j = 1; j < n; ++j)
        {
            t = s+j;
            cusum[cord_spec(i,j-1, p)] = sqrt(((double)(e-t))/((e-s)*(t-s))) *(cumsums[cord_spec(i,t,p)]-cumsums[cord_spec(i,s,p)])
            - sqrt(((double)(t-s))/((e-s)*(e-t)))*(cumsums[cord_spec(i,e,p)] - cumsums[cord_spec(i,t,p)]);

        }

    }

    return;
}

void internal_inspectOnSegment(double * cumsums, double * cusum, int * maxpos, double * maximum, int s, int e, int p, double lambda,
    double eps, int maxiter, double * mhat, double * mhatprod, double* v, double* v2,int debug){
    *maxpos  = e;
    *maximum = 0.0;
    if(e-s<2){
        return;
    }

    // compute CUSUM
    inspectCUSUM(cumsums, cusum, s, e, p);

    // find sparse SVD
    //internal_sparse_svd(double * Z, int r1, int c1, double lambda, double eps, int maxiter,
                        //double * mhat, double * mhatprod, double * v, double * v2)
    double * projvec = internal_sparse_svd(cusum, p, e-s-1, lambda*sqrt(log(p*(e-s-1))), eps, maxiter,
                        mhat, mhatprod, v, v2,debug);
    if(projvec==NULL){
      if(debug){
          printf("inspecting segment, s=%d, e=%d resulted in NULL projection. lambda = %f.\n", s,e,lambda*sqrt(log(  ((double)(p*(e-s-1))))));
      }

      return;
    }
    double * projected = v;
    if(projvec==v){
        projected = v2;
    }
    //double * v = REAL(retval);
    int n = e-s;
    double tmp;
    internal_matmult(projvec,cusum,  projected,1, p, p, n-1);

    for (int i = 0; i < n-1; ++i)
    {
        //t = s+i+1;
        tmp = fabs(projected[i]);
        if(tmp> *maximum){
            *maximum = tmp;
            *maxpos = s+i+1;
        }
    }
    if(debug){
        printf("inspecting segment, s=%d, e=%d, max_cusum = %f\n", s,e,*maximum);
    }


    return;
}

void myInspect_call(double * x, int s, int e, int n, int p, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter,
                double * maxval, double threshold,int adaptTresh, double *cumsums, int* lens, int lenLens, double lambda,
                double eps, int maxiter, int * segstarts, double * maxcusums, int* maxpos, int K, double * cusum, double * mhat,
                double * mhatprod, double * v, double * v2, int debug,int * coordchg){
    if(debug){
        printf("myInspectCall! s=%d, e=%d\n", s, e);
    }

    if(e-s<lens[0]){
        //printf("segment too short\n");
        return;
    }
    int argmax = s;
    double maximum = 0;

    //int tmpargmax;
    //double tmpmaximum;
    double tmp;
    int len;
    int jump;
    //int found=0;

    int i;

    double thresh = threshold;
    int j_max=0;
    int k_max = 0;
    for (int j = 0; j < lenLens; ++j)
    {
        len = lens[j];
        if(debug){
            printf("j=%d, len = %d\n", j, len);
        }

        jump = len /K;
        if(jump==0){
            jump = 1;
        }
        if(adaptTresh){
            thresh = threshold*sqrt(log(len*p));
        }
        if(e-s<len){
            break;
        }

        for (int k = 0; k < n; ++k)
        {
            i = segstarts[cord_spec(k,j,n)];

            if (i<s)
            {
                continue;
            }
            else if(i>e-len){
                break;
            }

            if(debug){
                printf("maxcusums[%d, %d] = %f\n", k , j , maxcusums[cord_spec(k,j,n)]);
            }

            if(maxcusums[cord_spec(k,j,n)]<=0.0){
                //this segment not computed!
                //inspectOnSegment(cumsums, cusum, &tmpargmax, &tmpmaximum, s, e, p, lambda,
                 //   eps, maxiter, projvec, cusum_proj);
                 //double * cumsums, double * cusum, int * maxpos, double * maximum,
                 internal_inspectOnSegment(cumsums, cusum, &(maxpos[cord_spec(k,j,n)]), &(maxcusums[cord_spec(k,j,n)]), i, i+len, p,
                 lambda,
                        eps, maxiter, mhat, mhatprod, v, v2,debug);
            }

            tmp = maxcusums[cord_spec(k,j,n)];
            if(tmp>maximum){
                maximum = tmp;
                argmax = maxpos[cord_spec(k,j,n)];
                j_max = j;
                k_max=k;
                //found=1;
            }


        }

        if(maximum>thresh){
            break;
        }
    }
    if(debug){
        printf("maximum=%f\n", maximum);
    }

    if(maximum > thresh){
        if(debug){
            printf("!!!!!! declared change-point in %d. val = %f, thresh =%f\n", argmax, maximum, thresh);
        }
        // identify in which coordinates the change happens:
        i = segstarts[cord_spec(k_max,j_max,n)];
        len = lens[j_max];
        int ss = i;
        int ee = i+len;
        inspectCUSUM(cumsums, cusum, ss, ee, p);
        double * projvec = internal_sparse_svd(cusum, p, ee-ss-1, lambda*sqrt(log(p*(ee-ss-1))), eps, maxiter,
                        mhat, mhatprod, v, v2,debug);
        for (int zz = 0; zz < p; ++zz)
        {
            if(fabs(projvec[zz])>1e-6){
                coordchg[cord_spec(zz,*changepoint_counter_ptr, p)]=1;
            }
        }
        
        changepoints[*changepoint_counter_ptr] = argmax;
        depthcounter[*changepoint_counter_ptr] = depth;
        maxval[*changepoint_counter_ptr] = maximum;
        (*changepoint_counter_ptr)++;
        //myInspect_call(x, s, argmax, n, p, depth+1, changepoints, changepoint_counter_ptr, depthcounter,  maxval, threshold,
        //        adaptTresh, cumsums, lens, lenLens, K, cusum,  mhat, mhatprod, v, v2);
        myInspect_call(x, s, argmax, n, p, depth+1, changepoints, changepoint_counter_ptr, depthcounter,
                maxval,threshold, adaptTresh, cumsums, lens, lenLens, lambda,
                eps, maxiter, segstarts, maxcusums, maxpos, K,  cusum, mhat,
                mhatprod, v, v2,debug,coordchg);
        //myInspect_call(x, argmax, e, n, p, depth+1, changepoints, changepoint_counter_ptr, depthcounter,  maxval, threshold,
        //        adaptTresh, cumsums, lens, lenLens, K, cusum, mhat, mhatprod, v, v2);
        myInspect_call(x, argmax, e, n, p, depth+1, changepoints, changepoint_counter_ptr, depthcounter,
                maxval,threshold, adaptTresh, cumsums, lens, lenLens, lambda,
                eps, maxiter, segstarts, maxcusums, maxpos, K,  cusum, mhat,
                mhatprod, v, v2,debug,coordchg);

    }

    return;
}

SEXP myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI,
    SEXP epsI, SEXP lambdaI, SEXP maxiterI, SEXP debugI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(lensI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(thresholdI);
    PROTECT(adaptTreshI);
    PROTECT(lenLensI);
    PROTECT(KI);
    PROTECT(epsI);
    PROTECT(lambdaI);
    PROTECT(maxiterI);
    PROTECT(debugI);

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double threshold = *(REAL(thresholdI));
    int adaptTresh = *(INTEGER(adaptTreshI));
    int *lens = INTEGER(lensI); /////// dobbeltsjekk at int er like stor som INTEGER!!!!!!!
    int lenLens = *(INTEGER(lenLensI));
    int K = *(INTEGER(KI));
    double eps = *REAL(epsI);
    double lambda = *REAL(lambdaI);
    int maxiter = *INTEGER(maxiterI);
    int debug = *INTEGER(debugI);
    UNPROTECT(11); // unprotecting all except X and lens
    if(debug){
      printf("p = %d\n", p);
      printf("lambda = %f\n", lambda);
    }

    SEXP out = PROTECT(allocVector(INTSXP, n));
    int * changepoints = INTEGER(out); //pointer to array
    memset(changepoints, 0, sizeof(int)*n);
    int changepoint_counter = 0;
    int* changepoint_counter_ptr = &changepoint_counter;
    SEXP maxvalSEXP = PROTECT(allocVector(REALSXP, n));
    double * maxval = REAL(maxvalSEXP); //pointer to array
    memset(maxval, 0, sizeof(double)*n);
    SEXP depthcounterSEXP= PROTECT(allocVector(INTSXP, n));
    int * depthcounter = INTEGER(depthcounterSEXP); //pointer to array
    memset(depthcounter, 0, sizeof(int)*n);
    SEXP coordschgSEXP = PROTECT(allocVector(INTSXP, n*p));
    int * coordchg = INTEGER(coordschgSEXP);
    memset(coordchg, 0, sizeof(int)*n*p);
    // first we compute all the cumulative sums of all
    // coordinates:

    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*sizeof(double));

    for (int j = 1; j <=n; ++j)
    {
        for (int i = 0; i < p; ++i)
        {
            //#define cord_spec(r,c, D) ((r) + (D)*(c))
            cumsums[cord_spec(i,j,p)] = X[cord_spec(i,j-1, p)] +cumsums[cord_spec(i,j-1, p)];
        }
    }

    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p*(n)));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p*(n));
    int maxlen = p;
    int minlen = n;
    if(n>p){
        maxlen = n;
        minlen = p;
    }
    SEXP vSEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP v2SEXP = PROTECT(allocVector(REALSXP, maxlen));
    SEXP mhatSEXP = PROTECT(allocVector(REALSXP, p*n));
    SEXP mhatprodSEXP = PROTECT(allocVector(REALSXP, minlen*minlen));

    double * v = REAL(vSEXP);
    memset(v, 0, sizeof(double)*maxlen);
    double * v2 = REAL(v2SEXP);
    memset(v2, 0, sizeof(double)*maxlen);
    double * mhat = REAL(mhatSEXP);
    memset(mhat, 0, sizeof(double)*p*n);
    double * mhatprod = REAL(mhatprodSEXP);
    memset(mhatprod, 0, sizeof(double)*minlen*minlen);



    // record all segment starts
    // we don't compute all cusums here but do it on the fly so we dont have to do anything
    // unecessary


    SEXP maxcusumsSEXP = PROTECT(allocVector(REALSXP, n*lenLens));
    double * maxcusums = REAL(maxcusumsSEXP);
    memset(maxcusums, 0, sizeof(double)*n*lenLens);
    SEXP maxposSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * maxpos = INTEGER(maxposSEXP);
    memset(maxpos, 0, sizeof(int)*n*lenLens);
    SEXP segstartsSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * segstarts = INTEGER(segstartsSEXP);
    memset(segstarts, 0, sizeof(int)*n*lenLens);

    int len;
    int jump;
    int counter = 0;
    for (int j = 0; j < lenLens; ++j)
    {
        counter=0;
        len = lens[j];
        jump = len/K;
        if(jump<1){
            jump=1;
        }

        for (int i = 0; i <= (n-len); i+=jump)
        {
            //cord_spec(r,c, D) ((r) + (D)*(c))
            //compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxcusums[cord_spec(i,j,n)]));
            segstarts[cord_spec(counter++,j,n)] = i;
            if(debug){
                printf("segstarts[%d, %d] = %d\n",counter-1, j, i );
            }


        }
    }











    //myInspect_call(X, 0, n, n, p, 1, changepoints, changepoint_counter_ptr, depthcounter, maxval, threshold,
    //            adaptTresh, cumsums, lens, lenLens, K, cusum, mhat, mhatprod,  v, v2);
    myInspect_call(X, 0, n, n, p, 1, changepoints, changepoint_counter_ptr, depthcounter,
                maxval, threshold, adaptTresh, cumsums, lens, lenLens, lambda,
                eps, maxiter, segstarts, maxcusums, maxpos, K, cusum, mhat,
                mhatprod, v, v2,debug,coordchg);

    // return:
    SEXP ret = PROTECT(allocVector(VECSXP, 4)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out);
    SET_VECTOR_ELT(ret, 1, maxvalSEXP);
    SET_VECTOR_ELT(ret, 2, depthcounterSEXP);
    SET_VECTOR_ELT(ret, 3, coordschgSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 4));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("changepoints"));
    SET_STRING_ELT(names, 1, mkChar("CUSUMval"));
    SET_STRING_ELT(names, 2, mkChar("depth"));
    SET_STRING_ELT(names, 3, mkChar("coordinate"));

    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(16);
    return ret;
}




    /*    int n = *(INTEGER(nI));
    double * X = REAL(XI);
    double threshold = *REAL(thresholdI);
    int *lens = (INTEGER(lensI));
    int lenLens = *INTEGER(lenLensI);

    //printf("lenLens:%d\n", lenLens);
    SEXP out1 = PROTECT(allocVector(INTSXP, n));
    int * changepoints = INTEGER(out1); //pointer to array
    memset(changepoints, 0, sizeof(int)*n);
    int changepoint_counter = 0;
    int* changepoint_counter_ptr = &changepoint_counter;
    //printf("size of int :%lu\n", sizeof(int));

    SEXP out4 = PROTECT(allocVector(REALSXP, n+1));
    double * cumsum = REAL(out4); //pointer to array
    SEXP out2 = PROTECT(allocVector(REALSXP, n));
    double * maxval = REAL(out2); //pointer to array
    SEXP out3= PROTECT(allocVector(INTSXP, n));
    memset(maxval, 0, sizeof(double)*n);
    int * depthcounter = INTEGER(out3); //pointer to array
    memset(depthcounter, 0, sizeof(int)*n);
    cumsum[0] = 0;
    for (int i = 0; i < n; ++i)
    {
        cumsum[i+1] = cumsum[i] + X[i];
        //printf("CUMSUM: %f\n", cumsum[i+1]);
    }

    SEXP out5 = PROTECT(allocVector(REALSXP, n*lenLens));
    double * maxcusums = REAL(out5);
    memset(maxcusums, -1, sizeof(double)*n*lenLens);
    SEXP out6 = PROTECT(allocVector(INTSXP, n*lenLens));
    int * maxpos = INTEGER(out6);

    int len;
    int jump;
    for (int j = 0; j < lenLens; ++j)
    {
        len = lens[j];
        jump = len/10;
        if(jump<1){
            jump=1;
        }

        for (int i = 0; i <= (n-len); i+=jump)
        {
            //cord_spec(r,c, D) ((r) + (D)*(c))
            compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxcusums[cord_spec(i,j,n)]));

        }
    }







    // call WBS here:
    //NOT_call(X, 0, n, n, 1, changepoints, changepoint_counter_ptr, depthcounter, maxval,threshold,cumsum, starts, stops, M);
    myWBS2_call(X, 0, n, n, 1, changepoints, changepoint_counter_ptr, depthcounter, maxval, threshold, maxcusums, maxpos, lens, lenLens);
    // creating a list to hold the two return values. this is the returned value.
    SEXP ret = PROTECT(allocVector(VECSXP, 3)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out1);
    SET_VECTOR_ELT(ret, 1, out2);
    SET_VECTOR_ELT(ret, 2, out3);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 3));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("changepoints"));
    SET_STRING_ELT(names, 1, mkChar("CUSUMval"));
    SET_STRING_ELT(names, 2, mkChar("depth"));

    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(10);
   return out;

}*/

/*void inspectOnSegment(double * cumsums, double * cusum, int * maxpos, double * maximum, int s, int e, int p, double lambda,
    double eps, int maxiter, double * projvec, double * cusum_proj){
    *maxpos  = e;
    *maximum = 0.0;
    if(e-s<2){
        return;
    }

    // compute CUSUM
    inspectCUSUM(cumsums, cusum, s, e, p)

    // find sparse SVD

    SEXP retval = PROTECT( internal_sparse_svd(cusum, p, e-s-1, lambda, eps, maxiter));

    double * v = REAL(retval);
    int n = e-s;
    int t;
    double tmp;
    internal_matmult(v,cusum,  cusum_proj,1, p, p, n-1);

    for (int i = 0; i < n-1; ++i)
    {
        t = s+i+1;
        tmp = fabs(cusum_proj[i]);
        if(tmp> *maximum){
            *maximum = tmp;
            *maxpos = s+i+1;
        }
    }

    return;
}*/
/*void myInspect_call(double * x, int s, int e, int n, int p, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter, double * maxval, double threshold,
                int adaptTresh, double *cumsums, int* lens, int lenLens, int K, double * temparr1, double * temparr2, double* tempvec1, double * tempvec2){

    if(e-s<lens[0]){
        return;
    }
    int argmax = s;
    double maximum = 0;

    int tmpargmax;
    double tmpmaximum;
    int tmp;
    int len;
    int jump;
    //int found=0;

    double * cusum = temparr1;
    double * cusum_proj = tempvec1;
    double * projvec = tempvec2;
    double tresh = threshold;

    for (int j = 0; j < lenLens; ++j)
    {
        len = lens[j];
        jump = len /K;
        if(jump==0){
            jump = 1;
        }
        if(adaptTresh){
            tresh = treshold*sqrt(log(len*p));
        }
        if(e-s<len){
            break;
        }
        for (int i = s; i < (e-len); i+=jump )
        {

            // run inspect on this part.
            inspectOnSegment(cumsums, cusum, &tmpargmax, &tmpmaximum, s, e, p, lambda,
                eps, maxiter, projvec, cusum_proj);

            tmp = maxcusums[cord_spec(i,j,n)];
            if(tmp>maximum){
                maximum = tmp;
                argmax = maxpos[cord_spec(i,j,n)];
                //found=1;
            }
        }
        if(maximum>threshold){
            break;
        }
    }


    if(maximum > threshold){
        changepoints[*changepoint_counter_ptr] = argmax;
        depthcounter[*changepoint_counter_ptr] = depth;
        maxval[*changepoint_counter_ptr] = maximum;
        (*changepoint_counter_ptr)++;
        myInspect_call(x, s, argmax, n, depth+1,changepoints,changepoint_counter_ptr, depthcounter,  maxval,threshold, maxcusums, maxpos, lens,lenLens);
        myInspect_call(x, argmax, e, n, depth+1,changepoints,changepoint_counter_ptr, depthcounter,  maxval,threshold, maxcusums, maxpos, lens,lenLens);
    }

    return;
}*/
/*
SEXP myInspect(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdI, SEXP adaptTreshI, SEXP lensI,SEXP lenLensI,SEXP KI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(lensI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(thresholdI);
    PROTECT(adaptTreshI);
    PROTECT(lenLensI);
    PROTECT(KI);

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double tresholdI = *(REAL(thresholdI));
    int adaptTresh = *(INTEGER(adaptTreshI));
    int *lens = INTEGER(lensI); /////// dobbeltsjekk at int er like stor som INTEGER!!!!!!!
    int lenLens = *(INTEGER(lenLensI));
    int K = *(INTEGER(KI));
    UNPROTECT(6); // unprotecting all except X and lens

    SEXP out = PROTECT(allocVector(INTSXP, n));
    int * changepoints = INTEGER(out); //pointer to array
    memset(changepoints, 0, sizeof(int)*n);
    int changepoint_counter = 0;
    int* changepoint_counter_ptr = &changepoint_counter;
    SEXP maxvalSEXP = PROTECT(allocVector(REALSXP, n));
    double * maxval = REAL(maxvalSEXP); //pointer to array
    memset(maxval, 0, sizeof(double)*n);
    SEXP depthcounterSEXP= PROTECT(allocVector(INTSXP, n));
    int * depthcounter = INTEGER(depthcounterSEXP); //pointer to array
    memset(depthcounter, 0, sizeof(int)*n);

    // first we compute all the cumulative sums of all
    // coordinates:

    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0

    memset(cumsums, 0, p*sizeof(double));

    for (int j = 0; j <=n; ++j)
    {
        for (int i = 0; i < p; ++i)
        {
            //#define cord_spec(r,c, D) ((r) + (D)*(c))
            cumsums[cord_spec(i,j,p)] = X[cord_spec(i,j-1, p)] +cumsums[cord_spec(i,j-1, p)];
        }
    }

    SEXP temparr1SEXP = PROTECT(allocVector(REALSXP, p*(n+1)));
    double * temparr1 = REAL(temparr1SEXP);
    memset(temparr1, 0, sizeof(double)*p*(n+1));
    SEXP temparr2SEXP = PROTECT(allocVector(REALSXP, p*(n+1)));
    double * temparr2 = REAL(temparr2SEXP);
    memset(temparr2, 0, sizeof(double)*p*(n+1));
    SEXP tempvec1SEXP = PROTECT(allocVector(REALSXP, p);
    double * tempvec1 = REAL(tempvec1SEXP);
    memset(tempvec1, 0, sizeof(double)*p);

    //void myInspect_call(double * x, int s, int e, int n, int p, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter, double * maxval, double threshold,
    //            int adaptTresh,double *cusums, int* lens, int lenLens, int K, double * temparr1, double * temparr2, double* tempvec1)
    void myInspect_call(X, 0, n, n, p, 1, changepoints, changepoint_counter_ptr, depthcounter, maxval, threshold, adaptTresh,
                cumsums, lens, lenLens, K, temparr1,  temparr2, tempvec1);

*/
    /*    int n = *(INTEGER(nI));
    double * X = REAL(XI);
    double threshold = *REAL(thresholdI);
    int *lens = (INTEGER(lensI));
    int lenLens = *INTEGER(lenLensI);

    //printf("lenLens:%d\n", lenLens);
    SEXP out1 = PROTECT(allocVector(INTSXP, n));
    int * changepoints = INTEGER(out1); //pointer to array
    memset(changepoints, 0, sizeof(int)*n);
    int changepoint_counter = 0;
    int* changepoint_counter_ptr = &changepoint_counter;
    //printf("size of int :%lu\n", sizeof(int));

    SEXP out4 = PROTECT(allocVector(REALSXP, n+1));
    double * cumsum = REAL(out4); //pointer to array
    SEXP out2 = PROTECT(allocVector(REALSXP, n));
    double * maxval = REAL(out2); //pointer to array
    SEXP out3= PROTECT(allocVector(INTSXP, n));
    memset(maxval, 0, sizeof(double)*n);
    int * depthcounter = INTEGER(out3); //pointer to array
    memset(depthcounter, 0, sizeof(int)*n);
    cumsum[0] = 0;
    for (int i = 0; i < n; ++i)
    {
        cumsum[i+1] = cumsum[i] + X[i];
        //printf("CUMSUM: %f\n", cumsum[i+1]);
    }

    SEXP out5 = PROTECT(allocVector(REALSXP, n*lenLens));
    double * maxcusums = REAL(out5);
    memset(maxcusums, -1, sizeof(double)*n*lenLens);
    SEXP out6 = PROTECT(allocVector(INTSXP, n*lenLens));
    int * maxpos = INTEGER(out6);

    int len;
    int jump;
    for (int j = 0; j < lenLens; ++j)
    {
        len = lens[j];
        jump = len/10;
        if(jump<1){
            jump=1;
        }

        for (int i = 0; i <= (n-len); i+=jump)
        {
            //cord_spec(r,c, D) ((r) + (D)*(c))
            compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxcusums[cord_spec(i,j,n)]));

        }
    }







    // call WBS here:
    //NOT_call(X, 0, n, n, 1, changepoints, changepoint_counter_ptr, depthcounter, maxval,threshold,cumsum, starts, stops, M);
    myWBS2_call(X, 0, n, n, 1, changepoints, changepoint_counter_ptr, depthcounter, maxval, threshold, maxcusums, maxpos, lens, lenLens);
    // creating a list to hold the two return values. this is the returned value.
    SEXP ret = PROTECT(allocVector(VECSXP, 3)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out1);
    SET_VECTOR_ELT(ret, 1, out2);
    SET_VECTOR_ELT(ret, 2, out3);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 3));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("changepoints"));
    SET_STRING_ELT(names, 1, mkChar("CUSUMval"));
    SET_STRING_ELT(names, 2, mkChar("depth"));

    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(10);*/
/*    return out;

}*/


void compute_cusum(double * cumsum, int s, int e, int * res, double * val ){
    *res  = e;
    *val = 0.0;
    if(e-s<2){
        return;
    }
    int n = e-s;
    int t;
    double cusumval;
    //printf("%f", cumsum[5]);
    for (int i = 1; i < n; ++i)
    {
        t = s+i;
        cusumval = sqrt(((double)(e-t))/((e-s)*(t-s))) *(cumsum[t]-cumsum[s]) - sqrt(((double)(t-s))/((e-s)*(e-t)))*(cumsum[e] - cumsum[t]);
        //printf("[%d, %d]\n", s,e);
        //printf("%f\n", cusumval);
        //printf("%f\n", (cumsum[t]-cumsum[s]) -(cumsum[e]-cumsum[t]));
        if(fabs(cusumval)>*val){
            *val = fabs(cusumval);
            *res = t;
        }
    }

    return;

}

void WBS_call(double * x, int s, int e, int n, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter, double * maxval, double threshold,
                double* cumsum, int * starts, int*stops, int M ){

    if(e-s<1){
        return;
    }
    int argmax = s;
    double maximum = 0;
    int s_m;
    int e_m;
    int tempargmax;
    double tempmax;
    for (int m = 0; m < M; ++m)
    {
        s_m = starts[m];
        e_m = stops[m];
        if(s_m >e || e_m < s){
            continue;
        }
        //s_m = max(s_m, s);
        if (s>s_m)
        {
            s_m=s;
        }
        //e_m = min(e_m, e);
        if(e<e_m){
            e_m = e;
        }
        if(e_m - s_m<1){
            continue;
        }

        compute_cusum(cumsum, s_m, e_m, &tempargmax, &tempmax);
        //printf("%f\n",tempmax);
        if (tempmax>maximum)
        {
            maximum = tempmax;
            argmax = tempargmax;
        }

    }
    //printf("largest value on interval [s,e] = [%d, %d] is %f\n", s,e,maximum);
    if (maximum >threshold)
    {
        changepoints[*changepoint_counter_ptr] = argmax;
        depthcounter[*changepoint_counter_ptr] = depth;
        maxval[*changepoint_counter_ptr] = maximum;
        (*changepoint_counter_ptr)++;
        WBS_call(x, s, argmax, n, depth+1,changepoints, changepoint_counter_ptr, depthcounter, maxval, threshold, cumsum, starts, stops, M);
        WBS_call(x, argmax, e, n, depth+1,changepoints, changepoint_counter_ptr, depthcounter, maxval, threshold, cumsum, starts, stops, M);
    }


}

SEXP WBS(SEXP nI, SEXP thresholdI, SEXP XI, SEXP MI, SEXP STARTS, SEXP STOPS){
    PROTECT(XI);
    int n = *(INTEGER(nI));
    int M = *(INTEGER(MI));
    double * X = REAL(XI);
    double threshold = *REAL(thresholdI);

    SEXP out1 = PROTECT(allocVector(INTSXP, n));
    int * changepoints = INTEGER(out1); //pointer to array
    memset(changepoints, 0, sizeof(int)*n);
    int changepoint_counter = 0;
    int* changepoint_counter_ptr = &changepoint_counter;
    //printf("size of int :%lu\n", sizeof(int));

    SEXP out4 = PROTECT(allocVector(REALSXP, n+1));
    double * cumsum = REAL(out4); //pointer to array
    SEXP out2 = PROTECT(allocVector(REALSXP, n));
    double * maxval = REAL(out2); //pointer to array
    SEXP out3= PROTECT(allocVector(INTSXP, n));
    memset(maxval, 0, sizeof(double)*n);
    int * depthcounter = INTEGER(out3); //pointer to array
    memset(depthcounter, 0, sizeof(int)*n);
    cumsum[0] = 0;
    for (int i = 0; i < n; ++i)
    {
        cumsum[i+1] = cumsum[i] + X[i];
        //printf("CUMSUM: %f\n", cumsum[i+1]);
    }

    protect(STARTS);
    int * starts = INTEGER(STARTS); //pointer to array

    protect(STOPS);
    int * stops = INTEGER(STOPS); //pointer to array



    // call WBS here:
    WBS_call(X, 0, n, n, 1, changepoints, changepoint_counter_ptr, depthcounter, maxval,threshold,cumsum, starts, stops, M);

    // creating a list to hold the two return values. this is the returned value.
    SEXP ret = PROTECT(allocVector(VECSXP, 3)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out1);
    SET_VECTOR_ELT(ret, 1, out2);
    SET_VECTOR_ELT(ret, 2, out3);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 3));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("changepoints"));
    SET_STRING_ELT(names, 1, mkChar("CUSUMval"));
    SET_STRING_ELT(names, 2, mkChar("depth"));

    /* assign names to list */
    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(9);
    return ret;



}


void NOT_call(double * x, int s, int e, int n, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter, double * maxval, double threshold,
                double* cumsum, int * starts, int*stops, int M ){

    if(e-s<2){
        return;
    }
    int argmax = s;
    double maximum = 0;
    int max_len = n+1;
    int s_m;
    int e_m;
    int tempargmax;
    double tempmax;
    int found = 0;
    int segment_len = 0;
    for (int m = 0; m < M; ++m)
    {
        s_m = starts[m];
        e_m = stops[m];
        if(s_m <s || e_m > e){
            continue;
        }
        segment_len = e_m - s_m;


        compute_cusum(cumsum, s_m, e_m, &tempargmax, &tempmax);

        //printf("%f\n",tempmax);
        if (tempmax>threshold && segment_len<max_len)
        {
            maximum = tempmax;
            argmax = tempargmax;
            max_len = segment_len;
            found = 1;
        }

    }
    //printf("largest value on interval [s,e] = [%d, %d] is %f\n", s,e,maximum);
    if (found)
    {
        changepoints[*changepoint_counter_ptr] = argmax;
        depthcounter[*changepoint_counter_ptr] = depth;
        maxval[*changepoint_counter_ptr] = maximum;
        (*changepoint_counter_ptr)++;
        WBS_call(x, s, argmax, n, depth+1,changepoints, changepoint_counter_ptr, depthcounter, maxval, threshold, cumsum, starts, stops, M);
        WBS_call(x, argmax, e, n, depth+1,changepoints, changepoint_counter_ptr, depthcounter, maxval, threshold, cumsum, starts, stops, M);
    }


}

SEXP NOT(SEXP nI, SEXP thresholdI, SEXP XI, SEXP MI, SEXP STARTS, SEXP STOPS){
    PROTECT(XI);
    int n = *(INTEGER(nI));
    int M = *(INTEGER(MI));
    double * X = REAL(XI);
    double threshold = *REAL(thresholdI);

    SEXP out1 = PROTECT(allocVector(INTSXP, n));
    int * changepoints = INTEGER(out1); //pointer to array
    memset(changepoints, 0, sizeof(int)*n);
    int changepoint_counter = 0;
    int* changepoint_counter_ptr = &changepoint_counter;
    //printf("size of int :%lu\n", sizeof(int));

    SEXP out4 = PROTECT(allocVector(REALSXP, n+1));
    double * cumsum = REAL(out4); //pointer to array
    SEXP out2 = PROTECT(allocVector(REALSXP, n));
    double * maxval = REAL(out2); //pointer to array
    SEXP out3= PROTECT(allocVector(INTSXP, n));
    memset(maxval, 0, sizeof(double)*n);
    int * depthcounter = INTEGER(out3); //pointer to array
    memset(depthcounter, 0, sizeof(int)*n);
    cumsum[0] = 0;
    for (int i = 0; i < n; ++i)
    {
        cumsum[i+1] = cumsum[i] + X[i];
        //printf("CUMSUM: %f\n", cumsum[i+1]);
    }

    protect(STARTS);
    int * starts = INTEGER(STARTS); //pointer to array

    protect(STOPS);
    int * stops = INTEGER(STOPS); //pointer to array



    // call WBS here:
    NOT_call(X, 0, n, n, 1, changepoints, changepoint_counter_ptr, depthcounter, maxval,threshold,cumsum, starts, stops, M);

    // creating a list to hold the two return values. this is the returned value.
    SEXP ret = PROTECT(allocVector(VECSXP, 3)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out1);
    SET_VECTOR_ELT(ret, 1, out2);
    SET_VECTOR_ELT(ret, 2, out3);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 3));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("changepoints"));
    SET_STRING_ELT(names, 1, mkChar("CUSUMval"));
    SET_STRING_ELT(names, 2, mkChar("depth"));

    /* assign names to list */
    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(9);
    return ret;



}

void BS_call(double * x, int s, int e, double n, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter, double * maxval, double threshold,
                double* cumsum ){

    if(e-s<2){
        return;
    }
    int argmax = s;
    double maximum = 0;

    compute_cusum(cumsum, s, e, &argmax, &maximum);


    //printf("largest value on interval [s,e] = [%d, %d] is %f\n", s,e,maximum);
    if (maximum >threshold)
    {
        changepoints[*changepoint_counter_ptr] = argmax;
        depthcounter[*changepoint_counter_ptr] = depth;
        maxval[*changepoint_counter_ptr] = maximum;
        (*changepoint_counter_ptr)++;
        BS_call(x, s, argmax, n, depth+1,changepoints, changepoint_counter_ptr, depthcounter, maxval, threshold, cumsum);
        BS_call(x, argmax, e, n, depth+1,changepoints, changepoint_counter_ptr, depthcounter, maxval, threshold, cumsum);
    }


}

SEXP BS(SEXP nI, SEXP thresholdI, SEXP XI){
    PROTECT(XI);
    int n = *(INTEGER(nI));
    double * X = REAL(XI);
    double threshold = *REAL(thresholdI);

    SEXP out1 = PROTECT(allocVector(INTSXP, n));
    int * changepoints = INTEGER(out1); //pointer to array
    memset(changepoints, 0, sizeof(int)*n);
    int changepoint_counter = 0;
    int* changepoint_counter_ptr = &changepoint_counter;
    //printf("size of int :%lu\n", sizeof(int));

    SEXP out4 = PROTECT(allocVector(REALSXP, n+1));
    double * cumsum = REAL(out4); //pointer to array
    SEXP out2 = PROTECT(allocVector(REALSXP, n));
    double * maxval = REAL(out2); //pointer to array
    SEXP out3= PROTECT(allocVector(INTSXP, n));
    memset(maxval, 0, sizeof(double)*n);
    int * depthcounter = INTEGER(out3); //pointer to array
    memset(depthcounter, 0, sizeof(int)*n);
    cumsum[0] = 0;
    for (int i = 0; i < n; ++i)
    {
        cumsum[i+1] = cumsum[i] + X[i];
        //printf("CUMSUM: %f\n", cumsum[i+1]);
    }


    // call WBS here:
    BS_call(X, 0, n, n, 1, changepoints, changepoint_counter_ptr, depthcounter, maxval,threshold,cumsum);

    // creating a list to hold the two return values. this is the returned value.
    SEXP ret = PROTECT(allocVector(VECSXP, 3)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out1);
    SET_VECTOR_ELT(ret, 1, out2);
    SET_VECTOR_ELT(ret, 2, out3);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 3));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("changepoints"));
    SET_STRING_ELT(names, 1, mkChar("CUSUMval"));
    SET_STRING_ELT(names, 2, mkChar("depth"));

    /* assign names to list */
    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(7);
    return ret;



}





SEXP CUSUM(SEXP n, SEXP p, SEXP A){
    // retrieving pointer for argument objects
    // here n is "int", p is "int" and A is flattened matrix
    // note that R flattens matrices to "by column".
    PROTECT(A);
    int nn = *(INTEGER(n));
    int pp = *(INTEGER(p));
    double * AA = REAL(A);

    // first return value, a vector
    // it is the flattened matrix of cusums
    SEXP out1 = PROTECT(allocVector(REALSXP, (nn) * (pp)));
    double * B = REAL(out1); //pointer to array

    // second return value, also vector
    //SEXP out2 = PROTECT(allocVector(REALSXP, (nn-1) ));
    //double * sumCUSUM= REAL(out2); //pointer to array

    // creating a list to hold the two return values. this is the returned value.
    SEXP ret = PROTECT(allocVector(VECSXP, 1)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out1);
    //SET_VECTOR_ELT(ret, 0, out2);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 1));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("CUSUM"));

    /* assign names to list */
    setAttrib(ret, R_NamesSymbol, names);


    //printf("HEI\n");
    int DEPTH = nn;
    //printf("n = %d\n", nn);
    //printf("p = %d\n", pp);

    // at first we let each col of B be a cumulative sum
    // as  Consistency results in multiple change-point problems (1993) page 10

    // first set top row of B to be the that of AA
    for (int j = 0; j < pp; ++j)
    {
        B[cord(0,j)] = AA[cord(0,j)];
    }

    // we don't want to keep the column-sumns in the resulting matrix, since
    // the cusum is not defined in the last element
    // so we keep this in a separate vector:
    SEXP temparr = PROTECT(allocVector(REALSXP, (pp) ));
    double * colsums= REAL(temparr); //pointer to array




   for (int j = 0; j < pp; ++j)
   {
        for (int i = 1; i < (nn-1); ++i)
        {
            B[cord(i,j)] = B[cord(i-1, j)] + AA[cord(i,j)];
        }
        colsums[j] = B[cord(nn-2, j)] + AA[cord(nn-1,j)];
   }



   double temp;

   for (int j = 0; j < pp; ++j)
   {
        for (int i = 0; i < (nn-1); ++i)
        {
            temp = ((double)(i+1) )/(nn) *colsums[j] -(B[cord(i,j)]);
            temp = temp *sqrt(nn) /sqrt(((double)(i+1))*(nn - i-1));
            B[cord(i,j)] = temp;
        }
   }

   UNPROTECT(5);
   return ret;

}

SEXP advancedsumCUSUM(SEXP n, SEXP p, SEXP A){
	// retrieving pointer for argument objects
	// here n is "int", p is "int" and A is flattened matrix
	// note that R flattens matrices to "by column".
	PROTECT(A);
	int nn = *(INTEGER(n));
	int pp = *(INTEGER(p));
	double * AA = REAL(A);

	// first return value, a vector
	// it is the flattened matrix of cusums
	SEXP out1 = PROTECT(allocVector(REALSXP, (nn) * (pp)));
	double * B = REAL(out1); //pointer to array

	// second return value, also vector
	SEXP out2 = PROTECT(allocVector(REALSXP, (nn-1) ));
	double * sumCUSUM= REAL(out2); //pointer to array

	// creating a list to hold the two return values. this is the returned value.
	SEXP ret = PROTECT(allocVector(VECSXP, 1)); //the list to be returned in R
  	//SET_VECTOR_ELT(ret, 0, out1);
  	SET_VECTOR_ELT(ret, 0, out2);

  	// creating list of names/titles to be returned in the output list
  	SEXP names = PROTECT(allocVector(STRSXP, 1));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("sumCUSUM"));

    /* assign names to list */
    setAttrib(ret, R_NamesSymbol, names);


    //printf("HEI\n");
    int DEPTH = nn;
    //printf("n = %d\n", nn);
    //printf("p = %d\n", pp);

    // at first we let each col of B be a cumulative sum
    // as  Consistency results in multiple change-point problems (1993) page 10

    // first set top row of B to be the that of AA
    for (int j = 0; j < pp; ++j)
    {
    	B[cord(0,j)] = AA[cord(0,j)];
    }

    // we don't want to keep the column-sumns in the resulting matrix, since
   	// the cusum is not defined in the last element
   	// so we keep this in a separate vector:
    SEXP temparr = PROTECT(allocVector(REALSXP, (pp) ));
    double * colsums= REAL(temparr); //pointer to array




   for (int j = 0; j < pp; ++j)
   {
   		for (int i = 1; i < (nn-1); ++i)
   		{
   			B[cord(i,j)] = B[cord(i-1, j)] + AA[cord(i,j)];
   		}
   		colsums[j] = B[cord(nn-2, j)] + AA[cord(nn-1,j)];
   }



   double temp;

   for (int j = 0; j < pp; ++j)
   {
   		for (int i = 0; i < (nn-1); ++i)
   		{
   			temp = ((double)(i+1) )/(nn) *colsums[j] -(B[cord(i,j)]);
   			temp = temp *sqrt(nn) /sqrt(((double)(i+1))*(nn - i-1));
   			B[cord(i,j)] = temp;
   		}
   }

   for (int i = 0; i < (nn -1); ++i)
   {
   	sumCUSUM[i] = 0;
   	for (int j = 0; j < (pp); ++j)
   	{
   		sumCUSUM[i] += fabs(B[cord(i,j)]);
   	}
   }

   UNPROTECT(6);
   return ret;

}

SEXP advancedsumCUSUMsquared(SEXP n, SEXP p, SEXP A){
    // retrieving pointer for argument objects
    // here n is "int", p is "int" and A is flattened matrix
    // note that R flattens matrices to "by column".
    PROTECT(A);
    int nn = *(INTEGER(n));
    int pp = *(INTEGER(p));
    double * AA = REAL(A);

    // first return value, a vector
    // it is the flattened matrix of cusums
    SEXP out1 = PROTECT(allocVector(REALSXP, (nn) * (pp)));
    double * B = REAL(out1); //pointer to array

    // second return value, also vector
    SEXP out2 = PROTECT(allocVector(REALSXP, (nn-1) ));
    double * sumCUSUM= REAL(out2); //pointer to array

    // creating a list to hold the two return values. this is the returned value.
    SEXP ret = PROTECT(allocVector(VECSXP, 1)); //the list to be returned in R
    //SET_VECTOR_ELT(ret, 0, out1);
    SET_VECTOR_ELT(ret, 0, out2);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 1));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("sumCUSUM"));

    /* assign names to list */
    setAttrib(ret, R_NamesSymbol, names);


    //printf("HEI\n");
    int DEPTH = nn;
    //printf("n = %d\n", nn);
    //printf("p = %d\n", pp);

    // at first we let each col of B be a cumulative sum
    // as  Consistency results in multiple change-point problems (1993) page 10

    // first set top row of B to be the that of AA
    for (int j = 0; j < pp; ++j)
    {
        B[cord(0,j)] = AA[cord(0,j)];
    }

    // we don't want to keep the column-sumns in the resulting matrix, since
    // the cusum is not defined in the last element
    // so we keep this in a separate vector:
    SEXP temparr = PROTECT(allocVector(REALSXP, (pp) ));
    double * colsums= REAL(temparr); //pointer to array




   for (int j = 0; j < pp; ++j)
   {
        for (int i = 1; i < (nn-1); ++i)
        {
            B[cord(i,j)] = B[cord(i-1, j)] + AA[cord(i,j)];
        }
        colsums[j] = B[cord(nn-2, j)] + AA[cord(nn-1,j)];
   }



   double temp;

   for (int j = 0; j < pp; ++j)
   {
        for (int i = 0; i < (nn-1); ++i)
        {
            temp = ((double)(i+1) )/(nn) *colsums[j] -(B[cord(i,j)]);
            temp = temp *sqrt(nn) /sqrt(((double)(i+1))*(nn - i-1));
            B[cord(i,j)] = temp;
        }
   }

   for (int i = 0; i < (nn -1); ++i)
   {
    sumCUSUM[i] = 0;
    for (int j = 0; j < (pp); ++j)
    {
        sumCUSUM[i] += (B[cord(i,j)])*(B[cord(i,j)]);
    }
   }

   UNPROTECT(6);
   return ret;

}


double internalsumCUSUM(int nn, int pp, double* AA, double *B, double *sumCUSUM, double * colsums){
    // retrieving pointer for argument objects
    // here n is "int", p is "int" and A is flattened matrix
    // note that R flattens matrices to "by column".


    // first return value, a vector
    // it is the flattened matrix of cusums
    /*SEXP out1 = PROTECT(allocVector(REALSXP, (nn) * (pp)));
    double * B = REAL(out1); //pointer to array

    // second return value, also vector
    SEXP out2 = PROTECT(allocVector(REALSXP, (nn-1) ));
    double * sumCUSUM= REAL(out2); //pointer to array*/

/*    // creating a list to hold the two return values. this is the returned value.
    SEXP ret = PROTECT(allocVector(VECSXP, 1)); //the list to be returned in R
    //SET_VECTOR_ELT(ret, 0, out1);
    SET_VECTOR_ELT(ret, 0, out2);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 1));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("sumCUSUM"));

    setAttrib(ret, R_NamesSymbol, names);*/


    //printf("HEI\n");
    int DEPTH = nn;
    //printf("n = %d\n", nn);
    //printf("p = %d\n", pp);

    // at first we let each col of B be a cumulative sum
    // as  Consistency results in multiple change-point problems (1993) page 10

    // first set top row of B to be the that of AA
    for (int j = 0; j < pp; ++j)
    {
        B[cord(0,j)] = AA[cord(0,j)];
    }

    // we don't want to keep the column-sumns in the resulting matrix, since
    // the cusum is not defined in the last element
    // so we keep this in a separate vector:
    /*SEXP temparr = PROTECT(allocVector(REALSXP, (pp) ));
    double * colsums= REAL(temparr); //pointer to array*/




   for (int j = 0; j < pp; ++j)
   {
        for (int i = 1; i < (nn-1); ++i)
        {
            B[cord(i,j)] = B[cord(i-1, j)] + AA[cord(i,j)];
        }
        colsums[j] = B[cord(nn-2, j)] + AA[cord(nn-1,j)];
   }



   double temp;

   for (int j = 0; j < pp; ++j)
   {
        for (int i = 0; i < (nn-1); ++i)
        {
            temp = ((double)(i+1) )/(nn) *colsums[j] -(B[cord(i,j)]);
            temp = temp *sqrt(nn) /sqrt(((double)(i+1))*(nn - i-1));
            B[cord(i,j)] = temp;
        }
   }

   for (int i = 0; i < (nn -1); ++i)
   {
    sumCUSUM[i] = 0;
    for (int j = 0; j < (pp); ++j)
    {
        sumCUSUM[i] += fabs(B[cord(i,j)]);
    }
   }

   double max = 0.0;
   double tempp = 0.0;
   for (int i = 0; i < (nn -1); ++i){
    tempp = fabs(sumCUSUM[i]);
    if (tempp>max)
    {
        max = tempp;
    }
   }
   //UNPROTECT(6);
   return max;

}

SEXP MCsumCUSUM(SEXP n, SEXP p, SEXP N,SEXP sd_ests) {
    // retrieving pointer for argument objects
    // here n is "int", p is "int" and A is flattened matrix
    // note that R flattens matrices to "by column".

    int nn = *(INTEGER(n));
    int pp = *(INTEGER(p));
    int NN = *(INTEGER(N));
    double * sds = REAL(sd_ests);

    // first return value, a vector
    // it is the flattened matrix of cusums
    SEXP out1 = PROTECT(allocVector(REALSXP, NN));
    double * maxes = REAL(out1); //pointer to array

    /*// second return value, also vector
    SEXP out2 = PROTECT(allocVector(REALSXP, (nn-1) ));
    double * sumCUSUM= REAL(out2); //pointer to array*/

    // creating a list to hold the two return values. this is the returned value.
    SEXP ret = PROTECT(allocVector(VECSXP, 1)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out1);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 1));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("sumCUSUM"));

    /* assign names to list */
    setAttrib(ret, R_NamesSymbol, names);


    //Mint DEPTH = nn;

    SEXP temp = PROTECT(allocVector(REALSXP, (nn) * (pp)));
    double * X = REAL(temp); //pointer to array

    SEXP temp2 = PROTECT(allocVector(REALSXP, (nn) * (pp)));
    double * B = REAL(temp2);

    SEXP temp3 = PROTECT(allocVector(REALSXP, (nn-1) ));
    double * sumCUSUM= REAL(temp3);

    SEXP temp4 = PROTECT(allocVector(REALSXP, (pp) ));
    double * colsums= REAL(temp4); //pointer to array



    GetRNGstate();

    for (int k = 0; k < NN; ++k)
    {
        // generate matrix
        for (int j = 0; j < pp; ++j){
            for (int i = 0; i < (nn); ++i)
            {
                X[cord_spec(i,j,nn)] = norm_rand()*sds[j];
            }
        }
        maxes[k] = internalsumCUSUM(nn,  pp, X, B, sumCUSUM, colsums);



    }
    PutRNGstate();
    UNPROTECT(7);
    return ret;

}





SEXP advancedmaxCUSUM(SEXP n, SEXP p, SEXP A) {
    // retrieving pointer for argument objects
    // here n is "int", p is "int" and A is flattened matrix
    // note that R flattens matrices to "by column".
    PROTECT(A);
    int nn = *(INTEGER(n));
    int pp = *(INTEGER(p));
    double * AA = REAL(A);

    // it is the flattened matrix of cusums
    SEXP out1 = PROTECT(allocVector(REALSXP, (nn) * (pp)));
    double * B = REAL(out1); //pointer to array

    // the return value, also vector
    SEXP out2 = PROTECT(allocVector(REALSXP, 1 ));
    double * maxCUSUM= REAL(out2); //pointer to array

    // creating a list to hold the two return values. this is the returned value.
    SEXP ret = PROTECT(allocVector(VECSXP, 1)); //the list to be returned in R
    //SET_VECTOR_ELT(ret, 0, out1);
    SET_VECTOR_ELT(ret, 0, out2);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(names, 0, mkChar("maxCUSUM"));
    //SET_STRING_ELT(names, 1, mkChar("sumCUSUM"));

    /* assign names to list */
    setAttrib(ret, R_NamesSymbol, names);


    //printf("HEI\n");
    int DEPTH = nn;
    //printf("n = %d\n", nn);
    //printf("p = %d\n", pp);

    // at first we let each col of B be a cumulative sum
    // as  Consistency results in multiple change-point problems (1993) page 10

    // first set top row of B to be the that of AA
    for (int j = 0; j < pp; ++j)
    {
        B[cord(0,j)] = AA[cord(0,j)];
    }

    // we don't want to keep the column-sumns in the resulting matrix, since
    // the cusum is not defined in the last element
    // so we keep this in a separate vector:
    SEXP temparr = PROTECT(allocVector(REALSXP, (pp) ));
    double * colsums= REAL(temparr); //pointer to array




   for (int j = 0; j < pp; ++j)
   {
        for (int i = 1; i < (nn-1); ++i)
        {
            B[cord(i,j)] = B[cord(i-1, j)] + AA[cord(i,j)];
        }
        colsums[j] = B[cord(nn-2, j)] + AA[cord(nn-1,j)];
   }



   double temp;

   for (int j = 0; j < pp; ++j)
   {
        for (int i = 0; i < (nn-1); ++i)
        {
            temp = ((double)(i+1) )/(nn) *colsums[j] -(B[cord(i,j)]);
            temp = temp *sqrt(nn) /sqrt(((double)(i+1))*(nn - i-1));
            B[cord(i,j)] = temp;
        }
   }
   temp = 0;
   for (int i = 0; i < (nn -1); ++i)
   {
    //sumCUSUM[i] = 0;
    for (int j = 0; j < (pp); ++j)
    {
        //sumCUSUM[i] += B[cord(i,j)];
        if(B[cord(i,j)]>temp){
            temp = B[cord(i,j)];
        }
    }
   }
   maxCUSUM[0] = temp;
   UNPROTECT(6);
   return ret;
}

double internalmaxCUSUM(int nn, int pp, double* AA, double * B, double * colsums) {
    // finds the MAXIMUM CUSUM value
    // AA, B and colsums are pre-allocated vectors
    // AA : nn \times pp
    // B : nn \times pp
    // AA is the data series. top row of A is one observed X_t in RR^p
    // returns the maximum cusum value (in abs)


    // DEPTH is necessary for matrix-indexing.
    // matrices are indexed column for column
    int DEPTH = nn;

    for (int j = 0; j < pp; ++j)
    {
        B[cord(0,j)] = AA[cord(0,j)];
    }





   for (int j = 0; j < pp; ++j)
   {
        for (int i = 1; i < (nn-1); ++i)
        {
            B[cord(i,j)] = B[cord(i-1, j)] + AA[cord(i,j)];
        }
        colsums[j] = B[cord(nn-2, j)] + AA[cord(nn-1,j)];
   }



   double temp;

   for (int j = 0; j < pp; ++j)
   {
        for (int i = 0; i < (nn-1); ++i)
        {
            temp = ((double)(i+1) )/(nn) *colsums[j] -(B[cord(i,j)]);
            temp = temp *sqrt(nn) /sqrt(((double)(i+1))*(nn - i-1));
            B[cord(i,j)] = temp;
        }
   }
   temp = 0;

    for (int j = 0; j < (pp); ++j)
   {
    //sumCUSUM[i] = 0;
    for (int i = 0; i < (nn -1); ++i)
    {
        //sumCUSUM[i] += B[cord(i,j)];
        if(B[cord(i,j)]>temp){
            temp = B[cord(i,j)];
        }
    }
   }
   return temp;
}




SEXP CUSUM_simulator(SEXP n, SEXP p, SEXP k) {

    int nn = *(INTEGER(n));
    int pp = *(INTEGER(p));
    int kk = *(INTEGER(k));

    printf("n: %d\n", nn);
    printf("p: %d\n", pp);
    printf("k: %d\n", kk);
    // it is the flattened matrix of max cusums
    SEXP out1 = PROTECT(allocVector(REALSXP, kk*(nn-1)));
    double * outmatrix = REAL(out1); //pointer to array


    SEXP internal1 = PROTECT(allocVector(REALSXP, pp*(nn)));
    double * AA = REAL(internal1); //pointer to array

    SEXP internal2 = PROTECT(allocVector(REALSXP, pp*(nn)));
    double * B = REAL(internal2); //pointer to array

    SEXP internal3 = PROTECT(allocVector(REALSXP, pp));
    double * colsums = REAL(internal3); //pointer to array




    SEXP ret = PROTECT(allocVector(VECSXP, 1)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out1);

    SEXP names = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(names, 0, mkChar("maxCUSUMs"));

    setAttrib(ret, R_NamesSymbol, names);


    //printf("HEI\n");
    int DEPTH = nn-1;

    GetRNGstate();
    int nnn;

    for (int i = 0; i < (nn-1); ++i)
    {
        //n = i+2
        nnn = i+2;
        for (int j = 0; j < kk; ++j)
        {
            // simulate series
            for (int jj = 0; jj < pp; ++jj)
            {
                for (int ii = 0; ii < nnn; ++ii)
                {
                    AA[cord_spec(ii,jj,nnn)] = norm_rand();
                }
            }
            outmatrix[cord(i,j)] = internalmaxCUSUM(nnn,pp, AA, B, colsums);
        }

    }



    PutRNGstate();
   UNPROTECT(6);
   return ret;
}








SEXP test2(SEXP a, SEXP b){
	//PROTECT(a = AS_NUMERIC(a));
	//PROTECT(a = AS_NUMERIC(a));
	int * ainput = INTEGER(a);
	//int k = asInteger(a);
	printf("%d",*ainput);
	double * bvec = REAL(b);
	SEXP out = PROTECT(allocVector(REALSXP, *ainput));
	double * pout = REAL(out);
	for (int i = 0; i < *ainput; i++) {
    	pout[i] = bvec[i];
  	}
  	UNPROTECT(1);

  return out;


/*	SEXP my_array = PROTECT(allocVector(REALSXP, 4));


	UNPROTECT(2);

	SEXP out = PROTECT(allocVector(REALSXP, n));

	px = REAL(x);
	pout = REAL(out);
	for (int i = 0; i < n; i++) {
	pout[i] = px[i] + 2;
	}
	UNPROTECT(1);

	return out;*/


}
