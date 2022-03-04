#include "header.h"

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

void CUSUM(double * cumsums, double * cusum, int s, int e, int p){
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
            cusum[cord_spec(i,j-1, p)] = sqrt(((double)(e-t))/((e-s)*(t-s))) *(cumsums[cord_spec(i,t+1,p)]-cumsums[cord_spec(i,s+1,p)])
            - sqrt(((double)(t-s))/((e-s)*(e-t)))*(cumsums[cord_spec(i,e+1,p)] - cumsums[cord_spec(i,t+1,p)]);

        }

    }

    return;
}

SEXP CUSUM_R(SEXP XI, SEXP sI, SEXP eI, SEXP pI, SEXP nI){
	PROTECT(XI);
	PROTECT(sI);
	PROTECT(eI);
	PROTECT(pI);
	PROTECT(nI);
	int s = *(INTEGER(sI));
	int p = *(INTEGER(pI));
	int e = *(INTEGER(eI));
	int n = *(INTEGER(nI));
	UNPROTECT(4);
	double * X = REAL(XI);

	SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*sizeof(double));

    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p * (n-1)));
    double *cusum = REAL(cusumSEXP); // p \times (n+1). first col is 0
    memset(cusum, 0, p*sizeof(double)*(n-1));

    for (int j = 1; j <=n; ++j)
    {
        for (int i = 0; i < p; ++i)
        {
            //#define cord_spec(r,c, D) ((r) + (D)*(c))
            cumsums[cord_spec(i,j,p)] = X[cord_spec(i,j-1, p)] +cumsums[cord_spec(i,j-1, p)];
        }
    }
    CUSUM(cumsums, cusum, s, e, p);
    UNPROTECT(3);
    return(cusumSEXP);
}

