/*#include "header.h"

void internal_threshold_matrix(double * matrix, int r1, int c1, double a, double nu_a, int previously_tresholded, 
                                double prev_nu_a ){

    // should be possible to optimize this code a bit
    double a_sq = a*a;
    double true_val = 0;
    if (previously_tresholded)
    {
        for (int i = 0; i < r1; ++i)
        {
            for (int j = 0; j < c1; ++j)
            {
                if (fabs(matrix[cord_spec(i, j, r1)])>1e-10)
                {
                    true_val = matrix[cord_spec(i, j, r1)]+prev_nu_a;

                    if (true_val>a_sq)
                    {
                        matrix[cord_spec(i, j, r1)] = true_val - nu_a;
                    }
                    else{
                        matrix[cord_spec(i, j, r1)] = 0.0;
                    }
                }

            }

        }
    }
    else{
        for (int i = 0; i < r1; ++i)
        {
            for (int j = 0; j < c1; ++j)
            {
                if (fabs(matrix[cord_spec(i, j, r1)])>a)
                {
                    matrix[cord_spec(i, j, r1)] = matrix[cord_spec(i, j, r1)]*matrix[cord_spec(i, j, r1)] - nu_a;
                    
                }
                else{
                    matrix[cord_spec(i, j, r1)] = 0.0;
                }

            }

        }
    }
    

}

void internal_colSum(double * matrix, int r1, int c1, double * vector){
    memset(vector, 0, c1);
    for (int i = 0; i < r1; ++i)
        {
            for (int j = 0; j < c1; ++j)
            {
                vector[j]+=matrix[cord_spec(i, j, r1)];

            }

        }
}
void internal_check_segment(double * cumsums, double * cusum, int * maxpos, double * maximum, int * maxa_pos,
                            int s, int e, int p, double * vector, double threshold_d, double threshold_s, int debug){
    // require as[0] = 0
    
    *maxpos  = e;
    *maximum = 0.0;
    *maxa_pos = 0;
    if(e-s<2){
        return;
    }

    // compute CUSUM
    CUSUM(cumsums, cusum, s, e, p); // dim is p \times (e-s-1)

    double prev_nu_a = 0.0;
    double threshold = threshold_d;
    //int n = e-s;
    for (int i = 0; i < len_as; ++i)
    {
        if (i==1)
        {
            threshold = threshold_s;
        }
        a = as[i];
        nu_a = nu_as[i];
        internal_threshold_matrix(cusum, p, e-s-1, a, nu_a, i>0, prev_nu_a );
        internal_colSum(cusum, p, e-s-1, vector);
        for (int i = 0; i < e-s-1; ++i)
        {
            //t = s+i+1;
            tmp = vector[i] / threshold;
            if(tmp> *maximum){
                *maximum = tmp;
                *maxpos = s+i+1;
                *maxa_pos = i;
            }
        }


    }


    return;
}


void cHDCD_call(double * x, int s, int e, int n, int p, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter,
                double * maxval, double threshold_d, double threshold_s , double *cumsums, int* lens, int lenLens, double * as, double * nu_as, int len_as,
                int * segstarts, double * maxvalues, int* maxpos, int K, double * cusum, double * matrix,
                double * vector,  int debug,int * coordchg){
    if(debug){
        printf("cHDCD_call! s=%d, e=%d\n", s, e);
    }

    if(e-s<lens[0]){
        //printf("segment too short\n");
        return;
    }
    int argmax = s;
    double maximum = 0;
    double maxa_pos = 0;

    //int tmpargmax;
    //double tmpmaximum;
    double tmp;
    int len;
    int jump;
    //int found=0;

    int i;

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

        if(e-s<len){
            break;
        }

        for (int k = 0; k < n; ++k)
        {
            i = segstarts[cord_spec(k,j,n)];

            if(i>e-len || i<-1){
                break;
            }
            else if (i<s)
            {
                continue;
            }

            if(debug){
                //printf("maxvalues[%d, %d] = %f\n", k , j , maxvalues[cord_spec(k,j,n)]);
            }

            if(maxvalues[cord_spec(k,j,n)]<=0.0){
                //this segment not computed!
                //inspectOnSegment(cumsums, cusum, &tmpargmax, &tmpmaximum, s, e, p, lambda,
                 //   eps, maxiter, projvec, cusum_proj);
                 //double * cumsums, double * cusum, int * maxpos, double * maximum,
                 internal_check_segment(cumsums, cusum, &(maxpos[cord_spec(k,j,n)]), &(maxvalues[cord_spec(k,j,n)]), &(maxas[cord_spec(k,j,n)]),
                            i, i+len, p, vector, threshold_d,threshold_s, debug);
                 //internal_inspectOnSegment(cumsums, cusum, &(maxpos[cord_spec(k,j,n)]), &(maxvalues[cord_spec(k,j,n)]), i, i+len, p,
                 //lambda,
                //      eps, maxiter, mhat, mhatprod, v, v2,debug);
            }

            tmp = maxvalues[cord_spec(k,j,n)];
            if(tmp>maximum){
                maximum = tmp;
                argmax = maxpos[cord_spec(k,j,n)];
                maxa_pos = maxas[cord_spec(k,j,n)];
                j_max = j;
                k_max=k;
                //found=1;
            }


        }

        if(maximum>1.0){
            break;
        }
    }
    if(debug){
        printf("maximum=%f\n", maximum);
    }

    if(maximum > 1.0){
        if(debug){
            printf("!!!!!! declared change-point in %d. val = %f", argmax, maximum);
        }
        // identify in which coordinates the change happens:
        i = segstarts[cord_spec(k_max,j_max,n)];
        len = lens[j_max];
        int ss = i;
        int ee = i+len;
        CUSUM(cumsums, cusum, , ee, p);
        internal_threshold_matrix(double * matrix, int r1, int c1, double a, double nu_a, int previously_tresholded, 
                                double prev_nu_a )

        double * projvec = internal_sparse_svd(cusum, p, ee-ss-1, lambda, eps, maxiter,
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
        //cInspect_call(x, s, argmax, n, p, depth+1, changepoints, changepoint_counter_ptr, depthcounter,  maxval, threshold,
        //        adaptTresh, cumsums, lens, lenLens, K, cusum,  mhat, mhatprod, v, v2);
        cInspect_call(x, s, argmax, n, p, depth+1, changepoints, changepoint_counter_ptr, depthcounter,
                maxval,xi, cumsums, lens, lenLens, lambda,
                eps, maxiter, segstarts, maxvalues, maxpos, K,  cusum, mhat,
                mhatprod, v, v2,debug,coordchg);
        //cInspect_call(x, argmax, e, n, p, depth+1, changepoints, changepoint_counter_ptr, depthcounter,  maxval, threshold,
        //        adaptTresh, cumsums, lens, lenLens, K, cusum, mhat, mhatprod, v, v2);
        cInspect_call(x, argmax, e, n, p, depth+1, changepoints, changepoint_counter_ptr, depthcounter,
                maxval,xi, cumsums, lens, lenLens, lambda,
                eps, maxiter, segstarts, maxvalues, maxpos, K,  cusum, mhat,
                mhatprod, v, v2,debug,coordchg);

    }

    return;
}

SEXP cHDCD(SEXP XI,SEXP nI, SEXP pI,SEXP xiI, SEXP lensI,SEXP lenLensI,SEXP KI,
    SEXP epsI, SEXP lambdaI, SEXP maxiterI, SEXP debugI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(lensI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(xiI);
    PROTECT(lenLensI);
    PROTECT(KI);
    PROTECT(epsI);
    PROTECT(lambdaI);
    PROTECT(maxiterI);
    PROTECT(debugI);

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double xi = *(REAL(xiI));
    int *lens = INTEGER(lensI); /////// dobbeltsjekk at int er like stor som INTEGER!!!!!!!
    int lenLens = *(INTEGER(lenLensI));
    int K = *(INTEGER(KI));
    double eps = *REAL(epsI);
    double lambda = *REAL(lambdaI);
    int maxiter = *INTEGER(maxiterI);
    int debug = *INTEGER(debugI);
    UNPROTECT(10); // unprotecting all except X and lens
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


    SEXP maxvaluesSEXP = PROTECT(allocVector(REALSXP, n*lenLens));
    double * maxvalues = REAL(maxvaluesSEXP);
    memset(maxvalues, 0, sizeof(double)*n*lenLens);
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

        len = lens[j]; //len here is -1 of the len in text
        jump = len/K;
        if(jump<1){
            jump=1;
        }

        for (int i = -1; i <= (n-len); i+=jump)
        {
            //cord_spec(r,c, D) ((r) + (D)*(c))
            //compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxvalues[cord_spec(i,j,n)]));
            segstarts[cord_spec(counter++,j,n)] = i;
            if(debug){
                printf("segstarts[%d, %d] = %d\n",counter-1, j, i );
            }


        }
    }



    cHDCD_call(X, 0, n, n, p, 1, changepoints, changepoint_counter_ptr, depthcounter,
                maxval, xi, cumsums, lens, lenLens, lambda,
                eps, maxiter, segstarts, maxvalues, maxpos, K, cusum, mhat,
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


*/