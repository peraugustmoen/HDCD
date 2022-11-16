#include "header.h"


int internal_check_segment_Pilliat(double * cumsums, double * cusum, int k, int len, int p, double * thresholds_partial, 
    int * thresholds_bj, double threshold_dense, int maxx,  int * detected,
    int * vec, double * vec2, int debug){
    // require as[0] = 0
    //if(debug){
     // Rprintf("checking segment (%d,%d]\n", s,e);

    //    }

    *detected = 0;

    if(len<1){
        return 0;
    }

    // compute CUSUM
    singleCUSUM(cumsums, cusum, k-len, k+len, p, k); 

    int z = 0;

    double cumsum = 0.0;
    int prev = -1;
    // double x; 
    // int N=0;
    int j = 0;

    //for (int j = 0; j < e-s-1; ++j){

    // first check dense:

    cumsum = -p;
    for (int i = 0; i < p; ++i)
    {
        cumsum+= cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
        vec2[i] = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
    }
    if(cumsum > threshold_dense){
        *detected = 1;
        return 0;
    }

    // then partial sum:
    R_qsort(vec2,1,p);

    cumsum = 0.0;
    prev = p;
    z = 1;
    int c = 0;
    int i = 0;
    while(1)
    {
        if(z>p){
            break;
        }
        for (i = prev-1; i >= p-z; --i)
        {
            cumsum +=vec2[i];
        }

        if (cumsum> thresholds_partial[c])
        {
            *detected = 1;
            return 1;
        }
        prev = p-z;
        c++;
        z = 2*z;


    }
    // then berk jones:
    //for (int h = 0; h <= maxx; ++h)
    //{
    int tmp = 0;
    memset(vec, 0, sizeof(int)*maxx );

        for (int i = 0; i < p; ++i)
        {
            tmp = (int)  fabs(cusum[cord_spec(i,j, p)]);
            //if(debug){
            //    Rprintf("val = %f, res = %d\n",fabs(cusum[cord_spec(i,j, p)]), tmp );
            //}
            if(tmp == 0){
                break;
            }else{

                for (int h = 0; h < tmp && h<maxx; ++h)
                {
                    vec[h]++;
                }
                
            }             
        }

        for (int h = 0; h < maxx; ++h)
        {
            if(vec[h] > thresholds_bj[h]){
                *detected = 1;
                if(debug){
                    Rprintf("Berk Jones detected at x = %d, seg [%d,%d). Count = %d, thresh = %d.\n", 
                        h, k-len, k+len, vec[h], thresholds_bj[h]);
                }
                return 2;
            }
        }
 
    //}
    return 0;

    //}


}


void cPilliat_call(double * x, int n, int p, int* changepoints,
                double * thresholds_partial, double threshold_dense, int * thresholds_bj, double *cumsums, int* lens, int lenLens,
                int * vec, double * vec2, int K, double * cusum, int* changepoint_counter, int * teststat,
                int * startpoints, int* endpoints, int maxx, int debug){
    if(debug){
        Rprintf("cPilliat_call! n= %d, p = %d", n, p);
    }


    //int intersects = 0;
    //int max=0;
    //int min=0;
    int detected = 0;
    //int stop = 0;

    //int l = 0;
    int jj=0;
    int jump = 0;
    int prevdetected = 0;
    int len = 0;
    int stop = 0;
    int k = 0;
    int res = 0;
    //int intersects = 0;
    int t1 = 0;
    int t2 = 0;
    int check = 1;
    for (int j = 0; j < lenLens; ++j)
    {
        len = lens[j];
        if(debug){
            Rprintf("j=%d, len = %d\n", j, len);
        }

        jump = len ;


        


        //int l = len -1;
        prevdetected=0;

        if(len ==1){
            for (k = 0; k < n-1; ++k)
            {
                // check segments
                if(debug){
                    Rprintf("checking k = %d, l = %d, interval [%d, %d)\n", k, len, k-1, k+1);
                }
                res = internal_check_segment_Pilliat(cumsums, cusum, k, len,  p, thresholds_partial, 
                thresholds_bj, threshold_dense, maxx, &detected,vec, vec2, debug);

                if(detected){
                    if(debug){
                        Rprintf("found new changepoint in [%d, %d). chgpt set to %d\n", k-1, k+1,k );
                    }
                    changepoints[*(changepoint_counter)] = k;
                    startpoints[(*(changepoint_counter))] = k-1;
                    endpoints[(*(changepoint_counter))] = k+1;
                    teststat[*(changepoint_counter)] = res;
                    (*(changepoint_counter))++;
                }
            }
        }
        else{
            stop=0;
            k = len-1;
            jj = 2;
            while(1){
                if(k + len>=n || k-len<-1){
                    break;
                }
                if(debug){
                    Rprintf("checking k = %d, l = %d, interval [%d, %d)\n", k, len, k-len, k+len);
                }
                for (int hh = 0; hh < *changepoint_counter; ++hh)
                {
                    check = 1;
                    detected = 0;
                    t1 = ( k- len >= startpoints[hh]) ? k-len : startpoints[hh];
                    t2 = (k+len <= endpoints[hh]) ? k+len : endpoints[hh];
                    if(t1 <= t2 - 2){
                        if(debug){
                            Rprintf("The interval [%d, %d) overlaps with [%d, %d) - ignored\n",
                                k - len, k+len, startpoints[hh],endpoints[hh]);
                        }

                        check = 0;
                        break;
                    }
                }
                if(check){
                    res = internal_check_segment_Pilliat(cumsums, cusum, k, len,  p, thresholds_partial, 
                            thresholds_bj, threshold_dense,  maxx, &detected,vec, vec2, debug);
                }
                

                if(detected && prevdetected){
                    changepoints[*(changepoint_counter)] = (k + startpoints[(*(changepoint_counter))])/2;
                    //startpoints[(*(changepoint_counter))] = k-1;
                    endpoints[(*(changepoint_counter))] = k+1;
                    prevdetected=1;
                    if(debug){
                        Rprintf("also found in [%d, %d). chgpt changed to %d\n", k-len, k+len,(
                            k + startpoints[(*(changepoint_counter))])/2  );
                    }
                    //(*(changepoint_counter))++;
                }
                else if(detected && !(prevdetected)){
                    //res = //check segment
                    if(debug){
                        Rprintf("found new changepoint in [%d, %d). chgpt set to %d\n", k-len, k+len,k );
                    }
                    changepoints[*(changepoint_counter)] = k;
                    startpoints[(*(changepoint_counter))] = k-len;
                    endpoints[(*(changepoint_counter))] = k+len;
                    prevdetected=1;
                    teststat[*(changepoint_counter)] = res;
                }
                else if(!detected && prevdetected){
                    (*(changepoint_counter))++;
                    prevdetected = 0;
                }

                if(stop){
                    if(detected && prevdetected){
                        (*(changepoint_counter))++;
                        prevdetected = 0;
                    }
                    break;
                }
                k = (jj++)*jump;
                if(k>= n - len-1){
                    k = n-len-1;

                    stop=1;
                }
            }
        }
    }


    return;
}

SEXP cPilliat(SEXP XI,SEXP nI, SEXP pI,SEXP thresholds_partialI, SEXP threshold_denseI, SEXP thresholds_bjI,
    SEXP lensI,SEXP lenLensI,SEXP KI, SEXP maxxI, SEXP rescale_variance_boolI, SEXP debugI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(thresholds_partialI);
    PROTECT(threshold_denseI);
    PROTECT(thresholds_bjI);
    PROTECT(lensI);
    PROTECT(lenLensI);
    //PROTECT(log2pI);
    //PROTECT(asI);
    //PROTECT(nu_asI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(KI);
    //PROTECT(len_asI);
    //PROTECT(twolognI);
    PROTECT(maxxI);
    PROTECT(debugI);
    PROTECT(rescale_variance_boolI);
    

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double *thresholds_partial = (REAL(thresholds_partialI));
    double threshold_dense = *(REAL(threshold_denseI));
    int *thresholds_bj = (INTEGER(thresholds_bjI));
    //int log2p = *(INTEGER(threshold_denseI));
    int *lens = INTEGER(lensI); /////// dobbeltsjekk at int er like stor som INTEGER!!!!!!!
    int lenLens = *(INTEGER(lenLensI));
    int maxx = *(INTEGER(maxxI));
    //int K = *(INTEGER(KI));
    //double * as = REAL(asI);
    //double * nu_as = REAL(nu_asI);
    //int len_as = *(INTEGER(len_asI));
    //int twologn = * (INTEGER(twolognI));
    int debug = *INTEGER(debugI);
    int rescale_variance_bool = *INTEGER(rescale_variance_boolI);
    int K = *INTEGER(KI);

    if(debug){
      Rprintf("p = %d\n", p);
      Rprintf("n = %d\n", n);
    }
    SEXP out = PROTECT(allocVector(INTSXP, n));
    int * changepoints = INTEGER(out); //pointer to array
    //memset(changepoints, -1, sizeof(int)*n);
    for (size_t i = 0; i < n; i++) {
      changepoints[i] = -1;
    }
    SEXP teststatSEXP = PROTECT(allocVector(INTSXP, n));
    int * teststat = INTEGER(teststatSEXP); //pointer to array
    //memset(changepoints, -1, sizeof(int)*n);
    for (size_t i = 0; i < n; i++) {
      teststat[i] = -1;
    }

    int changepoint_counter = 0;
    int* changepoint_counter_ptr = &changepoint_counter;
    SEXP startpointsSEXP = PROTECT(allocVector(INTSXP, n));
    int * startpoints = INTEGER(startpointsSEXP); //pointer to array
    //memset(startpoints, 0, sizeof(int)*n);
    SEXP endpointsSEXP = PROTECT(allocVector(INTSXP, n));
    int * endpoints = INTEGER(endpointsSEXP); //pointer to array
    //memset(endpoints, 0, sizeof(int)*n);
    for (int i = 0; i < n; ++i)
    {
        endpoints[i] = -2;
        startpoints[i] = -2;
    }
    // first we compute all the cumulative sums of all
    // coordinates:

    SEXP scalesSEXP = PROTECT(allocVector(REALSXP, p));
    double * scales = REAL(scalesSEXP);
    for (int i = 0; i < p; ++i)
    {
        scales[i] = 1;
    }

    if (rescale_variance_bool)
    {
        SEXP tmpvecSEXP = PROTECT(allocVector(REALSXP, n));
        double * tmpvec = REAL(tmpvecSEXP);
        memset(tmpvec, 0, sizeof(double)*n);
        rescale_variance(X, scales,  n, p, tmpvec);
        UNPROTECT(1);
    }


    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*(n+1)*sizeof(double));

    for (int j = 1; j <=n; ++j)
    {
        for (int i = 0; i < p; ++i)
        {
            //#define cord_spec(r,c, D) ((r) + (D)*(c))
            cumsums[cord_spec(i,j,p)] = X[cord_spec(i,j-1, p)] +cumsums[cord_spec(i,j-1, p)];
        }
    }

    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p);

    SEXP vecSEXP = PROTECT(allocVector(INTSXP, maxx));
    int * vec = INTEGER(vecSEXP);
    memset(vec, 0, sizeof(int)*maxx);
    SEXP vec2SEXP = PROTECT(allocVector(REALSXP, p));
    double * vec2 = REAL(vec2SEXP);
    memset(vec2, 0, sizeof(double)*p);




    // record all segment starts
    // we don't compute all cusums here but do it on the fly so we dont have to do anything
    // unecessary


    cPilliat_call(X, n, p, changepoints,
                thresholds_partial, threshold_dense, thresholds_bj,cumsums, lens, lenLens,
                vec, vec2, K, cusum, changepoint_counter_ptr, teststat,
                startpoints, endpoints, maxx, debug);

    /*cPilliat_call(X, -1, n-1, n, p, 1, changepoints,changepoint_counter_ptr,  depthcounter,
                thresholds , thresholds_test, cumsums, lens, lenLens,  as, nu_as, len_as,
                segstarts, maxvalues, maxpos, maxas, K, cusum, vector, coordchg,maxval, 
                startpoints, endpoints, maxaposes, tmpvec, twologn, ts,debug);*/

/*    cInspect_call(X, -1, n-1, n, p, 1, changepoints, changepoint_counter_ptr, depthcounter,
                maxval, xi, cumsums, lens, lenLens, lambda,
                eps, maxiter, segstarts, maxvalues, maxpos, K, cusum, mhat,
                mhatprod, v, v2,debug,coordchg);
*/
    // return:
    //Rprintf("found %d changepoints", *changepoint_counter_ptr);
    SEXP changepointnumSEXP = PROTECT(allocVector(INTSXP, 1));
    int * changepointnum = INTEGER(changepointnumSEXP);
    *changepointnum = changepoint_counter;

    SEXP ret = PROTECT(allocVector(VECSXP, 6)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out);
    SET_VECTOR_ELT(ret, 1, startpointsSEXP);
    SET_VECTOR_ELT(ret, 2, endpointsSEXP);
    SET_VECTOR_ELT(ret, 3, changepointnumSEXP);
    SET_VECTOR_ELT(ret, 4, teststatSEXP);
    SET_VECTOR_ELT(ret, 5, scalesSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 6));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("changepoints"));
    SET_STRING_ELT(names, 1, mkChar("startpoints"));
    SET_STRING_ELT(names, 2, mkChar("endpoints"));
    SET_STRING_ELT(names, 3, mkChar("number_of_changepoints"));
    SET_STRING_ELT(names, 4, mkChar("test_stat")); // 0 dense, 1 partial, 2 bj
    SET_STRING_ELT(names, 5, mkChar("scales")); 

    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(24);
    return ret;
}


/*
SEXP cHDCD_single(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdsI,
    SEXP asI, SEXP nu_asI, SEXP len_asI, SEXP debugI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(thresholdsI);
    
    PROTECT(asI);
    PROTECT(nu_asI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(len_asI);
    PROTECT(debugI);

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double *thresholds = (REAL(thresholdsI));
    
    double * as = REAL(asI);
    double * nu_as = REAL(nu_asI);
    int len_as = *(INTEGER(len_asI));
    int debug = *INTEGER(debugI);

    UNPROTECT(3); // unprotecting all except the arrays
    if(debug){
      Rprintf("p = %d\n", p);
      Rprintf("n = %d\n", n);
    }


    SEXP tmpvecSEXP = PROTECT(allocVector(REALSXP, len_as));
    double * tmpvec = REAL(tmpvecSEXP);
    memset(tmpvec, 0, sizeof(double)*len_as);

    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*(n+1)*sizeof(double));

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


    SEXP maxposSEXP = PROTECT(allocVector(INTSXP, 1));
    SEXP maxaposSEXP = PROTECT(allocVector(INTSXP, 1));
    int * maxpos = INTEGER(maxposSEXP);
    int * maxapos = INTEGER(maxaposSEXP);
    *maxpos = -10;
    double maximum = -1000000000000000000000.0;
    *maxapos = 0;
    int s = -1;
    int e = n-1;
    internal_check_segment( cumsums, cusum, maxpos,  &maximum, maxapos,
                             s, e, p, NULL, thresholds, thresholds,
                            as, nu_as, len_as,  tmpvec,0, NULL,debug);


    SEXP ret = PROTECT(allocVector(VECSXP, 2)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, maxposSEXP);
    SET_VECTOR_ELT(ret, 1, maxaposSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("pos"));
    SET_STRING_ELT(names, 1, mkChar("apos"));


    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(12);
    return ret;


}*/

SEXP cPilliat_calibrate(SEXP nI, SEXP pI,SEXP NI, SEXP tolnI, SEXP lensI,SEXP lenLensI,SEXP KI,
    SEXP maxxI, SEXP log2pI, SEXP rescale_variance_boolI, SEXP debugI){

    PROTECT(nI);
    PROTECT(pI);
    PROTECT(log2pI);
    PROTECT(lensI);
    PROTECT(lenLensI);
    PROTECT(KI);
    PROTECT(NI);
    PROTECT(tolnI);
    
    PROTECT(maxxI);
    PROTECT(rescale_variance_boolI);
    PROTECT(debugI);

    int n = *(INTEGER(nI));
    int toln = *(INTEGER(tolnI));
    int N = *(INTEGER(NI));
    int p = *(INTEGER(pI));
    int log2p = *(INTEGER(log2pI));
    int maxx = *INTEGER(maxxI);
    int lenLens = *(INTEGER(lenLensI));
    int * lens = INTEGER(lensI);
    int debug = *INTEGER(debugI);
    int rescale_variance_bool = *INTEGER(rescale_variance_boolI);

    if(debug){
        Rprintf("n = %d\n",n);
        Rprintf("p = %d\n",p);
        Rprintf("N = %d\n",N);
        Rprintf("toln = %d\n",toln);
        Rprintf("log2p = %d\n",log2p);
        Rprintf("maxx = %d\n",maxx);
        //return(NULL);
    }

    SEXP maxvals_partialSEXP = PROTECT(allocVector(REALSXP, log2p*N));
    double * maxvals_partial = REAL(maxvals_partialSEXP);
    memset(maxvals_partial, 0, sizeof(double)*log2p*N);

    SEXP maxvals_denseSEXP = PROTECT(allocVector(REALSXP, N));
    double * maxvals_dense = REAL(maxvals_denseSEXP);
    for (int i = 0; i < N; ++i)
    {
        maxvals_dense[i] = NEGINF;
    }

    SEXP maxvals_bjSEXP = PROTECT(allocVector(INTSXP, maxx*N));
    int * maxvals_bj = INTEGER(maxvals_bjSEXP);
    memset(maxvals_bj, 0, sizeof(int)*maxx*N);



    


    SEXP vecSEXP = PROTECT(allocVector(INTSXP, maxx));
    int * vec = INTEGER(vecSEXP);
    memset(vec, 0, sizeof(int)*maxx);

    SEXP vec2SEXP = PROTECT(allocVector(REALSXP, p));
    double * vec2 = REAL(vec2SEXP);
    memset(vec2, 0, sizeof(double)*p);


    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*(n+1)*sizeof(double));



    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p);

    SEXP vectorSEXP = PROTECT(allocVector(REALSXP, n));
    double * vector = REAL(vectorSEXP);
    memset(vector, 0, sizeof(double)*n);

    int jump;
    int jj;
    int tmp;
    int z;
    double cumsum;
    int prev;
    int j;
    int len = 0;
    int k = 0;
    int h;
    int i;
    int v;
    int c;
    int stop;

    GetRNGstate();

    SEXP XSEXP = PROTECT(allocVector(REALSXP, p * (n)));
    double *X = REAL(XSEXP); // p \times (n+1). first col is 0
    memset(X, 0, p*(n)*sizeof(double));

    for (k = 0; k < N; ++k)
    {
        
        // generate X

        for (i = 0; i < p; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                X[cord_spec(i,j,p)] = rnorm(0.0, 1.0);
            }
        }

        // normalize / rescale variance
        if (rescale_variance_bool)
        {
            rescale_variance(X, NULL,n, p, vector);
        }
        
        memset(cumsums, 0, p*sizeof(double));

        for (i = 0; i < p; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                cumsums[cord_spec(i,j+1,p)] = X[cord_spec(i,j,p)] +cumsums[cord_spec(i,j, p)];
            }
        }


        for (j = 0; j < lenLens; ++j)
        {
            len = lens[j];
            if(debug){
                Rprintf("j=%d, len = %d\n", j, len);
            }

            jump = len /2;

            if(len ==1){
                for (v = 0; v < n-1; ++v)
                {
                    // checv segments
                    if(debug){
                        Rprintf("checking v = %d, l = %d, interval [%d, %d)\n", v, len, v-1, v+1);
                    }
                    //res = internal_check_segment_Pilliat(cumsums, cusum, v, len,  p, thresholds_partial, 
                    //thresholds_bj, threshold_dense, maxx, &detected,vec, vec2, debug);

                

                    // compute CUSUM
                    singleCUSUM(cumsums, cusum, v-len, v+len, p, v); 

                    z = 0;

                    cumsum = 0.0;
                    prev = -1;
                    // double x; 
                    // int N=0;
                    


                    // first check dense:

                    cumsum = -p;
                    for (i = 0; i < p; ++i)
                    {
                        cumsum+= cusum[cord_spec(i,0, p)]*cusum[cord_spec(i,0, p)];
                        vec2[i] = cusum[cord_spec(i,0, p)]*cusum[cord_spec(i,0, p)];
                    }
                    if(cumsum > maxvals_dense[cord_spec(k,0,N)]){
                        maxvals_dense[cord_spec(k,0,N)] = cumsum;
                        //return 0;
                    }

                    // then partial sum:
                    R_qsort(vec2,1,p);

                    cumsum = 0.0;
                    prev = p;
                    z = 1;
                    int c = 0;
                    int i = 0;
                    while(1)
                    {
                        if(z>p){
                            break;
                        }
                        for (i = prev-1; i >= p-z; --i)
                        {
                            cumsum +=vec2[i];
                        }

                        if (cumsum> maxvals_partial[cord_spec(k,c,N)])
                        {
                            maxvals_partial[cord_spec(k,c,N)] = cumsum;
                        }
                        prev = p-z;
                        c++;
                        z = 2*z;


                    }
                    // then berk jones:
                    //for (int h = 0; h <= maxx; ++h)
                    //{
                    tmp = 0;
                    memset(vec, 0, sizeof(int)*maxx );

                    for (i = 0; i < p; ++i)
                    {
                        tmp = (int)  fabs(cusum[cord_spec(i,0, p)]);
                        //if(debug){
                        //    Rprintf("val = %f, res = %d\n",fabs(cusum[cord_spec(i,0, p)]), tmp );
                        //}
                        if(tmp == 0){
                            break;
                        }else{

                            for ( h = 0; h < tmp && h<maxx; ++h)
                            {
                                vec[h]++;
                            }
                            
                        }             
                    }

                    for (int h = 0; h < maxx; ++h)
                    {
                        if(vec[h] > maxvals_bj[cord_spec(k,h,N)]){
                            maxvals_bj[cord_spec(k,h,N)] = vec[h];
                            // if(debug){
                            //     Rprintf("Berv 0ones detected at x = %d, seg [%d,%d). Count = %d, thresh = %d.\n", 
                            //         h, v-len, v+len, vec[h], thresholds_b0[h]);
                            // }
                            //return 2;
                        }
                    }


                    
                }
            }
            else{
                //stop=0;
                v = len;
                jj = 2;
                stop = 0;
                while(1){
                    if(v + len>=n || v-len<-1){
                        break;
                    }
                    
                    //res = internal_checv_segment_Pilliat(cumsums, cusum, v, len,  p, thresholds_partial, 
                    //        thresholds_bj, threshold_dense,  maxx, &detected,vec, vec2, debug);
                    if(debug){
                        Rprintf("checking v = %d, l = %d, interval [%d, %d)\n", v, len, v-1, v+1);
                    }
                    //res = internal_check_segment_Pilliat(cumsums, cusum, v, len,  p, thresholds_partial, 
                    //thresholds_bj, threshold_dense, maxx, &detected,vec, vec2, debug);

                

                    // compute CUSUM
                    singleCUSUM(cumsums, cusum, v-len, v+len, p, v); 

                    z = 0;

                    cumsum = 0.0;
                    prev = -1;
                    // double x; 
                    // int N=0;
                    


                    // first check dense:

                    cumsum = -p;
                    for (i = 0; i < p; ++i)
                    {
                        cumsum+= cusum[cord_spec(i,0, p)]*cusum[cord_spec(i,0, p)];
                        vec2[i] = cusum[cord_spec(i,0, p)]*cusum[cord_spec(i,0, p)];
                    }
                    if(cumsum > maxvals_dense[cord_spec(k,0,N)]){
                        maxvals_dense[cord_spec(k,0,N)] = cumsum;
                        //return 0;
                    }

                    // then partial sum:
                    R_qsort(vec2,1,p);

                    cumsum = 0.0;
                    prev = p;
                    z = 1;
                    c = 0;
                    i = 0;
                    while(1)
                    {
                        if(z>p){
                            break;
                        }
                        for (i = prev-1; i >= p-z; --i)
                        {
                            cumsum +=vec2[i];
                        }

                        if (cumsum> maxvals_partial[cord_spec(k,c,N)])
                        {
                            maxvals_partial[cord_spec(k,c,N)] = cumsum;
                        }
                        prev = p-z;
                        c++;
                        z = 2*z;


                    }
                    // then berk jones:
                    //for (int h = 0; h <= maxx; ++h)
                    //{
                    tmp = 0;
                    memset(vec, 0, sizeof(int)*maxx );

                    for (i = 0; i < p; ++i)
                    {
                        tmp = (int)  fabs(cusum[cord_spec(i,0, p)]);
                        //if(debug){
                        //    Rprintf("val = %f, res = %d\n",fabs(cusum[cord_spec(i,0, p)]), tmp );
                        //}
                        if(tmp == 0){
                            break;
                        }else{

                            for ( h = 0; h < tmp && h<maxx; ++h)
                            {
                                vec[h]++;
                            }
                            
                        }             
                    }

                    for (h = 0; h < maxx; ++h)
                    {
                        if(vec[h] > maxvals_bj[cord_spec(k,h,N)]){
                            maxvals_bj[cord_spec(k,h,N)] = vec[h];
                            // if(debug){
                            //     Rprintf("Berv 0ones detected at x = %d, seg [%d,%d). Count = %d, thresh = %d.\n", 
                            //         h, v-len, v+len, vec[h], thresholds_b0[h]);
                            // }
                            //return 2;
                        }
                    }

                    if(stop){
                        break;
                    }
                    v = (jj++)*jump;
                    if(v>= n - len-1){
                        v = n-len-1;
                        stop=1;
                    }
                
                
                }
            }
        }

    }


 

    PutRNGstate();

    SEXP thresholds_partialSEXP = PROTECT(allocVector(REALSXP, log2p ));
    double * thresholds_partial = REAL(thresholds_partialSEXP);
    memset(thresholds_partial, 0, sizeof(double)*log2p);
    SEXP threshold_denseSEXP = PROTECT(allocVector(REALSXP, 1));
    double * threshold_dense = REAL(threshold_denseSEXP);
    *threshold_dense = 0.0;
    SEXP thresholds_bjSEXP = PROTECT(allocVector(INTSXP, maxx ));
    int * thresholds_bj = INTEGER(thresholds_bjSEXP);
    memset(thresholds_bj, 0, sizeof(int)*maxx);

    //int kkk = 1;
    for (int z = 0; z < log2p; ++z){
        sort_k_largest(maxvals_partial + cord_spec(0,z,N) , toln, 0, N);
        thresholds_partial[z] = maxvals_partial[cord_spec((toln-1), z, N)];
        //if(kkk >2 *sqrt(p*log(n))){
        //   thresholds_partial[z] = 100000000.0;
        //}

        //kkk *=2;

    }

    for (int z = 0; z < maxx; ++z){
        sort_k_largest_int(maxvals_bj + cord_spec(0,z,N) , toln, 0, N);
        thresholds_bj[z] = maxvals_bj[cord_spec((toln-1), z, N)];
        //thresholds_bj[z] = p;

    }
    sort_k_largest(maxvals_dense, toln, 0, N);
    threshold_dense[0] = maxvals_dense[toln-1];

    SEXP ret = PROTECT(allocVector(VECSXP, 3)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, thresholds_partialSEXP);
    SET_VECTOR_ELT(ret, 1, threshold_denseSEXP);
    SET_VECTOR_ELT(ret, 2, thresholds_bjSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 3));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("thresholds_partial"));
    SET_STRING_ELT(names, 1, mkChar("threshold_dense"));
    SET_STRING_ELT(names, 2, mkChar("thresholds_bj"));


    setAttrib(ret, R_NamesSymbol, names);


    UNPROTECT(25);
    return(ret);





}



SEXP cPilliat_test(SEXP XI,SEXP nI, SEXP pI,SEXP thresholds_partialI, SEXP threshold_denseI, SEXP thresholds_bjI,
    SEXP maxxI, SEXP debugI){

    PROTECT(XI);
    PROTECT(nI);
    PROTECT(pI);

    PROTECT(thresholds_partialI);
    PROTECT(threshold_denseI);
    PROTECT(thresholds_bjI);
    
    PROTECT(maxxI);
    PROTECT(debugI);

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double *thresholds_partial = (REAL(thresholds_partialI));
    double threshold_dense = *(REAL(threshold_denseI));
    int *thresholds_bj = (INTEGER(thresholds_bjI));
    int maxx = *INTEGER(maxxI);
    int debug = *INTEGER(debugI);


    SEXP out = PROTECT(allocVector(INTSXP, 1));
    int *result  = INTEGER(out);
    *result = 0;

    SEXP testposSEXP = PROTECT(allocVector(INTSXP, n));
    int *testpos  = INTEGER(testposSEXP);
    memset(testpos, 0, sizeof(int)*n);



    // find test poses
    int count = 0;
    int s = -1;
    int e = n-1;
    int middle = (s+e)/2;
    testpos[count++] = middle;
    int jump = 2;

 

    while(middle + jump < e || middle - jump > s){
        if(middle + jump <e){
            testpos[count++] = middle + jump;
        }
        if(middle - jump > s){
            testpos[count++] = middle - jump;
        }
        jump = jump*2;
    }

    R_qsort_int(testpos,1,count);

    SEXP vecSEXP = PROTECT(allocVector(INTSXP, maxx));
    int * vec = INTEGER(vecSEXP);
    memset(vec, 0, sizeof(int)*maxx);

    SEXP vec2SEXP = PROTECT(allocVector(REALSXP, p));
    double * vec2 = REAL(vec2SEXP);
    memset(vec2, 0, sizeof(double)*p);

    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*(n+1)*sizeof(double));

    for (int j = 1; j <=n; ++j)
    {
        for (int i = 0; i < p; ++i)
        {
            //#define cord_spec(r,c, D) ((r) + (D)*(c))
            cumsums[cord_spec(i,j,p)] = X[cord_spec(i,j-1, p)] +cumsums[cord_spec(i,j-1, p)];
        }
    }

    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p);

    int testpoint;
    //int detected = 0;
    int tmp;
    int z;
    double cumsum;
    int prev;
    //int NN;
    int j;
    for (int tt = 0; tt < count; ++tt)
    {
        testpoint = testpos[tt];
        singleCUSUM(cumsums, cusum, s, e, p, testpoint);
        //memset(tmpvec, 0, sizeof(double)*len_as);
        //singleCUSUM(cumsums, cusum, k-len,e, p, k); 
        z = 0;

        cumsum = 0.0;
        prev = -1;
        //NN=0;
        j = 0;


        // first check dense:

        cumsum = -p;
        for (int i = 0; i < p; ++i)
        {
            cumsum+= cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
            vec2[i] = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
        }
        if(cumsum > threshold_dense){
            if(debug){
                Rprintf("Detected dense in %d, val %f, thresh %f\n", testpoint, cumsum, threshold_dense);
            }
            *result = 1;
            break;
        }

        // then partial sum:
        R_qsort(vec2,1,p);

        cumsum = 0.0;
        prev = p;
        z = 1;
        int c = 0;
        int i = 0;
        while(1)
        {
            if(z>p){
                break;
            }
            for (i = prev-1; i >= p-z; --i)
            {
                cumsum +=vec2[i];
            }

            if (cumsum> thresholds_partial[c])
            {
                *result = 1;
                if(debug){
                    Rprintf("Partial stat detected at sparsity %d, pos %d. val = %f, thresh = %f.\n", 
                        z, testpoint, cumsum, thresholds_partial[c]);
                }
                break;
            }
            prev = p-z;
            c++;
            z = 2*z;


        }
        if(*result){
            break;
        }
        // then berk jones:
        //for (int h = 0; h <= maxx; ++h)
        //{
        tmp = 0;
        memset(vec, 0, sizeof(int)*maxx );

        for (int i = 0; i < p; ++i)
        {
            tmp = (int)  fabs(cusum[cord_spec(i,j, p)]);
            //if(debug){
            //    Rprintf("val = %f, res = %d\n",fabs(cusum[cord_spec(i,j, p)]), tmp );
            //}

            for (int h = 0; h < tmp; ++h)
            {
                if(h>=maxx){
                    break;
                }
                vec[h]++;
            }
                
        }

        for (int h = 0; h < maxx; ++h)
        {
            if(vec[h] > thresholds_bj[h]){
                *result = 1;
                if(debug){
                    Rprintf("Berk Jones detected at x = %d, pos %d. Count = %d, thresh = %d.\n", 
                        h, testpoint, vec[h], thresholds_bj[h]);
                }
                break;
            }
        }

        if(*result){
            break;
        }
    }


    UNPROTECT(14);
    return(out);

}


SEXP cPilliat_test_calibrate(SEXP nI, SEXP pI,SEXP NI, SEXP tolnI,
    SEXP maxxI, SEXP log2pI, SEXP debugI){

    PROTECT(nI);
    PROTECT(pI);
    PROTECT(log2pI);

    PROTECT(NI);
    PROTECT(tolnI);
    
    PROTECT(maxxI);
    PROTECT(debugI);

    int n = *(INTEGER(nI));
    int toln = *(INTEGER(tolnI));
    int N = *(INTEGER(NI));
    int p = *(INTEGER(pI));
    int log2p = *(INTEGER(log2pI));
    int maxx = *INTEGER(maxxI);
    int debug = *INTEGER(debugI);

    if(debug){
        Rprintf("n = %d\n",n);
        Rprintf("p = %d\n",p);
        Rprintf("N = %d\n",N);
        Rprintf("toln = %d\n",toln);
        Rprintf("log2p = %d\n",log2p);
        Rprintf("maxx = %d\n",maxx);
        //return(NULL);
    }

    SEXP maxvals_partialSEXP = PROTECT(allocVector(REALSXP, log2p*N));
    double * maxvals_partial = REAL(maxvals_partialSEXP);
    memset(maxvals_partial, 0, sizeof(double)*log2p*N);

    SEXP maxvals_denseSEXP = PROTECT(allocVector(REALSXP, N));
    double * maxvals_dense = REAL(maxvals_denseSEXP);
    for (int i = 0; i < N; ++i)
    {
        maxvals_dense[i] = NEGINF;
    }

    SEXP maxvals_bjSEXP = PROTECT(allocVector(INTSXP, maxx*N));
    int * maxvals_bj = INTEGER(maxvals_bjSEXP);
    memset(maxvals_bj, 0, sizeof(int)*maxx*N);



    

    SEXP testposSEXP = PROTECT(allocVector(INTSXP, n));
    int *testpos  = INTEGER(testposSEXP);
    memset(testpos, 0, sizeof(int)*n);

    // find test poses
    int count = 0;
    int s = -1;
    int e = n-1;
    int middle = (s+e)/2;
    testpos[count++] = middle;
    int jump = 2;

 

    while(middle + jump < e || middle - jump > s){
        if(middle + jump <e){
            testpos[count++] = middle + jump;
        }
        if(middle - jump > s){
            testpos[count++] = middle - jump;
        }
        jump = jump*2;
    }

    R_qsort_int(testpos,1,count);
    if(debug){
        Rprintf("Testposes:\n");
        for (int i = 0; i < count; ++i)
        {
            Rprintf("%d\n",testpos[i]);
        }
    }

    SEXP vecSEXP = PROTECT(allocVector(INTSXP, maxx));
    int * vec = INTEGER(vecSEXP);
    memset(vec, 0, sizeof(int)*maxx);

    SEXP vec2SEXP = PROTECT(allocVector(REALSXP, p));
    double * vec2 = REAL(vec2SEXP);
    memset(vec2, 0, sizeof(double)*p);


    SEXP cumsumsSEXP = PROTECT(allocVector(REALSXP, p * (n+1)));
    double *cumsums = REAL(cumsumsSEXP); // p \times (n+1). first col is 0
    memset(cumsums, 0, p*(n+1)*sizeof(double));



    SEXP cusumSEXP = PROTECT(allocVector(REALSXP, p));
    double * cusum = REAL(cusumSEXP);
    memset(cusum, 0, sizeof(double)*p);



    int testpoint;
    //int detected = 0;
    int tmp;
    int z;
    double cumsum;
    int prev;
    //int NN;
    int j;

    GetRNGstate();

    for (int k = 0; k < N; ++k)
    {
        if(debug){
            Rprintf("ITER %d\n", k);
        }
        // generate X
        memset(cumsums, 0, p*sizeof(double));

        for (int i = 0; i < p; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                cumsums[cord_spec(i,j+1,p)] = rnorm(0.0, 1.0) +cumsums[cord_spec(i,j, p)];
            }
        }




 
        for (int tt = 0; tt < count; ++tt)
        {
            if(debug && !(k)){
                Rprintf("pos %d\n",testpos[tt]);
            }
            testpoint = testpos[tt];
            singleCUSUM(cumsums, cusum, s, e, p, testpoint);
            //memset(tmpvec, 0, sizeof(double)*len_as);
            //singleCUSUM(cumsums, cusum, k-len,e, p, k); 
            z = 0;

            cumsum = 0.0;
            prev = -1;
            //NN=0;
            j = 0;


            // first check dense:

            cumsum = -p;
            for (int i = 0; i < p; ++i)
            {
                cumsum+= cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
                vec2[i] = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
            }

            if(cumsum > maxvals_dense[cord_spec(k,0, N)]){
                maxvals_dense[cord_spec(k,0, N)] = cumsum;
            }


            // then partial sum:
            R_qsort(vec2,1,p);

            cumsum = 0.0;
            prev = p;
            z = 1;
            int c = 0;
            int i = 0;
            while(1)
            {
                if(z>p){
                    break;
                }
                for (i = prev-1; i >= p-z; --i)
                {
                    cumsum +=vec2[i];
                }

                if (cumsum> maxvals_partial[cord_spec(k,c, N)])
                {
                    maxvals_partial[cord_spec(k,c, N)] = cumsum;
                    if(debug){
                        Rprintf("cumsum at partial %d = %f\n", z, cumsum);
                    }
                }
                prev = p-z;
                c++;
                z = 2*z;


            }
      
            // then berk jones:
            //for (int h = 0; h <= maxx; ++h)
            //{
            tmp = 0;
            memset(vec, 0, sizeof(int)*maxx );

            for (int i = 0; i < p; ++i)
            {
                tmp = (int)  fabs(cusum[cord_spec(i,j, p)]);
                //if(debug){
                //    Rprintf("val = %f, res = %d\n",fabs(cusum[cord_spec(i,j, p)]), tmp );
                //}
               

                for (int h = 0; h < tmp ; ++h)
                {
                    if(h>=maxx){
                        break;
                    }
                    vec[h]++;
                }
                    
                          
            }

            for (int h = 0; h < maxx; ++h)
            {
                if(vec[h] > maxvals_bj[cord_spec(k,h, N)]){
                    maxvals_bj[cord_spec(k,h, N)] = vec[h];
                    if(debug){
                        Rprintf("BJ: at h = %d we found %d counts\n", h, vec[h]);
                    }
                }
            }

        }
    }

    PutRNGstate();

    SEXP thresholds_partialSEXP = PROTECT(allocVector(REALSXP, log2p ));
    double * thresholds_partial = REAL(thresholds_partialSEXP);
    memset(thresholds_partial, 0, sizeof(double)*log2p);
    SEXP threshold_denseSEXP = PROTECT(allocVector(REALSXP, 1));
    double * threshold_dense = REAL(threshold_denseSEXP);
    *threshold_dense = 0.0;
    SEXP thresholds_bjSEXP = PROTECT(allocVector(INTSXP, maxx ));
    int * thresholds_bj = INTEGER(thresholds_bjSEXP);
    memset(thresholds_bj, 0, sizeof(int)*maxx);

    //int kkk = 1;
    for (int z = 0; z < log2p; ++z){
        sort_k_largest(maxvals_partial + cord_spec(0,z,N) , toln, 0, N);
        thresholds_partial[z] = maxvals_partial[cord_spec((toln-1), z, N)];
        //if(kkk >2 *sqrt(p*log(n))){
        //   thresholds_partial[z] = 100000000.0;
        //}

        //kkk *=2;

    }

    for (int z = 0; z < maxx; ++z){
        sort_k_largest_int(maxvals_bj + cord_spec(0,z,N) , toln, 0, N);
        thresholds_bj[z] = maxvals_bj[cord_spec((toln-1), z, N)];
        //thresholds_bj[z] = p;

    }
    sort_k_largest(maxvals_dense, toln, 0, N);
    threshold_dense[0] = maxvals_dense[toln-1];

    SEXP ret = PROTECT(allocVector(VECSXP, 3)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, thresholds_partialSEXP);
    SET_VECTOR_ELT(ret, 1, threshold_denseSEXP);
    SET_VECTOR_ELT(ret, 2, thresholds_bjSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 3));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("thresholds_partial"));
    SET_STRING_ELT(names, 1, mkChar("threshold_dense"));
    SET_STRING_ELT(names, 2, mkChar("thresholds_bj"));


    setAttrib(ret, R_NamesSymbol, names);


    UNPROTECT(20);
    return(ret);





}




