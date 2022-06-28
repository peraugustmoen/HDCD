#include "header.h"


int internal_check_segment_Pilliat(double * cumsums, double * cusum, int k, int len, int p, double * thresholds_partial, 
    int * thresholds_bj, double threshold_dense, int maxx,  int * detected,
    int * vec, int debug){
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
    double x; 
    int N=0;
    int j = 0;

    //for (int j = 0; j < e-s-1; ++j){

    // first check dense:

    cumsum = -p;
    for (int i = 0; i < p; ++i)
    {
        cumsum+= cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
    }
    if(cumsum > threshold_dense){
        *detected = 1;
        return 0;
    }

    // then partial sum:
    R_qsort(&(cusum[cord_spec(0,j, p)]),1,p);

    cumsum = 0.0;
    prev = -1;
    z = 1;
    int c = 0;
    int i = 0;
    while(1)
    {
        if(z>p){
            break;
        }
        for (i = prev+1; i <=z; ++i)
        {
            cumsum +=cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
        }

        if (cumsum> thresholds_partial[c])
        {
            *detected = 1;
            return 1;
        }
        prev = z;
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
                int * vec, int K, double * cusum, int* changepoint_counter, int * teststat,
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
    int intersects = 0;
    int t1 = 0;
    int t2 = 0;
    int check = 1;
    for (int j = 0; j < lenLens; ++j)
    {
        len = lens[j]-1;
        if(debug){
            Rprintf("j=%d, len = %d\n", j, len);
        }

        jump = len /2;


        


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
                thresholds_bj, threshold_dense, maxx, &detected,vec, debug);

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
                            thresholds_bj, threshold_dense,  maxx, &detected,vec, debug);
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
    SEXP lensI,SEXP lenLensI,SEXP KI, SEXP maxxI,SEXP debugI){
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




    // record all segment starts
    // we don't compute all cusums here but do it on the fly so we dont have to do anything
    // unecessary


    cPilliat_call(X, n, p, changepoints,
                thresholds_partial, threshold_dense, thresholds_bj,cumsums, lens, lenLens,
                vec, K, cusum, changepoint_counter_ptr, teststat,
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

    SEXP ret = PROTECT(allocVector(VECSXP, 5)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out);
    SET_VECTOR_ELT(ret, 1, startpointsSEXP);
    SET_VECTOR_ELT(ret, 2, endpointsSEXP);
    SET_VECTOR_ELT(ret, 3, changepointnumSEXP);
    SET_VECTOR_ELT(ret, 4, teststatSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 5));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("changepoints"));
    SET_STRING_ELT(names, 1, mkChar("startpoints"));
    SET_STRING_ELT(names, 2, mkChar("endpoints"));
    SET_STRING_ELT(names, 3, mkChar("number_of_changepoints"));
    SET_STRING_ELT(names, 4, mkChar("test_stat")); // 0 dense, 1 partial, 2 bj

    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(21);
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
