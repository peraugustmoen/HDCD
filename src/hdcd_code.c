#include "header.h"

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
    //memset(vector, 0, c1);
    for (int j = 0; j < c1; ++j)
    {
        vector[j]=0.0;
        for (int i = 0; i < r1; ++i)
        {
            vector[j]+=matrix[cord_spec(i, j, r1)];
        }
    }

}
void internal_check_segment(double * cumsums, double * cusum, int * maxpos, double * maximum, int * maxa_pos,
                            int s, int e, int p, double * vector, double * thresholds, double * thresholds_test,
                            double * as, double * nu_as, int len_as, double * tmpvec, int twologn, int * ts,int debug){
    // require as[0] = 0
    if(debug){
      Rprintf("checking segment (%d,%d]\n", s,e);

    }

    *maxpos  = s+(e-s)/2;
    *maximum = NEGINF+1;
    *maxa_pos = 0;
    if(e-s<2){
        return;
    }

    // compute CUSUM
    CUSUM(cumsums, cusum, s, e, p); // dim is p \times (e-s-1)

    //double prev_nu_a = as[0];
    //double a;
    //double nu_a;
    double tmp = 0;
    //double a_tmp=-100000;
    //int pos_a_tmp = 0;
    //double * val=0;
    //int n = e-s;

    //int detected = 0;
    double tmp2 =0;
    int localdetected = 0;
    for (int j = 0; j < e-s-1; ++j){
        memset(tmpvec, 0, sizeof(double)*len_as);
        
        // first aggregate thresholded CUSUMs
        for (int i = 0; i < p; ++i){
            tmp2 = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
            for (int z = 0; z < len_as; ++z){
                if (fabs(cusum[cord_spec(i,j, p)])>as[z])
                {
                    tmpvec[z] += tmp2 - nu_as[z];
                }
                else{
                    break;
                }
                
            }
        }

        // check if T >0
        localdetected = 0;
        for (int z = 0; z < len_as; ++z)
        {
            /*tmp = tmpvec[z] - thresholds[z];
            if (tmp> *maximum)
            {
                *maximum= tmp;
                *maxpos = s+j+1;
                *maxa_pos = z;
            }*/
            tmp = tmpvec[z] - thresholds_test[z];
            if(tmp>0 && ts[z] > twologn){
                localdetected = 1;
            }
        }

        // if T<= 0, check if P >0:
        if(twologn>0 && localdetected ==0){
            sort_k_largest(cusum + cord_spec(0,j,p) , twologn, 0, p);
            //partial_quicksort(cusum + cord_spec(0,j,p) , p , twologn);
            double cumsum = 0;
            int prev = 0;

            for (int z = len_as-1; z >= 0; --z)
            {
                if(ts[z] > twologn){
                    break;
                }
     
                if(z==0){
                    prev = -1;
                }
                for (int hh = prev+1; hh < ts[z]; ++hh)
                {
                    cumsum += cusum[cord_spec(hh,j,p)]*cusum[cord_spec(hh,j,p)];
                }
                prev = ts[z];
                tmp = cumsum;

                if(tmp>thresholds_test[z]){
                    localdetected = 1;
                }
            }
        }

        if (localdetected)
        {
            for (int z = 0; z < len_as; ++z)
            {
                if(fabs(tmpvec[z])<1e-10){
                    continue;
                }
                //tmp = tmpvec[z] - thresholds[z];
                //if(tmpvec[z] > thresholds_test[z]){
                    tmp = tmpvec[z] - thresholds[z];
                    if (tmp> *maximum)
                    {
                        *maximum= tmp;
                        *maxpos = s+j+1;
                        *maxa_pos = z;
                    }
                //}
                
            }
        }

        
        /*}else{
            for (int z = 0; z < len_as; ++z)
            {
                if(fabs(tmpvec[z])<1e-10){
                    continue;
                }
                //tmp = tmpvec[z] - thresholds[z];
                if(tmpvec[z] > thresholds_test[z]){
                    tmp = tmpvec[z] - thresholds[z];
                    if (tmp> *maximum)
                    {
                        *maximum= tmp;
                        *maxpos = s+j+1;
                        *maxa_pos = z;
                    }
                }
                
            }
        }*/
    }

/*    for (int i = 0; i < len_as; ++i)
    {

        a = as[i];
        nu_a = nu_as[i];
        internal_threshold_matrix(cusum, p, e-s-1, a, nu_a, i>0, prev_nu_a );
        internal_colSum(cusum, p, e-s-1, vector);
        prev_nu_a = nu_a;
        for (int j = 0; j < e-s-1; ++j)
        {
            //t = s+i+1;
            tmp = vector[j] -thresholds[i];
            if(tmp> *maximum){
                *maximum = tmp;
                *maxpos = s+j+1;
                *maxa_pos = i;
            }
            if(tmp>a_tmp){
                a_tmp = tmp;
                pos_a_tmp = s+j+1;

            }
        }

        if (debug)
        {

            Rprintf("for a=%f, max is %f (%f, thresh = %f) at pos %d\n", a, a_tmp, a_tmp*thresholds[i],thresholds[i], pos_a_tmp);
        }
        a_tmp = -100000;


    }*/

    if(debug){
      Rprintf("for segment (%d,%d] we found maximum in %d with val %f\n", s,e,*maxpos, *maximum);
    }

    return;
}


void cHDCD_call(double * x, int s, int e, int n, int p, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter,
                double * thresholds, double * thresholds_test, double *cumsums, int* lens, int lenLens, double * as, double * nu_as, int len_as,
                int * segstarts, double * maxvalues, int* maxpos,int * maxas, int K, double * cusum,
                double * vector, int * coordchg, double * maxval, int * startpoints, int* endpoints, int* maxaposes,
                double * tmpvec, int twologn, int * ts,int debug){
    if(debug){
        Rprintf("cHDCD_call! s=%d, e=%d\n", s, e);
    }

    if(e-s<lens[0]){
        //Rprintf("segment too short\n");
        return;
    }
    int argmax = s+1;
    double maximum = NEGINF;
    int maxa_pos = 0;

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
            Rprintf("j=%d, len = %d\n", j, len);
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
            if(debug){
                Rprintf("i= %d\n", i);
            }
            if(i>e-len || i<-1){
                if(debug){
                    Rprintf("i= %d is skipped\n", i);
                }
                break;
            }
            else if (i<s)
            {
                continue;
            }

            if(debug){
                //Rprintf("maxvalues[%d, %d] = %f\n", k , j , maxvalues[cord_spec(k,j,n)]);
            }

            if(maxvalues[cord_spec(k,j,n)]<=NEGINF){
                 if(debug){
                    Rprintf("segment (%d,%d] (k=%d, j=%d) not inspected, now checking!\n", i, i+len,k,j);
                 }
                //this segment not computed!
                //inspectOnSegment(cumsums, cusum, &tmpargmax, &tmpmaximum, s, e, p, lambda,
                 //   eps, maxiter, projvec, cusum_proj);
                 //double * cumsums, double * cusum, int * maxpos, double * maximum,
                 internal_check_segment(cumsums, cusum, &(maxpos[cord_spec(k,j,n)]), &(maxvalues[cord_spec(k,j,n)]), &(maxas[cord_spec(k,j,n)]),
                            i, i+len, p, vector, thresholds, thresholds_test, as, nu_as,len_as, tmpvec, twologn, ts,debug);
                 //internal_inspectOnSegment(cumsums, cusum, &(maxpos[cord_spec(k,j,n)]), &(maxvalues[cord_spec(k,j,n)]), i, i+len, p,
                 //lambda,
                //      eps, maxiter, mhat, mhatprod, v, v2,debug);
            }
            else if(maxvalues[cord_spec(k,j,n)] > NEGINF+1){
              if(debug){
                Rprintf("segment (%d,%d] (k=%d, j=%d) already inspected, with max val %f in %d\n", i, i+len,k,j,
                maxvalues[cord_spec(k,j,n)],maxpos[cord_spec(k,j,n)]);
              }
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

        if(maximum>NEGINF+1){
            break;
        }
/*        if(maximum>0.0){
            break;
        }
*/    }
    if(debug){
        Rprintf("maximum=%f\n", maximum);
    }

    if(maximum >NEGINF+1){
        if(debug){
          Rprintf("!!!!!! declared change-point in %d. val = %f. (s,e] = (%d,%d]\n", argmax, maximum,
          segstarts[cord_spec(k_max,j_max,n)], lens[j_max]+ segstarts[cord_spec(k_max,j_max,n)]);
          Rprintf("changeptcounter = %d\n", *changepoint_counter_ptr);
        }


        // identify in which coordinates the change happens:
        //i = segstarts[cord_spec(k_max,j_max,n)];
        //len = lens[j_max];
        //int ss = i;
        //int ee = i+len;
        if(maxa_pos==0){
            for (int zz = 0; zz < p; ++zz)
            {
                coordchg[cord_spec(zz,*changepoint_counter_ptr, p)]=1;
            }
        }
        else{
            i = segstarts[cord_spec(k_max,j_max,n)];
            len = lens[j_max];
            int ss = i;
            int ee = i+len;
            CUSUM(cumsums, cusum, ss, ee, p);
            internal_threshold_matrix(&(cusum[cord_spec(0,argmax-ss-1,p)]), p, 1, as[maxa_pos],  nu_as[maxa_pos], 0,
                                    0 );


            for (int zz = 0; zz < p; ++zz)
            {
                if(cusum[cord_spec(0,argmax-ss-1,p)+zz]>1e-10){
                    coordchg[cord_spec(zz,*changepoint_counter_ptr, p)]=1;
                }
            }

        }

        changepoints[*changepoint_counter_ptr] = argmax;
        depthcounter[*changepoint_counter_ptr] = depth;
        maxval[*changepoint_counter_ptr] = maximum;
        startpoints[*changepoint_counter_ptr] = segstarts[cord_spec(k_max,j_max,n)];
        endpoints[*changepoint_counter_ptr] = segstarts[cord_spec(k_max,j_max,n)] + lens[j_max];
        maxaposes[*changepoint_counter_ptr] = maxa_pos;
        (*changepoint_counter_ptr)++;
        if(*changepoint_counter_ptr >n){
          return;
        }
        //cHDCD_call(double * x, int s, int e, int n, int p, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter,
        //        double threshold_d, double threshold_s , double *cumsums, int* lens, int lenLens, double * as, double * nu_as, int len_as,
        //        int * segstarts, double * maxvalues, int* maxpos, int K, double * cusum, double * matrix,
        //        double * vector, int * coordchg, int debug)
        cHDCD_call(x, s, argmax, n, p, depth+1, changepoints,changepoint_counter_ptr,  depthcounter,
                thresholds, thresholds_test, cumsums, lens, lenLens,  as, nu_as, len_as,
                segstarts, maxvalues, maxpos,maxas, K, cusum, vector, coordchg, maxval, 
                startpoints, endpoints, maxaposes, tmpvec,twologn, ts,debug);
        cHDCD_call(x, argmax, e,n, p, depth+1, changepoints,changepoint_counter_ptr,  depthcounter,
                thresholds , thresholds_test, cumsums, lens, lenLens,  as, nu_as, len_as,
                segstarts, maxvalues, maxpos,maxas, K, cusum, vector, coordchg, maxval,
                startpoints, endpoints, maxaposes, tmpvec, twologn, ts,debug);

    }

    return;
}

SEXP cHDCD(SEXP XI,SEXP nI, SEXP pI,SEXP thresholdsI, SEXP thresholds_testI, SEXP lensI,SEXP lenLensI,SEXP KI,
    SEXP asI, SEXP nu_asI, SEXP len_asI, SEXP twolognI, SEXP tsI,SEXP debugI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(thresholdsI);
    PROTECT(thresholds_testI);
    PROTECT(lensI);
    PROTECT(asI);
    PROTECT(nu_asI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(lenLensI);
    PROTECT(KI);
    PROTECT(len_asI);
    PROTECT(twolognI);
    PROTECT(tsI);
    PROTECT(debugI);

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double *thresholds = (REAL(thresholdsI));
    double *thresholds_test = (REAL(thresholds_testI));
    int *lens = INTEGER(lensI); /////// dobbeltsjekk at int er like stor som INTEGER!!!!!!!
    int lenLens = *(INTEGER(lenLensI));
    int K = *(INTEGER(KI));
    double * as = REAL(asI);
    double * nu_as = REAL(nu_asI);
    int len_as = *(INTEGER(len_asI));
    int twologn = * (INTEGER(twolognI));
    int * ts = INTEGER(tsI);
    int debug = *INTEGER(debugI);

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
    int changepoint_counter = 0;
    int* changepoint_counter_ptr = &changepoint_counter;
    SEXP maxvalSEXP = PROTECT(allocVector(REALSXP, n));
    double * maxval = REAL(maxvalSEXP); //pointer to array
    memset(maxval, 0, sizeof(double)*n);
    SEXP startpointsSEXP = PROTECT(allocVector(INTSXP, n));
    int * startpoints = INTEGER(startpointsSEXP); //pointer to array
    memset(startpoints, 0, sizeof(int)*n);
    SEXP endpointsSEXP = PROTECT(allocVector(INTSXP, n));
    int * endpoints = INTEGER(endpointsSEXP); //pointer to array
    memset(endpoints, 0, sizeof(int)*n);
    SEXP maxaposesSEXP = PROTECT(allocVector(INTSXP, n));
    int * maxaposes = INTEGER(maxaposesSEXP); //pointer to array
    memset(maxaposes, 0, sizeof(int)*n);
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
    SEXP vectorSEXP = PROTECT(allocVector(REALSXP, maxlen));


    double * vector = REAL(vectorSEXP);
    memset(vector, 0, sizeof(double)*maxlen);




    // record all segment starts
    // we don't compute all cusums here but do it on the fly so we dont have to do anything
    // unecessary


    SEXP maxvaluesSEXP = PROTECT(allocVector(REALSXP, n*lenLens));
    double * maxvalues = REAL(maxvaluesSEXP);
    //memset(maxvalues, 0, sizeof(double)*n*lenLens);
    for (int i = 0; i < n*lenLens; ++i)
    {
        maxvalues[i] = NEGINF;
    }
    SEXP tmpvecSEXP = PROTECT(allocVector(REALSXP, len_as));
    double * tmpvec = REAL(tmpvecSEXP);
    memset(tmpvec, 0, sizeof(double)*len_as);
    SEXP maxposSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * maxpos = INTEGER(maxposSEXP);
    memset(maxpos, 0, sizeof(int)*n*lenLens);
    SEXP maxasSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * maxas = INTEGER(maxasSEXP);
    memset(maxas, 0, sizeof(int)*n*lenLens);
    SEXP segstartsSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * segstarts = INTEGER(segstartsSEXP);
    //memset(segstarts, -2, sizeof(int)*n*lenLens);
    for (size_t i = 0; i < n*lenLens; i++) {
      segstarts[i] = -2;
    }
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

        for (int i = -1; i < (n-len); i+=jump)
        {
            //cord_spec(r,c, D) ((r) + (D)*(c))
            //compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxvalues[cord_spec(i,j,n)]));
            segstarts[cord_spec(counter++,j,n)] = i;
            if(debug){
                //Rprintf("segstarts[%d, %d] = %d\n",counter-1, j, i );
            }


        }
    }

    if(debug){
      for (int j = 0; j < lenLens; ++j)
      {
          //counter=0;



          for (int i = 0; i < n; i++)
          {
              //cord_spec(r,c, D) ((r) + (D)*(c))
              //compute_cusum(cumsum, i, i+len, &(maxpos[cord_spec(i,j,n)]), &(maxvalues[cord_spec(i,j,n)]));
              //segstarts[cord_spec(counter++,j,n)] = i;
              if(debug){
                  Rprintf("segstarts[%d, %d] = %d\n",i, j, segstarts[cord_spec(i,j,n)] );
              }


          }
      }
    }




    cHDCD_call(X, -1, n-1, n, p, 1, changepoints,changepoint_counter_ptr,  depthcounter,
                thresholds , thresholds_test, cumsums, lens, lenLens,  as, nu_as, len_as,
                segstarts, maxvalues, maxpos, maxas, K, cusum, vector, coordchg,maxval, 
                startpoints, endpoints, maxaposes, tmpvec, twologn, ts,debug);

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

    SEXP ret = PROTECT(allocVector(VECSXP, 8)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, out);
    SET_VECTOR_ELT(ret, 1, maxvalSEXP);
    SET_VECTOR_ELT(ret, 2, depthcounterSEXP);
    SET_VECTOR_ELT(ret, 3, coordschgSEXP);
    SET_VECTOR_ELT(ret, 4, startpointsSEXP);
    SET_VECTOR_ELT(ret, 5, endpointsSEXP);
    SET_VECTOR_ELT(ret, 6, maxaposesSEXP);
    SET_VECTOR_ELT(ret, 7, changepointnumSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 8));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("changepoints"));
    SET_STRING_ELT(names, 1, mkChar("CUSUMval"));
    SET_STRING_ELT(names, 2, mkChar("depth"));
    SET_STRING_ELT(names, 3, mkChar("coordinate"));
    SET_STRING_ELT(names, 4, mkChar("startpoints"));
    SET_STRING_ELT(names, 5, mkChar("endpoints"));
    SET_STRING_ELT(names, 6, mkChar("maxaposes"));
    SET_STRING_ELT(names, 7, mkChar("changepointnumber"));

    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(32);
    return ret;
}



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


}
