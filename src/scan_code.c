#include "header.h"

/*void internal_colSum(double * matrix, int r1, int c1, double * vector){
    //memset(vector, 0, c1);
    for (int j = 0; j < c1; ++j)
    {
        vector[j]=0.0;
        for (int i = 0; i < r1; ++i)
        {
            vector[j]+=matrix[cord_spec(i, j, r1)];
        }
    }

}*/

void set_ind_vec(int * vector, int p){
    for (int i = 0; i < p; ++i)
    {
        vector[i]=i;
    }
}
void internal_scan_segment(double * cumsums, double * cusum, int * vector, double * vector2,int * maxpos, double * maximum, int * maxs,
                            int s, int e, int p, double * Ts, int debug){
    if(debug){
      Rprintf("checking segment (%d,%d]\n", s,e);

    }

    *maxpos = s+1;
    *maximum = NEGINF;
    *maxs = 1;


    CUSUM(cumsums, cusum, s, e, p);

    double cumsumm = 0;
    double tmp = 0;

    for (int j = 0; j < e-s-1; ++j)
    {
    	for (int i = 0; i < p; ++i)
    	{
    		cusum[cord_spec(i,j, p)] = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
    	}
    }

    for (int j = 0; j < e-s-1; ++j){
    	//R_rsort(&(cusum[cord_spec(0,j, p)]),p);
        //
        /*if(j<e-s-2){
            set_ind_vec(vector,p);
            memcpy(vector2, &(cusum[cord_spec(0,j+1, p)]), sizeof(double)*p );
        }*/
        R_qsort(&(cusum[cord_spec(0,j, p)]),1,p);
        //R_qsort_I(&(cusum[cord_spec(0,j, p)]), vector,1,p);
    	cumsumm = 0;

        /*if(j<e-s-2){
            for (int i = 0; i < p; ++i)
            {
                cumsumm +=cusum[cord_spec((p-1)-i,j, p)];
                tmp = (cumsumm - (i+1))/sqrt(2*(i+1))-Ts[i];
                if(tmp>*maximum){
                    *maximum = tmp;
                    *maxpos = s+j+1;
                    *maxs = i+1;
                }
                cusum[cord_spec(i,j+1, p)] = vector2[vector[i]];

            }
        }
        else{*/
        	for (int i = 0; i < p; ++i)
        	{
        		cumsumm +=cusum[cord_spec((p-1)-i,j, p)];
        		tmp = (cumsumm - (i+1))/sqrt(2*(i+1))-Ts[i];
        		if(tmp>*maximum){
        			*maximum = tmp;
        			*maxpos = s+j+1;
        			*maxs = i+1;
        		}
        	}

        //}

    }



    if(debug){
      Rprintf("for segment (%d,%d] we found maximum in %d with val %f\n", s,e,*maxpos, *maximum);
    }

    return;
}


void cScan_call(double * x, int s, int e, int n, int p, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter,
                double * Ts, double *cumsums, int* lens, int lenLens,
                int * segstarts, double * maxvalues, int* maxpos,int * maxss_all, int K, double * cusum, int* vector, double * vector2,
                int * coordchg, double * maxval, int * startpoints, int* endpoints, int* maxss,
                int debug){
    if(debug){
        Rprintf("cScan_call! s=%d, e=%d\n", s, e);
    }

    if(e-s<lens[0]){
        //Rprintf("segment too short\n");
        return;
    }
    int argmax = s+1;
    double maximum = NEGINF;
    int maxs = 0;

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
                 //internal_check_segment(cumsums, cusum, &(maxpos[cord_spec(k,j,n)]), &(maxvalues[cord_spec(k,j,n)]), &(maxss_all[cord_spec(k,j,n)]),
                  //          i, i+len, p, vector, Ts, as, nu_as,len_as, tmpvec, debug);
                 internal_scan_segment(cumsums, cusum, vector, vector2, &(maxpos[cord_spec(k,j,n)]), &(maxvalues[cord_spec(k,j,n)]), &(maxss_all[cord_spec(k,j,n)]),
                            i, i+len, p, Ts,  debug);
                 //internal_inspectOnSegment(cumsums, cusum, &(maxpos[cord_spec(k,j,n)]), &(maxvalues[cord_spec(k,j,n)]), i, i+len, p,
                 //lambda,
                //      eps, maxiter, mhat, mhatprod, v, v2,debug);
            }
            else{
              if(debug){
                Rprintf("segment (%d,%d] (k=%d, j=%d) already inspected, with max val %f in %d\n", i, i+len,k,j,
                maxvalues[cord_spec(k,j,n)],maxpos[cord_spec(k,j,n)]);
              }
            }

            tmp = maxvalues[cord_spec(k,j,n)];
            if(tmp>maximum){
                maximum = tmp;
                argmax = maxpos[cord_spec(k,j,n)];
                maxs = maxss_all[cord_spec(k,j,n)];
                j_max = j;
                k_max=k;
                //found=1;
            }


        }

        if(maximum>0.0){
            break;
        }
    }
    if(debug){
        Rprintf("maximum=%f\n", maximum);
    }

    if(maximum > 0.0){
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
        if(maxs==p){
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
            int jj = argmax-ss-1;
            for (int zz = 0; zz < p; ++zz)
            {
            	vector[zz] = zz;
            }

            //R_rsort(&(cusum[cord_spec(0,jj, p)]),p);
            rsort_with_index (&(cusum[cord_spec(0,jj, p)]), vector, p);

            int coord = 0;
            for (int zz = 0; zz < maxs; ++zz)
            {
            	//coord = vector[p-zz-1];
                coord = vector[zz];
            	coordchg[cord_spec(coord,*changepoint_counter_ptr, p)]=1;

            }
            

        }

        changepoints[*changepoint_counter_ptr] = argmax;
        depthcounter[*changepoint_counter_ptr] = depth;
        maxval[*changepoint_counter_ptr] = maximum;
        startpoints[*changepoint_counter_ptr] = segstarts[cord_spec(k_max,j_max,n)];
        endpoints[*changepoint_counter_ptr] = segstarts[cord_spec(k_max,j_max,n)] + lens[j_max];
        maxss[*changepoint_counter_ptr] = maxs;
        (*changepoint_counter_ptr)++;
        if(*changepoint_counter_ptr >n){
          return;
        }
        //cScan_call(double * x, int s, int e, int n, int p, int depth, int* changepoints,int* changepoint_counter_ptr, int * depthcounter,
        //        double threshold_d, double threshold_s , double *cumsums, int* lens, int lenLens, double * as, double * nu_as, int len_as,
        //        int * segstarts, double * maxvalues, int* maxpos, int K, double * cusum, double * matrix,
        //        double * vector, int * coordchg, int debug)
        cScan_call(x, s, argmax, n, p, depth+1, changepoints,changepoint_counter_ptr,  depthcounter,
                Ts , cumsums, lens, lenLens,  
                segstarts, maxvalues, maxpos,maxss_all, K, cusum, vector, vector2, coordchg, maxval, 
                startpoints, endpoints, maxss,debug);
        cScan_call(x, argmax, e,n, p, depth+1, changepoints,changepoint_counter_ptr,  depthcounter,
                Ts, cumsums, lens, lenLens,
                segstarts, maxvalues, maxpos,maxss_all, K, cusum, vector, vector2, coordchg, maxval,
                startpoints, endpoints, maxss, debug);

    }

    return;
}

SEXP cScan(SEXP XI,SEXP nI, SEXP pI,SEXP TsI, SEXP lensI,SEXP lenLensI,SEXP KI, SEXP debugI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(TsI);
    PROTECT(lensI);
    PROTECT(lenLensI);
    PROTECT(KI);
    PROTECT(debugI);

    double * X = REAL(XI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    double *Ts = (REAL(TsI));
    int *lens = INTEGER(lensI); /////// dobbeltsjekk at int er like stor som INTEGER!!!!!!!
    int lenLens = *(INTEGER(lenLensI));
    int K = *(INTEGER(KI));
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
    SEXP maxssSEXP = PROTECT(allocVector(INTSXP, n));
    int * maxss = INTEGER(maxssSEXP); //pointer to array
    memset(maxss, 0, sizeof(int)*n);
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

    SEXP vectorSEXP = PROTECT(allocVector(INTSXP, p));
    int * vector = INTEGER(vectorSEXP);
    memset(vector, 0, sizeof(int)*p);

    SEXP vector2SEXP = PROTECT(allocVector(REALSXP, p));
    double * vector2 = REAL(vector2SEXP);
    memset(vector2, 0, sizeof(double)*p);




    // record all segment starts
    // we don't compute all cusums here but do it on the fly so we dont have to do anything
    // unecessary


    SEXP maxvaluesSEXP = PROTECT(allocVector(REALSXP, n*lenLens));
    double * maxvalues = REAL(maxvaluesSEXP);
    for (int i = 0; i < n*lenLens; ++i)
    {
        maxvalues[i] = NEGINF;
    }
    SEXP maxposSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * maxpos = INTEGER(maxposSEXP);
    memset(maxpos, 0, sizeof(int)*n*lenLens);
    SEXP maxss_allSEXP = PROTECT(allocVector(INTSXP, n*lenLens));
    int * maxss_all = INTEGER(maxss_allSEXP);
    memset(maxss_all, 0, sizeof(int)*n*lenLens);
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




    cScan_call(X, -1, n-1, n, p, 1, changepoints,changepoint_counter_ptr,  depthcounter,
                Ts , cumsums, lens, lenLens,
                segstarts, maxvalues, maxpos, maxss_all, K, cusum, vector, vector2, coordchg,maxval, 
                startpoints, endpoints, maxss,  debug);

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
    SET_VECTOR_ELT(ret, 6, maxssSEXP);
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
    SET_STRING_ELT(names, 6, mkChar("maxss"));
    SET_STRING_ELT(names, 7, mkChar("changepointnumber"));

    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(26);
    return ret;
}





SEXP cScan_single(SEXP XI,SEXP nI, SEXP pI, SEXP debugI, SEXP TsI){
    // X : p \times n
    PROTECT(XI);
    PROTECT(nI);
    PROTECT(pI);
    PROTECT(debugI);
    PROTECT(TsI);


    double * X = REAL(XI);
    double * Ts = REAL(TsI);
    int n = *(INTEGER(nI));
    int p = *(INTEGER(pI));
    
    int debug = *INTEGER(debugI);


    if(debug){
      Rprintf("p = %d\n", p);
      Rprintf("n = %d\n", n);
    }


    /*SEXP tmpvecSEXP = PROTECT(allocVector(REALSXP, len_as));
    double * tmpvec = REAL(tmpvecSEXP);
    memset(tmpvec, 0, sizeof(double)*len_as);*/

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
    SEXP maxsSEXP = PROTECT(allocVector(INTSXP, 1));
    int * maxpos = INTEGER(maxposSEXP);
    int * maxs= INTEGER(maxsSEXP);

    *maxpos = -10;
    double maximum = -1000000000000000000000.0;
    *maxs = 0;
    int s = -1;
    int e = n-1;

    CUSUM(cumsums, cusum, s, e, p);

    double cumsumm = 0;
    double tmp = 0;

    for (int j = 0; j < e-s-1; ++j)
    {
    	for (int i = 0; i < p; ++i)
    	{
    		cusum[cord_spec(i,j, p)] = cusum[cord_spec(i,j, p)]*cusum[cord_spec(i,j, p)];
    	}
    }

    for (int j = 0; j < e-s-1; ++j){
    	R_rsort(&(cusum[cord_spec(0,j, p)]),p);
    	cumsumm = 0;
    	for (int i = 0; i < p; ++i)
    	{
    		/*if (i>0)
    		{
    			if(cusum[cord_spec(i-1,j, p)]>cusum[cord_spec(i,j, p)]){
    				Rprintf("ERROR");
    			}
    		}*/
    		cumsumm +=cusum[cord_spec((p-1)-i,j, p)];
    		tmp = (cumsumm - (i+1))/sqrt(2*(i+1))-Ts[i];
    		if(tmp>maximum){
    			maximum = tmp;
    			*maxpos = j;
    			*maxs = i+1;
    		}
    	}

    }



    SEXP ret = PROTECT(allocVector(VECSXP, 2)); //the list to be returned in R
    SET_VECTOR_ELT(ret, 0, maxposSEXP);
    SET_VECTOR_ELT(ret, 1, maxsSEXP);

    // creating list of names/titles to be returned in the output list
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    //SET_STRING_ELT(names, 0, mkChar("CUSUM"));
    SET_STRING_ELT(names, 0, mkChar("pos"));
    SET_STRING_ELT(names, 1, mkChar("s"));


    setAttrib(ret, R_NamesSymbol, names);

    UNPROTECT(11);
    return ret;


}