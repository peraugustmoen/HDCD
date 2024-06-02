#include "header.h"


	
void insertSort (double * a, int v, int h) {
	
	int i;
	double t;
	for ( int k = v ; k < h ; k++) {
		t = a [k+1] ;
		i=k;
		while ( i >= v && a [i]*a[i] < t*t ) {
			a [i+1] = a [i];
			i--; 
		}
		a[i+1] = t; 
	} 
	
} 


void insertSort_int(int * a, int v, int h) {
	
	int i;
	int t;
	for ( int k = v ; k < h ; k++) {
		t = a [k+1] ;
		i=k;
		while ( i >= v && a [i]*a[i] < t*t ) {
			a [i+1] = a [i];
			i--; 
		}
		a[i+1] = t; 
	} 
	
} 


/**
*Sorts the array a in accordance with algoritm A2
*
*@param a 	Array to be sorted
*@param k 	Number of largest elements to be found
**/

void sort_k_largest_abs(double * a, int k, int start, int stop){
	insertSort(a, start,start+ k-1);
	int i;
	double t ;

	for ( int j = k + start; j < stop ; j++) {
		if (fabs(a[j]) > fabs(a[k + start-1]) ) {
			
			t = a[j];
			a[j] = a[k+start-1];
			i=k-2 + start;
			while ( i >= start && fabs(a [i]) < fabs(t) ) {
				a [i+1] = a [i];
				i--; 
			}
			a[i+1] = t; 	
		}
	}
}

void sort_k_largest(double * a, int k, int start, int stop){
	insertSort(a, start,start+ k-1);
	int i;
	double t ;

	for ( int j = k + start; j < stop ; j++) {
		if (a[j] > a[k + start-1] ) {
			
			t = a[j];
			a[j] = a[k+start-1];
			i=k-2 + start;
			while ( i >= start && a [i] < t ) {
				a [i+1] = a [i];
				i--; 
			}
			a[i+1] = t; 	
		}
	}
}

void sort_k_largest_int(int* a, int k, int start, int stop){
	insertSort_int(a, start,start+ k-1);
	int i;
	int t ;

	for ( int j = k + start; j < stop ; j++) {
		if (a[j] > a[k + start-1] ) {
			
			t = a[j];
			a[j] = a[k+start-1];
			i=k-2 + start;
			while ( i >= start && a [i] < t ) {
				a [i+1] = a [i];
				i--; 
			}
			a[i+1] = t; 	
		}
	}
}


SEXP sort_k_largest_R(SEXP vecI, SEXP kI, SEXP startI, SEXP stopI){
	PROTECT(vecI);
	PROTECT(startI);
	PROTECT(kI);
	PROTECT(stopI);
	double * vec = REAL(vecI);
	int k = *(INTEGER(kI));
	int start = *(INTEGER(startI));
	int stop = *(INTEGER(stopI));
	sort_k_largest(vec,k, start, stop);
	UNPROTECT(4);
	return(vecI);
}





//function to swap variable
void swap(double* a, double* b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}

int partition (double * arr, int left, int right)
{
    double pivot = fabs(arr[right]); // pivot
    int i = (left - 1);
    int j = left;
    for (j = left; j <= (right - 1); j++)
    {
        
        if (fabs(arr[j]) > pivot)
        {
            
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[right]);
    return (i + 1);
}

int quickselect(double * arr, int left, int right, int k)
{
    
    if (k > 0 && k <= (right - left + 1))
    {
        
        int index = partition(arr, left, right);
        
        if (index - left == k - 1)
            return arr[index];
        
        if (index - left > k - 1)
            return quickselect(arr, left, index - 1, k);
        
        return quickselect(arr, index + 1, right,
                           k - index + left - 1);
    }
    return -1;
}

void rec_partial_quicksort(double * A, int i, int j, int m) {
	if (i<j){
		int index = partition(A, i, j);
		rec_partial_quicksort(A, i, index - 1, m);
		if (index < m - 1){
			rec_partial_quicksort(A, index + 1, j, m);
		}
	}
}


void partial_quicksort(double * A, int len, int k){
	rec_partial_quicksort(A, 0, len-1, k);
}


SEXP partial_quicksort_R(SEXP vecI, SEXP kI, SEXP lenI){
	PROTECT(vecI);
	PROTECT(kI);
	PROTECT(lenI);
	double * vec = REAL(vecI);
	int k = *(INTEGER(kI));
	int len = *(INTEGER(lenI));
	partial_quicksort(vec, len,k);
	UNPROTECT(3);
	return(vecI);
}

SEXP sort_test(SEXP xI, SEXP nI){
	PROTECT(xI);
	PROTECT(nI);

	double * x = REAL(xI);
	int n = *(INTEGER(nI));

	R_qsort(x,1,n);
	UNPROTECT(2);
	return(xI);

}

SEXP partialsum_test(SEXP xI, SEXP nI){
	PROTECT(xI);
	PROTECT(nI);

	double * x = REAL(xI);
	int n = *(INTEGER(nI));
	
	SEXP resSEXP= PROTECT(allocVector(REALSXP, n));
	double * res = REAL(resSEXP);
	memset(res, 0, sizeof(double)*n);
	

	R_qsort(x,1,n);

	int prev = n;
  int z = 1;
  int c = 0;
  int i = 0;
  int j=0;
  double cumsum=0;
  while(1)
  {
      if(z>n){
          break;
      }
      for (i = prev-1; i >= n-z; --i)
      {
          cumsum += x[cord_spec(i,j, n)]*x[cord_spec(i,j, n)];
      }
      
      res[c] = cumsum;
      prev = n-z;
      c++;
      z = 2*z;


  }

	UNPROTECT(3);
	return(resSEXP);

}



