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


/**
*Sorts the array a in accordance with algoritm A2
*
*@param a 	Array to be sorted
*@param k 	Number of largest elements to be found
**/

void sort_k_largest(double * a, int k, int start, int stop){
	insertSort(a, start,start+ k-1);
	int i;
	double t ;

	for ( int j = k + start; j < stop ; j++) {
		if (a[j]*a[j] > a[k + start-1]*a[k + start-1] ) {
			
			t = a[j];
			a[j] = a[k+start-1];
			i=k-2 + start;
			while ( i >= start && a [i]*a[i] < t*t ) {
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