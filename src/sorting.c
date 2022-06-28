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





//function to swap variable
void swap(double* a, double* b)
{
    double tmp = *a;
    *a = *b;
    *b = tmp;
}
/* partition function takes last element as pivot and places
   the pivot element at its correct position. It means all
   smaller element will be placed to left all greater elements
   to right of pivot
 */
int partition (double * arr, int left, int right)
{
    double pivot = fabs(arr[right]); // pivot
    int i = (left - 1);
    int j = left;
    for (j = left; j <= (right - 1); j++)
    {
        // If current element is smaller than the pivot
        if (fabs(arr[j]) > pivot)
        {
            i++; // increment index of smaller element
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[right]);
    return (i + 1);
}
// Function returns the k'th smallest
//element in the arr within `left…right`
// (i.e., `left <= k <= right`).
int quickselect(double * arr, int left, int right, int k)
{
    // If k is smaller than number of
    // elements in array
    if (k > 0 && k <= (right - left + 1))
    {
        // Partition the array around last
        // element and get position of pivot
        // element in sorted array
        int index = partition(arr, left, right);
        // If position is same as k
        if (index - left == k - 1)
            return arr[index];
        // If position is more, recur
        // for left subarray
        if (index - left > k - 1)
            return quickselect(arr, left, index - 1, k);
        // Else recur for right subarray
        return quickselect(arr, index + 1, right,
                           k - index + left - 1);
    }
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



