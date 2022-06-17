#define cord(r,c) ((r) + (DEPTH)*(c))
#define cord_spec(r,c, D) ((r) + (D)*(c))
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/RS.h>	// RK addition
#include <R_ext/Lapack.h> // RK addition
#include <R_ext/BLAS.h> // RK addition
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <float.h>

#define NEGINF -(FLT_MAX -10.0)

void internal_matmult(double *,
                      double *,
                      double *,
                      int, int, int, int);

void internal_matmultrightT(double *,
                            double *,
                            int , int );

void internal_matmultleftT(double * ,
                           double * ,
                           int , int );

double * internal_power_method(double *, int, double,
                               int , double * , double * ,int  );
  

void CUSUM(double *,double * , int , int , int );

// Sorting
void insertSort (double * , int , int );
void sort_k_largest(double *, int, int, int);