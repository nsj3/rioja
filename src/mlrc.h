#ifndef MLRCHPP
#define MLRCSHPP 1

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
/* #include <malloc.h> */
#include "mat.h"

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define ITMAX1 100
#define ITMAX2 200
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define TOL 2.0e-5
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)
#define FTOL 1.0e-11

double f1dim(double x, dMat &, dMat &);
void linmin(double *p, double *xi, int n, double *fret, dMat &params, dMat &SpecData, double(*func)(double *, dMat &, dMat &));
void powell(double *p, double **xi, int n, double ftol, int *iter, double *fret, dMat &, dMat &,  double(*func)(double *, dMat &, dMat &));
double brent(double ax, double bx, double cx, dMat &, dMat &, double (*f)(double, dMat &, dMat &), double tol, double *xmin);
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,	double *fc, dMat &params, dMat &, double (*func)(double *, dMat &, dMat &));

double calib_func(double *xt, dMat &params, dMat &SpecData);

double **matrix(int, int, int, int);
void free_matrix(double **m, int nrl, int nrh, int ncl, int nch);
double *vector(int, int);
void free_vector(double *v, int nl, int nh);

#endif
