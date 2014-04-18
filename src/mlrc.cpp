#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
/* #include <conio.h> */
#include <float.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "mat.h"
#include "mlrc.h"

#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

#define LOGITTOL 1.0E-12

double *pcom = 0, *xicom = 0;
int ncom = 0;
double (*nrfunc)(double *, dMat &, dMat &);

#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)
static double sqrarg;

/*
#ifdef _WIN32
#define finite _finite
#endiff
*/

//#define MAX(X,Y) ((X) > (Y) ? : (X) : (Y))

//int powell(double p[], double **xi, int n, double ftol, int *iter, double *fret, dMat & params, dMat &,  double(*func)(double *, dMat &, dMat &));

int logit(dMat &x, dMat &y, dMat &bhat, dMat &mmxinv, double tol, int maxiter, int verbose);

extern "C" {
/*
   __declspec(dllexport) 
*/   
   SEXP MLRC_regress(SEXP sexp_SpecData, SEXP sexp_Env, SEXP sexp_miter, SEXP sexp_verbose) 
{
   SEXP dims, R_Beta, R_IBeta, ret, retNames;
   dims = Rf_getAttrib(sexp_SpecData, R_DimSymbol);
   int nr = INTEGER(dims)[0];
   int nc = INTEGER(dims)[1];
   int verbose = INTEGER(sexp_verbose)[0];
   int maxiter = INTEGER(sexp_miter)[0];
   int i=0, j=0;

   dMat beta(nc, 3, 0.0);
   dMat Y(nr, nc, 0.0);
   dMat X(nr, 1, 0.0);
   PROTECT(sexp_SpecData);
   for (i=0;i<nr;i++) {
      for (int j=0;j<nc;j++) {
         Y(i,j) = REAL(sexp_SpecData)[i + nr*j];
      }
   }
   UNPROTECT(1);
   PROTECT(sexp_Env);
   for (i=0;i<nr;i++) {
      X(i,0) = REAL(sexp_Env)[i];
   }
   UNPROTECT(1);

   Index IBeta(nc, -1); 
 	 dMat bhat(3,1,1.0);
   dMat mmxinv;
 	 dMat ones(nr,1,1.0);
 	 dMat x2 = ones.concat(X.concat(X * X,ColWise), ColWise);
   dMat sp(nr, 1, 0.0);

   PROTECT(R_Beta = allocVector(REALSXP, nc*3));
   PROTECT(R_IBeta = allocVector(INTSXP, nc));
	 for (i=0;i<nc;i++) {
      for (j=0;j<nr;j++)
         sp(j,0) = Y(j,i);
      try {
         IBeta(i) = logit(x2, sp, bhat, mmxinv, LOGITTOL, maxiter, verbose);
      }
      catch (char *ex) {
         if (verbose) {
            REprintf("\n%s\n", ex);
         }
         IBeta(i) = -4;
      }
      if (IBeta(i) > -1 && IBeta(i) < maxiter) {
         REAL(R_Beta)[i] = -bhat(0,0);
         REAL(R_Beta)[i+nc] = -bhat(1,0);
         REAL(R_Beta)[i+nc*2] = -bhat(2,0);
      }
      else {
         REAL(R_Beta)[i] = NA_REAL;
         REAL(R_Beta)[i+nc] = NA_REAL;
         REAL(R_Beta)[i+nc*2] = NA_REAL;
      }
      INTEGER(R_IBeta)[i] = IBeta(i);
   }

   UNPROTECT(2);

   PROTECT(ret = allocVector(VECSXP, 2)); 
   PROTECT(retNames = allocVector(STRSXP, 2));
   SET_VECTOR_ELT(ret, 0, R_Beta);
   SET_VECTOR_ELT(ret, 1, R_IBeta);
   SET_STRING_ELT(retNames, 0, mkChar("Beta"));
   SET_STRING_ELT(retNames, 1, mkChar("IBeta"));
   SET_NAMES(ret, retNames);
   UNPROTECT(2);
   return(ret);
}
}

int logit(dMat &x, dMat &y, dMat &bhat, dMat &mmxinv, double tol, int maxiter, int verbose)
{
	dMat ys, r, p, delta;
	bhat = dMat(cols(x),1, 0.0);
	int iter=0;
	double maxii;
   char errorflag = 0;
	do {
    dMat t;
		ys = x.product(bhat);
    for (int j=0; j<rows(ys);j++) {
      if (ys(j, 0) < -300)
         ys(j,0) = -300;
      else if (ys(j, 0) > 300)
         ys(j,0) = 300;
    }
		r = 1.0 / (1.0 + exp(ys));
	  p = (r*(1.0-r));

//	   mmxinv = ((x*p).transpose().product(x)).inverse(errorflag);
    try {
      mmxinv = ((x).tproduct(x*p)).inverse(errorflag);
    }
    catch (const char *e) {
       if (verbose) {
          REprintf("\n%s\n", e);
       }
       errorflag = 3;
    }
    if (errorflag)
     	break;
    delta = mmxinv.product((x).tproduct(y-r));
		bhat = bhat - delta;
		double maxi, mini;
		delta /= bhat;
    maxmin(delta, mini, maxi);
		maxii = MAX(fabs(maxi), fabs(mini));
		iter++;
		if (iter == maxiter)
			break;
	} while (maxii >= tol);
  if (errorflag)
   	return -2;
	return iter;
}

extern "C" {
/*
__declspec(dllexport) 
*/
SEXP MLRC_predict(SEXP sexp_SpecData, SEXP sexp_Beta, SEXP sexp_meanX)
{
   SEXP dims, R_pred;
   dims = Rf_getAttrib(sexp_SpecData, R_DimSymbol);
   int nr = INTEGER(dims)[0];
   int nc = INTEGER(dims)[1];
   double meanX = REAL(sexp_meanX)[0];
   double **xi;
   double *p, fret;
   int iter;
   int nDimen = 1;
   int i=0, j=0;
   xi = matrix(1,nDimen+10,1,nDimen+10);
   p = vector(1,nDimen+5);
   dMat beta(nc, 3, 0.0);
   dMat Y(nr, nc, 0.0);
   PROTECT(sexp_SpecData);
   for (i=0;i<nr;i++) {
      for (int j=0;j<nc;j++) {
         Y(i,j) = REAL(sexp_SpecData)[i + nr*j];
      }
   }
   UNPROTECT(1);
   PROTECT(sexp_Beta);
   for (i=0;i<nc;i++) {
      for (j=0;j<3;j++) {
// check that this is -      
//         beta(i,j) = -REAL(sexp_Beta)[i + nc*j];
         beta(i,j) = REAL(sexp_Beta)[i + nc*j];
      }
   }
   UNPROTECT(1);
   PROTECT(R_pred = allocVector(REALSXP, nr));
   dMat sp(nc, 1, 0.0);
   for (i=0;i<nr;i++) {   
      for (j=0;j<nc;j++) {
         if (ISNA(beta(j, 0)))
            sp(j, 0) = -1;
         else
            sp(j, 0) = Y(i,j);
      }
      xi[1][1] = 1.0;
//         p[1] = -100.0;
      p[1] = meanX;
      p[2] = p[3] = 0.0;
      bool retval = FALSE;
      try {
          powell(p, xi, nDimen, FTOL, &iter, &fret, beta, sp, calib_func);
       }
       catch (const char *e) {
          REprintf("\n%s\n", e);
          retval = TRUE;
       }
       if (retval) {
          REAL(R_pred)[i] = NA_REAL;
//          REAL(R_pred)[i+nr] = NA_REAL;
       }
       else {
          REAL(R_pred)[i] = p[1];
//          REAL(R_pred)[i+nr] = fret;
       }
   }
   UNPROTECT(1);
   free_matrix(xi, 1, nDimen+10, 1, nDimen+10);
   free_vector(p, 1, nDimen+5);
   return(R_pred);
}
}


double calib_func(double *xt, dMat &params, dMat &SpecData)
{
	int j;
	double pred, prob, like;

	like = 0.0;
   int nsp = rows(params);

	for (j=0;j<nsp;j++) {
      if (SpecData(j, 0) < 0.0)
         continue;
      pred = (params(j, 0)+params(j, 1)*xt[1]+params(j, 2)*xt[1]*xt[1]);
      if (pred > 50) pred = 50;
	   if (pred < -50) pred = -50;
      prob = 1/(1+exp(-pred));
	   like += SpecData(j, 0)*log(prob)+(1.0-SpecData(j, 0))*log(1.0-prob+TINY);
   }
	return(-like);
}

double *vector(int nl, int nh)
{
   double *v;
   v = (double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
   if (!v) 
      throw("allocation failure in vector()");
   return v-nl;
}

void free_vector(double *v, int nl, int nh)
{
   free((char*) (v+nl));
}

double **matrix(int nrl, int nrh, int ncl, int nch)
{
   int i;
   double **m;
   m = (double **) malloc((unsigned)(nrh-nrl+1)*sizeof(double*));
   if (!m) 
      throw("allocation failure 1 in matrix");
   m -= nrl;
   for(i=nrl;i<=nrh;i++) {
      m[i]=(double *) malloc((unsigned)(nch-ncl+1)*sizeof(double));
      if (!m[i]) throw("allocation failure 2 in matrix");
      m[i] -= ncl;
   }
   return(m);
}

void free_matrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;
	for (i=nrh;i>=nrl;i--) free((char*)(m[i]+ncl));
	free((char*)(m+nrl));
}


void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
		double *fc, dMat &params, dMat &SpecData, double (*func)(double, dMat &, dMat &))
{
   double ulim,u,r,q,fu,dum;
   int iter = 0;
   *fa = (*func)(*ax, params, SpecData);
/*   printf("\nfa = %f", *fa);*/
   *fb = (*func)(*bx, params, SpecData);
/*   printf("\nfb = %f", *fb);*/
   if (*fb > *fa) {
      SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
   }
   *cx = (*bx)+GOLD*(*bx-*ax);
   *fc = (*func)(*cx, params, SpecData);
/*   printf("\nfc = %f", *fc);*/
   while (*fb > *fc) {
      iter++;
      r = (*bx-*ax)*(*fb-*fc);
      q = (*bx-*cx)*(*fb-*fa);
      u = (*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
            (2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
      ulim = (*bx)+GLIMIT*(*cx-*bx);
      if ((*bx-u)*(u-*cx) > 0.0) {
         fu = (*func)(u, params, SpecData);
/*	 printf("\nfu = %f", fu);*/
         if (fu < *fc) {
            *ax = (*bx);
            *bx = u;
            *fa = (*fb);
            *fb = fu;
            return;
         } else if (fu > *fb) {
            *cx = u;
            *fc = fu;
            return;
         }
         u = (*cx)+GOLD*(*cx-*bx);
         fu =(*func)(u, params, SpecData);
      } 
      else if ((*cx-u)*(u-ulim) > 0.0) {
         fu = (*func)(u, params, SpecData);
         if (fu < *fc) {
            SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
            SHFT(*fb,*fc,fu,(*func)(u, params, SpecData))
         }
      } 
      else if ((u-ulim)*(ulim-*cx) >= 0.0) {
         u = ulim;
         fu = (*func)(u, params, SpecData);
      } 
      else {
         u = (*cx)+GOLD*(*cx-*bx);
         fu = (*func)(u, params, SpecData);
      }
      SHFT(*ax,*bx,*cx,u)
      SHFT(*fa,*fb,*fc,fu)
      if (iter > 20000)
         throw("Too many iterations in mnbrak");
   }
}

void powell(double *p, double **xi, int n, double ftol, int *iterate,
            double *fret, dMat &params, dMat &SpecData, double (*func)(double *, dMat &, dMat &))
{
   int errorflag = 0;
   int i,ibig,j;
   double t,fptt,fp,del;
   double *pt,*ptt,*xit;

   pt = vector(1,n+5);
   ptt = vector(1,n+5);
   xit = vector(1,n+5);
   *fret = (*func)(p, params, SpecData);
   for (j = 1;j <= n; j++) pt[j] = p[j];
   *iterate = 0;
   while (1) {
      (*iterate)++;
/*
      printf("\niteration %d  ", *iterate);
      if (n == 1) 
         printf("b0 = %-.3f, likelihood = %-.9f", p[1], *fret);
      else  
         printf("b0 = %-.3f, b1 = %-.3f, b2 = %-.3f, likelihood = %-.9f",p[1],p[2],p[3],*fret);
*/
      fp = (*fret);
      ibig = 0;
      del = 0.0;
      for (i = 1; i <= n; i++) {
         for (j = 1; j <= n; j++) xit[j] = xi[j][i];
         fptt = (*fret);
         linmin(p,xit,n,fret,params, SpecData, func);
/*
         if (_finite(*fret) == FALSE) {
*/         
/*         if (finite(*fret) == FALSE) {  */
           if (R_FINITE(*fret) == FALSE) {
            errorflag = 1;
            throw ("NAN in routine brent");
         }
         if (fabs(fptt-(*fret)) >= del) {
            del = fabs(fptt-(*fret));
            ibig = i;
         }
      }
      if (ibig < 1) {
         errorflag = 1;
         throw("Error in routine POWELL");
      }
      if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
         free_vector(xit,1,n+5);
         free_vector(ptt,1,n+5);
         free_vector(pt,1,n+5);
         return;
//         return errorflag;
      }
      if (*iterate == ITMAX2) {
         errorflag = 1;
         throw("Too many iterations in routine POWELL");
//         nrerror("Too many iterations in routine POWELL");
      }
      for (j=1;j<=n;j++) {
         ptt[j] = 2.0*p[j]-pt[j];
         xit[j] = p[j]-pt[j];
         pt[j] = p[j];
      }
      fptt =(*func)(ptt, params, SpecData);
      if (fptt < fp) {
//         t = 2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
         t = 2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
         
         if (t < 0.0) {
            linmin(p,xit,n,fret,params, SpecData, func);
            for (j=1;j<=n;j++) {
               xi[j][ibig] = xit[j];
            }
         }
      }
   }
//   return errorflag;
}

void linmin(double *p, double *xi, int n, double *fret, dMat &params, dMat &SpecData, double(*func)(double *, dMat &, dMat &))
{
   int j;
   double xx, xmin,fx,fb,fa,bx,ax;

   ncom = n;
   pcom = vector(1,n+5);
   xicom = vector(1,n+5);
   nrfunc = func;
   for (j= 1; j<=n;j++) {
      pcom[j] = p[j];
      xicom[j] = xi[j];
   }
   ax = 0.0;
   xx = 1.0;
   bx = 2.0;
   mnbrak(&ax,&xx,&bx,&fa,&fx,&fb, params, SpecData, f1dim);
   *fret = brent(ax,xx,bx,params, SpecData, f1dim,TOL,&xmin);
   for (j=1;j<=n;j++) {
      xi[j] *= xmin;
      p[j] += xi[j];
   }
   free_vector(xicom,1,n+5);
   free_vector(pcom,1,n+5);
}

double f1dim(double x, dMat &params, dMat &SpecData)
{
   int j;
   double f,*xt;

   xt = vector(1,ncom);
   for (j=1;j<=ncom;j++) 
      xt[j]=pcom[j]+x*xicom[j];
   f = (*nrfunc)(xt, params, SpecData);
   free_vector(xt,1,ncom);
   return f;
}

double brent(double ax, double bx, double cx, dMat &params, dMat &SpecData, double (*f)(double, dMat&, dMat &), double tol, double *xmin)

{
   int iter;
   double a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
   double e=0.0;
   double d = 0.0;

   a = ((ax < cx) ? ax : cx);
   b = ((ax > cx) ? ax : cx);
   x = w = v = bx;
   fw = fv = fx = (*f)(x, params, SpecData);
   for (iter = 1; iter <= ITMAX1; iter++) {
      xm = 0.5*(a+b);
      tol2 = 2.0*(tol1 = tol*fabs(x)+ZEPS);
      if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
         *xmin = x;
         return fx;
      }
      if (fabs(e) > tol1) {
         r = (x-w)*(fx-fv);
         q = (x-v)*(fx-fw);
         p = (x-v)*q-(x-w)*r;
         q = 2.0*(q-r);
         if ( q > 0.0) p = -p;
         q = fabs(q);
         etemp = e;
         e = d;
         if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
            d = CGOLD*(e=(x >=xm ? a-x : b-x));
         else {
            d = p/q;
            u = x+d;
            if (u-a < tol2 || b-u < tol2)
               d = SIGN(tol1,xm-x);
         }
      } else {
         d = CGOLD*(e=(x >=xm ? a-x : b-x));
      }
      u = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      fu = (*f)(u, params, SpecData);
      if (fu <= fx) {
         if (u >= x) a = x; else b = x;
         SHFT(v,w,x,u)
         SHFT(fv,fw,fx,fu)
      } else {
         if (u < x) a = u; else b = u;
         if (fu <= fw || w == x) {
            v = w;
            w = u;
            fv = fw;
            fw = fu;
         } else if (fu <= fv || ((v == x) && (v == w))) {
            v = u;
            fv = fu;
         }
      }
   }
   throw("Too many iterations in BRENT");
   *xmin = x;
   return fx;
}

