#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "mat.h"

#ifndef max
   #define max(a,b) ((a)>(b)?(a):(b))
#endif

#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

extern "C" {
/*
__declspec(dllexport) 
*/

   SEXP WA_fit(SEXP sexp_SpecData, SEXP sexp_EnvData, SEXP sexp_MinTol, SEXP sexp_Lean)
{
   SEXP dims, retNames=R_NilValue;
   dims = Rf_getAttrib(sexp_SpecData, R_DimSymbol);
   int nr = INTEGER(dims)[0];
   int nc = INTEGER(dims)[1];
   double minTol = REAL(sexp_MinTol)[0];
//   double minTolReplace = REAL(sexp_MinTolReplace)[0]
//   bool bLean = bool(INTEGER(sexp_Lean)[0]);
   double smallValue = 1.0E-4;
   bool AdjustTolerancesN2 = true;

//   char *str = new char[1000]; 
   int i=0, j=0;

   dMat X(nr, nc, 0.0);
   dMat Y(nr, 1, 0.0);
   PROTECT(sexp_SpecData);
   PROTECT(sexp_EnvData);
   for (i=0;i<nr;i++) {
      Y(i,0) = REAL(sexp_EnvData)[i];
      for (int j=0;j<nc;j++) {
         X(i,j) = REAL(sexp_SpecData)[i + nr*j];
      }
   }
   UNPROTECT(2);

   dMat opt(nc, 2, 0.0);
   dMat sum3(nc, 1, 0.0);
   dMat est(nr, 1, 0.0);
   dMat estTol(nr, 1, 0.0);
   Index ICount(nc, 0);

   double **mm = dataptr(X);
   for (i=0;i<nr;i++) { 
      double env = Y(i, 0);
      for (j=0;j<nc; j++) {
         opt(j, 0) += mm[i][j] * env;
         sum3(j, 0) += mm[i][j];
         if (mm[i][j] > 1.0E-6)
            ICount(j) += 1;
      }
   }
   for (i=0;i<nc;i++) {
      if (ICount(i)<1) {
         sum3(i, 0) = 1.0;
         opt(i, 0) = -99.9;
      }
   }
   for (i=0;i<nc;i++) {
      opt(i,0) /= sum3(i,0);
   }
   for (i=0;i<nr;i++) { 
      double env = Y(i, 0);
      for (j=0;j<nc; j++) {
         if (ICount(j)>0) {
            double d = env - opt(j, 0);
            opt(j, 1) += mm[i][j] * d * d;
         }
      }
   }
   for (j=0;j<nc; j++) {
   	opt(j, 1) /= sum3(j,0);
   }
   dMat N2(nc, 1, 0.0);

   if (AdjustTolerancesN2) {
// Calc N2;
      for (i=0;i<nr;i++) { 
         for (j=0;j<nc; j++) {
				double d = mm[i][j] / sum3(j, 0);
				N2(j, 0) += d*d;
         }
      }
		for (i=0;i<nc;i++) {
         if (ICount(i) > 0) {
				N2(i, 0) = 1.0 / N2(i, 0);
   			opt(i, 1) /= max(1.0 - 1.0 / N2(i, 0), smallValue);
         }
		}
   }

   for (i=0;i<nc;i++) {
      opt(i,1) = sqrt(opt(i,1));
   }
   double meanTol = 0.0;
   double tolSum = 0.0;
   double tolN = 0.0;

   for (i=0;i<nc;i++) {
      if (ICount(i) > 0) {
         if (opt(i,1) > minTol) {
            tolSum += opt(i,1);
            tolN += 1.0;
         }
      }
   }
   meanTol = tolSum / tolN;

   for (i=0;i<nc;i++) {
      if (ICount(i) == 0) {
         opt(i, 0) = missingValue(opt);
         opt(i, 1) = missingValue(opt);
      }
      else if (opt(i, 1) <= minTol) {
         opt(i, 1) = meanTol;
      }
   }

   for (i=0;i<nr;i++) { 
      double sum2=0.0;
      double sum2Tol=0.0;
      for (j=0;j<nc; j++) {
         if (ICount(j)) {
            double t = opt(j, 1);
            t = t * t;
            sum2 += mm[i][j];
            sum2Tol += mm[i][j] / t;
            est(i,0) += mm[i][j] * opt(j, 0);
            estTol(i,0) += mm[i][j] * opt(j, 0) / t;
         }
      }
      est(i, 0) /= sum2;
      estTol(i, 0) /= sum2Tol;
   }

   double meanest = mean(est);
   dMat xx = est - meanest;
   double meany = mean(Y);
   dMat yy = Y - meany;
   double tb0inv=0.0, tb1inv = 0.0, tb0cla = 0.0, tb1cla = 0.0;

   for (i=0;i<nr;i++) {
      tb1inv += xx(i,0) * yy(i,0);                 
   }
   tb1cla = tb1inv;
   double ssqx = sumsq(xx);
   tb1inv = tb1inv / ssqx;
	tb0inv = meany - tb1inv * meanest;
  	dMat estinv = tb0inv + (tb1inv * est);
   double ssqy = sumsq(yy);
   tb1cla = tb1cla / ssqy;
	tb0cla = meanest - tb1cla * meany;
//  	dMat estcla = (est - tb0cla)  / tb1cla;
//   est = estinv.concat(estcla, ColWise);

   double meanestTol = mean(estTol);
   dMat xxTol = estTol - meanestTol;
   double tb0invTol=0.0, tb1invTol = 0.0, tb0claTol = 0.0, tb1claTol = 0.0;
   for (i=0;i<nr;i++) {
      tb1invTol += xxTol(i,0) * yy(i,0);                 
   }
   tb1claTol = tb1invTol;
   double ssqxTol = sumsq(xxTol);
   tb1invTol = tb1invTol / ssqxTol;
	tb0invTol = meany - tb1invTol * meanestTol;
  	dMat estinvTol = tb0invTol + (tb1invTol * estTol);
   tb1claTol = tb1claTol / ssqy;
	tb0claTol = meanestTol - tb1claTol * meany;
//  	dMat estclaTol = (estTol - tb0claTol)  / tb1claTol;
//   est.merge(estinvTol, ColWise);
//   est.merge(estclaTol, ColWise);

   SEXP R_Opt=R_NilValue, R_Coef=R_NilValue, ret;

   PROTECT(R_Opt = allocVector(REALSXP, nc*2));
   for (i=0;i<nc;i++) {
      if (ICount(i) == 0) {
         REAL(R_Opt)[i] = NA_REAL;
         REAL(R_Opt)[i + nc] = NA_REAL;
      }
      else {
         REAL(R_Opt)[i] = opt(i,0);
         REAL(R_Opt)[i + nc] = opt(i,1);
      }
   }
   PROTECT(R_Coef = allocVector(REALSXP, 8));
   REAL(R_Coef)[0] = tb0inv;
   REAL(R_Coef)[1] = tb0cla;
   REAL(R_Coef)[2] = tb0invTol;
   REAL(R_Coef)[3] = tb0claTol;
   REAL(R_Coef)[4] = tb1inv;
   REAL(R_Coef)[5] = tb1cla;
   REAL(R_Coef)[6] = tb1invTol;
   REAL(R_Coef)[7] = tb1claTol;
   UNPROTECT(2);

   PROTECT(ret = allocVector(VECSXP, 2)); 
   PROTECT(retNames = allocVector(STRSXP, 2));

   SET_VECTOR_ELT(ret, 0, R_Opt);
   SET_VECTOR_ELT(ret, 1, R_Coef);
   SET_STRING_ELT(retNames, 0, mkChar("Beta"));
   SET_STRING_ELT(retNames, 1, mkChar("DS_Coefficients"));

   SET_NAMES(ret, retNames);
   UNPROTECT(2);

//   delete str;
   return(ret);
}
}

extern "C" {
/*
__declspec(dllexport) 
*/
SEXP WA_predict(SEXP sexp_SpecData, SEXP sexp_Beta)
{
   SEXP dims, dims_beta, R_est = R_NilValue;
   dims = Rf_getAttrib(sexp_SpecData, R_DimSymbol);
   dims_beta = Rf_getAttrib(sexp_Beta, R_DimSymbol);
   int nr = INTEGER(dims)[0];
   int nc = INTEGER(dims)[1];

//   char *str = new char[1000]; 
   int i=0, j=0;
   
   dMat X(nr, nc, 0.0);
   dMat beta(nc, 2, 0.0);
   dMat est(nr, 4, -99.9);

   PROTECT(sexp_SpecData);
   PROTECT(sexp_Beta);
   for (i=0;i<nr;i++) {
      for (j=0;j<nc;j++) {
         X(i,j) = REAL(sexp_SpecData)[i + nr*j];
      }
   }
   for (i=0;i<nc;i++) {
      for (j=0;j<2;j++) {
         if (ISNA(REAL(sexp_Beta)[i + nc*j]))
            beta(i,j) = missingValue(beta);
         else
            beta(i,j) = REAL(sexp_Beta)[i + nc*j];
      }
   }
   UNPROTECT(2);

   PROTECT(R_est = allocVector(REALSXP, nr*4));
   for (i=0;i<nr;i++) {
       for (j=0;j<4;j++) {
           if (est.isMissing(i,j)) 
               REAL(R_est)[i + j*nr] = NA_REAL;
            else
               REAL(R_est)[i + j*nr] = est(i,j);
       }
   }
   UNPROTECT(1);

   return R_est;
}
}
