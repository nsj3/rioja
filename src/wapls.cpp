#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "mat.h"

using namespace std;

extern "C" {

#ifdef _MSC_VER
__declspec(dllexport) 
#endif

SEXP WAPLS_fit(SEXP sexp_SpecData, SEXP sexp_EnvData, SEXP sexpNPLS, SEXP sexpIsWAPLS, SEXP sexpStandX, SEXP sexpLean)
{
   SEXP dims, retNames=R_NilValue, R_meanY, R_meanT, R_sdX;
   dims = Rf_getAttrib(sexp_SpecData, R_DimSymbol);
   int nr = INTEGER(dims)[0];
   int nc = INTEGER(dims)[1];
   int nPLS = INTEGER(sexpNPLS)[0];
   bool bIsWAPLS = bool(INTEGER(sexpIsWAPLS)[0]);
   bool bStandX = bool(INTEGER(sexpStandX)[0]);
   bool bLean = bool(INTEGER(sexpLean)[0]);
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

   dMat SpecCount = count(X, ColWise);
   dMat R, C;
   double meanY;
	double Ytottot = 1.0;

	if (!bIsWAPLS) {

// we have linear method so centre y and calculate C, column weights

		meanY = mean(Y);
		Y -= meanY;
		R = dMat(nr, 1,1.0);
		if (bStandX) {
			C = transpose((sd(X, ColWise)));
			C *= sqrt(double(nr-1));
         for (j=0;j<nc;j++) {
            if (C(j,0)+1.0 <= 1.0) {
      			(C(j,0)=1.0);
            }
	      }
		}
		else {
			C = dMat(nc,1,1.0);
		}
	}
	else {

//  Weighted averaging, calculate column and row weights, centre y

		C = transpose(sum(X, ColWise));
		for (int j=0;j<nc;j++) {
			if (C(j,0)+1.0 <= 1.0)
				(C(j,0)=1.0);
		}
		R = sum(X, RowWise);
		Ytottot = sum(R);
		dMat wtr = R/Ytottot;
		meanY = sum(Y*wtr);
		Y -= meanY;
	}

   double meant=0.0, tau=0.0, alpha=0.0, gamma=0.0, gamma0=0.0;

   dMat b(nc,1,0.0), g, d, t, Wc = C / Ytottot, Wr = R / Ytottot, T, P;
	dMat gam(1,nPLS,0.0), meanT(1,nPLS,0.0), beta;

// initialisation: setup starting gradient g
// initialised and set to 0
// d is conjugate gradient, t is score vector

	g = X.tproduct(Y) / C;
	if (bIsWAPLS) {
		gamma0 = sum((g*g*Wc));
	}
	else {
		gamma0 = sumsq(g);
	}
	d = g;

// start of loop deriving pls components

	for (i=0;i<nPLS;i++) {
		if (bIsWAPLS) {
			t = X.product(d) / R;
			tau = sum((t*t*Wr));
		}
		else {
			gam(0,i) = sqrt(gamma0);
			t = X.product(d/C);
			if (i>0) {
				t -= (T / sumsq(T, ColWise) ).product( (transpose(T).product(t) ) );
			}
			meant = mean(t);
			t -= meant;
			tau = sumsq(t);
		}
		if ((tau+1.0) > 1.0) {
			alpha = gamma0 / tau;
			if (!bIsWAPLS) {
				meanT(0,i) = meant * alpha;    // save meant for future calibrations
			}

// update regression coefficients b and residuals r

			if (i==0) {                        // add jacknife later
				P = d * alpha;
			}
			else {
				P = P.concat(d * alpha, ColWise);
			}
			b += d * alpha;
			Y -= t * alpha;

// derive new gradient g from current residuals r

			g = X.tproduct(Y) / C;
			if (bIsWAPLS) {
				gamma = sum(g*g*Wc);
			}
			else {
				gamma = sumsq(g);
			}

			if ((1.0+gamma0) > 1.0) {
				d *= gamma/gamma0;
				d += g;
			}
			else {
				d.fill(0.0);
			}
			gamma0 = gamma;

		   if (i==0) {
				T = copy(t);
//				est = YOriginal - Y;
				if (bIsWAPLS) {
					beta = copy(b);
				}
				else {
					beta = b/C;
				}
			}
			else {
				T = T.concat(t,ColWise);
//				est = est.concat((YOriginal - Y),ColWise);
				if (bIsWAPLS) {
					beta = beta.concat(b,ColWise);
				}
				else {
					beta = beta.concat(b/C,ColWise);
				}
			}
      }
      else {
         break;
      }
   }
   if (i != nPLS) {
      Index II(i);
      gam = gam(II, ColWise);
//      sprintf(str, "Warning: Only %d components can be extracted", i);
//      SET_STRING_ELT(eMessage, 0, Rf_mkChar(str));
      nPLS = i;
   }
   if (bIsWAPLS) {
      for (j=0;j<nc;j++) {
         if (SpecCount(0,j) > 0.0) {
            for (int kk=0;kk<cols(beta);kk++)
   	         beta(j, kk) += meanY;  
         }
         else {
            for (int kk=0;kk<cols(beta);kk++)
   	         beta(j, kk) = missingValue(beta);  
         }
      }
   }
   else {
	   T /= gam;
      for (j=0;j<nc;j++) {
         if (SpecCount(0,j) < 1.0) {
            for (int kk=0;kk<cols(beta);kk++)
   	         beta(j, kk) = missingValue(beta);  
         }
      }
   }
   SEXP ret = R_NilValue, R_coef = R_NilValue;

   PROTECT(ret = Rf_allocVector(VECSXP, 6)); 
   PROTECT(retNames = Rf_allocVector(STRSXP, 6));

   PROTECT(R_coef = Rf_allocVector(REALSXP, nc*nPLS));
   for (i=0;i<nc;i++) {
       for (j=0;j<nPLS;j++) {
           if (beta.isMissing(i,j)) 
               REAL(R_coef)[i + j*nc] = NA_REAL;
            else
               REAL(R_coef)[i + j*nc] = (beta)(i,j);
       }
   }
   PROTECT(R_meanY = Rf_allocVector(REALSXP, 1));
   REAL(R_meanY)[0] = meanY;
   SET_VECTOR_ELT(ret, 0, R_coef);
   SET_VECTOR_ELT(ret, 1, R_meanY);
   SET_STRING_ELT(retNames, 0, Rf_mkChar("Beta"));
   SET_STRING_ELT(retNames, 1, Rf_mkChar("meanY"));
   UNPROTECT(2);

   if (!bLean) {
      SEXP R_T, R_P;
      PROTECT(R_T = Rf_allocVector(REALSXP, nr*nPLS));
      for (i=0;i<nr;i++) {
          for (j=0;j<nPLS;j++) {
              if (T.isMissing(i,j)) 
                  REAL(R_T)[i + j*nr] = NA_REAL;
               else
                  REAL(R_T)[i + j*nr] = T(i,j);
         }
      }
      PROTECT(R_P = Rf_allocVector(REALSXP, nc*nPLS));
      for (i=0;i<nc;i++) {
          for (j=0;j<nPLS;j++) {
              if (P.isMissing(i,j)) 
                  REAL(R_P)[i + j*nc] = NA_REAL;
               else
                  REAL(R_P)[i + j*nc] = P(i,j);
         }
      }
      SET_VECTOR_ELT(ret, 2, R_T);
      SET_VECTOR_ELT(ret, 3, R_P);
      UNPROTECT(2);
   }

   SET_STRING_ELT(retNames, 2, Rf_mkChar("T"));
   SET_STRING_ELT(retNames, 3, Rf_mkChar("P"));

   if (!bIsWAPLS) {
      PROTECT(R_meanT = Rf_allocVector(REALSXP, nPLS));
      for (j=0;j<nPLS;j++) {
         if (meanT.isMissing(0,j)) 
             REAL(R_meanT)[j] = NA_REAL;
          else
             REAL(R_meanT)[j] = meanT(0,j);
      }
      SET_VECTOR_ELT(ret, 4, R_meanT);
      PROTECT(R_sdX = Rf_allocVector(REALSXP, nc));
      for (j=0;j<nc;j++) {
         if (C.isMissing(0,j)) 
             REAL(R_sdX)[j] = NA_REAL;
          else
             REAL(R_sdX)[j] = C(j,0);
      }
      SET_VECTOR_ELT(ret, 5, R_sdX);
      UNPROTECT(2);
   }

   SET_STRING_ELT(retNames, 4, Rf_mkChar("meanT"));
   SET_STRING_ELT(retNames, 5, Rf_mkChar("sdX"));

   SET_NAMES(ret, retNames);
   UNPROTECT(2);

//   delete str;
   return(ret);
}
}

extern "C" {
#ifdef _MSC_VER
__declspec(dllexport) 
#endif
SEXP WAPLS_predict(SEXP sexp_SpecData, SEXP sexp_Beta, SEXP sexp_meanY, SEXP sexpIsWAPLS, SEXP sexpStandX, SEXP sexpSDX, SEXP sexp_meanT)
{
   SEXP dims, dims_beta, R_est = R_NilValue;
   dims = Rf_getAttrib(sexp_SpecData, R_DimSymbol);
   dims_beta = Rf_getAttrib(sexp_Beta, R_DimSymbol);
   int nr = INTEGER(dims)[0];
   int nc = INTEGER(dims)[1];
   int nPLS = INTEGER(dims_beta)[1];
   bool bIsWAPLS = bool(INTEGER(sexpIsWAPLS)[0]);
   double meanY = REAL(sexp_meanY)[0];
   dMat meanT;

//   char *str = new char[1000]; 
   int i=0, j=0, k=0;
   
   dMat X(nr, nc, 0.0);
   dMat beta(nc, nPLS, 0.0);
   PROTECT(sexp_SpecData);
   PROTECT(sexp_Beta);
   for (i=0;i<nr;i++) {
      for (j=0;j<nc;j++) {
         X(i,j) = REAL(sexp_SpecData)[i + nr*j];
      }
   }
   for (i=0;i<nc;i++) {
      for (j=0;j<nPLS;j++) {
         if (ISNA(REAL(sexp_Beta)[i + nc*j]))
            beta(i,j) = missingValue(beta);
         else
            beta(i,j) = REAL(sexp_Beta)[i + nc*j];
      }
   }
   UNPROTECT(2);

   if (!bIsWAPLS) {
//      bool bStandX = bool(INTEGER(sexpStandX)[0]);
      PROTECT(sexp_meanT);
      if (!bIsWAPLS) {
         meanT = dMat(1, nPLS);
         for (j=0;j<nPLS;j++) {
            meanT(0, j) = REAL(sexp_meanT)[j];
         }
      }
      UNPROTECT(1);
   }

   dMat est(nr, nPLS, 0.0);

 	for (j=0;j<nr;j++) {
   	double sum=0.0;
      for (k=0;k<nc;k++) {
         if (!beta.isMissing(k, 0)) {
            sum += X(j,k);
     		   for (i=0;i<cols(beta);i++) {
         		est(j,i) += X(j,k) * beta(k, i);
            }
         }
      }
	   for (i=0;i<cols(beta);i++) {
         if (bIsWAPLS) {
            if (sum > 1.0E-6) {
            	est(j,i) /= sum;
            }
            else {
               est(j,i) = missingValue(est);
            }
         }
         else {
            est(j,i) += meanY;
         }
      }
      if (!bIsWAPLS) {
         dMat mt = copy(meanT);
         for (i=1;i<nPLS;i++) {
			   mt(0,i) += mt(0,i-1);
         }
         for (i=0;i<nPLS;i++) {
            est(j, i) -= mt(0, i);
         }
      }
   }

   PROTECT(R_est = Rf_allocVector(REALSXP, nr*nPLS));
   for (i=0;i<nr;i++) {
       for (j=0;j<nPLS;j++) {
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
