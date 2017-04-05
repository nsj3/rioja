#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "mat.h"

using namespace std;

void NewCornellIn2(dataMat &S, char * fname, int etf, double missing_value, char &InFileType, long &nMissingValues, int &ColsWithNoData, int &RowsWithNoData, int &nCouplets, double impliedZero);

extern "C" {
/* __declspec(dllexport) */
SEXP ReadCornellFile(SEXP fN, SEXP mValue, SEXP impZero)
{
   dataMat dData;
   int ColsWithNoData=0;
   int RowsWithNoData=0;
   int nCouplets=0;
   long nMissingValues=0;
   char InFileType;
   const char *fName = CHAR(STRING_ELT(fN, 0));
   SEXP mat = R_NilValue, names = R_NilValue, ans = R_NilValue, rnames = R_NilValue, errorM = R_NilValue, retNames = R_NilValue, iSumm=R_NilValue, iSummNames = R_NilValue;
   PROTECT(ans = allocVector(VECSXP, 4)); /* create answer [0] = data, [1] = row names [2] = Erro message*/
   PROTECT(retNames = allocVector(STRSXP, 4));
   PROTECT(iSumm = allocVector(INTSXP, 2));
   SET_STRING_ELT(retNames, 0, mkChar("Data"));
   SET_STRING_ELT(retNames, 1, mkChar("RowNames"));
   SET_STRING_ELT(retNames, 2, mkChar("ErrorMessage"));
   SET_STRING_ELT(retNames, 3, mkChar("Summary"));
   bool bError = false;

   double *mV = REAL(PROTECT(mValue));
   double *iZ = REAL(PROTECT(impZero));
   try {
      NewCornellIn2(dData, (char *) fName, 1, (double) mV[0], InFileType, nMissingValues, ColsWithNoData, RowsWithNoData, nCouplets, iZ[0]);
   }
   catch (char *ErrorMessage) {
      bError = true;
      PROTECT(errorM = allocVector(STRSXP, 1));
      SET_STRING_ELT(errorM, 0, mkChar(ErrorMessage));
   }
   if (!bError) {
      int nr = rows(dData);
      int nc = cols(dData);
      if (nr * nc > 0) {
         PROTECT(mat = allocVector(VECSXP, nc));
         dMat *dm = getdMat(dData);
         PROTECT(names = allocVector(STRSXP, nc));
         for (int i=0;i<nc;i++) {
            SET_STRING_ELT(names, i, mkChar(dData.spName(i)));
            SET_VECTOR_ELT(mat, i, allocVector(REALSXP, nr));
            for (int j=0;j<nr;j++) {
               if (dm->isMissing(j,i)) 
                  REAL(VECTOR_ELT(mat, i))[j] = NA_REAL;
               else
                  REAL(VECTOR_ELT(mat, i))[j] = (*dm)(j,i);
            }
         }
         PROTECT(rnames = allocVector(STRSXP, nr));
         for (int j=0;j<nr;j++) {
            SET_STRING_ELT(rnames, j, mkChar(dData.samName(j)));
         }
      }
      SET_NAMES(mat, names);     
      PROTECT(iSummNames = allocVector(STRSXP, 2));
      SET_STRING_ELT(iSummNames, 0, mkChar("Number of missing values"));
      SET_STRING_ELT(iSummNames, 1, mkChar("Number of empty columns"));
      INTEGER(iSumm)[0] = nMissingValues;
      INTEGER(iSumm)[1] = ColsWithNoData;
      SET_NAMES(iSumm, iSummNames);     
   }
   SET_VECTOR_ELT(ans, 0, mat);
   SET_VECTOR_ELT(ans, 1, rnames);
   SET_VECTOR_ELT(ans, 2, errorM);
   SET_VECTOR_ELT(ans, 3, iSumm);
   SET_NAMES(ans, retNames);
   if (bError)
      UNPROTECT(6);
   else
      UNPROTECT(9);
   return(ans);
}
}
