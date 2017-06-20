#include <math.h>
#include "mat.h"
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

using namespace std;

extern "C" {

/* 
__declspec(dllexport) 
*/
SEXP WriteCornellFile(SEXP sexpData, SEXP sexpfName, SEXP sexpTitle, SEXP sexpSubFileType, SEXP sexpCornellFnPlaces, SEXP sexpCornellFnPerLine, SEXP sexpCornellCnPlaces, SEXP sexpCornellCnCouplets, SEXP sexpmV)
{
   int nFieldWidthWarnings = 0;
   char str[1000], str2[1000];
   int i, j;
   SEXP dims, eMessage = R_NilValue, dimnames;
   dims = Rf_getAttrib(sexpData, R_DimSymbol);
   int nr = INTEGER(dims)[0];
   int nc = INTEGER(dims)[1];
   double mV = REAL(sexpmV)[0];;
   dMat *pData = NULL;
   PROTECT(eMessage = allocVector(STRSXP, 1));
   PROTECT(dimnames = getAttrib(sexpData, R_DimNamesSymbol));
   
   try {
      pData = new dMat(nr, nc, 0.0);
      PROTECT(sexpData);
      for (int i=0;i<nr;i++) {
         for (int j=0;j<nc;j++) {
            if (ISNA(REAL(sexpData)[i + nr*j]))
               (*pData)(i,j) = mV;
            else
               (*pData)(i,j) = REAL(sexpData)[i + nr*j];
         }
      }
      UNPROTECT(1);
   }
   catch (char *eM) {
      delete pData;
      SET_STRING_ELT(eMessage, 0, mkChar(eM));
      UNPROTECT(2);
      return eMessage;
   }

   PROTECT(sexpfName);
   int nNumericFields = nc;
   int nRecords = nr;

   PROTECT(sexpTitle);
   int CornellFnPlaces = INTEGER(sexpCornellFnPlaces)[0];
   int CornellCnPlaces = INTEGER(sexpCornellCnPlaces)[0];
   int CornellCnCouplets = INTEGER(sexpCornellCnCouplets)[0];
   int CornellFnPerLine = INTEGER(sexpCornellFnPerLine)[0];

   int SubFileType = INTEGER(sexpSubFileType)[0];;  // 0 = Condensed, 1 = FULL;

//   FILE *fout; 

   const char * cc = CHAR(STRING_ELT(sexpfName, 0));

   FILE *fout = fopen(cc, "wt");
   if (!fout) {
      delete pData;
      sprintf(str, "Cannot open file %s", CHAR(STRING_ELT(sexpfName, 0)));
      SET_STRING_ELT(eMessage, 0, mkChar(str));
      UNPROTECT(5);
      return eMessage;
   }

//   fprintf(fout, "%5d%5d %-s", nNumericFields, nRecords, NULL); //, CHAR(STRING_ELT(sexpTitle, 0)));
   sprintf(str, "%5d%5d %-s", nNumericFields, nRecords, CHAR(STRING_ELT(sexpTitle, 0)));
   fwrite(str, sizeof(char), strlen(str), fout);

   int nItemsPerLine;
//   CString s;
   int nitems;
   int nAdjust = 7;
   int nMaxCharsPerLine = 0;
//  adjuster for correct field width output
   if (SubFileType == 1) {
      if (nNumericFields <= CornellFnPerLine) {
         sprintf(str, "\n(I5, 3X, %dF%d.0)", nNumericFields, CornellFnPlaces+nAdjust+1);
         fwrite(str, sizeof(char), strlen(str), fout);
      }
      else {
         sprintf(str, "\n(I5, 3X, %dF%d.0/(8X, %dF%d.0))", CornellFnPerLine, CornellFnPlaces+nAdjust+1, CornellFnPerLine, CornellFnPlaces+nAdjust+1);
         fwrite(str, sizeof(char), strlen(str), fout);
      }
      sprintf(str, "\n%d", nNumericFields);
      fwrite(str, sizeof(char), strlen(str), fout);
      nItemsPerLine = CornellFnPerLine;
      nMaxCharsPerLine = 8 + (CornellFnPerLine * (CornellFnPlaces+nAdjust+1));
   }
   else {
      sprintf(str2, "(I5,3X,%d(I4,F%d.0))", CornellCnCouplets, CornellCnPlaces+nAdjust+1);
      sprintf(str, "\n%-68s%2d", str2, CornellCnCouplets);
      fwrite(str, sizeof(char), strlen(str), fout);
      nItemsPerLine = CornellCnCouplets;
      nMaxCharsPerLine = 8 + (CornellCnCouplets * (CornellCnPlaces+nAdjust+1+4));
   }

   for (i=0;i<nRecords;i++) {
      nitems = 0;
      sprintf(str, "\n%5d   ", i+1);
      fwrite(str, sizeof(char), strlen(str), fout);
      for (j=0;j<nNumericFields;j++) {
         double x = (*pData)(i, j);
         if (SubFileType == 0 && (fabs(x) < 1.0E-20))
            continue;
         if (nitems == nItemsPerLine) {
            if (SubFileType == 0) {
               sprintf(str, "\n%5d   ", i+1);
               fwrite(str, sizeof(char), strlen(str), fout);
            }
            else {
               sprintf(str, "\n        ");
               fwrite(str, sizeof(char), strlen(str), fout);
            }
            nitems = 0;
         }
         if (SubFileType == 0) {
//            s.Format("%*.*g", CornellCnPlaces+nAdjust, CornellCnPlaces, x);
            sprintf(str2, "%*.*g", CornellCnPlaces+nAdjust, CornellCnPlaces, x);
            if (strlen(str2) > (size_t) (CornellCnPlaces+nAdjust))
               nFieldWidthWarnings++;
            sprintf(str, "%4d %s", j+1, str2);
            fwrite(str, sizeof(char), strlen(str), fout);
         }
         else {
//            s.Format("%*.*g", CornellFnPlaces+nAdjust, CornellFnPlaces, x);
            sprintf(str2, "%*.*g", CornellFnPlaces+nAdjust, CornellFnPlaces, x);
            if (strlen(str2) > (size_t) CornellFnPlaces+nAdjust)
               nFieldWidthWarnings++;
            sprintf(str, " %s", str2);
            fwrite(str, sizeof(char), strlen(str), fout);
         }
         nitems++;
      }
   }
   if (SubFileType == 0) {
      sprintf(str,"\n00000");
      fwrite(str, sizeof(char), strlen(str), fout);
   }
   else {
      nitems = 0;
      sprintf(str, "\n00000   ");
      fwrite(str, sizeof(char), strlen(str), fout);
      for (int j=0;j<nNumericFields;j++) {
         double x = 0;
         if (nitems == nItemsPerLine) {
            sprintf(str, "\n        ");
            fwrite(str, sizeof(char), strlen(str), fout);
            nitems = 0;
         }
         sprintf(str, " %*.*g", CornellFnPlaces+nAdjust, CornellFnPlaces, x);
         fwrite(str, sizeof(char), strlen(str), fout);
         nitems++;
      }
   }
   sprintf(str, "\n");
   fwrite(str, sizeof(char), strlen(str), fout);
   nitems = 0;
   for (j=0;j<nNumericFields;j++) {
      if (nitems == 10) {
         sprintf(str, "\n");
         fwrite(str, sizeof(char), strlen(str), fout);
         nitems = 0;
      }
      sprintf(str,"%8.8s", CHAR(STRING_ELT(VECTOR_ELT(dimnames, 1), j)));
      fwrite(str, sizeof(char), strlen(str), fout);
      nitems++;
   }
   sprintf(str, "\n");
   fwrite(str, sizeof(char), strlen(str), fout);
   nitems = 0;
   for (j=0;j<nRecords;j++) {
      if (nitems == 10) {
         sprintf(str, "\n");
         fwrite(str, sizeof(char), strlen(str), fout);
         nitems = 0;
      }
      sprintf(str,"%8.8s", CHAR(STRING_ELT(VECTOR_ELT(dimnames, 0), j)));
      fwrite(str, sizeof(char), strlen(str), fout);
      nitems++;
   }
   sprintf(str, "\n");
   fwrite(str, sizeof(char), strlen(str), fout);

   if (nFieldWidthWarnings) {
      sprintf(str, "\n%d Field width warnings", nFieldWidthWarnings);
      SET_STRING_ELT(eMessage, 0, mkChar(str));
   }
   if (nMaxCharsPerLine > 80) {
      sprintf(str, "Output contains lines upto %d characters long.\nThis may cause problems with some versions of CANOCO", nMaxCharsPerLine);
      SET_STRING_ELT(eMessage, 0, mkChar(str));
   }
   fclose(fout);
   delete pData;
   UNPROTECT(5);
   return eMessage;
}
}
