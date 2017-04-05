#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "tilfuncs.h"
#include "mat.h"

using namespace std;

typedef struct {
   char pollentypes[61];
   int num;
   char shortname[21];
   char sum;
} PTYPES;

bool TiliaBinIn(dataMat &dData, FILE *fin, char *title, PTYPES **pnames, double **pDepths, char *strError);
void RemoveSpaces(char *dest, char *src);
void RemoveTrailingSpaces(char *dest, char *src);

extern "C" {
/*
__declspec(dllexport) 
*/
SEXP ReadTiliaFile(SEXP fN)
{
   dataMat dData;
   const char *fName = CHAR(STRING_ELT(fN, 0));
   SEXP mat = R_NilValue, names=R_NilValue, ans=R_NilValue, rnames = R_NilValue, errorM = R_NilValue, retNames=R_NilValue, pFullNames=R_NilValue, pCodeNames=R_NilValue, pCodeNums=R_NilValue, pSums=R_NilValue, pDepths=R_NilValue;
   PROTECT(ans = allocVector(VECSXP, 8)); 
   PROTECT(retNames = allocVector(STRSXP, 8));
   SET_STRING_ELT(retNames, 0, mkChar("Data"));
   SET_STRING_ELT(retNames, 1, mkChar("RowNames"));
   SET_STRING_ELT(retNames, 2, mkChar("CodeNames"));
   SET_STRING_ELT(retNames, 3, mkChar("CodeNums"));
   SET_STRING_ELT(retNames, 4, mkChar("FullNames"));
   SET_STRING_ELT(retNames, 5, mkChar("Sum"));
   SET_STRING_ELT(retNames, 6, mkChar("Depths"));
   SET_STRING_ELT(retNames, 7, mkChar("ErrorMessage"));
   SET_NAMES(ans, retNames);
   PROTECT(errorM = allocVector(STRSXP, 1));
   bool bError = false;

   FILE *fin = fopen(fName, "rb");
   char title[100];
   char strError[100];
   PTYPES *pnames;
   double *pdDepths;
   if (fin!=NULL) {
      bError = TiliaBinIn(dData, fin, title, &pnames, &pdDepths, strError);
     fclose(fin);
   }
   else {
      bError = false;
      SET_STRING_ELT(errorM, 0, mkChar("Cannot open file"));
   }
   
   if (bError) {
     int nr = rows(dData);
     int nc = cols(dData);
      if (nr * nc > 0) {
         PROTECT(mat = allocVector(VECSXP, nc));
         PROTECT(names = allocVector(STRSXP, nc));
         PROTECT(pFullNames = allocVector(STRSXP, nc));
         PROTECT(pCodeNums = allocVector(INTSXP, nc));
         PROTECT(pSums = allocVector(STRSXP, nc));
         PROTECT(pCodeNames = allocVector(STRSXP, nc));
         dMat *dm = getdMat(dData);
         char str[5];
         for (int i=0;i<nc;i++) {
            SET_STRING_ELT(names, i, mkChar(dData.spName(i)));
            SET_STRING_ELT(pCodeNames, i, mkChar(pnames[i].shortname));
            SET_STRING_ELT(pFullNames, i, mkChar(pnames[i].pollentypes));
            INTEGER(pCodeNums)[i] = pnames[i].num;
            sprintf(str, "%c", pnames[i].sum);
            SET_STRING_ELT(pSums, i, mkChar(str));
            SET_VECTOR_ELT(mat, i, allocVector(REALSXP, nr));
            for (int j=0;j<nr;j++) {
              if (dm->isMissing(j,i)) 
                  REAL(VECTOR_ELT(mat, i))[j] = NA_REAL;
               else
                  REAL(VECTOR_ELT(mat, i))[j] = (*dm)(j,i);
            }
         }
         PROTECT(pDepths = allocVector(REALSXP, nr));
         PROTECT(rnames = allocVector(STRSXP, nr));
         for (int j=0;j<nr;j++) {
           SET_STRING_ELT(rnames, j, mkChar(dData.samName(j)));
            REAL(pDepths)[j] = (pdDepths)[j];
         }
//         for (int i=0;i<nc;i++) {
      }
      SET_NAMES(mat, names);
      SET_STRING_ELT(errorM, 0, mkChar(strError));
   }
   
   SET_VECTOR_ELT(ans, 0, mat);
   SET_VECTOR_ELT(ans, 1, rnames);
   SET_VECTOR_ELT(ans, 2, pCodeNames);
   SET_VECTOR_ELT(ans, 3, pCodeNums);
   SET_VECTOR_ELT(ans, 4, pFullNames);
   SET_VECTOR_ELT(ans, 5, pSums);
   SET_VECTOR_ELT(ans, 6, pDepths);
   SET_VECTOR_ELT(ans, 7, errorM);

   if (!bError)
      UNPROTECT(3);
   else
      UNPROTECT(11);
   return(ans);
}
}

bool TiliaBinIn(dataMat &S, FILE *fin, char *fname, PTYPES **pnames, double **pDepths, char *strError)
{
   char name[200];
   char code[20];
   int spnum, i;
   char shortcode[20];
   char version=1;
   size_t retval = fread(name,10,1,fin);
   if (retval == 0) {
		sprintf(strError,"Cannot read Tilia file");
		return false;
	}
	 name[10] = '\0';

  // first 6 characters must be: "tilia "

	if (strncmp(name,"tilia ",6)) {
		sprintf(strError,"This is not a TILIA file");
		return false;
	}
	double ver = atof(&name[5]);
	if (ver < 1.07) {
		sprintf(strError,"TILIA file must be at least version 1.07");
		return false;
	}
   if (ver > 1.99) {
      version = 2;
   }
   int n, m;
   if (version==1) {
      if (!Tilia1ReadFlags(fin, n, m))
         return false;
   }
   else {
      if (!Tilia2ReadFlags(fin, n, m))
         return false;
   }
   PTYPES *ppnames = new PTYPES[m];
   *pnames = ppnames;
   if (! ppnames) {
      sprintf(strError,"Out of memory allocating space for variable names");
   	return false;
   }
   char * tmp = new char[1000];
	for (i=0;i<m;i++) {
      char sum;
      if (version==1) {
         if (Tilia1ReadVar(fin, name, code, shortcode, spnum, sum)==FALSE) {
		      sprintf(strError, "Unexpected EOF while reading taxon names");
		   	return false;
   		}
         if (strlen(code) < 1)
            strcpy(code, shortcode);
         spnum = i+1;
      }
      else {
         if (Tilia2ReadVar(fin, name, code, spnum, sum)==FALSE) {
		      sprintf(strError,"Unexpected EOF while reading taxon names");
	   		return false;
		   }
      }
      strncpy(ppnames[i].pollentypes, name, 61);
		strncpy(ppnames[i].shortname, code, 19);
//    RemoveTrailingSpaces(ppnames[i].pollentypes, ppnames[i].pollentypes);
    RemoveTrailingSpaces(tmp, ppnames[i].pollentypes);
    strcpy(ppnames[i].pollentypes, tmp);
   	ppnames[i].pollentypes[60] = '\0';
   	ppnames[i].shortname[20] = '\0';
		ppnames[i].num = spnum;
    ppnames[i].sum = sum;
	}
	
	S.kill();

	S.setmatType(full);
	
	dMat *CC;
	CC = new dMat(n, m, 0.0);
	S.setdMat(CC);

	int *Dl = new int[n];
	char **samnam = new char*[n];
	char **spnam = new char*[m];
	char *sp = new char[m*9];
	char *sam = new char[n*9];
	if ((!Dl)||(!spnam)||(!samnam)||(!sp)||(!sam)) {
		sprintf(strError,"Out of memory allocating storage for sample names");
		return false;
	}
	for (i=0;i<m;i++)
		spnam[i] = &sp[i*9];
	for (i=0;i<n;i++)
		samnam[i] = &sam[i*9];

	S.setspName(spnam);
	S.setsamNo(Dl);
	S.setsamName(samnam);

	*pDepths = new double[n];

	for (i=0;i<n;i++) {
      float num;
      if (TiliaReadSample(fin, num, name)==FALSE) {
  			sprintf(strError,"Unexpected EOF while reading sample depths");
   		return false;
		}
      char str[30];
      strncpy(S.samName(i), name, 8);
      if (strlen(name) < 1) {
         sprintf(str, "%-g", num);
         strncpy(S.samName(i), str, 8);
      }
      S.samName(i)[8] = '\0';
		RemoveSpaces(S.samName(i),S.samName(i));
		S.samNo(i) = i+1;
      (*pDepths)[i] = num;
	}

  float x;
  char byte;
	for (int j=0;j<m;j++) {
		for (i=0;i<n;i++) {
  		   if (TiliaReadData(fin, byte, x)==FALSE) {
			   sprintf(strError,"Error reading data for taxon %4d", m+1);
			   return false;
		   }
         if (version==1) {
            if (byte==1) {
			      (*CC)(i,j) = -99.9;
            }
            else
			      (*CC)(i,j) = x;
         }
         else {
            if (byte==2) {
			      (*CC)(i,j) = x;
            }
            else
			      (*CC)(i,j) = -99.9;
         }
		}
	}
	

	for (i=0;i<m;i++) {
      if (strlen(ppnames[i].shortname) < 1) {
		   strncpy(S.spName(i), ppnames[i].pollentypes, 8);
         S.spName(i)[8] = '\0';
      }
      else {
		   strncpy(S.spName(i), ppnames[i].shortname, 8);
         S.spName(i)[8] = '\0';
      }
  	}
   strcpy(strError, "No error");

	return true;
}
