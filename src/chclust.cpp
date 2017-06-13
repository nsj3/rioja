#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "mat.h"
#include "nutil.h"

using namespace std;

bool Conslink(long nsam, double **DPtr, double **dend);
bool ConISS(long nsam, double **DPtr, double **dend);
bool CalcDissimilarity(dMat &Data, double ***DissimPtr, int coef);

extern "C" {
 /* __declspec(dllexport) */
SEXP chclust(SEXP sexpData, SEXP sexMethod)
{
   SEXP dims, eMessage = R_NilValue;
   dims = Rf_getAttrib(sexpData, R_DimSymbol);
   int method = INTEGER(sexMethod)[0];
   long nr = INTEGER(dims)[0];
   PROTECT(sexpData);
   double **DPtr = new double*[nr];
   double *diss;
   long i =0;
	for (i=1;i<nr;i++) {
		if ((DPtr[i]= new double[i])==NULL)
			return eMessage;
		diss=DPtr[i];
		for(long j=0;j<i;j++) {
            diss[j] = REAL(sexpData)[i + nr*j];
      }
	}
   UNPROTECT(1);
   double *dend = NULL;
   bool errorFlag = false;
   if (method==1) {
      if (!Conslink(nr, DPtr, &dend))  {
         PROTECT(eMessage = allocVector(STRSXP, 1));
         SET_STRING_ELT(eMessage, 0, mkChar("Error in Conslink C++ code"));
         errorFlag = true;
      }
   }
   else if (method==2) {
      if (!ConISS(nr, DPtr, &dend)) {
         PROTECT(eMessage = allocVector(STRSXP, 1));
         SET_STRING_ELT(eMessage, 0, mkChar("Error in ConISS C++ code"));
         errorFlag = true;
      }
   }
   else {
      PROTECT(eMessage = allocVector(STRSXP, 1));
      SET_STRING_ELT(eMessage, 0, mkChar("Unknown clustering method"));
         errorFlag = true;
   }
   SEXP sDend;
   PROTECT(sDend = allocVector(REALSXP, nr-1));
   for (i=1;i<nr;i++) {
      REAL(sDend)[i-1] = dend[i];
   }
  if (dend) 
    delete [] dend;
 	for (i=1;i<nr;i++)
      delete [] DPtr[i];
   delete [] DPtr;
   UNPROTECT(1);
   if (errorFlag) {
      UNPROTECT(1);
      return eMessage;
   }
   return sDend;
}

}

// I apologise for the fairly horrible code below.  Translated from the orginal FORTRAN ZONATION written program by John Birks & John Line when I was learning C/C++.

#define dc(a,b)  (((a) > (b)) ? (DPtr[(a)-1][(b)-1]) : (DPtr[(b)-1][(a)-1]))

double Update(double **DPtr, long j, long p, long q, long *nclus, long *name, double dshort, long np, long nq);
void Minim(double *diag, double *tiny, long *least, long *ncount, long nsam);
void Group(double **DPtr, double *diag, double tiny, double *prev, double *dend, long *least, long ncount, long *nclust, double large, long *nbit, long nsam, long *nlev, char *nsplur, FILE *fout);

double Update(double **DPtr, long j, long p, long q, long *nclus, long *name, double dshort, long np, long nq)
{
	long nr;
	nr = nclus[name[j-1]-1];
	return( ( ((double)(nr+np)) *dc(j,p)+((double)(nr+nq))*dc(j,q)-((double)nr)*dshort) / ((double)(nr+np+nq)));
}

bool ConISS(long nsam, double **DPtr, double **es)
{
	double *ess, dshort,e,de;
	long i=0 ,n=0, j=0, iter=0, p, q, np, nq, *name, *nclus, msiz,namp,namq;

	ess = new double[nsam];
	*es =  new double[nsam];
	nclus = new long[nsam];
	name = new long[nsam];
	if ((ess==NULL)||(es==NULL)||(nclus==NULL)||(name==NULL))
		return(false);

	for (i=0;i<nsam;i++) {
		nclus[i]=1;
		ess[i] =0.0;
		name[i] = i+1;
	}
	msiz = nsam-1;
	e = 0.0;

// FULL
/*
{
		fprintf(fout,"\n\n%65sMean\n%50sWithin-        within","","");
		fprintf(fout,"\n        Clusters    Increase in    Total          cluster        cluster");
		fprintf(fout,"\n Iter   merged      dispersion     dispersion     dispersion     dispersion\n");
}
*/
	for (iter=1;iter<=nsam-1;iter++) {
		dshort = dc(2,1);
		p = 1;
		for (n=2;n<=msiz;n++) {
			if (dc(n+1,n) < dshort) {
				dshort = dc(n+1,n);
				p = n;
			}
		}
		q = p + 1;
		namp = name[p-1];
		namq = name[q-1];
		np = nclus[namp-1];
		nq = nclus[namq-1];

		de = 0.5 * dshort;
		e += de;
		ess[namp-1] += ess[namq-1] + de;
		(*es)[namq-1] = e;

//FULL
//      fprintf(fout,"\n%4d %4d %4d %14.7f %14.7f %14.7f %14.7f", iter, namp, namq, de, e, ess[namp-1], (ess[namp-1]/ (double) (np+nq)));

		for (j=1;j<p;j++) {
			DPtr[p-1][j-1] = Update(DPtr,j,p,q,nclus,name,dshort,np,nq);
			for(i=q;i<=msiz;i++) {
				DPtr[i-1][j-1] = dc(i+1,j);
			}
		}
		for (i=q;i<=msiz;i++) {
			DPtr[i-1][p-1] = Update(DPtr,i+1,p,q,nclus,name,dshort,np,nq);
		}
		for (j=q;j<msiz;j++) {
         for (i=j+1;i<=msiz;i++) {
            DPtr[i-1][j-1] = dc(i+1,j+1);
         }
      }
      for (i=q;i<=msiz;i++)
			name[i-1] = name[i];
		nclus[namp-1] = np+nq;
		msiz--;
	}
//	fprintf(fout,"\n\nSample  Level  Total\nNumber         Dispersion\n");
//	DendrogramPlot(Data, es, fout);
//	delete ess;
//	delete nclus;
//	delete name;
	delete [] ess;
	delete [] nclus;
	delete [] name;
	return true;
}

void Minim(double *diag, double *tiny, long *least, long *ncount, long nsam)
{
   long i=0;
   double diff;

   *ncount=1;
   least[0]=1;
   *tiny = diag[1];
   for (i=2;i<nsam;i++) {
      diff = *tiny-diag[i];
      if (diff<0.0)  {
         continue;
      }
      else if (diff>1.0E-30) {
         *ncount = 1;
         *tiny = diag[i];
         least[*ncount-1] = i;
      }
      else if (diff<=1.0E-30) {
         (*ncount)++;
         least[*ncount-1] = i;
      }
   }
}

void Group(double **DPtr, double *diag, double tiny, double *prev, double *dend, long *least, long ncount, long *nclust, double large, long *nbit, long nsam, long *nlev, char *nsplur)
{
   long i=0,j=0,k=0,l=0,up=0, down=0, up1=0, up2=0, down1=0, down2=0;
   double diss;

   for (i=0;i<ncount;i++) {
      j = least[i];
      nsplur[j] = '*';
      diag[j] = large;
      nbit[j] = 1;
      
//      for(k=1;k<nsam;k++) 

      down=j;
      up = j;

/* find limits of new group, put in up & down */
      while ((down<nsam-1) && (nbit[down+1]==1)) {
         down++;
      }
      while ((up>1)&&(nbit[up-1]==1)) {
         up--;
      }
/* find limits of neighbouring groups */
      if (up>1) {
         down1=up-1;
         if (up>2) {
            up1=down1;
            while (up1>1) {
               if ((nbit[up1-1])==1) {  up1--; }
               else { break; }
            }
         }
         else
            down1=up1=1;
      }
      if (down<nsam-1) {
         up2=down+1;
         if (down<nsam-2) {
            down2=up2;
            while ((down2)<(nsam-1)) {
               if ((nbit[down2+1])==1) { down2++; }
               else { break; }
            }
         }
         else
            down2=up2=nsam-1;
      }

/* calc new dcs and update diag */
      diss = large;
      if (up>1) {
         for (k=up;k<=down+1;k++) {
            for (l=up1;l<=down1;l++) {
               if (dc(k,l) < diss) {
                  diss = dc(k,l);
               }
            }
         }
         diag[down1]=diss;
      }
      diss = large;
      if (down<nsam-1) {
         for (k=up;k<=down+1;k++) {
            for (l=up2+1;l<=down2+1;l++) {
               if (dc(k,l) < diss) {
                  diss = dc(k,l);
               }
            }
         }
         diag[up2]=diss;
      }
      if ((tiny-*prev)<0.0) {
         dend[j]=*prev;
      }
      else {
         dend[j]=tiny;
         *prev=tiny;
      }
      (*nlev)++;
      (*nclust)--;
      nsplur[j]=' ';
   }
}

bool Conslink(long nsam, double **DPtr, double **dend)
{
   double *diag, large, tiny;
   long nclust, *nbit, *least, nlev, ncount;
   char *nsplur;

   long i=0, nsam2;
   double prev;
   diag = new double[nsam+1];
   *dend = new double[nsam+1];
   nsplur = new char[nsam+1];
   nbit = new long[nsam+1];
   least = new long[nsam+1];
   if ((diag==NULL)||(*dend==NULL)||(nsplur==NULL)||(nbit==NULL)||(least==NULL))
      return false;

/*  Setup parameters   */
   large=0.0;
   nsam2 = unsigned(nsam);
   for (i=1;i<nsam2;i++) {
      diag[i] =  dc(i,i+1);
      nsplur[i]='\\';
      nbit[i]=0;
      large += diag[i];
      nbit[i] = 0;
   }
   prev = 0.0;
   nclust = nsam;
   nlev=0;
   while (nclust>1) {
      Minim(diag, &tiny, least, &ncount, nsam);
      Group(DPtr, diag, tiny, &prev, *dend, least, ncount, &nclust, large, nbit, nsam, &nlev, nsplur);
   }
   delete [] diag;
   delete [] nsplur;
   delete [] nbit;
   delete [] least;
   return true;
}

