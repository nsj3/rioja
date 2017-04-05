#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mat.h"
#include "nutil.h"

using namespace std;

#define MAX_SPECIES 5000
#define MAX_SAMPLES 50000
#define MIN(a, b) (a < b ? a : b) 

class sample {
public:
   double *m;
   Index *I;
   int samno;
   sample *next;
   sample(double *f, int i, int *ii, int no) { m = f; I = new Index(i,ii); samno = no; next = NULL; };
   ~sample() { if (I) delete I; };
};

void ccleanup(dataMat &S, void **memb, sample *first);
int TestForNumber(char *s, int *n, float *f);

#ifdef __cplusplus
extern "C" {
   int openf_(char *fname, char *title, char *format, int *ncoup, int *chan, int *tag, int *filenamelength);
   int getlin_(int *chan, char *format, int *ncoup, int *cursam, int *sp, double *abun, int *tag);
   int getl2_(int *chan, char *format, int *ncoup, int *cursam, double *abun, int *tag);
   int getnam_(int *chan, char *spnam, char *samnam, int *nsp, int *nsam, int *tag);
   int closef_(int chan);
}
#endif

void NewCornellIn2(dataMat &S, char * fname, int etf, double missing_value, char &type, long &nMissingValues, int &ColsWithNoData, int &RowsWithNoData, int &nCouplets, double impZero)
{
   int lsp, lsam, new_sample, ncoup, i, j, missing_species, missing_flag, nsam, nsp, current_sample;
   double zero_tolerance = 1.0E-08;
   type = full;
   nsam = 0;
   nsp = 0;
   lsp = 0;
   current_sample = 0;
   int channel = 7, tag = 0;
   int *iptr;
   sample *first = NULL;
   sample *next = NULL, *newS;
   unsigned nitem=0;
   ColsWithNoData = RowsWithNoData = 0;
   char *format = new char[201];
   char *buffer = new char[201];
   static char errorline[200];
   void *membuffers[10] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };

   membuffers[0] = (void *) format;
   membuffers[1] = (void *) buffer;

   int len2 = (int) strlen(fname);
   openf_(fname, buffer, format, &ncoup, &channel, &tag, &len2);
   if (tag) {
      ccleanup(S,membuffers, first);
      if (tag==1) sprintf(errorline, "Error opening %s - file does not exist", fname);
      else if (tag==2) sprintf(errorline, "Error reading header information from file %s", fname);
      else if (tag==3) sprintf(errorline, "Unexpected end of file while reading header");
      else if (tag==4) sprintf(errorline, "Integer for number of couplets not found");
      throw(errorline);
//		return ERROR;
   }
   char *chptr = format;
   int ft = 0;
   while (*chptr) { if ((*chptr == 'i' || *chptr == 'I')) ft++; chptr++; }
   if (ft > 1)
      type = sparse;

   char integ[11];
   strncpy(integ, buffer, 10);
   integ[10] = '\0';
   int nos;
   char *tit;
   if (ValidInt(integ, &nos) >= 1)
      tit = &buffer[10];
   else
      tit = buffer;
   char *tmp = new char[1000];
   RemoveLeadingSpaces(tmp,tit);
   RemoveTrailingSpaces(tit,tmp);
   delete [] tmp;
   int tlen = MIN((int)strlen(tit),80);
   char *title = new char[tlen+1];
   sprintf(title,"%-.80s",tit);
   lsam=0;
   missing_flag=0;
   if (nsam==0) { nsam = MAX_SAMPLES; }
   if (nsp==0) { nsp = MAX_SPECIES; }
   if (type==full) nsp = ncoup;
   int *sp2 = new int[ncoup];
   int *sp = new int[nsp];
   double *ab = new double[nsp];
   double *ab2 = new double[ncoup];
   double *m = new double[nsp];
   if ((sp==NULL)||(ab==NULL)||(sp2==NULL)||(ab2==NULL)||(m==NULL)) {
      if (sp2) delete[] sp2;
      if (ab2) delete[] ab2;
      if (sp) delete[] sp;
      if (ab) delete[] ab;
      if (m) delete[] m;
      ccleanup(S,membuffers, first);
      closef_(channel);
      throw("Out of memory allocating work space");
//      return ERROR;
	}
   membuffers[2] = (void *) sp;
   membuffers[3] = (void *) ab;
   membuffers[4] = (void *) m;
   membuffers[5] = (void *) sp2;
   membuffers[6] = (void *) ab2;
   do {
      if (type==sparse) {
         getlin_(&channel, format, &ncoup, &new_sample, sp2, ab2, &tag);
         if (tag > 0) {
            ccleanup(S,membuffers, first);
            if (tag==1) sprintf(errorline, "Error reading data for sample %d", new_sample);
            else if (tag==2) sprintf(errorline, "Unexpected end of file reading sample %d", new_sample);
            closef_(channel);
            throw(errorline);
//            return ERROR;
         }
         if (lsam==0) {
            current_sample=new_sample;
            lsam++;
         }
      }
      else {
         getl2_(&channel, format, &ncoup, &new_sample, ab2, &tag);
         if (tag > 0) {
            ccleanup(S,membuffers, first);
            if (tag==1) sprintf(errorline, "Error reading data for sample %d", new_sample);
            else if (tag==2) sprintf(errorline, "Unexpected end of file reading sample %d", new_sample);
            closef_(channel);
            throw(errorline);
//            return ERROR;
         }
         for (i=0;i<ncoup;i++) {
            sp[i] = i;
            ab[i] = ab2[i];
         }
         lsp = nsp;
         if (lsam==0) {
            lsam++;
         }
      }
      if ((new_sample==0)||(new_sample>current_sample)) {
         if (lsam>nsam) {
            ccleanup(S,membuffers, first);
            throw("Number of samples exceeds maximum");
            closef_(channel);
//            return ERROR;
         }
         double *mm;
         int ii=0;
         for (j=0;j<lsp;j++) {
      	   if (fabs(ab[j]) < zero_tolerance) {
               continue;
            }
            if (fabs(ab[j]-missing_value) < Mat::dTolerance) {
               missing_flag++;
            }
            m[ii] = ab[j];
            sp[ii] = sp[j];
            ii++;
            nitem++;
         }

// SJ added next 2 lines to fix memory leak 15/10/02

         if (type == full && new_sample == 0)
            break;
         mm = new double[ii];
         if (mm==NULL) {
            ccleanup(S,membuffers, first);
            closef_(channel);
            throw("Out of memory allocating storage space");
//            return ERROR;
         }
         memcpy(mm,m,(size_t)(sizeof(double)*ii));
         if (type == full)
            current_sample = new_sample;
         newS = new sample(mm,ii,sp, current_sample);
         if (lsam == 1)
               next = first = newS;
         else {
               next->next = newS;
				   next = next->next;
         }
         lsp=0;
         if (new_sample==0) {
            break;
         }
         else {
            lsam++;
            current_sample=new_sample;
         }
      }
      else if (((new_sample>0)&&(new_sample<current_sample))||(type==full && new_sample==current_sample)) {
         ccleanup(S,membuffers, first);
         sprintf(errorline, "Non sequential sample found after sample %d", current_sample);
         closef_(channel);
         throw(errorline);
//         return ERROR;
      }
      if (type == sparse) {
         for (i=0;i<ncoup;i++) {
            if (sp2[i] > 0) {
               sp[lsp] = sp2[i]-1;
               ab[lsp] = ab2[i];
            }
            else
               continue;
            if (sp2[i] >= nsp) {
               ccleanup(S,membuffers, first);
               sprintf(errorline, "Species number %d in sample %d\nis greater than maximum", sp2[i], new_sample);
               closef_(channel);
               throw(errorline);
//               return ERROR;
            }
            lsp++;
         }
      }
   } while (1);

/* ------- Find highest species number in condensed file
           & check for missing species                  ------ */

   if (type == full) 
      lsam--;
   nsam = lsam;
   int *occur;
   int maxsp=0;
   occur = new int[nsp];
   if (!occur) {
      ccleanup(S,membuffers, first);
      closef_(channel);
      throw("Out of memory allocating workspace");
//      return ERROR;
   }
   membuffers[7] = (void *) occur;
   for(i=0;i<nsp;i++) occur[i]=0;
   next=first;
   for (i=0;i<nsam;i++) {
//      iptr = next->I->dataptr();
      iptr = dataptr(*next->I);
      for (j=0;j<elementCount(*next->I);j++) {
         occur[iptr[j]]++;
	      if (iptr[j] > maxsp) {
            maxsp = iptr[j];
         }
      }
      next = next->next;
   }
//   if ((type==full)&&(maxsp+1!=nsp)) {
//      Message(1,"Warning: highest number variable (species) (%d)\ndoes not agree with number in header\npossibly after deletion species with zero vales values", maxsp+1);
//   }
   nsp = maxsp+1;
   int *DL = new int[nsam];
   if (!DL) {
      ccleanup(S,membuffers, first);
      closef_(channel);
      throw("Out of memory allocating storage for sample numbers");
//      return ERROR;
   }
   S.kill();
   S.setTitle(title);
   S.setsamNo(DL);
   S.setmatType((etf ? full : sparse));
   cMat *CC;
   if (!etf) {
      CC = new cMat(nsam);
      (*CC).setcols(nsp);
      S.setcMat(CC);
      missingValue(*CC) = missing_value;
   }
   else {
      S.setdMat(new dMat(nsam,nsp, impZero));
      missingValue(*getdMat(S)) = missing_value;
   }
   Index *DI = Indexptr(S);
   double **Dm = dataptr(S);
   if (!Dm) {
      ccleanup(S,membuffers, first);
      closef_(channel);
      throw("Memory error in CIN");
//		return ERROR;
   }

   newS = next = first;
   if (!etf) {
      for (i=0;i<nsam;i++) {
        newS = next;
        DI[i] = *newS->I;
        Dm[i] = newS->m;
        S.samNo(i) = newS->samno;
        next = next->next;
        delete newS;
      }
	}
   else {
      for (i=0;i<nsam;i++) {
         newS = next;
			S.samNo(i) = newS->samno;
         int *ii = dataptr(*newS->I);
         for (j=0;j<elementCount(*newS->I);j++) {
				Dm[i][ii[j]] = newS->m[j];
         }
         next = next->next;
         delete[] newS->m;
         delete newS;
      }
   }
   first = NULL;
   char **spName = new char*[nsp];
   char *sss = new char[nsp*9];
   if ((!spName)||(!sss)) {
      ccleanup(S,membuffers, first);
      closef_(channel);
      throw("Out of memory allocating storage for species names");
//      return ERROR;
	}
   for (i=0;i<nsp;i++) {
      spName[i] = sss;
      sss += 9;
   }

   char **samName = new char*[nsam];
	sss = new char[nsam*9];
   if ((!samName)||(!sss)) {
      ccleanup(S,membuffers, first);
      closef_(channel);
      throw("Out of memory allocating storage for sample names");
//      return ERROR;
   }
   for (i=0;i<nsam;i++) {
      samName[i] = sss;
      sss += 9;
   }
   S.setspName(spName);
   S.setsamName(samName);


/* -------------------  Finished reading abundance data
                        now get sample and species names ---------- */

	int nsam2 = int(S.samNo(nsam-1));
   char *samName2 = new char[nsam2*8];
   char *spName2 = new char[nsp*8];
   membuffers[8] = samName2;
   membuffers[9] = spName2;
   getnam_(&channel, spName2, samName2, &nsp, &nsam2, &tag);
   if (tag) {
   	ccleanup(S,membuffers, first);
      closef_(channel);
      if (tag==1) throw("Unexpected end of file reading species names");
      if (tag==2) throw("Error reading species names");
      if (tag==3) throw("Unexpected end of file reading sample names");
      if (tag==4) throw("Error reading sample names");
//      return ERROR;
   }
   char tsam[9];
   int bError = 0;
	for (j=0;j<nsp;j++) {
      strncpy((char *) tsam, &spName2[j*8],8);
      tsam[8] = '\0';
   	RemoveLeadingSpaces(spName[j],tsam);
      strcpy(tsam, spName[j]);
   	RemoveTrailingSpaces(spName[j],tsam);
      if (strlen(spName[j]) == 0)
         bError++;
   }
	for (j=0;j<nsam;j++) {
      int nn = (S.samNo(j)) - 1;
      strncpy((char *) tsam, &samName2[nn*8],8);
      tsam[8] = '\0';
   	RemoveLeadingSpaces(samName[j],tsam);
   	RemoveTrailingSpaces(tsam,samName[j]);
   	strcpy(samName[j], tsam);
      if (strlen(samName[j]) == 0)
         bError++;
   }
   if (bError) {
   	ccleanup(S,membuffers, first);
      sprintf(errorline, "Missing names for %d species / samples.", bError);
      closef_(channel);
      throw(errorline);
//      return ERROR;
   }
/* ---- Allocate and transfer data to permanent array ----- */

	missing_species=0;
	for(i=0;i<nsp;i++)  {
		if (occur[i]==0) {
			missing_species++;
		}
	}
	if (missing_species > 0) {
      if (missing_species == nsp) {
			ccleanup(S,membuffers, first);
      closef_(channel);
			throw("No data found in file");
//			return ERROR;
		}
		char *sptags = new char[nsp];
		if (!sptags) {
			ccleanup(S,membuffers, first);
      closef_(channel);
			throw("Out of memory deleting missing taxa");
//			return ERROR;
		}
		for(i=0;i<nsp;i++) sptags[i]= (occur[i] ? 0 : 1);
		S.deleteCols(sptags);
      ColsWithNoData = missing_species;
      delete [] sptags;
	}
/*
   if (nMissingValues != 0) {
	   if (missing_flag > 0) {
		   	Message(1,"Data contain %d missing value%s", missing_flag, (missing_flag == 1 ? "" : "s"));
	   }
   }
*/
	delete [] buffer;
	delete [] format;
	delete [] m;
	delete [] sp;
	delete [] ab;
	delete [] sp2;
	delete [] ab2;
  delete [] occur;
  delete [] spName2;
  delete [] samName2;
  closef_(channel);
  nMissingValues = missing_flag;
  nCouplets = ncoup;
  return;
}

void ccleanup(dataMat &S, void **memb, sample *first)
{
		if (memb[0])
        delete [] (char *) memb[0];
		if (memb[1])
        delete [] (char *) memb[1];
		if (memb[2])
        delete [] (char *) memb[2];
		if (memb[3])
        delete [] (double *) memb[3];
		if (memb[4])
        delete [] (double *) memb[4];
		if (memb[5])
        delete [] (int *) memb[5];
		if (memb[6])
        delete [] (double *) memb[6];
		if (memb[7])
        delete [] (int *) memb[7];
		if (memb[8])
        delete [] (char *) memb[8];
		if (memb[9])
        delete [] (char *) memb[9];

   sample *next;
   sample *newS = next = first;
   while (newS) {
      next = next->next;
      delete [] newS->m;
      delete newS;
      newS = next;
   }
   S.kill();
}


