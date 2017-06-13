#include <stdlib.h>
#include <stdarg.h>
#include <math.h>     
#include <string.h>
#include "mat.h"

using namespace std;

bool Mat::BoundsCheck = true;              
double Mat::dMissingValue = -99.9;
double Mat::dTolerance = 1.0E-06;
int Mat::maxCols = 10000;
int Mat::maxRows = 20000;
int Mat::maxIndex = 20000;
char Mat::MATVER[10] = "1.00";

void Mat::BoundsError(void)
{
   throw("Index out of range");
}

bool dMat::isExEConformable(const dMat &f) const
{
	if (p->r != f.p->r) {
		if ((p->r != 1) && (f.p->r != 1))
			return false;
	}
	if (p->c != f.p->c) {
		if ((p->c != 1) && (f.p->c != 1))
         return false;
   }
   return true;
}

char dMat::isVectorOrScalar(void) const    // returns 1 = scalar
{                                         //         2 = row vec
   if (p->r==1) {                         //         3 = col vec
		if (p->c==1)
         return 1;
      else
         return 2;
   }
   else if (p->c==1)
      return 3;
   else
      return 0;
}

bool cMat::isExEConformable(const dMat &f) const
{
   if (p->r != f.p->r) {
      if ((p->r != 1) && (f.p->r != 1))
         return false;
   }
   if (p->c != f.p->c) {
		if ((p->c != 1) && (f.p->c != 1))
         return false;
   }
   return true;
}

char cMat::isVectorOrScalar(void) const    // returns 1 = scalar
{                                         //         2 = row vec
   if (p->r==1) {                         //         3 = col vec
      if (p->c==1)
			return 1;
      else
         return 2;
   }
   else if (p->c==1)
      return 3;
   else
      return 0;
}

void sort (double *item, int count)
{
   int i,j,h;
   double v;
   for (h=1;h<=count/9;h=3*h+1);
   for ( ; h>0; h/=3) {
      for (i=h+1; i<=count;i+=1) {
         v = item[i-1];
         j = i;
         while (j>h && item[j-h-1]>v) {
            item[j-1] = item[j-h-1];
            j -= h;
         }
         item[j-1] = v;
      }
   }
}

void sort (char *item, int count)
{
   int i,j,h;
   char v;
   for (h=1;h<=count/9;h=3*h+1);
   for ( ; h>0; h/=3) {
      for (i=h+1; i<=count;i+=1) {
         v = item[i-1];
         j = i;
         while (j>h && item[j-h-1]>v) {
            item[j-1] = item[j-h-1];
            j -= h;
         }
         item[j-1] = v;
      }
   }
}
void sort (int *item, int count)
{
   int i,j,h;
   int v;
   for (h=1;h<=count/9;h=3*h+1);
   for ( ; h>0; h/=3) {
      for (i=h+1; i<=count;i+=1) {
         v = item[i-1];
         j = i;
         while (j>h && item[j-h-1]>v) {
            item[j-1] = item[j-h-1];
            j -= h;
         }
         item[j-1] = v;
      }
   }
}

Index::Index()
{
   p = new IndexRep;
   if (!p)
      throw("Error: Out of memory in Index");
   p->I= NULL;
   p->n=0;
   p->refs=1;
#ifdef MATDEBUG
   nindexes++;
#endif
}

Index::Index(int n)
{
   p = new IndexRep;
   if (!p)
      throw("Out of memory in Index(int)");
   p->I = new int[n];
   if (!p->I)
      throw("Out of memory in Index(int)");
   for (int i=0;i<n;i++)
      p->I[i] = i;
   p->n=n;
   p->refs=1;
#ifdef MATDEBUG
   nindexes++;
#endif
}

Index::Index(int n, int initvalue)
{
   p = new IndexRep;
   if (!p)
      throw("Out of memory in Index(int)");
   p->I = new int[n];
   if (!p->I)
      throw("Out of memory in Index(int)");
   for (int i=0;i<n;i++)
      p->I[i] = initvalue;
   p->n=n;
   p->refs=1;
#ifdef MATDEBUG
   nindexes++;
#endif
}

Index::Index(int n, int *initvalues)
{
   p = new IndexRep;
   if (!p)
      throw("Error: Out of memory in Index(int, int*)");
   p->I= new int[n];
   if (!p->I)
      throw("Error: Out of memory in Index(int, int *)");
   memcpy(p->I,initvalues,(size_t)(sizeof(int)*n));
   p->n=n;
   p->refs=1;
#ifdef MATDEBUG
   nindexes++;
#endif
}

Index::~Index()
{
   if (--p->refs==0) {
      if (p->I) {
//         delete p->I;
         delete[] p->I;
      }
      delete p;
#ifdef MATDEBUG
      nindexes--;
#endif
   }
}


Index Index::operator=(const Index &I)
{
   if (--p->refs==0) {
      if (p->I)
         delete p->I;
      delete p;
#ifdef MATDEBUG
      nindexes--;
#endif
   }
   p = I.p;
   I.p->refs++;
   return *this;
}

Index copy(const Index &I)
{
   Index target(elementCount(I),dataptr(I));
   return target;
}

Index fsort(dMat &f)
{
   Index I(rows(f));
   Index T(rows(f));
   int i;
   for (i=0;i<rows(f);i++)
	T(i) = 1;
   for (i=0;i<rows(f);i++) {
      double mini = 1.0E30;
      int cur=0;
      for (int j=0;j<rows(f);j++) {
	      if (T(j) && (f(j, 0) < mini)) {
	         mini = f(j, 0);
	         cur = j;
	      }
      }
		T(cur) = 0;
      I(i) = cur;
   }
   return I;
}

int rows(const dataMat &D)
{
   dataMatRep *DMP = D.p;
   if (DMP->C)
      return rows(*(DMP->C));
   else if (DMP->F)
      return rows(*(DMP->F));
   else
      return 0;
}

int cols(const dataMat &D)
{
   dataMatRep *DMP = D.p;
   if (DMP->C)
      return cols(*(DMP->C));
   else if (DMP->F)
      return cols(*(DMP->F));
   else
      return 0;
}

double ** dataptr(const dataMat &D)
{
   dataMatRep *DMP = D.p;
   if (DMP->C)
      return dataptr(*(DMP->C));
   else if (DMP->F)
      return dataptr(*(DMP->F));
   else
      return NULL;
}

Index * Indexptr(const dataMat &D)
{
   if ((mattype(D.p->mType) == sparse)&&(D.p->C!=NULL))
      return D.p->C->p->I;
   else
      return NULL;
}

dataMat::dataMat()
{
   p = new dataMatRep;
   if (!p)
      throw("Out of memory in dataMat");
   p->spNam=NULL;
   p->samNam=NULL;
   p->samNum=NULL;
   p->title=NULL;
   p->C=NULL;
   p->F=NULL;
   p->mType=0;
   p->refs=1;
}

dataMat::~dataMat()
{
   if (--p->refs==0)
      kill();
   delete p;
}

dataMat::dataMat(const dataMat &D)
{
   p->refs++;
   p = D.p;
}

void dataMat::kill()
{
	if (p->spNam) {
		if (p->spNam[0])
//			delete p->spNam[0];
			delete[] p->spNam[0];
//		delete p->spNam;
		delete[] p->spNam;
		p->spNam=NULL;
	}
	if (p->samNam) {
		if (p->samNam[0])
//			delete p->samNam[0];
			delete[] p->samNam[0];
//		delete p->samNam;
		delete[] p->samNam;
      p->samNam=NULL;
   }
   if (p->samNum) {
//      delete p->samNum;
      delete[] p->samNum;
      p->samNum=NULL;
   }
   if (p->C) {
      delete p->C;
      p->C=NULL;
   }
   if (p->F) {
      delete p->F;
      p->F=NULL;
   }
   if (p->title) {
//      delete p->title;
      delete[] p->title;
      p->title=NULL;
   }
   p->mType = mattype(undefined);
}

double & missingValue(dataMat &d)
{
	if (d.p->C)
		return missingValue(*d.p->C);
	else if (d.p->F)
		return missingValue(*d.p->F);
	else
		throw("Trying to set missingvalue in empty dMat");
	return missingValue(*d.p->F);	
}

bool dataMat::deleteRows(char *ii)
{
   int r = rows(*this);
   int *n = new int[r];
   if (!n)
      throw("Out of memory in deleteRows dataMat");
   if (mattype(p->mType)==sparse) {
      if (p->C)
         p->C->deleteRows(ii);
      else
         return false;
   }
   else if (mattype(p->mType)==full) {
      if (p->F)
         p->F->deleteRows(ii);
      else
        return false;
   }
   else
      return false;
   int i, j=0;
   for(i=0;i<r;i++) {
      if (!ii[i]) {
         n[i] = i-j;
      }
      else
         j++;
   }
   j = r - j;
   if (p->samNam) {
      char **s = new char*[j];
      if (s==NULL)
         return 1;
      for (i=0;i<r;i++) {
         if (!ii[i])
           s[n[i]] = p->samNam[i];
      }
//      delete p->samNam;
      delete[] p->samNam;
      p->samNam = s;
   }
   if (p->samNum) {
      int *l = new int[j];
      if (l==NULL)
         return 1;
      for (i=0;i<r;i++) {
         if (!ii[i])
            l[n[i]] = p->samNum[i];
      }
//      delete p->samNum;
      delete[] p->samNum;
      p->samNum = l;
   }
   delete n;
   return true;
}

bool dataMat::deleteCols(char *ii)
{
   int c = cols(*this);
   int *n = new int[c];
   if (!n)
      throw("Out of memory in deleteCols dataMat");
   if (mattype(p->mType)==sparse) {
      if (p->C)
          p->C->deleteCols(ii);
      else
	 return false;
   }
   else if (mattype(p->mType)==full) {
      if (p->F)
         p->F->deleteCols(ii);
      else
         return 1;
   }
   else
      return 1;
   int i, j=0;
   for(i=0;i<c;i++) {
      if (!ii[i]) {
         n[i] = i-j;
      }
      else
         j++;
   }
   j = c - j;
   if (p->spNam) {
      char **s = new char*[j];
      if (s==NULL)
         return 1;
      for (i=0;i<c;i++) {
         if (!ii[i])
            s[n[i]] = p->spNam[i];
      }
//      delete p->spNam;
      delete[] p->spNam;
      p->spNam = s;
   }
   delete n;
   return true;
}

#define MAX(a, b) (a > b ? a : b) 

// constructors

// create an 'empty' fMat

dMat::dMat()
{
   p = new dMatRep;
   if (!p)
      throw("Error: Out of memory in fMat");
   p->m = NULL;
   p->r = 0;
   p->c = 0;
   p->refs = 1;
   parent=NULL;
	p->missingValue = dMissingValue;
#ifdef MATDEBUG
   nmats++;
   mats++;
   p->matNo = mats;
#endif
}

// create a new fMat and fill it with initvals

dMat::dMat(int r, int c, double initval)
{
   if (r > maxRows) 
      throw("nRows too high in dMat::dMat(r, c)");
   if (c > maxCols) 
      throw("nCols too high in dMat::dMat(r, c)");
   p = new dMatRep;
   if (!p)
      throw("Error: Out of memory in dMat");
   p->r = r;
   p->c = c;
   p->m = new double*[r];
   if (!p->m)
      throw("Out of memory in dMat");
   if (c==1) {
      double *m = new double[r];
		if (!m)
         throw("Out of memory in dMat");
      p->m[0] = m;
      for (int i=0;i<r;i++) {
         m[i] = initval;
         p->m[i] = &m[i];
      }
   }
   else {
      for (int i=0;i<r;i++) {
         double *m = new double[c];
         if (!m)
            throw("Out of memory in dMat");
         p->m[i] = m;
         if (!i)                                  // if nrow > 1, then just
            for (int j=0;j<c;j++)                 // memcpy row 0 - slightly
               *m++ = initval;                    // faster
         else
				memcpy(m,p->m[0],(size_t)(sizeof(double)*c));
      }
   }
   p->missingValue = dMissingValue;
   p->refs=1;
   parent=NULL;
#ifdef MATDEBUG
   nmats++;
   mats++;
   p->matNo = mats;
#endif
}

// create a new fMat and fill it with values pointed to by initvalues
// values in initvalues stored columnwise

dMat::dMat(int r, int c, double *initvalues)
{
   if (r > maxRows) 
      throw("nRows too high in dMat::dMat(r, c)");
   if (c > maxCols) 
      throw("nCols too high in dMat::dMat(r, c)");
   p = new dMatRep;
   if (!p)
      throw("Out of memory in dMat");
   p->r = r;
	p->c = c;
   p->m = new double*[r];
   if (!p->m)
      throw("Out of memory in dMat");
   if (c==1) {
      double *m = new double[r];
      if (!m)
         throw("Out of memory in dMat");
      p->m[0] = m;
      for (int i=0;i<r;i++) {
         m[i] = initvalues[i];
         p->m[i] = &m[i];
      }
   }
   else {
      for (int i=0;i<r;i++) {
         p->m[i] = new double[c];
         if (!p->m[i])
				throw("Out of memory in dMat");
         memcpy(p->m[i],&initvalues[i*c],(size_t)(sizeof(double)*c));
      }
   }
   p->missingValue = dMissingValue;
   p->refs=1;
   parent=NULL;
#ifdef MATDEBUG
   nmats++;
   mats++;
   p->matNo = mats;
#endif
}

// destructor

dMat::~dMat()
{
   if (--p->refs==0) {
      if (parent==NULL) {
         if (p->c==1) {
//            delete p->m[0];
            delete[] p->m[0];
         }
			else {
            for(int i=0;i<p->r;i++) {
//               delete p->m[i];
               delete[] p->m[i];
            }
         }
      }
      else if (--parent->refs==0) {
         if (parent->m) {
            if (parent->c==1) {
//               delete parent->m[0];
               delete[] parent->m[0];
            }
            else {
               for(int i=0;i<parent->r;i++) {
//                  delete parent->m[i];
                  delete[] parent->m[i];
               }
            }
            if (parent->m)
//               delete parent->m;
               delete[] parent->m;
			}
         delete parent;
#ifdef MATDEBUG
         nmats--;
#endif
      }
      if (p->m)
//         delete p->m;
         delete[] p->m;
      delete p;
#ifdef MATDEBUG
      nmats--;
#endif
   }
   else if (parent)
      --parent->refs;
}

dMat dMat::operator=(const dMat &f)     // makes a reference to f
{                                         // not a new copy
   if (--p->refs==0) {
      if (parent==NULL) {
         if (p->c==1) {
//				delete p->m[0];
				delete[] p->m[0];
         }
         else {
            for(int i=0;i<p->r;i++) {
//               delete p->m[i];
               delete[] p->m[i];
            }
         }
      }
      else if (--parent->refs==0) {
         if (parent->m) {
            if (parent->c==1) {
//               delete parent->m[0];
               delete[] parent->m[0];
            }
            else {
               for(int i=0;i<parent->r;i++) {
//                  delete parent->m[i];
                  delete[] parent->m[i];
               }
            }
//				delete parent->m;
				delete[] parent->m;
#ifdef MATDEBUG
            nmats--;
#endif
         }
         delete parent;
      }
      if (p->m)
//         delete p->m;
         delete [] p->m;
      delete p;
#ifdef MATDEBUG
      nmats--;
#endif
      p = NULL;
   }
   else if (parent)
      --parent->refs;
   p = f.p;
   parent=f.parent;
   f.p->refs++;
   if (f.parent)
		f.parent->refs++;
	return *this;
}

// method to return an fMat which is contains a subset of rows
// it references the original data, so the parent pointer is set
// to *this and this->refs is incremented
// but if *this is a col vec then a new copy is made, since the fast
// vector access would not work

dMat dMat::operator()(const Index &I, int dim) const
{
	int n = elementCount(I);
	int *ii = dataptr(I);
	if (!dim) {
		if (p->c==1) {            // we have a col vec so make a new copy
			dMat target(n);
			double *m1 = target.p->m[0], *m2 = p->m[0];
			for (int i=0;i<n;i++)
				m1[i] = m2[*ii++];
			missingValue(target) = p->missingValue;
			return target;
		}
		else {
			dMat target;
			target.p->r = n;
			target.p->c = p->c;
			target.p->m = new double*[n];
			if (!target.p->m)
				throw("Out of memory in fMat");
			for (int i=0;i<n;i++) {
				target.p->m[i] = p->m[*ii++];
			}
			target.parent=p;
			target.parent->refs++;
			missingValue(target) = p->missingValue;
			return target;
		}
	}
	else {
      if (elementCount(I) > p->c) {
         throw("Column index out of bounds in dMat(Index, dir)");
         dMat d;
         return d;
      }
		dMat target(p->r, elementCount(I));
		double **m = p->m;
		double **m1 = target.p->m;
		ii = dataptr(I);
		for (int j=0;j<elementCount(I);j++) {
			 for (int i=0;i<p->r;i++) 
				 m1[i][j] = m[i][*ii];
			 ii++;
		}
		missingValue(target) = p->missingValue;
		return target;
	}
}

// binary operators *(), +(), -(), and /()
// all are a bit large, but needed to cope with the various
// combinations of array types (mat*mat, colvec*mat, rowvec*mat etc.
// they should behave in the same way as the equivalent M++ operators
// comments for -, / & + are as for *, except - and / have to have
// separate parts for eg. f1 - f2 and f2 - f1



// performs as M++ concat ie:
// dir = 0 - concatenated to rows
// dir = 1 - concatenated to cols

dMat dMat::concat(const dMat &f, int dir)
{
   int i;
   if (dir==0) {
      if (p->c+f.p->c > maxCols) 
         throw("nCols too high in dMat::concat");
      if (p->c == f.p->c) {
         dMat target(p->r+f.p->r,p->c);
         double **m = target.p->m;
         double **m1 = p->m;
         for(i=0;i<p->r;i++)
            memcpy(m[i],m1[i],(size_t) (sizeof(double) * p->c));
         m1 = f.p->m;
         for(i=0;i<f.p->r;i++)
            memcpy(m[i+p->r],m1[i],(size_t) (sizeof(double) * p->c));
         return target;
      }
      else
         throw("Arrays are not conformable in function concat");
   }
   else if (dir==1) {
      if (p->r+f.p->r > maxRows) 
         throw("nRows too high in dMat::concat");
      if (p->r == f.p->r) {
         dMat target(p->r,p->c+f.p->c);
         double **m = target.p->m;
         double **m1 = p->m;
         for(i=0;i<p->r;i++)
            memcpy(m[i],m1[i],(size_t) (sizeof(double) * p->c));
         m1 = f.p->m;
         for(i=0;i<f.p->r;i++)
            memcpy(&m[i][p->c],m1[i],(size_t) (sizeof(double) * f.p->c));
         return target;
      }
      else
         throw("Arrays are not conformable in function concat");
   }
   else
      throw("Integer out of range in concat (must be 0 or 1)");         
   return *this;   
}

void dMat::merge(const dMat &f, int dir)
{
	int i;
	if (dir==0) {
      if (p->r+f.p->r > maxRows) 
         throw("nRows too high in dMat::merge");
		if (p->c == f.p->c) {
			double **m1 = p->m;
			double **m = new double*[p->r+f.p->r];
			if (!m)
				throw("Out of memory in merge(dMat &)");
			if (p->c==1) {
				double *ff = new double[p->r+f.p->r];
				if (!ff)
					throw("Out of memory in merge(dMat &)");
				m[0] = ff;
				memcpy(ff, m1[0], (size_t) (sizeof(double) * p->r));
//				delete m1[0];
				delete[] m1[0];
				memcpy(&ff[p->r], f.p->m[0], (size_t) (sizeof(double) * f.p->r));
				for(i=0;i<p->r+f.p->r;i++)
					m[i] = &ff[i];
			}
			else {
				for(i=0;i<p->r;i++)
					m[i] = m1[i];
				m1 = f.p->m;
				for(i=0;i<f.p->r;i++) {
					double *ff = new double[p->c];
					if (!ff)
						throw("Out of memory in merge(dMat &)");
					m[i+p->r] = ff;
					memcpy(m[i+p->r],m1[i],(size_t) (sizeof(double) * p->c));
				}
			}
//			delete p->m;
			delete[] p->m;
			p->m = m;
			p->r += f.p->r;
		}
		else
         throw("Arrays are not conformable in function merge");
   }
   else if (dir==1) {
      if (p->c+f.p->c > maxCols) 
         throw("nRows too high in dMat::merge");
      if (p->r == f.p->r) {
         double **m = p->m;
         double **m1 = f.p->m;
         double *m2 = p->m[0];
         for(i=0;i<p->r;i++) {
            double *ff = new double[p->c+f.p->c];
            if (!ff)
               throw("Out of memory in merge(fMat &)");
            if (p->c>1) {
               memcpy(ff,m[i],(size_t) (sizeof(double) * p->c));
               delete m[i];
            }
            else
               ff[0] = m[i][0];
            m[i]=ff;
            if (f.p->c>1)
               memcpy(&m[i][p->c],m1[i],(size_t) (sizeof(double) * f.p->c));
            else
               m[i][p->c] = m1[i][0];
         }
         if (p->c==1)
            delete m2;
         p->c += f.p->c;
      }
      else
         throw("Arrays are not conformable in function merge");
   }
   else
      throw("Integer out of range in merge (must be 0 or 1)");
}

dMat copy(const dMat &f)     // makes a direct copy of f
{
   dMat target(f.p->r,f.p->c);
   double **m = target.p->m;
   double **m1 = f.p->m;
   if (f.p->c==1)
      memcpy(m[0],m1[0],(size_t) sizeof(double)*f.p->r);
   else {
      for (int i=0;i<f.p->r;i++) {
         memcpy(m[i],m1[i],(size_t) sizeof(double)*f.p->c);
      }
   }
   target.p->missingValue = f.p->missingValue;
   return target;
}

// friend function operators

void dMat::fill(double f)
{
   double **m = p->m;
   for (int i=0;i<p->r;i++) {
      for (int j=0;j<p->c;j++) {
         m[i][j] = f;
      }
   }
}

const char mes1[] = "\nError: Cannot deleteRows - fMat has references";
const char mes2[] = "\nError: Out of memory in fMat";

int dMat::deleteRows(char *ii)
{
   if (p->refs > 1)
      throw(mes1);
   long *n = new long[p->r];
   long i;
   if (!n)
      throw(mes2);
   long j=0;
   for(i=0;i<p->r;i++) {
      if (!ii[i]) {
         n[i] = i-j;
      }
      else
         j++;
   }
   j = p->r-j;
   double **m = new double*[j];
   if (m==NULL) {
      throw(mes2);
   }
   if (p->c == 1) {
      double *mm = new double[j];
		if (!mm)
         throw(mes2);
      for (i=0;i<j;i++)
         m[i] = &mm[i];
      for (i=0;i<p->r;i++) {
         if(!ii[i])  {
            mm[n[i]] = p->m[i][0];
         }
      }
//      delete[] p->m[0];
      delete p->m[0];
   }
   else {
      for(i=0;i<p->r;i++) {
         if(ii[i]) {
//            delete p->m[i];
            delete[] p->m[i];
         }
         else {
            m[n[i]] = p->m[i];
         }
      }
   }
   delete n;
//   delete[] p->m;
   delete p->m;
   p->m = m;
   p->r = j;
   return 0;
}

int dMat::deleteCols(char *ii)
{
   if (p->refs > 1)
      throw(mes1);
   long *n = new long[p->c];
   if (!n)
      return 1;
   long i, j=0;
   for(i=0;i<p->c;i++) {
      if (!ii[i]) {
	 n[i] = i-j;
      }
      else
	 j++;
   }
   j = p->c-j;
   if (j==1) {
      double *m = new double[j];
      for (long k=0;k<p->c;k++) {
	      if (!ii[k]) {
            for(i=0;i<p->r;i++)
	            m[n[k]] = p->m[i][k];
         }
      }
//      delete p->m[i];
      delete[] p->m[i];
      for(i=0;i<p->r;i++)
         p->m[i] = &m[i];
   }
   else {
      for(i=0;i<p->r;i++) {
         double *m = new double[j];
         if (!m)
   	      throw(mes2);
         for (long k=0;k<p->c;k++) {
   	      if (!ii[k])
	         m[n[k]] = p->m[i][k];
         }
//         delete p->m[i];
         delete[] p->m[i];
         p->m[i] = m;
      }
   }
   p->c = j;
   delete n;
   return 0;
}


dMat operator/(const dMat &f1, const dMat &f2)
{
   int t1, t2;
   if (!f1.isExEConformable(f2))
      throw("Arrays are not binary conformable in operator/()");
   t1=f1.isVectorOrScalar();
   t2=f2.isVectorOrScalar();
   if (t1==1)                                   // f1 is a scalar
      return (f1.p->m[0][0] / f2);
   if (t2==1)                                   // f2 is a scalar
      return (f1 / f2.p->m[0][0]);
   int r = MAX(f1.p->r, f2.p->r);
   int c = MAX(f1.p->c, f2.p->c);
   dMat target(r,c);
   double **m = target.p->m;
   double **m1 = f1.p->m;
   double **m2 = f2.p->m;                    // f1 and f2 have same
   if (f1.isBinaryConformable(f2)) {        // dimensions - simple divide
      if (c==1) {
         for (int i=0;i<r;i++) {
            m[0][i] = m1[0][i] / m2[0][i];
         }
      }
      else {
         for (int i=0;i<r;i++) {
            for (int j=0;j<c;j++) {
               m[i][j] = m1[i][j] / m2[i][j];
            }
         }
      }
      return target;
   }
   if (t1==2) {                      // f1 or f2 = row vector
      if (t2==3) {
         for (int i=0;i<r;i++) {
            for(int j=0;j<c;j++) {
               m[i][j] = m1[0][j] / m2[i][0];
            }
         }
         return target;
      }
      else {
         for (int i=0;i<r;i++) {
            for(int j=0;j<c;j++) {
               m[i][j] = m1[0][j] / m2[i][j];
            }
         }
         return target;
      }
   }
   if (t2==2) {                      // f2 = row vector
      if (t1==3) {
         for (int i=0;i<r;i++) {
            for(int j=0;j<c;j++) {
               m[i][j] = m1[i][0] / m2[0][j];
            }
         }
         return target;
      }
      else {
         for (int i=0;i<r;i++) {
            for(int j=0;j<c;j++) {
               m[i][j] = m1[i][j] / m2[0][j];
            }
         }
         return target;
      }
   }
   if (t1==3) {                      // f1 = col vector
      for (int i=0;i<r;i++) {
         for(int j=0;j<c;j++) {
            m[i][j] = m1[i][0] / m2[i][j];
         }
      }
      return target;
   }
   if (t2==3) {                      // f2 = col vector
      for (int i=0;i<r;i++) {
         for(int j=0;j<c;j++) {
            m[i][j] = m1[i][j] / m2[i][0];
         }
      }
      return target;
   }                                      
   return target;
}


dMat operator/(const dMat &f, double fl)
{
   dMat target(f.p->r,f.p->c);
   double **m = target.p->m;
   double **m1 = f.p->m;
   if (f.isVectorOrScalar()) {
      int n = MAX(f.p->r,f.p->c);
      for (int i=0;i<n;i++)
         m[0][i] = m1[0][i] / fl;
   }
   else {
      for (int i=0;i<f.p->r;i++) {
         for (int j=0;j<f.p->c;j++) {
            m[i][j] = m1[i][j] / fl;
         }
      }
   }
   return target;
}

dMat operator/(double fl, const dMat &f)
{
   dMat target(f.p->r,f.p->c);
   double **m = target.p->m;
   double **m1 = f.p->m;
   if (f.isVectorOrScalar()) {
      int n = MAX(f.p->r,f.p->c);
      for (int i=0;i<n;i++)
         m[0][i] = fl / m1[0][i];
   }
   else {
      for (int i=0;i<f.p->r;i++) {
         for (int j=0;j<f.p->c;j++) {
            m[i][j] = fl / m1[i][j];
         }
      }
   }
   return target;
}

dMat & dMat::operator/=(double fl)
{
   double **m = p->m;
   for (int i=0;i<p->r;i++) {
      for (int j=0;j<p->c;j++) {
         m[i][j] = m[i][j] / fl;
      }
   }
   return *this;
}

dMat & dMat::operator/=(const dMat &f)
{
   int t1, t2;
   if (!isExEConformable(f))
		throw("Arrays are not conformable in fMat::operator/=()");
	t1=isVectorOrScalar();
	t2=f.isVectorOrScalar();
	if (t2==1)                                   // f is a scalar
		return (*this /= f.p->m[0][0]);
	if (t2==1)                                   // this is a scalar
		throw("lhs is a scalar, rhs is not, in operator/=()");
	double **m = p->m;
	double **m1 = f.p->m;
	if (isBinaryConformable(f)) {        // this and f have same
		for (int i=0;i<p->r;i++) {        // dimensions - simple subtraction
			for (int j=0;j<p->c;j++) {
				m[i][j] = m[i][j] / m1[i][j];
			}
		}
		return *this;
	}
	if ((t1==2)||(t1==3)) {        // lhs is a row or col vector, rhs is not
		throw("lhs is a row or col vector, rhs is not, in operator/=()");
	}
	if (t2==2) {                   // rhs is a row vec, lhs a mat
		for (int i=0;i<p->r;i++) {
         for(int j=0;j<p->c;j++) {
            m[i][j] = m[i][j] / m1[0][j];
         }
      }
		return *this;
   }
   if (t2==3) {                   // rhs is a col vec, lhs a mat
		for (int i=0;i<p->r;i++) {
         for(int j=0;j<p->c;j++) {
            m[i][j] = m[i][j] / m1[i][0];
         }
		}
		return *this;
   }            
   return *this;
}


dMat operator*(const dMat &f1, const dMat &f2)
{
   int t1, t2;
   if (!f1.isExEConformable(f2))
      throw("Arrays are not binary conformable in operator*()");
   t1=f1.isVectorOrScalar();
   t2=f2.isVectorOrScalar();
   if (t1==1)                                   // f1 is a scalar
      return (f2 * f1.p->m[0][0]);
   if (t2==1)                                   // f2 is a scalar
      return (f1 * f2.p->m[0][0]);
   int r = MAX(f1.p->r, f2.p->r);
   int c = MAX(f1.p->c, f2.p->c);
   dMat target(r,c);
   double **m = target.p->m;
   double **m1 = f1.p->m;
   double **m2 = f2.p->m;                    // f1 and f2 have same
   if (f1.isBinaryConformable(f2)) {        // dimensions - simple multiply
      if (c==1) {                           // col vectors
         for (int i=0;i<r;i++) {
            m[0][i] = m1[0][i] * m2[0][i];
         }
      }
      else {
         for (int i=0;i<r;i++) {
            for (int j=0;j<c;j++) {
               m[i][j] = m1[i][j] * m2[i][j];
            }
         }
      }
      return target;
   }
   if ((t1==2)||(t2==2)) {                   // f1 or f2 is a row vector
      if (t2==2) {                           // f1 = row vec, f2 = mat
         m2 = f1.p->m;
         m1 = f2.p->m;
      }
      if ((t1==3)||(t2==3)) {
         for (int i=0;i<r;i++) {             // f1 or f2 is a col vector
            for(int j=0;j<c;j++) {
               m[i][j] = m1[0][j] * m2[i][0];
            }
         }
         return target;
      }
      else {                                 // vec * mat, the swap of
         for (int i=0;i<r;i++) {             // m1 & m2 above put them in
            for(int j=0;j<c;j++) {
               m[i][j] = m1[0][j] * m2[i][j];
            }
         }
         return target;
      }
   }
   if ((t1==3)||(t2==3)) {                   // f1 or f2 = col vector
      if (t2==3) {                           // if f2 col vec swap
         m2 = f1.p->m;
         m1 = f2.p->m;
      }
      for (int i=0;i<r;i++) {                // col vec * mat
         for(int j=0;j<c;j++) {
            m[i][j] = m1[i][0] * m2[i][j];
         }
      }
      return target;
   }                                  
   dMat d;
   return d;
}

dMat operator*(const dMat &f, double fl)
{
   dMat target(f.p->r,f.p->c);
   double **m = target.p->m;
   double **m1 = f.p->m;
   if (f.isVectorOrScalar()) {
      int n = MAX(f.p->r,f.p->c);
      for (int i=0;i<n;i++)
         m[0][i] = m1[0][i] * fl;
   }
   else {
      for (int i=0;i<f.p->r;i++) {
         for (int j=0;j<f.p->c;j++) {
            m[i][j] = m1[i][j] * fl;
         }
      }
   }
   return target;
}

dMat & dMat::operator*=(double fl)
{
   double **m = p->m;
   for (int i=0;i<p->r;i++) {
      for (int j=0;j<p->c;j++) {
         m[i][j] = m[i][j] * fl;
      }
   }
   return *this;
}

dMat & dMat::operator*=(const dMat &f)
{
   int t1, t2;
   if (!isExEConformable(f))
      throw("Arrays are not conformable in fMat::operator*=()");
   t1=isVectorOrScalar();
   t2=f.isVectorOrScalar();
   if (t2==1)                                   // f is a scalar
      return (*this *= f.p->m[0][0]);
   if (t2==1)                                   // this is a scalar
      throw("lhs is a scalar, rhs is not, in operator*=()");
   double **m = p->m;
   double **m1 = f.p->m;
   if (isBinaryConformable(f)) {        // this and f have same
      for (int i=0;i<p->r;i++) {           // dimensions - simple multiply
         for (int j=0;j<p->c;j++) {
            m[i][j] = m[i][j] * m1[i][j];
         }
      }
      return *this;
   }
   if ((t1==2)||(t1==3)) {        // lhs is a row or col vector, rhs is not
      throw("lhs is a row or col vector, rhs is not, in operator*=()");
   }
   if (t2==2) {                   // rhs is a row vec, lhs a mat
      for (int i=0;i<p->r;i++) {
         for(int j=0;j<p->c;j++) {
            m[i][j] = m[i][j] * m1[0][j];
         }
      }
      return *this;
   }
   if (t2==3) {                   // rhs is a col vec, lhs a mat
      for (int i=0;i<p->r;i++) {
         for(int j=0;j<p->c;j++) {
            m[i][j] = m[i][j] * m1[i][0];
         }
      }
      return *this;
   }
   return *this;
}                     

dMat operator+(const dMat &f1, const dMat &f2)
{
   int t1, t2;
   if (!f1.isExEConformable(f2))
      throw("Arrays are not binary conformable in operator+()");
   t1=f1.isVectorOrScalar();
   t2=f2.isVectorOrScalar();
   if (t1==1)                                   // f1 is a scalar
      return (f2 + f1.p->m[0][0]);
   if (t2==1)                                   // f2 is a scalar
      return (f1 + f2.p->m[0][0]);
   int r = MAX(f1.p->r, f2.p->r);
   int c = MAX(f1.p->c, f2.p->c);
   dMat target(r,c);
   double **m = target.p->m;
   double **m1 = f1.p->m;
   double **m2 = f2.p->m;                     // f1 and f2 have same
   if (f1.isBinaryConformable(f2)) {         // dimensions - simple add
      if (c==1) {
         for (int i=0;i<r;i++) {
            m[0][i] = m1[0][i] + m2[0][i];
         }
      }
      else {
         for (int i=0;i<r;i++) {
            for (int j=0;j<c;j++) {
               m[i][j] = m1[i][j] + m2[i][j];
            }
         }
      }
      return target;
   }
   if ((t1==2)||(t2==2)) {                      // f1 or f2 = row vector
      if (t2==2) {                              // f1 = row vec, f2 = mat
         m2 = f1.p->m;
         m1 = f2.p->m;
      }
      if ((t1==3)||(t2==3)) {
         for (int i=0;i<r;i++) {
            for(int j=0;j<c;j++) {
               m[i][j] = m1[0][j] + m2[i][0];
            }
         }
         return target;
      }
      else {
         for (int i=0;i<r;i++) {
            for(int j=0;j<c;j++) {
               m[i][j] = m1[0][j] + m2[i][j];
            }
         }
         return target;
      }
   }
   if ((t1==3)||(t2==3)) {                      // f1 or f2 = col vector
      if (t2==3) {
         m2 = f1.p->m;
         m1 = f2.p->m;
      }
      for (int i=0;i<r;i++) {
         for(int j=0;j<c;j++) {
            m[i][j] = m1[i][0] + m2[i][j];
         }
      }
      return target;
   }
   dMat d;
   return d;
}

dMat operator+(const dMat &f, double fl)
{
   dMat target(f.p->r,f.p->c);
   double **m = target.p->m;
   double **m1 = f.p->m;
   if (f.isVectorOrScalar()) {
      int n = MAX(f.p->r,f.p->c);
      for (int i=0;i<n;i++)
         m[0][i] = m1[0][i] + fl;
   }
   else {
      for (int i=0;i<f.p->r;i++) {
         for (int j=0;j<f.p->c;j++) {
            m[i][j] = m1[i][j] + fl;
         }
      }
   }
   return target;
}

dMat & dMat::operator+=(const dMat &f)
{
   int t1, t2;
	if (!isExEConformable(f))
      throw("Arrays are not conformable in fMat::operator+=()");
   t1=isVectorOrScalar();
   t2=f.isVectorOrScalar();
   if (t2==1)                                   // f is a scalar
		return (*this += f.p->m[0][0]);
   if (t2==1)                                   // this is a scalar
      throw("lhs is a scalar, rhs is not, in operator+=()");
   double **m = p->m;
   double **m1 = f.p->m;
   if (isBinaryConformable(f)) {        // this and f have same
      for (int i=0;i<p->r;i++) {        // dimensions - simple addition
         for (int j=0;j<p->c;j++) {
            m[i][j] = m[i][j] + m1[i][j];
         }
      }
		return *this;
   }
   if ((t1==2)||(t1==3)) {        // lhs is a row or col vector, rhs is not
      throw("lhs is a row or col vector, rhs is not, in operator+=()");
   }
   if (t2==2) {                   // rhs is a row vec, lhs a mat
      for (int i=0;i<p->r;i++) {
         for(int j=0;j<p->c;j++) {
            m[i][j] = m[i][j] + m1[0][j];
         }
      }
		return *this;
   }
   if (t2==3) {                   // rhs is a col vec, lhs a mat
      for (int i=0;i<p->r;i++) {
			for(int j=0;j<p->c;j++) {
            m[i][j] = m[i][j] + m1[i][0];
         }
      }
		return *this;
   }
   return *this;
}


dMat & dMat::operator+=(double fl)
{
   double **m = p->m;
   for (int i=0;i<p->r;i++) {
      for (int j=0;j<p->c;j++) {
         m[i][j] = m[i][j] + fl;
      }
   }
	return *this;
}

dMat dMat::product(const dMat &f)
{
   if (isVectorOrScalar()==1)
      return (f * p->m[0][0]);
   else if (f.isVectorOrScalar()==1)
      return (*this * f.p->m[0][0]);
   else if (!innerConformingDimensions(f))
      throw("Dimensions are not conformable in function product");
   else {
      dMat target(p->r,f.p->c);
      if (f.p->c==1) {
         double *m = target.p->m[0];
         double *m2 = f.p->m[0];
         for (int i=0;i<p->r;i++) {
            double sum = 0.0;
            double *m1 = p->m[i];
            for (int j=0;j<p->c;j++)
               sum += m1[j] * m2[j];
            m[i] = sum;
         }
      }
      else if (p->c==1) {
         double **m = target.p->m;
         double *m1 = p->m[0];
         double *m2 = f.p->m[0];
         for (int i=0;i<p->r;i++) {
            for (int j=0;j<f.p->c;j++)
               m[i][j] = m1[i] * m2[j];
         }
      }
      else {
         double **m = target.p->m;
         double **m1 = p->m;
         double **m2 = f.p->m;
         for (int i=0;i<p->r;i++) {
            for (int j=0;j<f.p->c;j++) {
               double sum = 0.0;
               double *m3 = m1[i];
               double **m4 = m2;
               for (int k=0;k<p->c;k++)
                  sum += *m3++ * *(*m4++ + j);
               m[i][j] = sum;
            }
         }
      }
      return target;
   }                        
   dMat d;
   return d;
}

dMat dMat::tproduct(const dMat &f)
{
   if (isVectorOrScalar()==1)
      return (f * p->m[0][0]);
   else if (f.isVectorOrScalar()==1)
      return (*this * f.p->m[0][0]);
   else if (p->r != f.p->r)
      throw("Dimensions are not conformable in function tproduct");
   else {
      dMat target(p->c,f.p->c);
      if (f.p->c==1) {
         double *m = target.p->m[0];
         double *m2 = f.p->m[0];
         double **m1 = p->m;
         for (int i=0;i<p->r;i++) {
            for (int j=0;j<p->c;j++)
               m[j] += m2[i] * m1[i][j];
         }
      }
      else {
         double **m = target.p->m;
         double **m1 = p->m;
         double **m2 = f.p->m;
         for (int i=0;i<p->c;i++) {
            for (int j=0;j<f.p->c;j++) {
               double sum = 0.0;
               for (int k=0;k<p->r;k++)
                  sum += m1[k][i] * m2[k][j];
               m[i][j] = sum;
            }
         }
      }
      return target;
   }                          
   dMat d;
   return d;
}

dMat operator-(const dMat &f, double fl)
{
   dMat target(f.p->r,f.p->c);
   double **m = target.p->m;
   double **m1 = f.p->m;
   if (f.isVectorOrScalar()) {
      int n = MAX(f.p->r,f.p->c);
      for (int i=0;i<n;i++)
         m[0][i] = m1[0][i] - fl;
   }
   else {
      for (int i=0;i<f.p->r;i++) {
         for (int j=0;j<f.p->c;j++) {
            m[i][j] = m1[i][j] - fl;
         }
      }
   }
   return target;
}

dMat operator-(double fl, const dMat &f)
{
   dMat target(f.p->r,f.p->c);
   double **m = target.p->m;
   double **m1 = f.p->m;
   if (f.isVectorOrScalar()) {
      int n = MAX(f.p->r,f.p->c);
      for (int i=0;i<n;i++)
         m[0][i] = fl - m1[0][i];
   }
   else {
      for (int i=0;i<f.p->r;i++) {
         for (int j=0;j<f.p->c;j++) {
            m[i][j] = fl - m1[i][j];
         }
      }
   }
   return target;
}

dMat operator-(const dMat &f1, const dMat &f2)
{
   int t1, t2;
   if (!f1.isExEConformable(f2))
      throw("Arrays are not binary conformable in operator-()");
   t1=f1.isVectorOrScalar();
   t2=f2.isVectorOrScalar();
   if (t1==1)                                   // f1 is a scalar
      return (f1.p->m[0][0] - f2);
   if (t2==1)                                   // f2 is a scalar
      return (f1 - f2.p->m[0][0]);
   int r = MAX(f1.p->r, f2.p->r);
   int c = MAX(f1.p->c, f2.p->c);
   dMat target(r,c);
   double **m = target.p->m;
   double **m1 = f1.p->m;
   double **m2 = f2.p->m;                   // f1 and f2 have same
   if (f1.isBinaryConformable(f2)) {       // dimensions - simple subtract
      if (c==1) {
         for (int i=0;i<r;i++) {
               m[0][i] = m1[0][i] - m2[0][i];
         }
      }
      else {
         for (int i=0;i<r;i++) {
            for (int j=0;j<c;j++) {
               m[i][j] = m1[i][j] - m2[i][j];
            }
         }
      }
      return target;
   }
   if (t1==2) {                      // f1 or f2 = row vector
      if (t2==3) {
         for (int i=0;i<r;i++) {
            for(int j=0;j<c;j++) {
               m[i][j] = m1[0][j] - m2[i][0];
            }
         }
         return target;
      }
      else {
         for (int i=0;i<r;i++) {
            for(int j=0;j<c;j++) {
               m[i][j] = m1[0][j] - m2[i][j];
            }
         }
         return target;
      }
   }
   if (t2==2) {                      // f2 = row vector
      if (t1==3) {
         for (int i=0;i<r;i++) {
            for(int j=0;j<c;j++) {
               m[i][j] = m1[i][0] - m2[0][j];
            }
         }
         return target;
      }
      else {
         for (int i=0;i<r;i++) {
            for(int j=0;j<c;j++) {
               m[i][j] = m1[i][j] - m2[0][j];
            }
         }
         return target;
      }
   }
   if (t1==3) {                      // f1 = col vector
      for (int i=0;i<r;i++) {
         for(int j=0;j<c;j++) {
            m[i][j] = m1[i][0] - m2[i][j];
         }
      }
      return target;
   }
   if (t2==3) {                      // f2 = col vector
      for (int i=0;i<r;i++) {
         for(int j=0;j<c;j++) {
            m[i][j] = m1[i][j] - m2[i][0];
         }
      }
      return target;
   }        
   dMat d;
   return d;
}

dMat & dMat::operator-=(const dMat &f)
{
   int t1, t2;
   if (!isExEConformable(f))
      throw("Arrays are not conformable in fMat::operator-=()");
   t1=isVectorOrScalar();
   t2=f.isVectorOrScalar();
	if (t2==1)                                   // f is a scalar
		return (*this -= f.p->m[0][0]);
   if (t2==1)                                   // this is a scalar
      throw("lhs is a scalar, rhs is not, in operator-=()");
   double **m = p->m;
   double **m1 = f.p->m;
   if (isBinaryConformable(f)) {        // this and f have same
      for (int i=0;i<p->r;i++) {           // dimensions - simple subtraction
         for (int j=0;j<p->c;j++) {
            m[i][j] = m[i][j] - m1[i][j];
         }
      }
		return *this;
   }
   if ((t1==2)||(t1==3)) {        // lhs is a row or col vector, rhs is not
      throw("lhs is a row or col vector, rhs is not, in operator-=()");
	}
   if (t2==2) {                   // rhs is a row vec, lhs a mat
      for (int i=0;i<p->r;i++) {
         for(int j=0;j<p->c;j++) {
            m[i][j] = m[i][j] - m1[0][j];
         }
      }
		return *this;
   }
   if (t2==3) {                   // rhs is a col vec, lhs a mat
      for (int i=0;i<p->r;i++) {
         for(int j=0;j<p->c;j++) {
            m[i][j] = m[i][j] - m1[i][0];
         }
      }
		return *this;
   }                            
   return *this;
}

dMat & dMat::operator-=(double fl)
{
   double **m = p->m;
   for (int i=0;i<p->r;i++) {
      for (int j=0;j<p->c;j++) {
         m[i][j] = m[i][j] - fl;
      }
   }
	return *this;
}

cMat::cMat() 
{
   p = new cMatRep;
   if (!p)
      throw("Out of memory in cMat");
   p->m = NULL;
   p->I = NULL;
   p->r = 0;
   p->c = 0;
   p->refs = 1;
   p->missingValue = dMissingValue;
   parent = NULL;
#ifdef MATDEBUG
   nmats++;
   mats++;
   p->matNo = mats;
#endif
}

void dMat::MinMax(double &min, double &max)
{
   double **m = p->m;
   max=min=m[0][0];
   for (int i=0;i<rows(*this);i++) {
      for (int j=0;j<cols(*this);j++) {
         if(isMissing(i, j))
            continue;
			if (max < m[i][j]) max = m[i][j];
         if (min > m[i][j]) min = m[i][j];
      }
   }
}

cMat::cMat(int r)
{
   p = new cMatRep;
   if (!p)
      throw("Out of memory in cMat");
   p->m = new double*[r];
   if (!p->m)
      throw("Out of memory in cMat");
   for (int i=0;i<r;i++)
      p->m[i] = NULL;
   p->I = new Index[r];
   p->r = r;
   p->c = 0;
   p->refs = 1;
   p->missingValue = dMissingValue;
   parent = NULL;
#ifdef MATDEBUG
   nmats++;
   mats++;
   p->matNo = mats;
#endif
}

cMat::~cMat()
{
   int i=0;
	if (--p->refs==0) {
		if (parent==NULL) {
         if (p->m) {
			   for(i=0;i<p->r;i++) {
               if (p->m[i])
//				      delete p->m[i];
				      delete[] p->m[i];
            }
			}
		}
		else if (--parent->refs==0) {
			if (parent->m) {
				for(int i=0;i<parent->r;i++) {
               if (parent->m[i])
//					   delete parent->m[i];
					   delete[] parent->m[i];
				}
//				delete parent->m;
				delete[] parent->m;
			}
			if (parent->I)
				delete [] parent->I;
			delete parent;
#ifdef MATDEBUG
         nmats--;
#endif
      }
      if (p->m)
//         delete p->m;
         delete[] p->m;
      if (p->I)
			delete [] p->I;
      delete p;
#ifdef MATDEBUG
      nmats--;
#endif
   }
   else if (parent)
      --parent->refs;
}


// Operators

cMat cMat::operator=(const cMat &c)
{
   int i = 0;
   if (--p->refs==0) {
		if (parent==NULL) {
         if (p->m) {
            for(i=0;i<p->r;i++) {
               if (p->m)
//                  delete p->m[i];
                  delete[] p->m[i];
            }
         }
      }
      else if (--parent->refs==0) {
         if (parent->m) {
            for(i=0;i<parent->r;i++) {
               if (parent->m[i])
//                  delete parent->m[i];
                  delete[] parent->m[i];
            }
//            delete parent->m;
            delete[] parent->m;
         }
         if (parent->I)
            delete [] parent->I;
         delete parent;
      }
//      delete p->m;
      delete[] p->m;
      if (p->I)
         delete [] p->I;
      delete p;
   }
   else if (parent)
      --parent->refs;
   p = c.p;
   parent = c.parent;
   if (parent)
      c.parent->refs++;
   c.p->refs++;
   return *this;
}

cMat cMat::operator()(const Index &I, int dim)
{
   int n = elementCount(I);
   cMat target(n);
   int *ii = dataptr(I);
   target.p->c = p->c;
   for (int i=0;i<n;i++) {
      target.p->m[i] = p->m[*ii];
      target.p->I[i] = p->I[*ii++];
   }
   target.parent=p;
   target.parent->refs++;
   missingValue(target) = p->missingValue;
   return target;
}

double cMat::operator()(int r, int c)
{
   int *I = dataptr(p->I[r]);
   int n = elementCount(p->I[r]);
   for (int i=0;i<n;i++) {
      if (I[i]==c)
    return p->m[r][i];
   }
   return 0.0;
}

cMat copy(const cMat &c)
{
   int col = c.p->c;
   cMat target(c.p->r);
   cMatRep *tp = target.p, *cp = c.p;
   double **m = tp->m;
   double **m1 = cp->m;
   for (int i=0;i<cp->r;i++) {
      int n = elementCount(cp->I[i]);
      tp->I[i] = copy(cp->I[i]);
      m[i] = new double[n];
      if (!m[i])
         throw("Error: Out of memory in copy(const cMat &)");
      memcpy(m[i],m1[i],n*sizeof(double));
   }
   target.setcols(col);
   return target;
}

bool cMat::deleteRows(char *ii)
{
   if (p->refs > 1)
      throw("Cannot deleteRows - cMat has reference");
   int *n = new int[p->r];
   if (!n)
      return false;
   int i, j=0;
   for(i=0;i<p->r;i++) {
      if (!ii[i]) {
	      n[i] = i-j;
      }
      else
	 j++;
   }
   j = p->r-j;
   double **m = new double*[j];
   if (m==NULL) {
      delete n;
      return 1;
   }
   Index *I = new Index[j];
   if (I==NULL) {
//      delete m;
      delete[] m;
      delete n;
      return 1;
   }
   for(i=0;i<p->r;i++) {
      if(ii[i]) {
//         delete p->m[i];
         delete[] p->m[i];
      }
      else {
         m[n[i]] = p->m[i];
         I[n[i]] = p->I[i];
      }
   }
   delete n;
//   delete [p->r] p->I;
//	delete p->I;
//   delete p->m;
	delete [] p->I;
   delete [] p->m;
   p->I = I;
   p->m = m;
   p->r = j;
   return true;
}


bool cMat::deleteCols(char *ii)
{
   if (p->refs > 1)
      throw("Cannot deleteCols - cMat has reference");
   int *n = new int[p->c];
   if (!n)
      return false;
   int j=0, i, k;
   for(i=0;i<p->c;i++) {
      if (!ii[i]) {
         n[i] = i-j;
      }
      else
         j++;
   }
   j = p->c-j;
   int *newsp = new int[p->c];
   if (!newsp)
      return 1;
   for(i=0;i<p->r;i++) {
      int nsp = elementCount(p->I[i]);
      int *sp = dataptr(p->I[i]);
      int count = 0;
      for (k=0;k<nsp;k++) {
         if (!ii[sp[k]])
            count++;
      }
      if (count==nsp) {
         int *jj = dataptr(p->I[i]);
         for (k=0;k<nsp;k++)
            jj[k] = n[jj[k]];
         continue;
      }
      double *m = new double[count];
      if (!m)
         return 1;
      count = 0;
      for (k=0;k<nsp;k++) {
         if (!ii[sp[k]]) {
            m[count] = p->m[i][k];
            newsp[count] = n[sp[k]];
            count++;
         }
      }
//      delete p->m[i];
      delete [] p->m[i];
      p->m[i] = m;
      Index III(count,newsp);
      p->I[i] = III;
   }
   p->c = j;
   delete newsp;
   delete n;
   return true;
}

cMat & cMat::operator/=(double fl)
{
   double **m = p->m;
   for (int i=0;i<p->r;i++) {
      for (int j=0;j<elementCount(p->I[i]);j++) {
         m[i][j] = m[i][j] / fl;
      }
   }
   return *this;
}

cMat & cMat::operator/=(const dMat &f)
{
   int t1, t2;
   if (!isExEConformable(f))
      throw("Arrays are not conformable in cMat::operator/=()");
   t1=isVectorOrScalar();
   t2=f.isVectorOrScalar();
   if (t2==1)                                   // f is a scalar
      return (*this /= f.p->m[0][0]);
   if (t2==1)                                   // *this is a scalar
      throw("lhs is a scalar, rhs is not, in operator/=()");
   double **m = p->m;
   double **m1 = f.p->m;
   if (isBinaryConformable(f)) {        // *this and f have same
      for (int i=0;i<p->r;i++) {        // dimensions - simple division
         int *ii = p->I[i].p->I;
         for (int j=0;j<elementCount(p->I[i]);j++) {
            m[i][j] = m[i][j] / m1[i][*ii++];
         }
      }
      return *this;
   }
	if ((t1==2)||(t1==3)) {        // lhs is a row or col vector, rhs is not
      throw("\nlhs is a row or col vector, rhs is not, in operator*=()");
   }
   if (t2==2) {                   // rhs is a row vec, lhs a mat
      for (int i=0;i<p->r;i++) {
         int *ii = p->I[i].p->I;
         for (int j=0;j<elementCount(p->I[i]);j++) {
            m[i][j] = m[i][j] / m1[0][*ii++];
         }
      }
      return *this;
   }
   if (t2==3) {                   // rhs is a col vec, lhs a mat
      for (int i=0;i<p->r;i++) {
         for (int j=0;j<elementCount(p->I[i]);j++) {
            m[i][j] = m[i][j] / m1[i][0];
         }
		}
      return *this;
   }
   return *this;
}

cMat & cMat::operator*=(double fl)
{
   double **m = p->m;
   for (int i=0;i<p->r;i++) {
      for (int j=0;j<elementCount(p->I[i]);j++) {
         m[i][j] = m[i][j] * fl;
      }
   }
   return *this;
}

cMat & cMat::operator*=(const dMat &f)
{
   int t1, t2;
   if (!isExEConformable(f))
      throw("Arrays are not conformable in cMat::operator*=()");
   t1=isVectorOrScalar();
   t2=f.isVectorOrScalar();
   if (t2==1)                                   // f is a scalar
      return (*this *= f.p->m[0][0]);
   if (t2==1)                                   // *this is a scalar
      throw("lhs is a scalar, rhs is not, in operator*=()");
   double **m = p->m;
   double **m1 = f.p->m;
   if (isBinaryConformable(f)) {        // *this and f have same
      for (int i=0;i<p->r;i++) {        // dimensions - simple multiply
         int *ii = p->I[i].p->I;
         for (int j=0;j<elementCount(p->I[i]);j++) {
            m[i][j] = m[i][j] * m1[i][*ii++];
         }
      }
      return *this;
   }
   if ((t1==2)||(t1==3)) {        // lhs is a row or col vector, rhs is not
      throw("lhs is a row or col vector, rhs is not, in operator*=()");
   }
   if (t2==2) {                   // rhs is a row vec, lhs a mat
      for (int i=0;i<p->r;i++) {
         int *ii = p->I[i].p->I;
         for (int j=0;j<elementCount(p->I[i]);j++) {
            m[i][j] = m[i][j] * m1[0][*ii++];
         }
      }
      return *this;
   }
   if (t2==3) {                   // rhs is a col vec, lhs a mat
      for (int i=0;i<p->r;i++) {
         for (int j=0;j<elementCount(p->I[i]);j++) {
            m[i][j] = m[i][j] * m1[i][0];
         }
      }
      return *this;
   }                                     
   return *this;
}


dMat cMat::product(const dMat &f)
{
   if (isVectorOrScalar()==1)
      return (f * p->m[0][0]);
   if (f.isVectorOrScalar()==1)
      return (cMat2dMat(*this,0.0) * f.p->m[0][0]);
   if (!innerConformingDimensions(f))
      throw("Dimensions are not conformable in function cMat::product");
   dMat target(p->r,f.p->c);
   Index *I = p->I;
   if (f.p->c==1) {
      double *m2 = f.p->m[0];
      double *m = target.p->m[0];
      for (int i=0;i<p->r;i++) {
         int *ii = dataptr(I[i]);
         double *m1 = p->m[i];
         double sum = 0.0;
         for (int k=0;k<elementCount(I[i]);k++)
            sum += *m1++ * m2[*ii++];
         m[i] = sum;
      }
   }
   else {
      double **m2 = f.p->m;
      for (int i=0;i<p->r;i++) {
         double *m = target.p->m[i];
         for (int j=0;j<f.p->c;j++) {
            double sum = 0.0;
            int *ii = dataptr(I[i]);
            double *m1 = p->m[i];
            for (int k=0;k<elementCount(I[i]);k++)
               sum += *m1++ * m2[*ii++][j];
            m[j] = sum;
         }
      }
   }
   return target;
}

dMat cMat::tproduct(const dMat &f)
{
   double *dd;
   int *ii, i;
   if (isVectorOrScalar()==1)
      return (f * p->m[0][0]);
   if (f.isVectorOrScalar()==1)
      return (cMat2dMat(*this) * f.p->m[0][0]);
   if (p->r != f.p->r)
      throw("Dimensions are not conformable in function cMat::tproduct");
   dMat target(p->c,f.p->c);
   if (f.p->c==1) {
      dd = new double[p->c];
      if (!dd)
         throw("Out of memory in tproduct(cMat)");
      for (int k=0;k<p->c;k++)
         dd[k] = 0.0;
      double *m = target.p->m[0];
      for (i=0;i<p->r;i++) {
         ii = dataptr(p->I[i]);
         double *m1 = p->m[i];
         double m2 = f.p->m[0][i];
         for (int j=0;j<elementCount(p->I[i]);j++) {
            dd[*ii++] += m1[j] * m2;
         }
      }
      for (i=0;i<p->c;i++) {
         m[i] = dd[i];
      }
//      delete dd;
      delete[] dd;
   }
   else {
      double **m = target.p->m;
      dd = new double[p->c];
      if (!dd)
         throw("Out of memory in tproduct(cMat)");
      for (int k=0;k<f.p->c;k++) {
         for (i=0;i<p->c;i++)
            dd[i] = 0.0;
         for (i=0;i<p->r;i++) {
            ii = dataptr(p->I[i]);
            double *m2 = f.p->m[i];
            double *m1 = p->m[i];
            for (int j=0;j<elementCount(p->I[i]);j++) {
               dd[*ii++] += m1[j] * m2[k];
            }
         }
         for (i=0;i<p->c;i++) {
            m[i][k] = dd[i];
         }
      }
//      delete dd;
      delete[] dd;
   }
   return target;
}

cMat dMat2cMat(const dMat &f, double missing_value)
{
   int r = rows(f);
   int c = cols(f);
   cMat target;
   double **mm = new double*[r];
   Index *II = new Index[r];
   if ((!mm)||(!II))
      throw("Out of memory in dMat2cMat");
   target.setdataptr(mm);
   target.setIndexptr(II);
   target.setrows(r);
   target.setcols(c);
   double *dd = new double[c];
   int *ii = new int[c];
   if ((!dd)||(!ii))
      throw("Error: Out of memory in fMat2cMat");
   for (int i=0;i<r;i++) {
      int n=0;
      double *d;
      double *ddd = rowptr(f,i);
      for (int j=0;j<c;j++) {
         if (fabs(ddd[j]-missing_value) < f.dTolerance)
            continue;
         dd[n] = ddd[j];
         ii[n++] = j;
      }
      Index I(n,ii);
      Indexptr(target,i) = I;
      d = new double[n];
      if (!d)
         throw("Error: Out of memory in dMat2cMat");
      target.setrowptr(d,i);
      memcpy(d,dd,(size_t)(sizeof(double)*n));
   }
//   delete dd;
   delete [] dd;
//   delete ii;
   delete [] ii;
   return target;
}

dMat cMat2dMat(const cMat &c, double missing_value)
{
   dMat target(rows(c),cols(c),missing_value);
   for (int i=0;i<rows(c);i++) {
      int *I = dataptr(Indexptr(c,i));
      double *tm = rowptr(target,i);
      double *m = rowptr(c,i);
      for (int j=0;j<elementCount(Indexptr(c,i));j++)
      tm[I[j]] = m[j];
   }
   return target;
}

// Math stuff

dMat count(const dMat &f, enumDirection dir)
{
   int col = cols(f);
   int row = rows(f);
   double **mm = dataptr(f);
   if (dir==ColWise) {                   // column sums
		dMat target(1,col);
      double *m1 = dataptr(target)[0];
		for (int i=0;i<row;i++) {
         for (int j=0;j<col;j++) {
            if (fabs(mm[i][j]) > 0.0)
               m1[j]++;
         }
      }
      return target;
   }
   else if (dir==RowWise) {                      // row sums
      dMat target(row,1);
      double *ss = dataptr(target)[0];
		for (int i=0;i<row;i++) {
         for (int j=0;j<col;j++) {
            if (fabs(mm[i][j]) > 0.0)
               ss[i]++;
         }
      }
      return target;
   }
   else
      throw("Direction out of range in dMat::count(dir)");
   return NULL;
}


dMat sum(const dMat &f, enumDirection dir)
{
   if (dir==ColWise) {                   // column sums
      dMat target(1,f.p->c);
      double **m = f.p->m;
      for (int i=0;i<f.p->c;i++) {
         double s=0.0;
         for (int j=0;j<f.p->r;j++)
            s += m[j][i];
         target(0,i) = s;
      }
      return target;
   }
   else if (dir==RowWise) {                      // row sums
      dMat target(f.p->r,1);
      double **m = f.p->m;
      for (int i=0;i<f.p->r;i++) {
         double s=0.0;
         for (int j=0;j<f.p->c;j++)
            s += m[i][j];
         target(i, 0) = s;
      }
      return target;
   }
   else
      throw("Integer out of range in sum (must be 0 or 1)");
   return NULL;
}


double sum(const dMat &f)
{
   double s=0.0;
   double **m = f.p->m;
   for (int i=0;i<f.p->r;i++) {
      for (int j=0;j<f.p->c;j++)
         s += m[i][j];
   }
   return s;
}

dMat sumsq(const dMat &f, enumDirection dir)
{
   double d;
   if (dir==ColWise) {                   // column sumsq
      dMat target(1,f.p->c);
      double **m = f.p->m;
      for (int i=0;i<f.p->c;i++) {
         double s=0.0;
         for (int j=0;j<f.p->r;j++) {
            d = m[j][i];
            s += d * d;
         }
         target(0,i) = s;
      }
      return target;
   }
   else if (dir==RowWise) {                      // row sumsq
      dMat target(f.p->r,1);
      double **m = f.p->m;
      for (int i=0;i<f.p->r;i++) {
         double s=0.0;
         for (int j=0;j<f.p->c;j++) {
            d = m[i][j];
            s += d * d;
         }
         target(i, 0) = s;
      }
      return target;
   }
   else
      throw("Integer out of range in sumsq (must be 0 or 1)");
   return NULL;
}

double sumsq(const dMat &f)
{
   double s = 0.0;
   if (f.isVectorOrScalar()) {
      double *m = f.p->m[0];
      int n = MAX(f.p->r,f.p->c);
      for (int i=0;i<n;i++)
         s += m[i]*m[i];
   }
   else {
      double **m = f.p->m;
      for (int i=0;i<f.p->r;i++) {
         double d;
         for (int j=0;j<f.p->c;j++) {
            d = m[i][j];
            s += d * d;
         }
      }
   }
   return s;
}

dMat sqrt(const dMat &f)
{
   dMat target(rows(f),cols(f));
   double **m = dataptr(target);
   double **m1 = dataptr(f);
   for (int i=0;i<rows(f);i++) {
      for (int j=0;j<cols(f);j++) {
			m[i][j] = sqrt(m1[i][j]);
      }
   }
   return target;
}

dMat exp(const dMat &f)
{
   dMat target(rows(f), cols(f));
   for (int i=0;i<rows(f);i++) {
      for (int j=0;j<cols(f);j++) {
         if ((f)(i,j) > 70.0)
            target(i,j) = exp(double(70));
         else
            target(i,j) = exp(f(i,j));
      }
   }
   return target;
}

double mean(const dMat &f)
{
   return (sum(f) / (f.p->r * f.p->c));
}


dMat sd(const dMat &f, enumDirection dir)
{
	if (dir==RowWise) {
		dMat target = sumsq(f-mean(f,RowWise),RowWise);
		target /= (double) (f.p->c - 1);
		return sqrt(target);
	}
	else if (dir==ColWise) {
		dMat target = sumsq(f-mean(f,ColWise),ColWise);
		target /= (double) (f.p->r - 1);
		return sqrt(target);
	}
	else
		throw("Integer out of range in sum (must be 0 or 1)");
   return NULL;
}


dMat mean(const dMat &f, enumDirection dir)
{
   if (dir==RowWise) {
		dMat target = sum(f,RowWise);
      target /= (double) f.p->c;
      return target;
   }
   else if (dir==ColWise) {
		dMat target = sum(f,ColWise);
      target /= (double) f.p->r;
      return target;
   }
   else
      throw("Integer out of range in sum (must be 0 or 1)");
   return NULL;
}

void maxmin(const dMat &f, double &min, double &max)
{
   double **m = f.p->m;
   max=min=m[0][0];
   for (int i=0;i<rows(f);i++) {
      for (int j=0;j<cols(f);j++) {
			if (max < m[i][j]) max = m[i][j];
         if (min > m[i][j]) min = m[i][j];
      }
   }
}

dMat transpose(const dMat &f)
{
   dMat target(f.p->c,f.p->r);
   if (f.isVectorOrScalar()) {
      int n = MAX(f.p->c,f.p->r);
      memcpy(target.p->m[0],f.p->m[0],(size_t) sizeof(double)*n);
   }
   else {
      double **m = target.p->m;
      double **m1 = f.p->m;
      for (int i=0;i<f.p->r;i++) {
         for (int j=0;j<f.p->c;j++) {
				m[j][i] = m1[i][j];
         }
      }
   }
   return target;
}

#undef MAX

