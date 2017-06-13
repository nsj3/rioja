#ifndef _MAT_HPP
#define _MAT_HPP 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

enum mattype { undefined = 0, sparse = 1, full = 2 };
enum enumDirection { RowWise = 0, ColWise = 1 };

void sort(char *, long);
void sort(int *, long);
void sort(double *, long);

class cMat;
class dMat;
class dataMat;

// #define MATDEBUG 1

struct IndexRep {
   int *I;
   int n;
   int refs;
};

struct dMatRep {
   double **m;
   int r;
   int c;
   int refs;
   double missingValue;
};

class Mat {
public:
//   static void error(char *s,...);
   static char MATVER[10];
   static double dTolerance;
   static bool BoundsCheck;
   static int maxCols;
   static int maxRows;
   static int maxIndex;
   static double dMissingValue;
   static void BoundsCheckOn(void) { BoundsCheck = 1; };
   static void BoundsCheckOff(void) { BoundsCheck = 0; };
   static void BoundsError(void);
   static void SetMaxRows(int nRows) { maxRows = nRows; };
   static void SetMaxCols(int nCols) { maxCols = nCols; };
   static void SetIndex(int nIndex) { maxIndex = nIndex; };
   static void SetMissingValue(double dMV) { dMissingValue = dMV; };
};

class Index : Mat {
protected:
   IndexRep *p;
public:
   Index();
   Index(int n);
   Index(int n, int initvalue);
   Index(int n, int *initvalues);
inline Index(const Index &I) { I.p->refs++; p=I.p; };
   ~Index();
#ifdef MATDEBUG
   static int nindexes;
#endif
   Index operator=(const Index &I);
   inline int & operator()(int n) {   if ((BoundsCheck) && ((n >= p->n)||(n < 0))) BoundsError(); return p->I[n]; }
   void setelementCount(int n) { p->n = n; };
   void setptr(int *I) { p->I = I; };
   Index reference(const Index &I);
   friend class cMat;
   friend class dMat;
   friend class dataMat;
   friend Index copy(const Index &);
//   friend cMat copy(const cMat &c);
   friend int * dataptr(const Index &I) { return I.p->I; };
   friend int elementCount(const Index &I) { return I.p->n; };
};


struct cMatRep : dMatRep {
   Index *I;
};

struct dataMatRep {
   char **spNam;
   char **samNam;
   int *samNum;
   int mType;
   char *title;
   cMat *C;
   dMat *F;
   int refs;
};

class dMat : public Mat {
public:
   dMatRep *p;
   dMatRep *parent;
   inline bool innerConformingDimensions(const dMat &d) const
      { if (p->c==d.p->r) return true; else return false; }
   char isVectorOrScalar(void) const;
   bool isExEConformable(const dMat &d) const;
   inline bool isBinaryConformable(const dMat &d) const
      { if ((p->r==d.p->r)&&(p->c==d.p->c)) return true; else return false; }

public:

   friend class cMat;

// Constructors

   dMat();
   dMat(int r, int c=1, double initval=0.0);
   dMat(int r, int c, double* initvalues);
   inline dMat(const dMat &d) { d.p->refs++; p = d.p; parent = d.parent; if (parent) d.parent->refs++; };
   ~dMat();

// Operators

   dMat operator=(const dMat &d);
   double & operator()(int r, int c) const { if ((BoundsCheck)&&((r >= p->r)||(r < 0)||(c >= p->c)||(c < 0))) BoundsError(); return p->m[r][c]; }
   dMat operator()(const Index &I, int dim) const;
   dMat & operator*=(double d);
   dMat & operator*=(const dMat &d);
   dMat & operator+=(double d);
   dMat & operator+=(const dMat &d);
   dMat & operator-=(double d);
   dMat & operator-=(const dMat &d);
   dMat & operator/=(double d);
   dMat & operator/=(const dMat &d);
   dMat concat(const dMat &d, int dir = 0);
   void merge(const dMat &d, int dir = 0);
   dMat product(const dMat &d);
   dMat tproduct(const dMat &d);
   void fill(double d);
   void rand(void);
   int deleteRows(char *i);
   int deleteCols(char *i);
   int deleteRows(Index &I);
   int deleteCols(Index &I);
   dMat scale(char &errorflag);
   bool isMissing(int r, int c) { return (fabs(p->m[r][c] - p->missingValue) < 1.0E-06); }

// friend functions

   friend double & missingValue(const dMat &d) { return d.p->missingValue; };
   friend dMat operator*(const dMat &d, double m);
   friend dMat operator*(double m, const dMat &d) { return d * m; };
   friend dMat operator*(const dMat &d1, const dMat &d2);
   friend dMat operator/(const dMat &d, double m);
   friend dMat operator/(double m, const dMat &d);
   friend dMat operator/(const dMat &d1, const dMat &d2);
   friend dMat operator-(const dMat &d, double m);
   friend dMat operator-(double m, const dMat &d);
   friend dMat operator-(const dMat &d1, const dMat &d2);
   friend dMat operator+(const dMat &d, double m);
   friend dMat operator+(double m, const dMat &d) {return d + m; };
   friend dMat operator+(const dMat &d1, const dMat &d2);
   friend dMat copy(const dMat &d);
   friend dMat copy(const dMat &d, const Index &Ir);
   friend dMat copy(const dMat &d, const Index &Ir, const Index &Ic);
   friend dMat sum(const dMat &d, enumDirection dir);
   friend dMat sum(const cMat &c, enumDirection dir);
   friend double sum(const dMat &d);
   friend dMat sumsq(const dMat &d, enumDirection dir);
   friend double sumsq(const dMat &d);
   friend dMat sqrt(const dMat &d);
   friend dMat exp(const dMat &f);
   friend dMat mean(const dMat &d, enumDirection dir);
   friend double mean(const dMat &d);
   friend double sd(const dMat &d);
   friend dMat sd(const dMat &d, enumDirection dir);
   friend void maxmin(const dMat &d, double &max, double &min);
   friend dMat count(const dMat &d, enumDirection dir);
   friend dMat transpose(const dMat &d);
   friend int rows(const dMat &m) { return m.p->r; };
   friend int cols(const dMat &m) { return m.p->c; };
   friend double ** dataptr(const dMat &d) { return d.p->m; };
   friend double * rowptr(const dMat &d, int row) { return d.p->m[row]; };

// INDEX friends from matmath

//   friend double IndexMinimum(const dMat &d, Index *I=NULL);
//   friend double IndexMaximum(const dMat &d, Index *I=NULL);
//   friend double IndexSd(const dMat &d, Index *I=NULL);
//   friend double IndexMean(const dMat &d, Index *I=NULL);
//   friend double IndexMedian(const dMat &d, Index *I=NULL);
//   friend double IndexQuartile(const dMat &d, double perc, Index *I=NULL);

   dMat diag();
   dMat inverse(char &errorflag);
   double determinant(char &errorflag);
//   dMat scale(char &errorflag);
   void switch_columns(int col1, int col2);
   dMat lu_decompose(dMat &indx, double &d, char &errorflag);
   dMat lu_dcmp(dMat &indx, double &d, char &errorflag);
   void lu_back_subst(dMat &indx, dMat &b);
   void copy_column(dMat &, int, int);
   void MinMax(double &min, double &max);
};

class cMat : public Mat {

protected:
   cMatRep *p;
   cMatRep * parent;
   inline bool innerConformingDimensions(const dMat &d) const
      { if (p->c==d.p->r) return true; else return false; };
   char isVectorOrScalar(void) const;
   bool isExEConformable(const dMat &d) const;
   inline bool isBinaryConformable(const dMat &d) const
      { if ((p->r==d.p->r)&&(p->c==d.p->c)) return true; else return false; };

public:

   friend class dataMat;
   friend class dMat;

// constructors

   cMat() ;
   cMat(int r);
   inline cMat(const cMat &c) { c.p->refs++; p = c.p; parent = c.parent; if (parent) c.parent->refs++; };
   ~cMat();

// Operators

   cMat operator=(const cMat &c);
   double operator()(int r, int c);
   bool deleteRows(char *i);
   bool deleteCols(char *i);

   cMat & operator*=(double m);
   cMat & operator*=(const dMat &d);
   cMat & operator/=(double m);
   cMat & operator/=(const dMat &d);
   cMat operator()(const Index &I, int dim);
   dMat product(const dMat &d);
   dMat tproduct(const dMat &d);
   dMat product(const dMat &d, const dMat &d1);
   dMat tproduct(const dMat &d, const dMat &d1);

   void setcols(int c) { p->c = c; };
   void setrows(int r) { p->r = r; };
   void setrowptr(double *f, int i) { p->m[i] = f; };
   void setIndexptr(Index *I) { p->I = I; };
   void setdataptr(double **m) { p->m = m; };


// friend functions

   friend double & missingValue(cMat &c) { return c.p->missingValue; };
   friend dMat sum(const cMat &c, int n);
   friend double sum(const cMat &c);
   friend dMat sumsq(const cMat &c, int n);
   friend double sumsq(const cMat &c);
   friend dMat mean(const cMat &c, int n);
   friend double mean(const cMat &c);
   friend dMat sd(const cMat &c, int n);
   friend double sd(const cMat &c);
   friend cMat copy(const cMat &c);
   friend Index * Indexptr(const dataMat &D);
   friend Index * Indexptr(const cMat &c) { return c.p->I; };
   friend int rows(const cMat &m) { return m.p->r; };
   friend int cols(const cMat &m) { return m.p->c; };
   friend double ** dataptr(const cMat &c) { return c.p->m; };
   friend double * rowptr(const cMat &c, int row) { return c.p->m[row]; };
   friend Index & Indexptr(const cMat &c, int row) { return c.p->I[row]; };
   friend int * Indexdataptr(const cMat &c, int row) { return dataptr(c.p->I[row]); };
};

class dataMat: Mat {
public:

   dataMatRep *p;
   dataMat();
   ~dataMat();
   dataMat(const dataMat &S);

   mattype matType(void) { return (mattype) p->mType; };
   void setspName(char **s) { p->spNam = s; };
   void setsamName(char **s) { p->samNam = s; };
   char ** getspName(void) { return p->spNam; };
   char ** getsamName(void) { return p->samNam; };
   void setsamNo(int *l) { p->samNum = l; };
   int * getsamNo(void) { return p->samNum; };
   void setmatType(mattype i) { p->mType = i; };
   void setcMat(cMat *C) { p->C = C; };
   void setdMat(dMat *F) { p->F = F; };
   void setTitle(char *t) { p->title = t; };
   char * spName(int i) { return p->spNam[i]; };
   char * samName(int i) { return p->samNam[i]; };
   int & samNo(int i) { return p->samNum[i]; };
   void kill(void);
   bool deleteRows(char *i);
   bool deleteCols(char *i);
friend char * gettitle(const dataMat &d) {return d.p->title; };
friend cMat * getcMat(const dataMat &d) { return d.p->C; };
friend dMat * getdMat(const dataMat &d) { return d.p->F; };
friend double & missingValue(dataMat &d);
friend Index * Indexptr(const dataMat &D);
friend double ** dataptr(const dataMat &);
friend int rows(const dataMat &);
friend int cols(const dataMat &);

};

cMat dMat2cMat(const dMat &d, double missing_value = 0.0);
dMat cMat2dMat(const cMat &c, double missing_value = 0.0);

Index fsort(dMat &f);

#endif
