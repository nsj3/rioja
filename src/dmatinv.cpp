#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mat.h"

using namespace std;

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#define TINY 1.0e-20;

dMat dMat::diag()
{
	if (p->r != p->c)
		throw("Matrix must be square for diag");
   dMat target(p->r, 1);
   for (int i=0;i<p->r;i++)
      target(i, 0) = p->m[i][i];
   return target;
}

double dMat::determinant(char &errorflag)
{
	if (p->r != p->c)
		throw("Matrix must be square for determinant");
	dMat indx(p->c);
	double d;
	dMat decomp = lu_decompose(indx,d,errorflag);
   if (errorflag)
      return 0;
	double determinant = d;
	for (int i=0;i<p->c;i++)
		determinant *= decomp(i,i);
	return determinant;
}

dMat dMat::inverse(char &errorflag)
{
   errorflag = 0;
	if (p->r != p->c)
		throw("Matrix must be square for inverse");
	dMat Y(p->r, p->r, 0.0);
	for (int i=0;i<p->r;i++)
		Y(i,i) = 1.0;
	dMat indx(p->r);
	dMat B(p->r);
	double d;
	dMat decomp = lu_decompose(indx,d,errorflag);
   if (errorflag)
      return Y;
	for (int col=0;col<p->c;col++) {
		B.copy_column(Y,col,0);
		decomp.lu_back_subst(indx,B);
		Y.copy_column(B,0,col);
	}
	return transpose(Y);
}

void dMat::copy_column(dMat &mm, int from, int to)
{
	for (int i=0;i<p->r;i++)
		p->m[i][to] = mm.p->m[i][from];
}

void dMat::switch_columns(int col1, int col2)
{
   int i;
   dMat temp(p->r);
	for (i=0;i<p->r;i++)
		temp.p->m[i][0] = p->m[i][col1];
	for (i=0;i<p->r;i++)
		p->m[i][col1] = p->m[i][col2];
	for (i=0;i<p->r;i++)
		p->m[i][col2] = temp.p->m[i][0];
}

dMat dMat::scale(char &errorflag)
{
	double temp;
	if (p->r != p->c)
		throw("Matrix must be square for scale");
	dMat scale_vector(p->r);
	for (int i=0;i<p->c;i++) {
		double maxi = 0.0;
		for (int j=0;j<p->r;j++) {
			if ((temp = fabs(p->m[j][i])) > maxi)
				maxi = temp;
		}
		if (maxi == 0.0) {
         errorflag = 1;
			throw("Singular matrix in scale");
         return scale_vector;
      }
		scale_vector(i, 0) = 1.0 / maxi;
	}
	return scale_vector;
}


dMat dMat::lu_dcmp(dMat &indx, double &d, char &errorflag)
{
	if (p->r != p->c)
		throw("Matrix must be square for lu_decompose");
	d = 1.0;
	int i, j, k, imax=0, n;
   n = p->r;
	double dum, sum, big, temp;
	dMat lu_decomp = copy(*this);
	dMat vv(n, 1);
   for (i=0;i<n;i++) {
      big = 0.0;
      for (j=0;j<n;j++)
         if ((temp=fabs(lu_decomp.p->m[i][j])) > big) big = temp;
      if (big == 0.0) {
         errorflag = 1;
			throw("Singular matrix in scale");
         return lu_decomp;
      }
      vv(i, 0) = 1.0 / big;
   }
	for (j=0;j<n;j++) {
   	for (i=0;i<j;i++) {
			sum = lu_decomp.p->m[i][j];
			for (k=0;k<i;k++)
				sum -= lu_decomp.p->m[i][k] * lu_decomp.p->m[k][j];
			lu_decomp.p->m[i][j] = sum;
		}
		big = 0.0;
		for (i=j;i<n;i++) {
			sum = lu_decomp.p->m[i][j];
			for (k=0;k<j;k++)
				sum -= lu_decomp.p->m[i][k] * lu_decomp.p->m[k][j];
			lu_decomp.p->m[i][j] = sum;
         if ( (dum = vv.p->m[i][0] * fabs(sum) ) >= big) {
				imax = i;
				big = dum;
			}
		}
		if (j != imax) {
         for (k=0;k<n;k++) {
            dum = lu_decomp.p->m[imax][k];
            lu_decomp.p->m[imax][k] = lu_decomp.p->m[j][k];
            lu_decomp.p->m[j][k] = dum;
         }
			d = -d;
			dum = vv.p->m[imax][0];
			vv.p->m[imax][0] = vv.p->m[j][0];
			vv.p->m[j][0] = dum;
		}
		indx.p->m[j][0] = imax;

		if(lu_decomp.p->m[j][j] == 0.0)
			lu_decomp.p->m[j][j] = TINY;
      if (j != n-1) {
   		dum = 1.0 / lu_decomp.p->m[j][j];
   		for (i=j+1;i<n;i++)
	   		lu_decomp.p->m[i][j] *= dum;
      }
	}
	return lu_decomp;
}


dMat dMat::lu_decompose(dMat &indx, double &d, char &errorflag)
{
	if (p->r != p->c)
		throw("Matrix must be square for lu_decompose");
	d = 1;
	int row, col, k, col_max=0, n;
  n = p->c;
	double dum, sum, maxi;
	dMat lu_decomp = copy(*this);
	dMat scale_vector = scale(errorflag);
   if (errorflag)
      return lu_decomp;
	for (row=0;row<p->r;row++) {
		for (col=0;col<row;col++) {
			sum = lu_decomp.p->m[row][col];
			if (col>0) {
				for (k=0;k<=col-1;k++)
					sum -= lu_decomp.p->m[row][k] * lu_decomp.p->m[k][col];
				lu_decomp.p->m[row][col] = sum;
			}
		}
		maxi = 0.0;

    for (col=row;col<n;col++) {
			sum = lu_decomp.p->m[row][col];
			if (row>0) {
				for (k=0;k<row;k++)
					sum -= lu_decomp.p->m[k][col] * lu_decomp.p->m[row][k];
				lu_decomp.p->m[row][col] = sum;
			}
			dum = scale_vector.p->m[col][0] * fabs(sum);
			if (dum >= maxi) {
				col_max = col;
				maxi = dum;
			}
		}
		if (row != col_max) {
			lu_decomp.switch_columns(col_max, row);
			d *= -1;
			dum = scale_vector.p->m[col_max][0];
			scale_vector.p->m[col_max][0] = scale_vector.p->m[row][0];
			scale_vector.p->m[row][0] = dum;
		}
		indx.p->m[row][0] = col_max;
      if(lu_decomp.p->m[row][row] == 0.0)
         throw ("Matrix singular in lu_decompose");
//				lu_decomp.p->m[row][row] = TINY;
		if (row != n-1) {
			dum = 1.0 / lu_decomp.p->m[row][row];
			for (col=row+1;col<n;col++)
				lu_decomp.p->m[row][col] *= dum;
		}
	}
//	if (lu_decomp.p->m[p->r-1][p->c-1] == 0.0)
//      throw ("Matrix singular in lu_decompose");
//		lu_decomp.p->m[p->r-1][p->c-1] = TINY;
	return lu_decomp;
}

void dMat::lu_back_subst(dMat &indx, dMat &b)
{
	if (p->r != p->c)
		throw("Matrix must be square for lu_back_subst");
	int row, col, ll, ii=0;
	double sum;
	for (col=0;col<p->c;col++) {
		ll = (int) indx.p->m[col][0];
		sum = b.p->m[ll][0];
		b.p->m[ll][0] = b.p->m[col][0];
		if (ii >= 0) {
			for (row=ii;row<=col-1;row++)
				sum -= p->m[row][col] * b.p->m[row][0];
		}
		else if (sum != 0)
		  ii = col;
		b.p->m[col][0] = sum;
	}
	for (col=p->c-1;col>=0;col--) {
		sum = b.p->m[col][0];
		if (col < p->c-1) {
			for (row=col+1;row<=p->r-1;row++)
				sum -= p->m[row][col] * b.p->m[row][0];
		}
		b.p->m[col][0] = sum / p->m[col][col];
	}
}


