#include <math.h>
#include <stdlib.h> /* for min */

#ifndef min
   #define min(a,b) ((a)<(b)?(a):(b))
#endif

double Dissimilarity(double *data1, double *data2, int len, int coef);
/*
__declspec(dllexport) void Dissim(double *x, double *res, long *nrow, long *ncol, long *coef) {
*/

void Dissim(double *x, double *res, int *nrow, int *ncol, int *coef) {
   int i, j;
   int nn, nnn;
   for (i=0;i<*ncol;i++) {
      nn = i * (*nrow);
      for (j=i+1;j<*ncol;j++) {
         nnn = j * (*nrow);
         res[(i * *ncol) + j] = res[(j * *ncol) + i] = Dissimilarity(x+nn, x+nnn, *nrow, *coef);
      }
   }
}

/*
__declspec(dllexport) void Dissim2(double *x1, double *x2, double *res, long *nrow, long *ncol1, long *ncol2, long *coef) {
*/
void Dissim2(double *x1, double *x2, double *res, int *nrow, int *ncol1, int *ncol2, int *coef) {
   int i, j;
   int nn1, nn2;
   for (i=0;i<*ncol1;i++) {
      nn1 = i * (*nrow);
      for (j=0;j<*ncol2;j++) {
         nn2 = j * (*nrow);
         res[i + (j * *ncol1) ] = Dissimilarity(x1+nn1, x2+nn2, *nrow, *coef);
      }
   }
}

double Dissimilarity(double *data1, double *data2, int len, int coef)
{
/*
 1    " Euclidean Distance          ",
 2    " Squared Euclidean Distance  ",
 3    " Mean Euclidean distance     ",
 4    " Absolute Distance           ",
 5    " Mean Absolute Distance      ",
 6    " Percent Dissimilarity (B&C) ",
 7    " Canberra Metric             ",
 8    " Chi-squared Distance        ",
 9    " Squared Chi-squared distance",
 10   " Relative Euclidean Distance ",
 11   " Relative Absolute Distance  ",
 12   " Chord Distance              ",
 13   " Squared Chord Distance      ",
 14   " Geodesic Distance           ",
 15   " Dice Coefficient            ",
 16   " Jaccard Coefficient         ",
 17   " Ochiai Coefficient          ",
 18   " Edwards Chord               ",
 19   " Overpeck Sq Chord Dist      ",
 20   " Overpeck Chord Dist         ",
*/
   int i;
   double s1=0.0, s2=0.0, ss1=0.0, ss2=0.0;
   double f =0.0, d=0.0, dd=0.0;
   double *d1, *d2;
   int a, b, c;
	d1=data1;
	d2=data2;
	switch (coef) {
		case 1:  /*  ED */
			for(i=0;i<len;i++) {
				f = (*d1++ - (*d2++));
				d += f*f;
			}
			d = sqrt(d);
         break;
      case 2:   /* SED */
         for(i=0;i<len;i++) {
            f = (*d1++ - (*d2++));
            d += f*f;
         }
         break;
      case 3:   /* MED */
         for(i=0;i<len;i++) {
            f = (*d1++ - (*d2++));
            d += f*f;
         }
         d = sqrt(d)/ (float) len;
         break;
      case 4:   /* AD  */
         for(i=0;i<len;i++) d += fabs(*d1++ - (*d2++));
         break;
      case 5:   /* MAD */
         for(i=0;i<len;i++) d += fabs(*d1++ - (*d2++));
         d = d/ (float) len;
         break;
      case 6:   /* PD or Bray Curtiss */
         for (i=0;i<len;i++) {
            f += min(*d1, *d2);
            s1 += *d1++;
            s2 += *d2++;
         }
         d = 1 - (2*f/(s1+s2));
         break;
      case 7:   /* Canberra metric  */
			for (i=0;i<len;i++) {
				dd = *d1 + *d2;
				if (dd < 1.0E-10) {
					d1++;
					d2++;
					continue;
				}
				f += (fabs(*d1++ - *d2++) / dd);
			}
			d = f / (float) len;
	 break;
		case 8:  /* Chi - squared */
			for (i=0;i<len;i++) {
				dd = *d1 + *d2;
				if (dd < 1.0E-10) { d1++; d2++; continue; }
				f = *d1++ - *d2++;
				f = f*f;
				d += f / dd;
			}
			d = sqrt(d);
	 break;
		case 9:  /* Squared Chi - squared */
			for (i=0;i<len;i++) {
				dd = *d1 + *d2;
				if (dd<1.0E-10) { d1++; d2++; continue; }
				f = *d1++ - *d2++;
				f = f*f;
				d += f / dd;
			}
	 break;
		case 10:  /* RED */
			for (i=0;i<len;i++) {
				s1 += *d1++;
				s2 += *d2++;
			}
			d1=data1;
			d2=data2;
			for (i=0;i<len;i++) {
				f = (((*d1++)/s1) - ((*d2++)/s2));
				d += f*f;
			}
			d = sqrt(d);
			break;
		case 11:   /* RAD */
			for (i=0;i<len;i++) {
				s1 += *d1++;
				s2 += *d2++;
			}
			d1=data1;
			d2=data2;
			for (i=0;i<len;i++) {
				d += fabs(((*d1++)/s1) - ((*d2++)/s2));
			}
			break;
		case 12:   /* chord distance  */
			for (i=0;i<len;i++) {
				ss1 += *d1 * (*d1);
				ss2 += *d2 * (*d2);
				f += (*d1++ * (*d2++));
			}
			d = f / (sqrt(ss1*ss2));
			d = 2.0*(1.0-d);
			if (d<0)
				d = 0.0;
			else
				d = sqrt(d);
			break;
		case 13:   /* squared chord distance  */
			for (i=0;i<len;i++) {
				ss1 += *d1 * (*d1);
				ss2 += *d2 * (*d2);
				f += (*d1++ * (*d2++));
			}
			d = f / (sqrt(ss1*ss2));
			d = 2.0*(1.0-d);
			break;
		case 14:  /* Geodesic distance */
			for (i=0;i<len;i++) {
				ss1 += *d1 * (*d1);
				ss2 += *d2 * (*d2);
				f += (*d1++ * (*d2++));
			}
			f = f / (sqrt(ss1*ss2));
			f = (d>0.01) ? d : 0.01 ;
			d = sqrt(1-(f*f));
			d = d / f;
			d = (d<1000.0) ? d : 1000.0 ;
			d = atan(d);
			break;
		case 15:  /* Dice coefficient */
			d = 0;
			break;
		case 16:  /* Jaccard */
			a=0,b=0,c=0;
			for (i=0;i<len;i++) {
				if (*d1 > 1.0E-08)
					if (*d2 > 1.0E-08)
						a++;
					else
						b++;
				else
					if (*d2 > 1.0E-08)
						c++;
				d1++;
				d2++;
			}
			d = (float) a / ((float) a + (float) b + (float) c);
			break;
		case 17:  /* Ochiai  */
			d = 0;
			break;
		case 18:   /* Edwards and Cavalli-Svorza */
			for(i=0;i<len;i++) {
				f = (*d1++) - *d2++;
				d += f*f;
			}
			break;
		case 19:   /* Overpeck Chord dist */
			for(i=0;i<len;i++) {
				f = sqrt(*d1++) - sqrt(*d2++);
				d += f*f;
			}
			break;
		case 20:   /* Overpeck Chord dist */
			for(i=0;i<len;i++) {
				f = sqrt(*d1++) - sqrt(*d2++);
				d += f*f;
			}
         d = sqrt(d);
			break;
		default :
			d = -99.9;
			break;
	}
	return d;
}
