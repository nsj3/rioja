#include <math.h>
#include <Rinternals.h>
#include <Rdefines.h>

/*
__declspec(dllexport) void GetSetRand(double *x) {
*/
void GetSetRand(double *x) {

  static int IX;
  static int IY;
  static int IZ;
  double val = 0.0;
  long seed = 0;
  if (*x < 0.0) {
    seed = (long) fabs(*x);
    IX = seed;
    IY = seed;
    IZ = seed;
/*    
    for (i=0;i<10;i++)
       GetRand(&d);
*/       
   }
   else {
     IX=171*(IX % 177)-2*(IX/177);
     IY=172*(IY % 176)-35*(IY/176);
     IZ=170*(IZ % 178)-63*(IZ/178);
     if (IX < 0)
        IX=IX+30269;
     if (IY < 0)
        IY=IY+30307;
     if (IZ < 0)
        IZ=IZ+30323;
     val = ( ((float)IX) / 30269.0 + ((float)IY) / 30307.0 + ((float)IZ) / 30323.0);
     *x = val - ((double) ((long) val));
   }
}


/*
__declspec(dllexport) void SetSeed(long *seed) {

   double d = 0.0;
   int i;

   IX = seed;
   IY = seed;
   IZ = seed;
   
   for (i=0;i<10;i++)
      GetRand(&d);
}

*/

