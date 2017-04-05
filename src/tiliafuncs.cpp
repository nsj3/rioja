#include <stdio.h>
#include <string.h>
#include "tilfuncs.h"
#include "mat.h"
#include <R_ext/Print.h>

/* Must be compiled with 1 Byte Structure Member Alignment under MSVC */

using namespace std;

bool TiliaReadHeader(FILE *fin, char *header, int size)
{
   if (fread(header,size,1,fin) == 0)
      return false;
   return true;
}

bool Tilia1ReadFlags(FILE *fin, int &n, int &m)
{
  TILIAFLAGS flags;
   if (fread(&flags,sizeof(flags),1,fin) != 1)
      return false;
   n = flags.levs;
   m = flags.vars;
   return true;
}

bool Tilia2ReadFlags(FILE *fin, int &n, int &m)
{
   TILIA2FLAGS flags;
  if (fread(&flags,sizeof(flags),1,fin) != 1)
      return false;
   n = flags.levs;
   m = flags.vars;
   return true;
}

void TiliaWriteFlags(FILE *fout, int n, int m)
{
	TILIA2FLAGS flags2;
  memset(&flags2,0,sizeof(flags2));
   flags2.vars            = m;
	flags2.levs            = n;
   flags2.showVarNums     = 0;
   flags2.dict_loaded     = 1;
   flags2.data_loaded     = 1;
	flags2.data_smoothed   = 0;
   flags2.sums_calculated = 0;
   flags2.conc_calculated = 0;
   flags2.sample_toggle   = 0;
   flags2.dates           = 0;
   flags2.sums            = 0;
   flags2.maxdates        = 0;
   flags2.dummy2          = 0;
   flags2.dummy3          = 0;
	flags2.dummy4          = 0;
	flags2.dummy5          = 0;
	flags2.dummy6          = 0;
	flags2.dummy7          = 0;
   fwrite(&flags2,sizeof(flags2),1,fout);
}

bool Tilia1ReadVar(FILE *fin, char *name, char *longcode, char *shortcode, int &spnum, char &sum)
{
   TILIAVARS var;
  if (fread(&var,sizeof(var),1,fin) != 1) {
      return false;
   }
   spnum = (short) var.cam_code;
   strcpy(name, (char *) var.name);
   strcpy(longcode, (char *) var.codename);
   strcpy(shortcode, (char *) var.code);
   sum = var.sum;
   return true;
}

bool Tilia2ReadVar(FILE *fin, char *name, char *longcode, int &spnum, char &sum)
{
   TILIA2VARS var;
  if (fread(&var,sizeof(var),1,fin) != 1)
      return false;
   strcpy(name, (char *) var.name);
   strcpy(longcode, (char *) var.VarCode);
   if (var.VarNum.null == 0)
      spnum = var.VarNum.n;
   else
      spnum = 0;
   sum = var.sum;
   return true;
}

void TiliaWriteVar(FILE *fout, const char *code, const char *name, int num, char sum)
{	
   TILIA2VARS var;
  strncpy((char *) var.VarCode, code, 8);
   var.VarCode[8] = '\0';
   if (name == NULL) {
      strncpy((char *) var.name, code, 8);
      var.name[8] = '\0';
   }
   else {
      strncpy((char *) var.name, name, 60);
      var.name[60] = '\0';
   }
   var.VarNum.null = 0;
   var.VarNum.n = num;
   var.sum = sum;
	fwrite(&var,sizeof(var),1,fout);
}

void TiliaWriteSample(FILE *fout, float num, const char *s)
{
   TILIA2SAMPLES level;
  level.num = (float) num;
	strcpy((char *) level.name, s);
   fwrite(&level,sizeof(level),1,fout);
}

bool TiliaReadSample(FILE *fin, float &num, char *s)
{
   TILIASAMPLES level;
  if (fread(&level,sizeof(level),1,fin) == 0)
      return false;
	num = level.num;
	strcpy(s, (char *) level.name);
   return true;
}

bool TiliaReadData(FILE *fin, char &byte, float &x)
{
   TILIA2DATA data;
   if (fread(&data,sizeof(data),1,fin) == 0)
      return false;
   byte = data.flag.byte;
   x = data.x;
   return true;
}

void TiliaWriteData(FILE *fout, char byte, float x)
{
 	TILIA2DATA data;
   data.flag.byte = byte;
   if (byte == 2)
      data.x = x;
	fwrite(&data,sizeof(data),1,fout);
}

