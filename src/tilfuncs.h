#ifndef TILFUNCS_HPP
#define TILFUNCS_HPP 1

//#ifndef _WIN32
   #define WORD unsigned short
//#endif

//these are Eric's Tilia 2.0 structs:

typedef struct
{
  WORD vars;
  WORD levs;
  WORD showVarNums;
  WORD dict_loaded;
  WORD data_loaded;
  WORD data_smoothed;
  WORD sums_calculated;
  WORD conc_calculated;
  WORD sample_toggle;
  WORD dates;
  WORD sums;
  WORD maxdates;
  WORD dummy2;
  WORD dummy3;
  WORD dummy4;
  WORD dummy5;
  WORD dummy6;
  WORD dummy7;
} TILIA2FLAGS ;

typedef struct
{
  WORD vars;
  WORD levs;
  WORD standard_dict;
  WORD dict_loaded;
  WORD data_loaded;
  WORD data_smoothed;
  WORD sums_calculated;
  WORD conc_calculated;
  WORD sample_toggle;
  WORD dates;
  WORD sums;
  WORD maxdates;
  WORD dummy2;
  WORD dummy3;
  WORD dummy4;
  WORD dummy5;
  WORD dummy6;
  WORD dummy7;
} TILIAFLAGS;

typedef struct
{
  char null;
#ifndef _MSC_VER
  WORD  n __attribute__ ((packed));
#else
  WORD  n;
#endif
} NULLINT;

typedef struct
{
  NULLINT VarNum;
  unsigned char sum;
  unsigned char name[61];
  unsigned char VarCode[9];
} TILIA2VARS;

typedef struct {
  unsigned char code[3];
#ifndef _MSC_VER
  WORD          cam_code __attribute__ ((packed));
#else
  WORD          cam_code ;
#endif
  unsigned char sum;
  unsigned char name[41];
  unsigned char codename[10];
} TILIAVARS;

typedef struct
{
#ifndef _MSC_VER
  float         num __attribute__ ((packed));
#else
  float         num ;
#endif
  unsigned char name[11];
} TILIA2SAMPLES;

typedef struct
{
#ifndef _MSC_VER
  float         num __attribute__ ((packed));
#else
  float         num ;
#endif
  unsigned char name[11];
} TILIASAMPLES;

typedef struct
{
#ifndef _MSC_VER
  float x __attribute__ ((packed));
#else
  float x ;
#endif 
  union
  {
    unsigned char byte;
    struct
    {
      unsigned char null   : 1;
      unsigned char num    : 1; //set for numeric data, ie data.flag.byte=2;
      unsigned char unused : 6;
    } attr;
  } flag;
} TILIA2DATA;

typedef struct
{
#ifndef _MSC_VER
  float x __attribute__ ((packed));
#else
  float x ;
#endif
  union
  {
	 unsigned char byte;
	 struct
	 {
		unsigned char null   : 1;
		unsigned char unused : 7;
	 } attr;
  } flag;
} TILIADATA;

/*
typedef struct
{
  double depth;
  double age;
  double sd;
} TILIAAGES;
*/
void TiliaWriteHeader(FILE *fout, char *header);
void TiliaWriteFlags(FILE *fout, int n, int m);
void TiliaWriteVar(FILE *fout, const char *s, const char *code = NULL, int num = 0, char sum = 'A');
void TiliaWriteSample(FILE *fout, float i, const char *s);
void TiliaWriteData(FILE *fout, char byte, float x);
bool TiliaReadHeader(FILE *fin, char *header, int size);
bool Tilia1ReadFlags(FILE *fin, int &n, int &m);
bool Tilia2ReadFlags(FILE *fin, int &n, int &m);
bool Tilia1ReadVar(FILE *fout, char *name, char *longcode, char *shortcode, int &spnum, char &sum);
bool Tilia2ReadVar(FILE *fout, char *name, char *code, int &spnum, char &sum);
bool TiliaReadSample(FILE *fin, float &num, char *s);
bool TiliaReadData(FILE *fin, char &byte, float &x);
#endif

/*

//the schematic for writing a Tilia 2.0 .TIL file is:

static char szTILversion[] = "tilia 2.00";
FILE*       fpTIL;
N           n;

fwrite(szTILversion,strlen(szTILversion),1,fpTIL);
memset(&n,0,sizeof(n));
n.vars        = nvars;      //var count
n.levs        = nlevs;      //sample count
n.dict_loaded = 1;
n.data_loaded = 1;
fwrite(&n,sizeof(n),1,fpTIL);
for (i=0; i<n.vars; i++)
{
  VARS v;
  ...
  fwrite(&v,sizeof(v),1,fpTIL);
}
for (i=0; i<n.levs; i++)
{
  SAMPLES levels;
  ...
  fwrite(&levels,sizeof(levels),1,fpTIL);
}
for (i=0; i<n.vars * n.levs; i++)
{
  DATA    data;
  ...
  fwrite(&data,sizeof(data),1,fpTIL);
}

//I think that is about it - have fun
//john

*/
