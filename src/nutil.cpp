#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include "nutil.h"

using namespace std;

void RLS(char *s)
{
	char *s1 = s;
	while (*s1) {
      if (*s1==' ') {
         char *s2 = s;
         char *s3 = s;
         while (*s3) {
            s3++;
            *s2 = *s3;
				s2++;
         }
         *s2 = '\0';
      }
      else return;
   }
   return;
}

void RemoveEOL(char *str)
{
   while (*str) {
		if (*str=='\n') {
         *str='\0';
         return;
      }
      str++;
   }
}

void RemoveLeadingSpaces(char *dest, char *src)
{
   while (((*src==' ')||(*src=='\t'))&&(*src!='\0')) src++;
   strcpy(dest,src);
}

void RemoveTrailingSpaces(char *dest, char *src)
{
   int i;
   char *c;
   strcpy(dest,src);
   i = (int) strlen(src);
   c = &dest[i-1];
   while (((*c==' ')||(*c=='\t'))&&(--i>1)) c--;
   dest[i] = '\0';
}

void RemoveSpaces(char *dest, char *src)
{
   while (*src) {
      if ((*src==' ')||(*src=='\t')) {
         src++;
      }
		else {
         *dest=*src;
         src++;
         dest++;
      }
   }
   *dest='\0';
}


int ValidInt(char *c, int *n)
{
	char *s;
   int digit=0;

   s=c;
   while(*s) {
      if (*s =='\n')
         break;
      if(isdigit(*s))
			digit++;
      else if ((*s!=' ')&&(*s!='\t'))
         return(0);
      s++;
   }
   if (digit==0)
      return(2);
   *n = atoi(c);
   return(1);
}
int ValidLong(char *c, long *n)
{
	char *s;
   int digit=0;

   s=c;
   while(*s) {
      if (*s =='\n')
         break;
      if(isdigit(*s))
			digit++;
      else if ((*s!=' ')&&(*s!='\t'))
         return(0);
      s++;
   }
   if (digit==0)
      return(2);
   *n = atol(c);
   return(1);
}

int ValidFloat(char *c, float *f)
{
   char *s;
   int digit=0;

   s=c;
   while(*s) {
      if (*s =='\n')
         break;
		if((isdigit(*s))||(*s=='.')||(*s=='-')||(*s=='+')||(*s=='E')||(*s=='e'))
         digit++;
      else if ((*s!=' ')&&(*s!='\t'))
         return(0);
      s++;
   }
   if (digit==0)
      return(2);
   *f = (float) atof(c);
   return(1);
}

void SubStr(char *dest, char *src, int start, int len)
{
   int i;
   char *s=dest;

   for (i=start;i<start+len;i++) {
      if (src[i]=='\n') break;
      *s=src[i];
		s++;
   }
   *s='\0';
}

void BubbleSort (int *item, int count)
{
   int i=0,j=0,n;

   for (i=1;i<count;i++) {
      for (j=count-1;j>=i;--j) {
         if (item[j-1] > item[j]) {
				n = item[j-1];
				item[j-1]=item[j];
				item[j]=n;
			}
		}
	}
}

void PadWithZeros(char *s)
{
   while(*s) {
      if (*s==' ')
         *s='0';
      s++;
   }
}
