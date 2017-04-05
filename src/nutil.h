#ifndef NUTILHPP
#define NUTILHPP 1

void RLS(char *s);
int ValidInt(char *c, int *n);
int ValidLong(char *c, long *n);
int ValidFloat(char *c, float *f);
void SubStr(char *dest, char *src, int start, int len);
void PadWithZeros(char *s);
void RemoveEOL(char *str);
void RemoveSpaces(char *dest, char *src);
void RemoveLeadingSpaces(char *dest, char *src);
void RemoveTrailingSpaces(char *dest, char *src);

#endif
