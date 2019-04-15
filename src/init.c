#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> /* for NULL */
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void Dissim(double *, double *, int *, int *, int *);
extern void Dissim2(double *, double *, double *, int *, int *, int *, int *);
extern void GetSetRand(double *x);

/* .Call calls */
extern SEXP WAPLS_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP WAPLS_predict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MLRC_regress(SEXP, SEXP, SEXP, SEXP);
extern SEXP MLRC_predict(SEXP, SEXP, SEXP, SEXP);
extern SEXP chclust_c(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"Dissim", (DL_FUNC) &Dissim, 5},
  {"Dissim2",  (DL_FUNC) &Dissim2,  7},
  {"GetSetRand",  (DL_FUNC) &GetSetRand,  1},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"WAPLS_fit",   (DL_FUNC) &WAPLS_fit,   6},
  {"WAPLS_predict",(DL_FUNC) &WAPLS_predict, 7},
  {"MLRC_regress",   (DL_FUNC) &MLRC_regress,   4},
  {"MLRC_predict",(DL_FUNC) &MLRC_predict, 4},
  {"chclust_c",(DL_FUNC) &chclust_c, 2},
  {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
  {NULL, NULL, 0}
};

void R_init_rioja(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}