#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void Dissim(void *, void *, void *, void *, void *);
extern void Dissim2(void *, void *, void *, void *, void *, void *, void *);
extern void GetSetRand(void *);

/* .Call calls */
extern SEXP chclust_c(SEXP, SEXP);
extern SEXP MLRC_predict(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP MLRC_regress(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ReadCornellFile(SEXP, SEXP, SEXP);
extern SEXP ReadTiliaFile(SEXP);
extern SEXP WAPLS_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP WAPLS_predict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP WriteCornellFile(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"Dissim",     (DL_FUNC) &Dissim,     5},
    {"Dissim2",    (DL_FUNC) &Dissim2,    7},
    {"GetSetRand", (DL_FUNC) &GetSetRand, 1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"chclust_c",        (DL_FUNC) &chclust_c,        2},
    {"MLRC_predict",     (DL_FUNC) &MLRC_predict,     5},
    {"MLRC_regress",     (DL_FUNC) &MLRC_regress,     5},
    {"ReadCornellFile",  (DL_FUNC) &ReadCornellFile,  3},
    {"ReadTiliaFile",    (DL_FUNC) &ReadTiliaFile,    1},
    {"WAPLS_fit",        (DL_FUNC) &WAPLS_fit,        6},
    {"WAPLS_predict",    (DL_FUNC) &WAPLS_predict,    7},
    {"WriteCornellFile", (DL_FUNC) &WriteCornellFile, 9},
    {NULL, NULL, 0}
};

void R_init_rioja(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
