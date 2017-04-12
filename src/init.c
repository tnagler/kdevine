#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP kdevine_eval_kde1d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kdevine_eval_pkde1d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kdevine_eval_qkde1d(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP kdevine_ikern_gauss(SEXP);
extern SEXP kdevine_kern_gauss(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"kdevine_eval_kde1d",  (DL_FUNC) &kdevine_eval_kde1d,  5},
    {"kdevine_eval_pkde1d", (DL_FUNC) &kdevine_eval_pkde1d, 5},
    {"kdevine_eval_qkde1d", (DL_FUNC) &kdevine_eval_qkde1d, 5},
    {"kdevine_ikern_gauss", (DL_FUNC) &kdevine_ikern_gauss, 1},
    {"kdevine_kern_gauss",  (DL_FUNC) &kdevine_kern_gauss,  1},
    {NULL, NULL, 0}
};

void R_init_kdevine(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
