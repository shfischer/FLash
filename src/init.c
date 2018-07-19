#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

/* .Call calls */
extern SEXP CalcF(SEXP, SEXP, SEXP);
extern SEXP fwd_adolc_FLStock(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"CalcF", (DL_FUNC) &CalcF, 3},
  {"fwd_adolc_FLStock", (DL_FUNC) &fwd_adolc_FLStock, 10},
  {NULL, NULL, 0}
};

void R_init_FLash(DllInfo *info)
{
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_RegisterCCallable("FLash", "CalcF", (DL_FUNC) &CalcF);
  R_RegisterCCallable("FLash", "fwd_adolc_FLStock", (DL_FUNC) &CalcF);
}
