#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP W_plus_y_k(SEXP R_W, SEXP R_y, SEXP R_nrow, SEXP R_ncol, SEXP R_i_k);
SEXP W_plus_y(SEXP R_W, SEXP R_y, SEXP R_nrow, SEXP R_ncol);
SEXP W_divided_by_y(SEXP R_W, SEXP R_y, SEXP R_nrow, SEXP R_ncol);

static const R_CallMethodDef callMethods[] = {
	{"W_plus_y_k", (DL_FUNC) &W_plus_y_k, 5},
	{"W_plus_y", (DL_FUNC) &W_plus_y, 4},
	{"W_divided_by_y", (DL_FUNC) &W_divided_by_y, 4},

	/* Finish R_CallMethodDef. */
	{NULL, NULL, 0}
};
/* End of the callMethods[]. */


void R_init_MixfMRI(DllInfo *info){
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
} /* End of R_init_MixfMRI(). */
