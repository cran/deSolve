#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP get_deSolve_gparms(void);

void
R_init_deSolve(DllInfo *info) {
  R_RegisterCCallable("deSolve","get_deSolve_gparms", get_deSolve_gparms);
}
