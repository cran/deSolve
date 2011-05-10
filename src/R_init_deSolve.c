#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "deSolve.h"

SEXP get_deSolve_gparms(void);

double glob_timesteps[] = {0, 0};

void
R_init_deSolve(DllInfo *info) {
  R_RegisterCCallable("deSolve", "get_deSolve_gparms", (DL_FUNC) get_deSolve_gparms);
  /* initialize global variables */
  timesteps = glob_timesteps;
}
