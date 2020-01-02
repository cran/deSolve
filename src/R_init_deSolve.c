#ifndef R_R_H
# include <R.h>
#endif

#ifndef R_EXT_DYNLOAD_H_
# include <R_ext/Rdynload.h>
#endif

#define EXTERN
#include "deSolve.h"
#undef EXTERN

#include <Rinternals.h>
#include <stdlib.h> // for NULL

/* register native routines ------------------------------------------------ */

/*
   ToDo:
   - consider replacing SEXP with REALSXP, INTSXP, STRSXP (character), VEXSXP (lists) etc.
   - unlock
*/

/* .C calls */
extern void unlock_solver();

/* Examples (manually added) */
extern void initccl4(void (* odeparms)(int *, double *));
extern void eventfun(int *n, double *t, double *y);
extern void derivsccl4(int *neq, double *t, double *y, double *ydot, double *out, int *ip);

extern void initparms(void (* daspkparms)(int *, double *));
extern void initforcs(void (* daspkforcs)(int *, double *));
extern void chemres (double *t, double *y, double *ydot, double *cj, double *delta, int *ires, double *out, int *ip);

extern void scocpar(void (* odeparms)(int *, double *));
extern void scocforc(void (* odeforcs)(int *, double *));
extern void scocder (int *neq, double *t, double *y, double *ydot, double *out, int *ip);

extern void iniaqua(void (* odeparms)(int *, double *));
extern void initaqforc(void (* odeforc)(int *, double *));
extern void aquaphy (int *neq, double *t, double *y, double *ydot, double *out, int *ip);
extern void aquaphyforc (int *neq, double *t, double *y, double *ydot, double *out, int *ip);


/* .Call calls */
extern SEXP call_daspk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_DLL(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_euler(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_iteration(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_lsoda(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_radau(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_rk4(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_rkAuto(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_rkFixed(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_rkImplicit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP call_zvode(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getLagDeriv(SEXP, SEXP);
extern SEXP getLagValue(SEXP, SEXP);
extern SEXP getTimestep();


static const R_CMethodDef CEntries[] = {
    {"unlock_solver", (DL_FUNC) &unlock_solver, 0},
    {"initccl4",     (DL_FUNC) &initccl4,    1},
    {"initparms",    (DL_FUNC) &initparms,   1},
    {"initforcs",    (DL_FUNC) &initforcs,   1},
    {"eventfun",     (DL_FUNC) &eventfun,    3},
    {"derivsccl4",   (DL_FUNC) &derivsccl4,  6},
    {"chemres",      (DL_FUNC) &chemres,     8},
    {"scocpar",      (DL_FUNC) &scocpar,     1},
    {"scocforc",     (DL_FUNC) &scocforc,    1},
    {"scocder",      (DL_FUNC) &scocder,     6},
    {"iniaqua",      (DL_FUNC) &iniaqua,     1},
    {"initaqforc",   (DL_FUNC) &initaqforc,  1},
    {"aquaphy",      (DL_FUNC) &aquaphy,     6},
    {"aquaphyforc",  (DL_FUNC) &aquaphy,     6},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"call_daspk",      (DL_FUNC) &call_daspk,      28},
    {"call_DLL",        (DL_FUNC) &call_DLL,        11},
    {"call_euler",      (DL_FUNC) &call_euler,      11},
    {"call_iteration",  (DL_FUNC) &call_iteration,  12},
    {"call_lsoda",      (DL_FUNC) &call_lsoda,      28},
    {"call_radau",      (DL_FUNC) &call_radau,      26},
    {"call_rk4",        (DL_FUNC) &call_rk4,        11},
    {"call_rkAuto",     (DL_FUNC) &call_rkAuto,     21},
    {"call_rkFixed",    (DL_FUNC) &call_rkFixed,    17},
    {"call_rkImplicit", (DL_FUNC) &call_rkImplicit, 17},
    {"call_zvode",      (DL_FUNC) &call_zvode,      21},
    {"getLagDeriv",     (DL_FUNC) &getLagDeriv,      2},
    {"getLagValue",     (DL_FUNC) &getLagValue,      2},
    {"getTimestep",     (DL_FUNC) &getTimestep,      0},
    {NULL, NULL, 0}
};



/* C callable functions ---------------------------------------------------- */
SEXP get_deSolve_gparms(void);

void lagvalue(double T, int* nr, int N, double* ytau);

void lagderiv(double T, int* nr, int N, double* ytau);

double glob_timesteps[] = {0, 0};

/* Initialization ---------------------------------------------------------- */
void R_init_deSolve(DllInfo *dll) {

  // thpe 2017-03-22, register entry points
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);

  // the following two lines protect against accidentially finding entry points
  R_useDynamicSymbols(dll, FALSE); // disable dynamic searching
  //R_forceSymbols(dll, TRUE);       // entry points as R objects, not as strings

  /*
    thpe: register C callable to support compiled dede functions
    The direct way would be:
      R_RegisterCCallable("deSolve", "get_deSolve_gparms", (DL_FUNC) get_deSolve_gparms);
    while the following macro (taken from package Matrix) makes this more compact.
  */
  #define RREGDEF(name)  R_RegisterCCallable("deSolve", #name, (DL_FUNC) name)
  RREGDEF(get_deSolve_gparms);
  RREGDEF(lagvalue);
  RREGDEF(lagderiv);

  /* initialize global variables */
  timesteps = glob_timesteps;
}
