/*==========================================================================*/
/* Runge-Kutta Solvers, (C) Th. Petzoldt, License: GPL >=2                  */
/* Definitions and Utilities needed by Runge-Kutta Solvers                  */
/*==========================================================================*/


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>

#include <R_ext/Applic.h>
#include <R_ext/Boolean.h>

#ifdef HAVE_LONG_DOUBLE
# define LDOUBLE long double
#else
# define LDOUBLE double
#endif

#include "deSolve.h"

void R_test_call(DllInfo *info);

void R_unload_test_call(DllInfo *info);

SEXP getvar(SEXP name, SEXP Rho);

SEXP getInputs(SEXP symbol, SEXP Rho);

void blas_matprod1(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z);

void matprod(int m, int n, int o, double* a, double* b, double* c);

double maxdiff(double *x, double *y, int n);

double maxerr(double *y1, double *y2, double* Atol, double* Rtol, int n);

void derivs(SEXP Func, double t, double* y, SEXP Parms, SEXP Rho,
	    double *ydot, double *yout, int j, int neq, int *ipar, 
            int isDll, int isForcing);
	    
void denspar(double *FF, double *y0, double *y1, double dt, double *d,
  int neq, int stage, double *r);

void densout(double *r, double t0, double t, double dt, double* res, int neq);

void neville(double *xx, double *y, double tnew, double *ynew, int n, int ksig);

void shiftBuffer (double *x, int n, int k);

void setIstate(SEXP R_yout, SEXP R_istate, int *istate,
  int it_tot, int stage, int fsal, int qerr);
  
