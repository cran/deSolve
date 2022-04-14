/*==========================================================================*/
/* Runge-Kutta Solvers, (C) Th. Petzoldt, License: GPL >=2                  */
/* Definitions and Utilities needed by Runge-Kutta Solvers                  */
/*==========================================================================*/

/* USE_FC_LEN_T to ensure compatibility with Fortran BLAS/LAPACK */
#ifndef USE_FC_LEN_T
# define USE_FC_LEN_T
#endif

/* Load headers needed by the R interface */

#include <R_ext/Rdynload.h>
#include <R_ext/Applic.h>
#include <R_ext/Boolean.h>
#include "deSolve.h"

#ifdef HAVE_LONG_DOUBLE
# define LDOUBLE long double
#else
# define LDOUBLE double
#endif



/* sign of a number */
#define sign(x) (( x > 0 ) - ( x < 0 ))


/*==========================================================================*/
/* general utilies and interpolation                                        */
/*==========================================================================*/

void R_test_call(DllInfo *info);

void R_unload_test_call(DllInfo *info);

SEXP getvar(SEXP name, SEXP Rho);

SEXP getInputs(SEXP symbol, SEXP Rho);

void blas_matprod1(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z);

void matprod(int m, int n, int o, double* a, double* b, double* c);

double maxdiff(double *x, double *y, int n);

double maxerr(double *y0, double *y1, double *y2, double* Atol, double* Rtol, int n);

void derivs(SEXP Func, double t, double* y, SEXP Parms, SEXP Rho,
	    double *ydot, double *yout, int j, int neq, int *ipar, 
            int isDll, int isForcing);
	    
void denspar(double *FF, double *y0, double *y1, double dt, double *d,
  int neq, int stage, double *r);

void densout(double *r, double t0, double t, double dt, double* res, int neq);

void densoutck(double t0, double t, double dt, double * y0,   
  double* FF, double* dy, double* res, int neq);

void neville(double *xx, double *y, double tnew, double *ynew, int n, int ksig);

void shiftBuffer (double *x, int n, int k);

void setIstate(SEXP R_yout, SEXP R_istate, int *istate,
  int it_tot, int stage, int fsal, int qerr, int nrej);
  
/*==========================================================================*/
/* core functions (main loop) for solvers with variable / fixed step size   */
/*==========================================================================*/

void rk_auto(
  /* integers */
  int fsal, int neq, int stage,
  int isDll, int isForcing, int verbose,
  int nknots, int interpolate, int densetype, int maxsteps, int nt,
  /* int pointers */
  int* _iknots, int* _it, int* _it_ext, int* _it_tot, int *_it_rej,
  int* istate,  int* ipar,
  /* double */
  double t, double tmax, double hmin, double hmax, 
  double alpha, double beta,
  /* double pointers */
  double* _dt, double* _errold,
  /* arrays */
  double* tt, double* y0, double* y1, double* y2, double* dy1, double* dy2,
  double* f, double* y, double* Fj, double* tmp,
  double* FF, double* rr, double* A, double* out, 
  double* bb1, double* bb2, double* cc, double* dd, 
  double* atol, double* rtol, double* yknots, double* yout,
  /* SEXPs */
  SEXP Func, SEXP Parms, SEXP Rho
 );


void rk_fixed(
  /* integers */
  int fsal, int neq, int stage,
  int isDll, int isForcing, int verbose,
  int nknots, int interpolate, int maxsteps, int nt,
  /* int pointers */
  int* _iknots, int* _it, int* _it_ext, int* _it_tot,  
  int* istate,  int* ipar,
  /* double */
  double t, double tmax, double hini,
  /* double pointers */
  double* _dt,
  /* arrays */
  double* tt, double* y0, double* y1,double* dy1, 
  double* f, double* y, double* Fj, double* tmp,
  double* FF, double* rr, double* A, double* out, 
  double* bb1, double* cc, 
  double* yknots, double* yout,
  /* SEXPs */
  SEXP Func, SEXP Parms, SEXP Rho
);

 
void rk_implicit(double * alfa, int *index, 
       /* integers */
       int fsal, int neq, int stage,
       int isDll, int isForcing, int verbose,
       int nknots, int interpolate, int maxsteps, int nt,
       /* int pointers */
       int* _iknots, int* _it, int* _it_ext, int* _it_tot, 
       int* istate,  int* ipar,
       /* double */
        double t, double tmax, double hini,
       /* double pointers */
       double* _dt,
       /* arrays */
       double* tt, double* y0, double* y1, double* dy1, 
       double* f, double* y, double* Fj, 
       double* tmp, double* tmp2, double *tmp3,
       double* FF, double* rr, double* A, double* out, 
       double* bb1, double* cc, 
       double* yknots, double* yout,
       /* SEXPs */
       SEXP Func, SEXP Parms, SEXP Rho
); 
