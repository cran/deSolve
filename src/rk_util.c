/*==========================================================================*/
/* Runge-Kutta Solvers, (C) Th. Petzoldt, License: GPL >=2                  */
/* Definitions and Utilities needed by Runge-Kutta Solvers                  */
/*==========================================================================*/

/* Load headers needed by the R interface */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>     /* AS_NUMERIC ... */

#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Boolean.h>

#ifdef HAVE_LONG_DOUBLE
# define LDOUBLE long double
#else
# define LDOUBLE double
#endif

#include "deSolve.h"

/*============================================================================*/
/*   Interface to functions written in compiled languages                     */
/*============================================================================*/

/* give name to data types */
typedef void deriv_func(int *, double *, double *,double *, double *, int *);

typedef void init_func (void (*)(int *, double *));

/*============================================================================*/
/*   DLL specific functions                                                   */
/*============================================================================*/

void R_test_call(DllInfo *info) {
  /* Register routines, allocate resources. */
  Rprintf("test_call DLL loaded\n");
}

void R_unload_test_call(DllInfo *info) {
  /* Release resources. */
  Rprintf("test_call DLL unloaded\n");
}

/*============================================================================*/
/*   Functions for processing complex R arguments                             */
/*============================================================================*/

/* -------- getvar from environment ------------------------------------------*/
SEXP getvar(SEXP name, SEXP Rho) {
  SEXP ans;
  if(!isString(name) || length(name) != 1)
    error("name is not a single string");
  if(!isEnvironment(Rho))
    error("Rho should be an environment");
  ans = findVar(install(CHAR(STRING_ELT(name, 0))), Rho);
  return(ans);
}

SEXP getInputs(SEXP symbol, SEXP Rho) {
  if(!isEnvironment(Rho)) error("Rho should be an environment");
  return(getvar(symbol, Rho));
}

/* -- getvar from list ------------------------------------------------------ */
SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < length(list); i++)
   if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
     elmt = VECTOR_ELT(list, i);
     break;
   }
  return elmt;
}

/*============================================================================*/
/*   Arithmetic utilities                                                     */
/*============================================================================*/

/*----------------------------------------------------------------------------*/
/* Matrix Multiplikation                                                      */
/* a reduced version without NA checking, this is ensured otherwise           */
/*----------------------------------------------------------------------------*/

void blas_matprod1(double *x, int nrx, int ncx,
		    double *y, int nry, int ncy, double *z)
{
    const char *transa = "N", *transb = "N";
    int i;
    double one = 1.0, zero = 0.0;

    if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0) {
	    F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
			    x, &nrx, y, &nry, &zero, z, &nrx);
    } else /* zero-extent operations should return zeroes */
    	for(i = 0; i < nrx*ncy; i++) z[i] = 0;
}


/* -- Simple Matrix Multiplikation ------------------------------------------ */
void matprod(int m, int n, int o, double* a, double* b, double* c) {
  int i, j, k;
  for (i = 0; i < m; i++) {
    for (j = 0; j < o; j++) {
    	c[i + m * j] = 0;
    	for (k = 0; k < n; k++) {
    	  c[i + m * j] += a[i + m * k] * b[k + n * j];
    	}
    }
  }
}

double maxdiff(double *x, double *y, int n) {
  double d = 0.0;
  for (int i = 0; i < n; i++) d = fmax(d, fabs(x[i] - y[i]));
  return(d);
}

double maxerr(double *y1, double *y2, double* Atol, double* Rtol, int n) {
  double err = 0, serr = 0, scal, delta;
  for (int i = 0; i < n; i++) {
    scal  = Atol[i] +  fmax(fabs(y1[i]), fabs(y2[i])) * Rtol[i]; /* min?? */
    delta = fabs(y2[i] - y1[i]);
    /* there are also other too */
    if (scal > 0) {
      err   = fmax(err, delta / scal);
      serr  = err + pow(delta/scal, 2.0);
    }
  }
  err = sqrt(serr); /* euclidean norm */
  return(err);
}

/*==========================================================================*/
/*   CALL TO THE MODEL FUNCTION                                             */
/*==========================================================================*/
void derivs(SEXP Func, double t, double* y, SEXP Parms, SEXP Rho,
	    double *ydot, double *yout, int j, int neq, int *ipar, int isDll) {
  SEXP Val, R_fcall;
  SEXP R_t;
  SEXP R_y;
  int i = 0;
  int nout = ipar[0];
  double *yy;
  double ytmp[neq];

  // Rprintf("i0, i1, i2, %d  %d  %d\n", ipar[0], ipar[1], ipar[2]);

  if (isDll) {
    /*------------------------------------------------------------------------*/
    /*   Function is a DLL function                                           */
    /*------------------------------------------------------------------------*/
    deriv_func *cderivs;
    cderivs = (deriv_func *) R_ExternalPtrAddr(Func);
    cderivs(&neq, &t, y, ytmp, yout, ipar);
    if (j >= 0)
      for (i = 0; i < neq; i++)  ydot[i + neq * j] = ytmp[i];
  } else {
    /*------------------------------------------------------------------------*/
    /* Function is an R function                                              */
    /*------------------------------------------------------------------------*/
    PROTECT(R_t = ScalarReal(t)); incr_N_Protect();
    PROTECT(R_y = allocVector(REALSXP, neq)); incr_N_Protect();
    yy = REAL(R_y);
    for (i=0; i< neq; i++) yy[i] = y[i];

    PROTECT(R_fcall = lang4(Func, R_t, R_y, Parms)); incr_N_Protect();
    PROTECT(Val = eval(R_fcall, Rho)); incr_N_Protect();
    /* extract the states of list "val" */
    if (j >= 0)
      for (i = 0; i < neq; i++)  ydot[i + neq * j] = REAL(VECTOR_ELT(Val, 0))[i];
    /* extract outputs from second list element */
    if (j < 0)
      for (i = 0; i < nout; i++)  yout[i] = REAL(VECTOR_ELT(Val, 1))[i];
    my_unprotect(4);
  }
}

/*============================================================================*/
/*   Interpolation functions                                                  */
/*============================================================================*/

/*----------------------------------------------------------------------------*/
/* "dense output"                                                             */
/* is a specific polynomial interpolation that uses intermediate rk steps     */
/*----------------------------------------------------------------------------*/
void denspar(double *FF, double *y0, double *y1, double dt, double *d,
  int neq, int stage, double *r) {
  double ydiff, bspl;
  int i, j;
  for (i=0; i< neq; i++) {
   r[i]           = y0[i];
   ydiff          = y1[i] - y0[i];
   r[i + neq]     = ydiff;
   bspl           = dt * FF[i] - ydiff;
   r[i + 2 * neq] = bspl;
   r[i + 3 * neq] = ydiff - dt * FF[i + (stage - 1) * neq] - bspl;
   r[i + 4 * neq] = 0;
   for (j=0; j < stage; j++)
     r[i + 4 * neq] = r[i + 4 * neq] + d[j] * FF[i + j * neq];
     r[i + 4 * neq] = r[i + 4 * neq] * dt;
  }
}

void densout(double *r, double t0, double t, double dt, double* res, int neq) {
  double s  = (t - t0) / dt;
  double s1 = 1.0 - s;
  for (int i = 0; i < neq; i++)
    res[i] = r[i] + s * (r[i +     neq] + s1 * (r[i + 2 * neq]
                  + s * (r[i + 3 * neq] + s1 * (r[i + 4 * neq]))));
}

/*----------------------------------------------------------------------------*/
/* Polynomial interpolation                                                   */
/*    ksig: number of signals                                                 */
/*    n:    number of knots per signal                                        */
/*    x[0 .. n-1]:          vector of x values                                */
/*    y[0 .. n-1, 0 .. ksig] array  of y values                               */
/*                                                                            */
/*    ToDo: check if ringbuffer is faster; rewrite eventually                 */
/*----------------------------------------------------------------------------*/
void neville(double *xx, double *y, double tnew, double *ynew, int n, int ksig) {
  int i, j, k;
  double x[n];
  double yy[n * ksig]; /* temporary workspace */
  double tscal = xx[n-1] - xx[0];
  double t = tnew / tscal;
  for (i = 0; i < n; i++) x[i] = xx[i]/tscal;

  for (i=0; i < n * ksig; i++) yy[i] = y[i];

  for (k = 0; k < ksig; k++) {
    for (j = 1; j < n; j++)
      for (i = n - 1; i >= j; i--) {
        yy[i + k * n] = ((t - x[i - j]) * yy[i + k * n]
          - (t - x[i]) * yy[i - 1 + k * n]) / (x[i] - x[i - j]);
      }
    ynew[k] = yy[n - 1 + k * n];
  }
}

/*============================================================================*/
/*   Specific utility functions                                               */
/*============================================================================*/

void shiftBuffer (double *x, int n, int k) {
  /* n = rows, k=columns */
  for (int i = 0; i < (n - 1); i++)
    for (int j = 0; j < k; j++)
      x[i + j * n] = x[i + 1 + j * n];
}

void initParms(SEXP Initfunc, SEXP Parms) {
  if (inherits(Initfunc, "NativeSymbol"))  {
    PROTECT(de_gparms = Parms); incr_N_Protect();
    init_func *initializer;
    initializer = (init_func *) R_ExternalPtrAddr(Initfunc);
    initializer(Initdeparms);
  }
}

void setIstate(SEXP R_yout, SEXP R_istate, int *istate,
  int it_tot, int stage, int fsal, int qerr) {
  /* note that indices are 1 smaller in C than in R  */
  istate[11] = it_tot;                  /* number of steps */
  istate[12] = it_tot * (stage - fsal); /* number of function evaluations */
  istate[14] = qerr;                    /* order of the method */
  setAttrib(R_yout, install("istate"), R_istate);
}

