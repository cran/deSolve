/* Functions to test compiled code implementation of ODE and DAE */

#include <time.h>
#include <string.h>
#include "deSolve.h"
#include "externalptr.h"

SEXP call_DLL(SEXP y, SEXP dY, SEXP time, SEXP func, SEXP initfunc, SEXP parms,
		          SEXP nOut, SEXP Rpar, SEXP Ipar, SEXP Type, SEXP flist)
{
  SEXP   yout;

  double *ytmp, *dy, tin, *delta, cj;
  int    ny, j,  type, ires, isDll, isForcing, nout=0, ntot=0;

  C_deriv_func_type *derivs;
  C_res_func_type *res;

  int nprot = 0;

  ny   = LENGTH(y);
  type = INTEGER(Type)[0];

/* function is a dll ?*/
  if (inherits(func, "NativeSymbol")) {
   isDll = 1;
  } else {
   isDll = 0;
  }

/* initialise output, parameters, forcings ... */
  initOutR(isDll,  &nout, &ntot, ny, nOut, Rpar, Ipar);

  //initParms(initfunc, parms);
  if (initfunc != NA_STRING) {
    if (inherits(initfunc, "NativeSymbol")) {
      init_func_type *initializer;
      PROTECT(de_gparms = parms); nprot++;
      initializer = (init_func_type *) R_ExternalPtrAddrFn_(initfunc);
      initializer(Initdeparms);
    }
  }
  // end inline initParms


  isForcing = initForcings(flist);

  PROTECT(yout = allocVector(REALSXP,ntot)); nprot++;

  tin = REAL(time)[0];

  ytmp = (double *) R_alloc(ny, sizeof(double));
    for (j = 0; j < ny; j++) ytmp[j] = REAL(y)[j];

  dy   = (double *) R_alloc(ny, sizeof(double));
    for (j = 0; j < ny; j++) dy[j] = REAL(dY)[j];

  if(isForcing == 1)  updatedeforc(&tin);

  if (type == 1)   {
    derivs = (C_deriv_func_type *) R_ExternalPtrAddrFn_(func);

    derivs (&ny, &tin, ytmp, dy, out, ipar) ;
    for (j = 0; j < ny; j++)  REAL(yout)[j] = dy[j];

  } else {

    res = (C_res_func_type  *) R_ExternalPtrAddrFn_(func);
    delta = (double *) R_alloc(ny, sizeof(double));
    for (j = 0; j < ny; j++) delta[j] = 0.;

    res    (&tin, ytmp, dy, &cj, delta, &ires, out, ipar) ;
    for (j = 0; j < ny; j++)  REAL(yout)[j] = delta[j];

  }

  if (nout > 0)   {

	   for (j = 0; j < nout; j++)
	       REAL(yout)[j + ny] = out[j];
  }

  UNPROTECT(nprot);
  return(yout);
}
