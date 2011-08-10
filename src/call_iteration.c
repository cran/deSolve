/*==========================================================================*/
/* Fixed Step time stepping routine - NO Integration                        */
/*==========================================================================*/

#include "rk_util.h"

SEXP call_iteration(SEXP Xstart, SEXP Times, SEXP Nsteps, SEXP Func, SEXP Initfunc,
	        SEXP Parms, SEXP Nout, SEXP Rho, SEXP Verbose, SEXP Rpar, SEXP Ipar,
          SEXP Flist) {

  /* Initialization */
  long int old_N_Protect = save_N_Protected();

  double *tt = NULL, *xs = NULL;
  double *ytmp, *out;

  SEXP R_y0, R_yout, R_t=NULL, R_y=NULL;
  SEXP Val, R_fcall;

  double *y0, *yout, *yy;

  double t, dt;
  int i = 0, j = 0, it = 0, nt = 0, nst = 0, neq = 0;
  int isForcing;
  C_deriv_func_type *cderivs = NULL;

  /*------------------------------------------------------------------------*/
  /* Processing of Arguments                                                */
  /*------------------------------------------------------------------------*/
  int nsteps = INTEGER(Nsteps)[0];

  PROTECT(Times = AS_NUMERIC(Times));         incr_N_Protect();
  tt = NUMERIC_POINTER(Times);
  nt = length(Times);

  PROTECT(Xstart = AS_NUMERIC(Xstart));       incr_N_Protect();
  xs  = NUMERIC_POINTER(Xstart);
  neq = length(Xstart);

  ytmp =  (double *) R_alloc(neq, sizeof(double));

  int nout  = INTEGER(Nout)[0]; /* n of global outputs if func is in a DLL */
  int verbose = INTEGER(Verbose)[0];

  /*------------------------------------------------------------------------*/
  /* timesteps (e.g. for advection computation in ReacTran)                 */
  /*------------------------------------------------------------------------*/
  for (i = 0; i < 2; i++) timesteps[i] = (tt[1] - tt[0])/nsteps;

  /*------------------------------------------------------------------------*/
  /* DLL, ipar, rpar (for compatibility with lsoda)                         */
  /*------------------------------------------------------------------------*/
  int isDll = FALSE;
  int lrpar= 0, lipar = 0;
  int *ipar = NULL;

  if (inherits(Func, "NativeSymbol")) { /* function is a dll */
    isDll = TRUE;
    if (nout > 0) isOut = TRUE;
    lrpar = nout + LENGTH(Rpar);  /* length of rpar; LENGTH(Rpar) is always >0 */
    lipar = 3    + LENGTH(Ipar);  /* length of ipar */
    cderivs = (C_deriv_func_type *) R_ExternalPtrAddr(Func);

  } else {                        /* function is not a dll */
    isDll = FALSE;
    isOut = FALSE;
    lipar = 3;
    lrpar = nout;
    PROTECT(R_y = allocVector(REALSXP, neq)); incr_N_Protect();
  }
  out   = (double *) R_alloc(lrpar, sizeof(double));
  ipar  = (int *) R_alloc(lipar, sizeof(int));


  ipar[0] = nout;              /* first 3 elements of ipar are special */
  ipar[1] = lrpar;
  ipar[2] = lipar;
  if (isDll == 1) {
    /* other elements of ipar are set in R-function via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar); j++) ipar[j+3] = INTEGER(Ipar)[j];
    /* rpar is passed via "out"; first nout elements of out are reserved for
       output variables; other elements are set via argument rpar */
    for (j = 0; j < nout; j++)         out[j] = 0.0;
    for (j = 0; j < LENGTH(Rpar); j++) out[nout+j] = REAL(Rpar)[j];
  }

  /*------------------------------------------------------------------------*/
  /* Allocation of Workspace                                                */
  /*------------------------------------------------------------------------*/
  PROTECT(R_y0 = allocVector(REALSXP, neq));  incr_N_Protect();
  y0 = REAL(R_y0);

  /* matrix for holding the outputs */
  PROTECT(R_yout = allocMatrix(REALSXP, nt, neq + nout + 1)); incr_N_Protect();
  yout = REAL(R_yout);

  /* attribute that stores state information, similar to lsoda */
  SEXP R_istate;
  int *istate;
  PROTECT(R_istate = allocVector(INTSXP, 22)); incr_N_Protect();
  istate = INTEGER(R_istate);
  istate[0] = 0; /* assume succesful return */
  for (i = 0; i < 22; i++) istate[i] = 0;

  /*------------------------------------------------------------------------*/
  /* Initialization of Parameters (for DLL functions)                       */
  /*------------------------------------------------------------------------*/
  initParms(Initfunc, Parms);
  isForcing = initForcings(Flist);

  /*------------------------------------------------------------------------*/
  /* Initialization of Loop                                                 */
  /*------------------------------------------------------------------------*/
  yout[0] = tt[0]; /* initial time */
  for (i = 0; i < neq; i++) {
    y0[i]              = xs[i];
    yout[(i + 1) * nt] = y0[i];
  }

  /*------------------------------------------------------------------------*/
  /* Main Loop                                                              */
  /*------------------------------------------------------------------------*/
  t = tt[0];
  for (it = 0; it < nt ; it++) {
    dt = (tt[it + 1] - t)/nsteps;
    timesteps[0] = timesteps[1];
    timesteps[1] = dt;
    if (verbose)
        Rprintf("Time steps = %d / %d time = %e\n", it + 1, nt, t);

    if (it == (nt-1)) nsteps = 1;  /* to make sure last step is saved */
    
    for (nst = 0; nst < nsteps; nst++) {

      if (nst == 0) {
        yout[it] = t;
        for (i = 0; i < neq; i++)  yout[it + nt * (1 + i)] = y0[i];
      }
      
      if (isDll) {
        if (isForcing) updatedeforc(&t);
        cderivs(&neq, &t, y0, ytmp, out, ipar);
        for (i = 0; i < neq; i++)  y0[i] = ytmp[i];

      } else {
        yy = REAL(R_y);
        PROTECT(R_t = ScalarReal(t));                     incr_N_Protect();

        for (i = 0; i < neq; i++) yy[i] = y0[i];

        PROTECT(R_fcall = lang4(Func, R_t, R_y, Parms)); incr_N_Protect();
        PROTECT(Val = eval(R_fcall, Rho));               incr_N_Protect();

        for (i = 0; i < neq; i++)  y0[i] = REAL(VECTOR_ELT(Val, 0))[i];

       /* extract outputs from second and following list elements */
        if (nst == (nsteps - 1)) {
          int elt = 1, ii = 0, l;
          for (i = 0; i < nout; i++)  {
            l = LENGTH(VECTOR_ELT(Val, elt));
            if (ii == l) {
	            ii = 0; elt++;
	          }
            out[i] = REAL(VECTOR_ELT(Val, elt))[ii];
            ii++;
          }
        }
        my_unprotect(3);
      }  /* isDLL*/
      t = t + dt;

      if (nst == 0)
        for (i = 0; i < nout; i++) yout[it + nt * (1 + neq + i)] = out[i];
    } /* nsteps*/
  } /* end of main loop */

  /* attach essential internal information (codes are compatible to lsoda) */

  setIstate(R_yout, R_istate, istate, it, 1, 0, 1, 0);

  /* reset timesteps pointer to saved state, release R resources */
  timesteps[0] = 0;
  timesteps[1] = 0;
  restore_N_Protected(old_N_Protect);
  return(R_yout);
}
