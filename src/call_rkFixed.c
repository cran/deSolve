/*==========================================================================*/
/* Runge-Kutta Solvers, (C) Th. Petzoldt, License: GPL >=2                  */
/* General RK Solver for methods with fixed step size                       */
/*==========================================================================*/

#include "rk_util.h"

SEXP call_rkFixed(SEXP Xstart, SEXP Times, SEXP Func, SEXP Initfunc,
  SEXP Parms, SEXP eventfunc, SEXP elist, SEXP Nout, SEXP Rho,
  SEXP Tcrit, SEXP Verbose, SEXP Hini, SEXP Rpar, SEXP Ipar,
      SEXP Method, SEXP Maxsteps, SEXP Flist) {

  /**  Initialization **/
  long int old_N_Protect = save_N_Protected();

  double *tt = NULL, *xs = NULL;

  double *y,  *f,  *Fj, *tmp, *FF, *rr;
  SEXP  R_yout;
  double *y0,  *y1, *dy1, *out, *yout;

  double t, dt, tmax;

  int fsal = FALSE;       /* fixed step methods have no FSAL */
  int interpolate = TRUE; /* polynomial interpolation is done by default */

  int i = 0, j=0, it=0, it_tot=0, it_ext=0, nt = 0, neq=0;
  int isForcing, isEvent;

  /**************************************************************************/
  /****** Processing of Arguments                                      ******/
  /**************************************************************************/
  double  tcrit = REAL(Tcrit)[0];
  double  hini  = REAL(Hini)[0];
  int  maxsteps = INTEGER(Maxsteps)[0];
  int  nout     = INTEGER(Nout)[0]; /* number of global outputs if func is in a DLL */
  int  verbose  = INTEGER(Verbose)[0];

  int stage     = (int)REAL(getListElement(Method, "stage"))[0];

  SEXP R_A, R_B1, R_C;
  double  *A, *bb1, *cc=NULL;

  PROTECT(R_A = getListElement(Method, "A")); incr_N_Protect();
  A = REAL(R_A);

  PROTECT(R_B1 = getListElement(Method, "b1")); incr_N_Protect();
  bb1 = REAL(R_B1);

  PROTECT(R_C = getListElement(Method, "c")); incr_N_Protect();
  if (length(R_C)) cc = REAL(R_C);
  
  double  qerr  = REAL(getListElement(Method, "Qerr"))[0];

  PROTECT(Times = AS_NUMERIC(Times)); incr_N_Protect();
  tt = NUMERIC_POINTER(Times);
  nt = length(Times);

  PROTECT(Xstart = AS_NUMERIC(Xstart)); incr_N_Protect();
  xs  = NUMERIC_POINTER(Xstart);
  neq = length(Xstart);

  /*------------------------------------------------------------------------*/
  /* timesteps (for advection computation in ReacTran)                      */
  /*------------------------------------------------------------------------*/
  if (hini > 0)
    for (i = 0; i < 2; i++) timesteps[i] = fmin(hini, tt[1] - tt[0]);
  else
    for (i = 0; i < 2; i++) timesteps[i] = tt[1] - tt[0];
  
  /**************************************************************************/
  /****** DLL, ipar, rpar (to be compatible with lsoda)                ******/
  /**************************************************************************/
  int isDll = FALSE;
  //int ntot  = 0;
  int lrpar= 0, lipar = 0;
  int *ipar = NULL;

  if (inherits(Func, "NativeSymbol")) { /* function is a dll */
    isDll = TRUE;
    if (nout > 0) isOut = TRUE;
    //ntot  = neq + nout;         /* length of yout */
    lrpar = nout + LENGTH(Rpar);  /* length of rpar; LENGTH(Rpar) is always >0 */
    lipar = 3    + LENGTH(Ipar);  /* length of ipar */

  } else {                              /* function is not a dll */
    isDll = FALSE;
    isOut = FALSE;
    //ntot = neq;
    lipar = 3;    /* in lsoda = 1 */
    lrpar = nout; /* in lsoda = 1 */
  }
  out   = (double *) R_alloc(lrpar, sizeof(double)); 
  ipar  = (int *) R_alloc(lipar, sizeof(int));

  ipar[0] = nout;              /* first 3 elements of ipar are special */
  ipar[1] = lrpar;
  ipar[2] = lipar;
  if (isDll == 1) {
    /* other elements of ipar are set in R-function lsodx via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar); j++) ipar[j+3] = INTEGER(Ipar)[j];
    /* out:  first nout elements of out are reserved for output variables
       other elements are set via argument *rpar*  */
    for (j = 0; j < nout; j++)         out[j] = 0.0;                
    for (j = 0; j < LENGTH(Rpar); j++) out[nout+j] = REAL(Rpar)[j];
  }

  /*------------------------------------------------------------------------*/
  /* Allocation of Workspace                                                */
  /*------------------------------------------------------------------------*/
  y0  =  (double *) R_alloc(neq, sizeof(double));
  y1  =  (double *) R_alloc(neq, sizeof(double));
  dy1 =  (double *) R_alloc(neq, sizeof(double));
  f   =  (double *) R_alloc(neq, sizeof(double));
  y   =  (double *) R_alloc(neq, sizeof(double));
  Fj  =  (double *) R_alloc(neq, sizeof(double));
  tmp =  (double *) R_alloc(neq, sizeof(double));
  FF  =  (double *) R_alloc(neq * stage, sizeof(double));
  rr  =  (double *) R_alloc(neq * 5, sizeof(double));

  /* matrix for polynomial interpolation */
  SEXP R_nknots;
  int nknots = 6;  /* 6 = 5th order polynomials by default*/
  int iknots = 0;  /* counter for knots buffer */
  double *yknots;

  PROTECT(R_nknots = getListElement(Method, "nknots")); incr_N_Protect();
  if (length(R_nknots)) nknots = INTEGER(R_nknots)[0] + 1;
  if (nknots < 2) {nknots=1; interpolate = FALSE;}
  yknots = (double *) R_alloc((neq + 1) * (nknots + 1), sizeof(double));

  /* matrix for holding states and global outputs */
  PROTECT(R_yout = allocMatrix(REALSXP, nt, neq + nout + 1)); incr_N_Protect();
  yout = REAL(R_yout);
  /* initialize outputs with NA first */
  for (i = 0; i < nt * (neq + nout + 1); i++) yout[i] = NA_REAL;

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
  PROTECT(Y = allocVector(REALSXP,(neq)));        incr_N_Protect(); 
  
  initParms(Initfunc, Parms);
  isForcing = initForcings(Flist);
  isEvent = initEvents(elist, eventfunc, 0);
  if (isEvent) interpolate = FALSE;
  
  /*------------------------------------------------------------------------*/
  /* Initialization of Integration Loop                                     */
  /*------------------------------------------------------------------------*/
  yout[0]   = tt[0];              /* initial time                 */
  yknots[0] = tt[0];              /* for polynomial interpolation */
  for (i = 0; i < neq; i++) {
    y0[i]        = xs[i];         /* initial values               */
    yout[(i + 1) * nt] = y0[i];   /* output array                 */
    yknots[iknots + nknots * (i + 1)] = xs[i]; /* for polynomials */
  }
  iknots++;

  t = tt[0];                   
  tmax = fmax(tt[nt - 1], tcrit);

  /* Initialization of work arrays (to be on the safe side, remove this later) */
  for (i = 0; i < neq; i++)  {
    y1[i] = 0;
    Fj[i] = 0;
    for (j= 0; j < stage; j++)  {
      FF[i + j * neq] = 0;
    }
  }

  /*------------------------------------------------------------------------*/
  /* Main Loop                                                              */
  /*------------------------------------------------------------------------*/
  it     = 1; /* step counter; zero element is initial state   */
  it_ext = 0; /* counter for external time step (dense output) */
  it_tot = 0; /* total number of time steps                    */

  if (interpolate) {
  /* integrate over the whole time step and interpolate internally */
    rk_fixed(
         fsal, neq, stage, isDll, isForcing, verbose, nknots, interpolate, 
         maxsteps, nt,
         &iknots, &it, &it_ext, &it_tot,
         istate, ipar,
         t, tmax, hini,
         &dt,
         tt, y0, y1, dy1, f, y, Fj, tmp, FF, rr, A,
         out, bb1, cc, yknots,  yout,
         Func, Parms, Rho
    );
  } else {
  /* integrate until next time step and return */
   for (int j = 0; j < nt - 1; j++) {
       t = tt[j];
       tmax = fmin(tt[j + 1], tcrit);
       dt = tmax - t;
       if (isEvent) {
         updateevent(&t, y0, istate);
       }
       rk_fixed(
         fsal, neq, stage, isDll, isForcing, verbose, nknots, interpolate, 
         maxsteps, nt,
         &iknots, &it, &it_ext, &it_tot,
         istate, ipar,
         t, tmax, fmin(hini, fabs(dt)) * sign(dt),      // <----- hini for backward steps (still experimental)
         &dt,
         tt, y0, y1, dy1, f, y, Fj, tmp, FF, rr, A,
         out, bb1, cc, yknots,  yout,
         Func, Parms, Rho
      );
      /* in this mode, internal interpolation is skipped,
         so we can simply store the results at the end of each call */
      yout[j + 1] = tmax;
      for (i = 0; i < neq; i++) yout[j + 1 + nt * (1 + i)] = y1[i];
    }
  }
  
  /*====================================================================*/
  /* call derivs again to get global outputs                            */
  /* j = -1 suppresses unnecessary internal copying                     */
  /*====================================================================*/
  if(nout > 0) {
    for (int j = 0; j < nt; j++) {
      t = yout[j];
      for (i = 0; i < neq; i++) tmp[i] = yout[j + nt * (1 + i)];
      derivs(Func, t, tmp, Parms, Rho, FF, out, -1, neq, ipar, isDll, isForcing);
      for (i = 0; i < nout; i++) {
        yout[j + nt * (1 + neq + i)] = out[i];
      }
    }
  }

  /* attach diagnostic information (codes are compatible to lsoda) */
  setIstate(R_yout, R_istate, istate, it_tot, stage, fsal, qerr, 0);

  /* verbose printing in debugging mode*/
  if (verbose) {
    Rprintf("Number of time steps it = %d, it_ext = %d, it_tot = %d\n", it, it_ext, it_tot);
    Rprintf("Maxsteps %d\n", maxsteps);
  }
  /* release R resources */
  timesteps[0] = 0;
  timesteps[1] = 0;
  restore_N_Protected(old_N_Protect);
  return(R_yout);
}
 
