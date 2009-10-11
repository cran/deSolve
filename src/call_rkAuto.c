/*==========================================================================*/
/* Runge-Kutta Solvers, (C) Th. Petzoldt, License: GPL >=2                  */
/* General RK Solver for methods with adaptive step size                    */
/*==========================================================================*/

#include "rk_util.h"

SEXP call_rkAuto(SEXP Xstart, SEXP Times, SEXP Func, SEXP Initfunc,
  SEXP Parms, SEXP Nout, SEXP Rho,
  SEXP Rtol, SEXP Atol, SEXP Tcrit, SEXP Verbose,
  SEXP Hmin, SEXP Hmax, SEXP Hini, SEXP Rpar, SEXP Ipar,
  SEXP Method, SEXP Maxsteps, SEXP Flist) {

  /**  Initialization **/
  init_N_Protect();

  double *tt = NULL, *xs = NULL;

  double *y,  *f,  *Fj, *tmp, *FF, *rr;
  SEXP  R_yout;
  double *y0,  *y1,  *y2,  *dy1,  *dy2, *out, *yout;

  double err=0, errold=0, dtnew=0, t, dt, t_ext, tmax;

  SEXP R_FSAL, Alpha, Beta;
  int fsal = FALSE;       /* assume no FSAL */
  int interpolate = TRUE; /* polynomial interpolation is done by default */

  int i = 0, j=0, j1=0, k, it=0, it_tot=0, it_ext=0, nt = 0, neq=0;
  int accept = 0;
  int one=1, isForcing;

  /*------------------------------------------------------------------------*/
  /* Processing of Arguments                                                */
  /*------------------------------------------------------------------------*/
  int lAtol = LENGTH(Atol);
  double *atol = (double *) R_alloc((int) lAtol, sizeof(double));

  int lRtol = LENGTH(Rtol);
  double *rtol = (double *) R_alloc((int) lRtol, sizeof(double));

  for (j = 0; j < lRtol; j++) rtol[j] = REAL(Rtol)[j];
  for (j = 0; j < lAtol; j++) atol[j] = REAL(Atol)[j];

  double  tcrit = REAL(Tcrit)[0];
  double  hmin  = REAL(Hmin)[0];
  double  hmax  = REAL(Hmax)[0];
  double  hini  = REAL(Hini)[0];
  int  maxsteps = INTEGER(Maxsteps)[0];
  int  nout     = INTEGER(Nout)[0]; /* number of global outputs is func is in a DLL */
  int  verbose  = INTEGER(Verbose)[0];

  int stage     = (int)REAL(getListElement(Method, "stage"))[0];

  SEXP R_A, R_B1, R_B2, R_C, R_D;
  double  *A, *bb1, *bb2=NULL, *cc=NULL, *dd=NULL;

  PROTECT(R_A = getListElement(Method, "A")); incr_N_Protect();
  A = REAL(R_A);

  PROTECT(R_B1 = getListElement(Method, "b1")); incr_N_Protect();
  bb1 = REAL(R_B1);

  PROTECT(R_B2 = getListElement(Method, "b2")); incr_N_Protect();
  if (length(R_B2)) bb2 = REAL(R_B2);

  PROTECT(R_C = getListElement(Method, "c")); incr_N_Protect();
  if (length(R_C)) cc = REAL(R_C);

  PROTECT(R_D = getListElement(Method, "d")); incr_N_Protect();
  if (length(R_D)) dd = REAL(R_D);

  double  qerr  = REAL(getListElement(Method, "Qerr"))[0];
  double  beta  = 0;      //0.4/qerr;

  PROTECT(Beta = getListElement(Method, "beta")); incr_N_Protect();
  if (length(Beta)) beta = REAL(Beta)[0];

  double  alpha = 1/qerr - 0.75 * beta;
  PROTECT(Alpha = getListElement(Method, "alpha")); incr_N_Protect();
  if (length(Alpha)) alpha = REAL(Alpha)[0];

  PROTECT(R_FSAL = getListElement(Method, "FSAL")); incr_N_Protect();
  if (length(R_FSAL)) fsal = INTEGER(R_FSAL)[0];

  PROTECT(Times = AS_NUMERIC(Times)); incr_N_Protect();
  tt = NUMERIC_POINTER(Times);
  nt = length(Times);

  PROTECT(Xstart = AS_NUMERIC(Xstart)); incr_N_Protect();
  xs  = NUMERIC_POINTER(Xstart);
  neq = length(Xstart);

  /*------------------------------------------------------------------------*/
  /* DLL, ipar, rpar (for compatibility with lsoda)                         */
  /*------------------------------------------------------------------------*/
  int isDll = FALSE;
  int ntot  =  0;
  int lrpar= 0, lipar = 0;
  int *ipar = NULL;

  /* code adapted from lsoda to improve compatibility */
  if (inherits(Func, "NativeSymbol")) { /* function is a dll */
    isDll = TRUE;
    if (nout > 0) isOut = TRUE;
    ntot  = neq + nout;           /* length of yout */
    lrpar = nout + LENGTH(Rpar);  /* length of rpar; LENGTH(Rpar) is always >0 */
    lipar = 3    + LENGTH(Ipar);  /* length of ipar */

  } else {                              /* function is not a dll */
    isDll = FALSE;
    isOut = FALSE;
    ntot = neq;
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
    /* out: first nout elements of out are reserved for output variables
       other elements are set via argument *rpar* */
    for (j = 0; j < nout; j++)         out[j] = 0.0;                
    for (j = 0; j < LENGTH(Rpar); j++) out[nout+j] = REAL(Rpar)[j];
  }

  /*------------------------------------------------------------------------*/
  /* Allocation of Workspace                                                */
  /*------------------------------------------------------------------------*/
  y0  =  (double *) R_alloc(neq, sizeof(double));
  y1  =  (double *) R_alloc(neq, sizeof(double));
  y2  =  (double *) R_alloc(neq, sizeof(double));
  dy1 =  (double *) R_alloc(neq, sizeof(double));
  dy2 =  (double *) R_alloc(neq, sizeof(double));
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
  initParms(Initfunc, Parms);
  isForcing = initForcings(Flist);

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
  tmax = fmax(tt[nt], tcrit);
  dt = fmin(hmax, hini);
  hmax = fmin(hmax, tmax - t);

  /* Initialize work arrays (to be on the safe side, remove this later) */
  for (i = 0; i < neq; i++)  {
    y1[i] = 0;
    y2[i] = 0;
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

  do { 
    /******  save former results of last step if the method allows this
            (first same as last)                                       ******/
    if (fsal && accept){
      j1 = 1;
      for (i = 0; i < neq; i++) FF[i] = FF[i + neq * (stage - 1)];
    } else {
      j1 = 0;
    }
    /******  Prepare Coefficients from Butcher table ******/
    for (j = j1; j < stage; j++) {
      for(i = 0; i < neq; i++) Fj[i] = 0;
        k = 0;
        while(k < j) {
          for(i = 0; i < neq; i++)
            Fj[i] = Fj[i] + A[j + stage * k] * FF[i + neq * k] * dt;
          k++;
        }
        for (int i = 0; i < neq; i++) {
          tmp[i] = Fj[i] + y0[i];
        }
        /******  Compute Derivatives ******/
        /* pass option to avoid unnecessary copying in derivs */
        derivs(Func, t + dt * cc[j], tmp, Parms, Rho, FF, out, j, neq, 
               ipar, isDll, isForcing);
    }

    /*====================================================================*/
    /* Estimation of new values                                           */
    /*====================================================================*/

    /* use BLAS wrapper with reduced error checking */
    blas_matprod1(FF, neq, stage, bb1, stage, one, dy1);
    blas_matprod1(FF, neq, stage, bb2, stage, one, dy2);

    it_tot++; /* count total number of time steps */
    for (i = 0; i < neq; i++) {
      y1[i] = y0[i] +  dt * dy1[i];
      y2[i] = y0[i] +  dt * dy2[i];
    }

    /*====================================================================*/
    /*      stepsize adjustment                                           */
    /*====================================================================*/
    
    err = maxerr(y1, y2, atol, rtol, neq);
    dtnew = dt;
    //accept = TRUE;
    if (err == 0) {  /* use max scale if all tolerances are zero */
      dtnew = fmin(dt * 10, hmax);
      errold = fmax(err, 1e-4); // 1e-4 taken from Press et al.
      accept = TRUE;
    } else if (err < 1.0) {
      //dtnew = fmin(hmax, dt * 0.9 * pow(err, -1.0/qerr));
      // increase step size only if last one was accepted
      if (accept) dtnew = fmin(hmax, dt * 0.9 * pow(err, -alpha) * pow(errold, beta));
      errold = fmax(err, 1e-4); // 1e-4 taken from Press et al.
      accept = TRUE;
    } else if (err > 1.0) {
      accept = FALSE;
      //dtnew = dt * fmax(0.9 * pow(err, -1.0/qerr), 0.2);
      dtnew = dt * fmax(0.9 * pow(err, -alpha), 0.2);
    }

    if (dtnew < hmin) {
      accept=TRUE;
      if (verbose) Rprintf("warning, h < Hmin\n"); // remove this later ...
      istate[0] = -2;
      dtnew = hmin;
    }
    /*====================================================================*/
    /*      Interpolation and Data Storage                                */
    /*====================================================================*/
    if (accept) {
      if (interpolate) {
      /*--------------------------------------------------------------------*/
      /* case A) "Dense Output": built-in polynomial interpolation          */
      /* available for certain rk formulae, e.g. for rk45dp7                */
      /*--------------------------------------------------------------------*/
	if (dd) { /* i.e. if dd is not Zero */
        denspar(FF, y0, y1, dt, dd, neq, stage, rr);
        t_ext = tt[it_ext];
        while (t_ext <= t + dt) {
          densout(rr, t, t_ext, dt, tmp, neq);
          /* store outputs */
          if (it_ext < nt) {
            yout[it_ext] = t_ext;
            for (i = 0; i < neq; i++)
              yout[it_ext + nt * (1 + i)] = tmp[i];
          }
          if(it_ext < nt) t_ext = tt[++it_ext]; else break;
        }
        /*--------------------------------------------------------------------*/
        /* case B) "Neville-Aitken-Interpolation" for integrators             */
        /* without dense output                                               */
        /*--------------------------------------------------------------------*/
        } else {
          /* (1) collect number "nknots" of knots in advanve */
          yknots[iknots] = t + dt;   /* time is first column */
          for (i = 0; i < neq; i++) yknots[iknots + nknots * (1 + i)] = y2[i];
          if (iknots < (nknots - 1)) {
            iknots++;
          } else {
	   /* (2) do polynomial interpolation */
           t_ext = tt[it_ext];
           while (t_ext <= t + dt) {
            neville(yknots, &yknots[nknots], t_ext, tmp, nknots, neq);
            /* (3) store outputs */
            if (it_ext < nt) {
              yout[it_ext] = t_ext;
              for (i = 0; i < neq; i++)
                yout[it_ext + nt * (1 + i)] = tmp[i];
            }
            if(it_ext < nt) t_ext = tt[++it_ext]; else break;
           }
           shiftBuffer(yknots, nknots, neq + 1);
          }
        }
      } else {
        /*--------------------------------------------------------------------*/
        /* Case C) no interpolation at all (for testing purposes)             */
        /* Note that this works only if time steps match exactly              */
        /* i.e. not in all cases because of floating point inaccuracy         */
        /* but that's the price for "no interpolation"                        */
        /*--------------------------------------------------------------------*/
        t_ext = tt[it_ext];
        while (t_ext <= t + dt) {
          if (it_ext < nt) {
            yout[it_ext] = t_ext;
            if (it < nt) {
              /* all time steps */
              yout[it_ext] = t_ext;
              /* only matching time steps */
              if (t_ext == t + dt)
                for (i = 0; i < neq; i++) yout[it_ext + nt * (1 + i)] = y2[i];
            }
          }
          if(it_ext < nt) t_ext = tt[++it_ext]; else break;
        }
      }
      /*--------------------------------------------------------------------*/
      /* next time step                                                     */
      /*--------------------------------------------------------------------*/
      t = t + dt;
      it++;
      for (i=0; i < neq; i++) y0[i] = y2[i];
    } /* else rejected time step */
    dt = fmin(dtnew, tmax - t);
    if (it_ext > nt) {
      Rprintf("error in rk_solvers.c - call_rkauto: output buffer overflow\n");
      break;
    }
    if (it_tot > maxsteps) {
      if (verbose) Rprintf("Max. number of steps exceeded\n");
      istate[0] = -1;
      break;
    }
  } while (t < tmax); /* end of rk main loop */

  /*====================================================================*/
  /* call derivs again to get global outputs                            */
  /* j = -1 suppresses unnecessary internal copying                     */
  /*====================================================================*/

  for (int j = 0; j < nt; j++) {
    t = yout[j];
    for (i = 0; i < neq; i++) tmp[i] = yout[j + nt * (1 + i)];
    derivs(Func, t, tmp, Parms, Rho, FF, out, -1, neq, ipar, isDll, isForcing);
    for (i = 0; i < nout; i++) {
      yout[j + nt * (1 + neq + i)] = out[i];
    }
  }
  /* attach essential internal information (codes are compatible to lsoda) */
  /* ToDo: respect function evaluations due to global outputs              */
  setIstate(R_yout, R_istate, istate, it_tot, stage, fsal, qerr);

  /* release R resources */
  if (verbose) 
    Rprintf("\nNumber of time steps it = %d, it_ext = %d, it_tot = %d\n", 
      it, it_ext, it_tot);
  unprotect_all();
  return(R_yout);
}
