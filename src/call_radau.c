#include <time.h>
#include <string.h>
#include "deSolve.h"
#include "externalptr.h"

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   RADAU: Implicit runge-Kutta of order 5
   due to Hairer and Wanner, with stepsize control and dense output

   The C-wrappers that provide the interface between FORTRAN codes and R-code
   are: C_deriv_func_rad: interface with R-code "func", passes derivatives
        C_deriv_out_rad : interface with R-code "func", passes derivatives +
                                                               output variables

   C_deriv_func_forc_rad provides the interface between the function specified in
   a DLL and the integrator, in case there are forcing functions.

   version 1.9.1: added time lags -> delay differential equations
                  added root function
                  added events
   version 1.10: mass matrix for func in a DLL
   karline soetaert
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* globals for radau */

  int    maxt, it, nout, isDll, ntot;
  double *xdytmp, *ytmp, *tt, *rwork, *root, *oldroot;
  int    *iwork, *jroot;
  int    iroot, nroot, nr_root, islag, isroot, isEvent, endsim;
  double tin, tprevroot;

  typedef void C_root_func_type (int *, double *, double *,int *, double *);
  C_root_func_type      *root_func = NULL;
  C_deriv_func_type     *deriv_func;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 definition of the calls to the FORTRAN subroutines in file radau.f           */

void F77_NAME(radau5)( int *,
         void (*)(int *, double *, double *, double *, double *, int *), // func
		     double *, double *, double *, double *,
		     double *, double *, int *,
 	       void (*)(int *, double *, double *, int *, int *,
                    double *, int *, double *, int *),                   // jac
		     int *, int *, int *,
 	       void (*)(int *, double *, int *, double *, int *),              // mas
		     int *, int *, int *,
         void (*)(int *, double *, double *, double *, double *,
			            int *, int *, double *, int *, int *, double *),   // soloutrad
		     int *, double *, int *, int *, int*, double *, int*, int*);

/* continuous output formula for radau (used in radau.c and lags.c) */
void F77_NAME (contr5) (int *, double *, double *, int *, double *);

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   interface R with FORTRAN functions                                         */

/* wrapper above the derivate function in a dll that first estimates the
values of the forcing functions                                               */

static void C_deriv_func_forc_rad (int *neq, double *t, double *y,
                         double *ydot, double *yout, int *iout)
{
  updatedeforc(t);
  DLL_deriv_func(neq, t, y, ydot, yout, iout);
}

/* Fortran code calls C_deriv_func_rad(N, t, y, ydot, yout, iout)
   R code called as R_deriv_func(time, y) and returns ydot                    */

static void C_deriv_func_rad (int *neq, double *t, double *y,
                          double *ydot, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, Time, ans;

  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(Time = ScalarReal(*t));
  PROTECT(R_fcall = lang3(R_deriv_func,Time,Y));
  PROTECT(ans = eval(R_fcall, R_envir));

  for (i = 0; i < *neq; i++)   ydot[i] = REAL(ans)[i];

  UNPROTECT(3);
}

/* mass matrix function                                                       */

static void C_mas_func_rad (int *neq, double *am, int *lmas,
                             double *yout, int *iout)
{
  int i;
  SEXP NEQ, LM, R_fcall, ans;

  PROTECT(NEQ = NEW_INTEGER(1));
  PROTECT(LM = NEW_INTEGER(1));

  INTEGER(NEQ)[0] = *neq;
  INTEGER(LM) [0] = *lmas;
  PROTECT(R_fcall = lang3(R_mas_func,NEQ,LM));
  PROTECT(ans = eval(R_fcall, R_envir));

  for (i = 0; i <*lmas * *neq; i++)   am[i] = REAL(ans)[i];

  UNPROTECT(4);
}

/* deriv output function - for ordinary output variables                      */

static void C_deriv_out_rad (int *nOut, double *t, double *y,
                       double *ydot, double *yout)
{
  int i;
  SEXP R_fcall, Time, ans;

  for (i = 0; i < n_eq; i++)
      REAL(Y)[i] = y[i];

  PROTECT(Time = ScalarReal(*t));
  PROTECT(R_fcall = lang3(R_deriv_func,Time, Y));
  PROTECT(ans = eval(R_fcall, R_envir));

  for (i = 0; i < *nOut; i++) yout[i] = REAL(ans)[i + n_eq];

  UNPROTECT(3);
}

/* save output in R-variables                                                 */

static void saveOut (double t, double *y) {
  int j;

    REAL(YOUT)[(it)*(ntot+1)] = t;
	  for (j = 0; j < n_eq; j++)
	    REAL(YOUT)[(it)*(ntot + 1) + j + 1] = y[j];

    /* if ordinary output variables: call function again */
    if (nout>0)   {
      if (isDll == 1)   /* output function in DLL */
        deriv_func (&n_eq, &t, y, xdytmp, out, ipar) ;
      else
        C_deriv_out_rad(&nout, &t, y, xdytmp, out);
      for (j = 0; j < nout; j++)
        REAL(YOUT)[(it)*(ntot + 1) + j + n_eq + 1] = out[j];
    }
}

/* save lagged variables                                                      */

static void C_saveLag(int ini, double *t, double *y, double *con, int *lrc,
                      double *rpar, int *ipar) {
   /* estimate dy (xdytmp) */
   if (isDll == 1)
      deriv_func (&n_eq, t, y, xdytmp, rpar, ipar) ;
   else
      C_deriv_func_rad (&n_eq, t, y, xdytmp, rpar, ipar) ;

   if (ini == 1)
    updatehistini(*t, y, xdytmp, rpar, ipar);
   else
    updatehist(*t, y, xdytmp, con, lrc);
}

/* root function                                                              */

static void C_root_radau (int *neq, double *t, double *y, int *ng, double *gout)
{
  int i;
  SEXP R_fcall, Time, ans;

  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(Time = ScalarReal(*t));
  PROTECT(R_fcall = lang3(R_root_func,Time,Y));
  PROTECT(ans = eval(R_fcall, R_envir));

  for (i = 0; i < *ng; i++)   gout[i] = REAL(ans)[i];

  UNPROTECT(3);
}
/* function for brent's root finding algorithm                                */

double f (double t, double *Con, int *Lrc) {
   F77_CALL(contr5) (&n_eq, &t, Con, Lrc, ytmp);    /* ytmp = value of y at t */
   if (isDll == 1)
     root_func (&n_eq, &t, ytmp, &nroot, root);    /* root at t, ytmp */
   else
     C_root_radau (&n_eq, &t, ytmp, &nroot, root);
   return root[iroot] ;
}

/* function called by Fortran to check for output, lags, events, roots        */

static void C_soloutrad(int * nr, double * told, double * t, double * y,
  double * con, int * lrc, int * neq, double * rpar, int * ipar,
  int * irtrn, double * xout)

{
  int i, j;
  int istate, iterm;
  double tr, tmin;
  double tol = 1e-9;			/* Acceptable tolerance		*/
  int maxit = 100;				/* Max # of iterations */
  extern double brent(double, double,	double, double,
                      double (double, double *, int *), double *, int *,
                      double, int);

  if (*told == *t) return;
  timesteps[0] = *told-*t;
  timesteps[1] = *told-*t;

  if (islag == 1) C_saveLag(0, t, y, con, lrc, rpar, ipar);
  *irtrn = 0;

  if (isEvent && ! rootevent) {
    if (*told <= tEvent && tEvent < *t) {
      tin = tEvent;
      F77_CALL(contr5) (&n_eq, &tEvent, con, lrc, y);
      updateevent(&tin, y, &istate);
      *irtrn = -1;
    }
  }
  tmin = *t;
  iroot = -1;
  if (isroot & (fabs(*t - tprevroot) > tol)) {
    if (isDll == 1)
     root_func (&n_eq, t, y, &nroot, root);    /* root at t, ytmp */
    else
     C_root_radau (&n_eq, t, y, &nroot, root);

    for (i = 0; i < nroot; i++)
     if (fabs(root[i]) <  tol) {
       iroot = i;
       jroot[i] = 1;
       *irtrn = -1;
       endsim = 1;
       tprevroot = *t;
     } else if (fabs(oldroot[i]) >= tol && root[i] * oldroot[i] < 0) {
       iroot = i;
       jroot[i] = 1;
       tr = brent(*told, *t, oldroot[i], root[i], f, con, lrc, tol, maxit);
       if (fabs(tprevroot - tr) > tol) {
       F77_CALL(contr5) (&n_eq, &tr, con, lrc, ytmp);
       *irtrn = -1;
        endsim = 1;
        if (tr < tmin) {
          tmin = tr;
          tprevroot = tmin;
          for (j = 0; j < n_eq; j++) y[j] = ytmp[j];
        }
       }
     } else jroot[i] = 0;
    for (i = 0; i < nroot; i++) oldroot[i] = root[i];
  }

  while (*told <= tt[it] && tt[it] < tmin) {
    F77_CALL(contr5) (neq, &tt[it], con, lrc, ytmp);
    saveOut(tt[it], ytmp);
    it++;
    if ( it >= maxt) break;
  }
   if ((*irtrn == -1) && rootevent) {
     *t = tmin;
     tin = *t;
     tEvent = tin;
     if (nr_root < Rootsave) {
       troot[nr_root] = tin;
       for (j = 0; j < nroot; j++)
         if (jroot[j] == 1) nrroot[nr_root] = j+1;
       for (j = 0; j < n_eq; j++)
         valroot[nr_root* n_eq + j] = y[j];
     }
     iterm = 0;      /* check if simulation should be terminated */
     for (j = 0; j < nroot; j++)
       if (jroot[j] == 1 && termroot[j] == 1) iterm = 1;

     if (iterm == 0) {
       nr_root++;
       updateevent(&tin, y, &istate);
       endsim = 0;
     } else {
       endsim = 1;
     }
   }
}

/* interface to jacobian function                                             */

static void C_jac_func_rad(int *neq, double *t, double *y, int *ml,
		    int *mu, double *pd, int *nrowpd, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, Time, ans;

  for (i = 0; i < *neq; i++) REAL(Y)[i] = y[i];

  PROTECT(Time = ScalarReal(*t));
  PROTECT(R_fcall = lang3(R_jac_func,Time,Y));
  PROTECT(ans = eval(R_fcall, R_envir));

  for (i = 0; i < *neq * *nrowpd; i++)  pd[i] = REAL(ans)[i];

  UNPROTECT(3);
}

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 give name to data types                                                    */
typedef void C_solout_type (int *, double *, double *, double *,
  double *, int *, int *, double *, int *, int *, double *) ;

typedef void C_mas_type (int *, double *, int *, double *, int *);

// to be changed...
typedef void C_jac_func_type_rad(int *, double *, double *, int *,
                     int *, double *, int*, double *, int *);


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  MAIN C-FUNCTION, CALLED FROM R-code                       */

SEXP call_radau(SEXP y, SEXP times, SEXP derivfunc, SEXP masfunc, SEXP jacfunc,
    SEXP parms, SEXP rtol, SEXP atol,
    SEXP Nrjac, SEXP Nrmas,
		SEXP rho, SEXP initfunc, SEXP rWork, SEXP iWork,
    SEXP nOut, SEXP lRw, SEXP lIw,
    SEXP Rpar, SEXP Ipar, SEXP Hini, SEXP flist, SEXP elag,
    SEXP rootfunc, SEXP nRoot, SEXP eventfunc, SEXP elist )

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

  int  j, nt, latol, lrtol, lrw, liw,
       ijac, mljac, mujac, imas, mlmas, mumas;
  int  isForcing;
  double *xytmp, tout, *Atol, *Rtol, hini=0;
  int itol, iout, idid;
  int nprot = 0;

  SEXP TROOT, NROOT, VROOT, IROOT;

  /* pointers to functions passed to FORTRAN */
  C_solout_type         *solout = NULL;
  C_jac_func_type_rad   *jac_func = NULL;
  C_mas_type            *mas_func = NULL;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/
/*                      #### initialisation ####                              */

  lock_solver(); /* prevent nested call of solvers that have global variables */

  n_eq = LENGTH(y);             /* number of equations */
  nt   = LENGTH(times);         /* number of output times */
  maxt = nt;
  nroot  = INTEGER(nRoot)[0];   /* number of roots  */
  isroot = 0; nr_root = 0;
  if (nroot > 0) isroot = 1;

  tt = (double *) R_alloc(nt, sizeof(double));
  for (j = 0; j < nt; j++) tt[j] = REAL(times)[j];

  ijac  = INTEGER(Nrjac)[0];
  mljac = INTEGER(Nrjac)[1];
  mujac = INTEGER(Nrjac)[2];
  imas  = INTEGER(Nrmas)[0];
  mlmas = INTEGER(Nrmas)[1];
  mumas = INTEGER(Nrmas)[2];
  /* is function a dll ?*/
  isDll = inherits(derivfunc, "NativeSymbol");

  /* initialise output ... */
  initOutC(isDll, &nout, &ntot, n_eq, nOut, Rpar, Ipar);

  /* copies of variables that will be changed in the FORTRAN subroutine */
  xytmp = (double *) R_alloc(n_eq, sizeof(double));
  for (j = 0; j < n_eq; j++) xytmp[j] = REAL(y)[j];

  ytmp = (double *) R_alloc(n_eq, sizeof(double));

  latol = LENGTH(atol);
  Atol = (double *) R_alloc((int) latol, sizeof(double));
  for (j = 0; j < latol; j++) Atol[j] = REAL(atol)[j];

  lrtol = LENGTH(rtol);
  Rtol = (double *) R_alloc((int) lrtol, sizeof(double));
  for (j = 0; j < lrtol; j++) Rtol[j] = REAL(rtol)[j];

  /* tolerance specifications */
  if (latol == 1 ) itol = 0;
  else             itol = 1;

  hini = REAL(Hini)[0];

  /* work vectors */
  liw = INTEGER (lIw)[0];
  iwork = (int *) R_alloc(liw, sizeof(int));
  for (j=0; j<LENGTH(iWork); j++) iwork[j] = INTEGER(iWork)[j];
  for (j=LENGTH(iWork); j<liw; j++) iwork[j] = 0;

  lrw = INTEGER(lRw)[0];
  rwork = (double *) R_alloc(lrw, sizeof(double));
  for (j=0; j<length(rWork); j++) rwork[j] = REAL(rWork)[j];
  for (j=length(rWork); j<lrw; j++) rwork[j] = 0.;

  /* initialise global R-variables...  */
  //initglobals (nt, ntot);
  PROTECT(Y = allocVector(REALSXP, (n_eq))); nprot++;
  PROTECT(YOUT = allocMatrix(REALSXP, ntot+1, nt)); nprot++;

  //timesteps = (double *) R_alloc(2, sizeof(double));
  for (j=0; j<2; j++) timesteps[j] = 0.;

  /* Initialization of Parameters, Forcings (DLL), lags */
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
  isEvent = initEvents(elist, eventfunc, nroot);
  islag = initLags(elag, 10, nroot);

  if (nout > 0 || islag) {
     xdytmp= (double *) R_alloc(n_eq, sizeof(double));
     for (j = 0; j < n_eq; j++) xdytmp[j] = 0.;
  }

 /* pointers to functions deriv_func, jac_func, passed to FORTRAN */
  if (isDll)  { /* DLL address passed to FORTRAN */
      deriv_func = (C_deriv_func_type *) R_ExternalPtrAddrFn_(derivfunc);

 	   /* overruling deriv_func if forcing */
      if (isForcing) {
        DLL_deriv_func = deriv_func;
        deriv_func = (C_deriv_func_type *) C_deriv_func_forc_rad;
      }
  } else {
      /* interface function between FORTRAN and C/R passed to FORTRAN */
      deriv_func = (C_deriv_func_type *) C_deriv_func_rad;
      /* needed to communicate with R */
      R_deriv_func = derivfunc;
  }
  R_envir = rho;           /* karline: this to allow merging compiled and R-code (e.g. events)*/

  if (!isNull(jacfunc))   {
      if (isDll)
	      jac_func = (C_jac_func_type_rad *) R_ExternalPtrAddrFn_(jacfunc);
	    else  {
	      R_jac_func = jacfunc;
	      jac_func= C_jac_func_rad;
	    }
    }
  if (!isNull(masfunc))   {
	   R_mas_func = masfunc;
	   mas_func= C_mas_func_rad;
     if (isDll)
       R_envir = rho;

  }

 	solout = C_soloutrad;

  iout = 2;                           /* solout called after each step OR 1???*/
  idid = 0;

/*                   ####      integration     ####                           */
  it   = 0;
  tin  = REAL(times)[0];
  tout = REAL(times)[nt-1];
  saveOut (tin, xytmp);               /* save initial condition */
  it++;

  if (nroot > 0)  {      /* also must find a root */
    jroot = (int *) R_alloc(nroot, sizeof(int));
    for (j = 0; j < nroot; j++) jroot[j] = 0;

    root = (double *) R_alloc(nroot, sizeof(double));
    oldroot = (double *) R_alloc(nroot, sizeof(double));

    if (isDll) {
      root_func = (C_root_func_type *) R_ExternalPtrAddrFn_(rootfunc);
    } else {
      root_func = (C_root_func_type *) C_root_radau;
      R_root_func = rootfunc;
    }

    /* value of oldroot */
    if (isDll == 1)
      root_func (&n_eq, &tin, xytmp, &nroot, oldroot);    /* root at t, ytmp */
    else
      C_root_radau (&n_eq, &tin, xytmp, &nroot, oldroot);

    tprevroot = tin; /* to make sure that roots are not too close */
  }
  endsim = 0;
  do {
    if (islag == 1) C_saveLag(1, &tin, xytmp, out, ipar, out, ipar);

    F77_CALL(radau5) ( &n_eq, deriv_func, &tin, xytmp, &tout, &hini,
		     Rtol, Atol, &itol, jac_func, &ijac, &mljac, &mujac,
         mas_func, &imas, &mlmas, &mumas, solout, &iout,
		     rwork, &lrw, iwork, &liw, out, ipar, &idid);
	} while (tin < tout && idid >= 0 && endsim == 0);

  if (idid == -1)
     warning("input is not consistent");
  else if (idid == -2)
     warning("larger maxsteps needed");
  else if (idid == -3)
     warning("step size becomes too small");
  else if (idid == -4)
     warning("problem is probably stiff - interrupted");

/*                   ####  an error occurred   ####                           */
  if(it <= nt-1) saveOut (tin, xytmp);              /* save final condition */
  if (idid < 0) {
    it = it-1;
    PROTECT(YOUT2 = allocMatrix(REALSXP,ntot+1,(it+2))); nprot++;
    returnearly (1, it, ntot);
  } else if (idid == 2) {
    it = it-1;
	PROTECT(YOUT2 = allocMatrix(REALSXP,ntot+1,(it+2))); nprot++;   
    returnearly (0, it, ntot);
    idid = -2;
  }
/*                   ####   returning output   ####                           */
  rwork[0] = hini;
  rwork[1] = tin ;

  PROTECT(ISTATE = allocVector(INTSXP, 7)); nprot++;
  PROTECT(RWORK = allocVector(REALSXP, 5)); nprot++;
  terminate(idid,iwork,7,13,rwork,5,0);

  if (iroot >= 0 || nr_root > 0)  {
    PROTECT(IROOT = allocVector(INTSXP, nroot)); nprot++;
    for (j = 0; j < nroot; j++) INTEGER(IROOT)[j] = jroot[j];
    PROTECT(NROOT = allocVector(INTSXP, 1)); nprot++;
    INTEGER(NROOT)[0] = nr_root;

    if (nr_root == 0) {
      PROTECT(TROOT = allocVector(REALSXP, 1)); nprot++;
      REAL(TROOT)[0] = tin;
    } else {
      if (nr_root > Rootsave) nr_root = Rootsave;

      PROTECT(TROOT = allocVector(REALSXP, nr_root)); nprot++;
      for (j = 0; j < nr_root; j++) REAL(TROOT)[j] = troot[j];

      PROTECT(VROOT = allocVector(REALSXP, nr_root*n_eq)); nprot++;
      for (j = 0; j < nr_root*n_eq; j++) REAL(VROOT)[j] = valroot[j];

      PROTECT(IROOT = allocVector(INTSXP, nr_root)); nprot++;
      for (j = 0; j < nr_root; j++) INTEGER(IROOT)[j] = nrroot[j];

      if (idid == 1) {
        setAttrib(YOUT, install("valroot"), VROOT);
        setAttrib(YOUT, install("indroot"), IROOT);
      }
      else  {
        setAttrib(YOUT2, install("valroot"), VROOT);
        setAttrib(YOUT2, install("indroot"), IROOT);
      }
    }

    if (idid == 1 ) {
      setAttrib(YOUT, install("troot"), TROOT);
      setAttrib(YOUT, install("nroot"), NROOT);
    } else  {
      setAttrib(YOUT2, install("iroot"), IROOT);
      setAttrib(YOUT2, install("troot"), TROOT);
      setAttrib(YOUT2, install("nroot"), NROOT);
    }
  }

/*                   ####     termination      ####                           */
  unlock_solver();
  UNPROTECT(nprot);
  						
// thpe: after reworking PROTECT/UNPROTECT, I checked how YOUT, YOUT2 is handled
// and see that the following is not consistent, because YOUT is only set when idid==1
// Is this a (still) hidden bug?
//

// original version
//  if (idid > 0)
//    return(YOUT);
//  else
//    return(YOUT2);

// thpe: test version (currently disabled)
  if (idid > 0)
    return(YOUT);
  else
    return(YOUT2);	
}

