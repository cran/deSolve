#include <time.h>
#include <string.h>
#include "deSolve.h"

/* definition of the calls to the FORTRAN functions - in file opkdmain.f
and in file dvode.f**/

void F77_NAME(dlsoda)(void (*)(int *, double *, double *, double *, double *, int *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *,int *,int *, int *,
		     void (*)(int *, double *, double *, int *,
			      int *, double *, int *, double *, int *),
		     int *, double *, int *);

void F77_NAME(dlsode)(void (*)(int *, double *, double *, double *, double *, int *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *,int *,int *, int *,
		     void (*)(int *, double *, double *, int *,
			      int *, double *, int *, double *, int *),
		     int *, double *, int *);

void F77_NAME(dlsodes)(void (*)(int *, double *, double *, double *, double *, int *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *, int *, int *, int *,
		     void (*)(int *, double *, double *, int *,
			      int *, int *, double *, double *, int *),   /* jacvec */
		     int *, double *, int *);

void F77_NAME(dlsodar)(void (*)(int *, double *, double *, double *, double *, int *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *,int *,int *, int *,
		     void (*)(int *, double *, double *, int *,
			      int *, double *, int *, double *, int *), int *, 
		     void (*)(int *, double *, double *, int *, double *),  /* rootfunc */
         int *, int *, double *, int *);

void F77_NAME(dvode)(void (*)(int *, double *, double *, double *,
                              double *, int *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *,int *,int *, int *,
		     void (*)(int *, double *, double *, int *,
			            int *, double *, int *, double*, int*),
		     int *, double *, int *);

/* KS: wrapper above the derivate function that first estimates the
values of the forcing functions */

static void forc_lsoda (int *neq, double *t, double *y,
                         double *ydot, double *yout, int *iout)
{
  updatedeforc(t);
  derfun(neq, t, y, ydot, yout, iout);
}

/* interface between FORTRAN function call and R function
   Fortran code calls lsoda_derivs(N, t, y, ydot, yout, iout) 
   R code called as odesolve_deriv_func(time, y) and returns ydot 
   Note: passing of parameter values and "..." is done in R-function lsodx*/

static void lsoda_derivs (int *neq, double *t, double *y,
                          double *ydot, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;

                              REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(odesolve_deriv_func,Time,Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, odesolve_envir));           incr_N_Protect();

  for (i = 0; i < *neq; i++)   ydot[i] = REAL(VECTOR_ELT(ans,0))[i];

  my_unprotect(2);
}

/* only if lsodar: 
   interface between FORTRAN call to root and corresponding R function */

static void lsoda_root (int *neq, double *t, double *y, int *ng, double *gout)
{
  int i;
  SEXP R_fcall, ans;
                              REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(odesolve_root_func,Time,Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, odesolve_envir));          incr_N_Protect();

  for (i = 0; i < *ng; i++)   gout[i] = REAL(ans)[i];

  my_unprotect(2);
}

/* interface between FORTRAN call to jacobian and R function */

static void lsoda_jac (int *neq, double *t, double *y, int *ml,
		    int *mu, double *pd, int *nrowpd, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;

                             REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++) REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(odesolve_jac_func,Time,Y));    incr_N_Protect();
  PROTECT(ans = eval(R_fcall, odesolve_envir));          incr_N_Protect();

  for (i = 0; i < *neq * *nrowpd; i++)  pd[i] = REAL(ans)[i];

  my_unprotect(2);
}

/* only if lsodes: 
   interface between FORTRAN call to jacvec and corresponding R function */

static void lsoda_jacvec (int *neq, double *t, double *y, int *j,
		    int *ian, int *jan, double *pdj, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans, J;
  PROTECT(J = NEW_INTEGER(1));                  incr_N_Protect();
                             INTEGER(J)[0] = *j;
                             REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++) REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang4(odesolve_jac_vec,Time,Y,J));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, odesolve_envir));          incr_N_Protect();

  for (i = 0; i < *neq ; i++)  pdj[i] = REAL(ans)[i];

  my_unprotect(3);
}


/* give name to data types */
typedef void root_func (int *, double *, double *,int *, double *);
typedef void jac_func  (int *, double *, double *, int *,
		                    int *, double *, int *, double *, int *);
typedef void jac_vec   (int *, double *, double *, int *,
		                    int *, int *, double *, double *, int *);

/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_lsoda(SEXP y, SEXP times, SEXP func, SEXP parms, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc,
		SEXP verbose, SEXP iTask, SEXP rWork, SEXP iWork, SEXP jT, 
    SEXP nOut, SEXP lRw, SEXP lIw, SEXP Solver, SEXP rootfunc, 
    SEXP nRoot, SEXP Rpar, SEXP Ipar, SEXP Type, SEXP flist)

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP yout, yout2=NULL, ISTATE, RWORK, IROOT=NULL;    

  int  i, j, k, nt, repcount, latol, lrtol, lrw, liw;
  int  maxit, solver, isForcing;
  double *xytmp, *rwork, tin, tout, *Atol, *Rtol, *dy=NULL, ss;
  int neq, itol, itask, istate, iopt, *iwork, jt, mflag,  is;
  int nroot, *jroot=NULL, isroot,  isDll, type;

  deriv_func *derivs;
  jac_func   *jac=NULL;
  jac_vec    *jacvec;
  root_func  *root=NULL;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    
  init_N_Protect();


  jt  = INTEGER(jT)[0];         /* method flag */
  neq = LENGTH(y);              /* number of equations */ 
  nt  = LENGTH(times);
  
  maxit = 10;                   /* number of iterations */ 
  mflag = INTEGER(verbose)[0];
 
  nroot  = INTEGER(nRoot)[0];   /* number of roots (lsodar) */
  solver = INTEGER(Solver)[0];  /* 1=lsoda,2=lsode,3=lsodeS,4=lsodar,5=vode */

/* is function a dll ?*/
  if (inherits(func, "NativeSymbol")) {
   isDll = 1;
  } else {
   isDll = 0;
  }

/* initialise output ... */
  initOut(isDll, neq, nOut, Rpar, Ipar);

/* copies of all variables that will be changed in the FORTRAN subroutine */

  xytmp = (double *) R_alloc(neq, sizeof(double));
  for (j = 0; j < neq; j++) xytmp[j] = REAL(y)[j];
 
  latol = LENGTH(atol);
  Atol = (double *) R_alloc((int) latol, sizeof(double));

  lrtol = LENGTH(rtol);
  Rtol = (double *) R_alloc((int) lrtol, sizeof(double));

  liw = INTEGER (lIw)[0];
  iwork = (int *) R_alloc(liw, sizeof(int));
     for (j=0; j<LENGTH(iWork); j++) iwork[j] = INTEGER(iWork)[j];

  lrw = INTEGER(lRw)[0];
  rwork = (double *) R_alloc(lrw, sizeof(double));
     for (j=0; j<length(rWork); j++) rwork[j] = REAL(rWork)[j];

/* if a 1-D or 2-D special-purpose problem (lsodes)
   iwork will contain the sparsity structure */

  if (solver ==3)
  {
    type   = INTEGER(Type)[0];
    if (type == 2)        /* 1-D problem ; Type contains further information */
       sparsity1D( Type, iwork, neq, liw) ;
    else if (type == 3)  /* 2-D problem */
       sparsity2D( Type, iwork, neq, liw);
    else if (type == 4)  /* 3-D problem */
     sparsity3D (Type, iwork, neq, liw);
  }

/* initialise global R-variables... */

  PROTECT(Time = NEW_NUMERIC(1));                  incr_N_Protect();
  PROTECT(Y = allocVector(REALSXP,(neq)));         incr_N_Protect();
  PROTECT(yout = allocMatrix(REALSXP,ntot+1,nt));  incr_N_Protect();

  /**************************************************************************/
  /****** Initialization of Parameters and Forcings (DLL functions)    ******/
  /**************************************************************************/
  initParms(initfunc, parms);
  isForcing = initForcings(flist);

/* pointers to functions derivs, jac, jacvec and root, passed to FORTRAN */

  if (isDll) 
    { /* DLL address passed to FORTRAN */
      derivs = (deriv_func *) R_ExternalPtrAddr(func);  
      /* no need to communicate with R - but output variables set here */
      if (isOut) {dy = (double *) R_alloc(neq, sizeof(double));
                  for (j = 0; j < neq; j++) dy[j] = 0.; }
	  
	  /* here overruling derivs if forcing */
      if (isForcing) {
        derfun = (deriv_func *) R_ExternalPtrAddr(func);
        derivs = (deriv_func *) forc_lsoda;
      }
    } else {
      /* interface function between FORTRAN and R passed to FORTRAN */
      derivs = (deriv_func *) lsoda_derivs; 
      /* needed to communicate with R */
      odesolve_deriv_func = func;
      odesolve_envir = rho;
    }

  if (!isNull(jacfunc) && solver !=3)  /* lsodes uses jacvec */
    {
      if (isDll)
	    {
	     jac = (jac_func *) R_ExternalPtrAddr(jacfunc);
	    } else  {
	     odesolve_jac_func = jacfunc;
	     jac = lsoda_jac;
	    }
    }

  if (!isNull(jacfunc) && solver ==3)   /*lsodes*/
    {
      if (isDll)
	    {
	     jacvec = (jac_vec *) R_ExternalPtrAddr(jacfunc);
	    } else  {
	     odesolve_jac_vec = jacfunc;
	     jacvec = lsoda_jacvec;
	    }
    }

  if (solver == 4 && nroot > 0)        /* lsodar */
  { jroot = (int *) R_alloc(nroot, sizeof(int));
     for (j=0; j<nroot; j++) jroot[j] = 0;
  
    if (isDll) 
    {
      root = (root_func *) R_ExternalPtrAddr(rootfunc);
    } else {
      root = (root_func *) lsoda_root;
      odesolve_root_func = rootfunc; 
    }
  }

/* tolerance specifications */
  if (latol == 1 && lrtol == 1 ) itol = 1;
  if (latol  > 1 && lrtol == 1 ) itol = 2;
  if (latol == 1 && lrtol  > 1 ) itol = 3;
  if (latol  > 1 && lrtol  > 1 ) itol = 4;

  for (j = 0; j < lrtol; j++) Rtol[j] = REAL(rtol)[j];
  for (j = 0; j < latol; j++) Atol[j] = REAL(atol)[j];

  itask = INTEGER(iTask)[0];   
  istate = 1;

  iopt = 0;
  ss = 0.;
  is = 0 ;
  for (i = 5; i < 8 ; i++) ss = ss+rwork[i];
  for (i = 5; i < 10; i++) is = is+iwork[i];
  if (ss >0 || is > 0) iopt = 1; /* non-standard input */

/*                      #### initial time step ####                           */    

  REAL(yout)[0] = REAL(times)[0];
  for (j = 0; j < neq; j++) REAL(yout)[j+1] = REAL(y)[j];

  if (isOut == 1) {  /* function in DLL and output */
    tin = REAL(times)[0];
    derivs (&neq, &tin, xytmp, dy, out, ipar) ;
    for (j = 0; j < nout; j++) REAL(yout)[j + neq + 1] = out[j]; 
                  }

/*                     ####   main time loop   ####                           */    
  for (i = 0; i < nt-1; i++) {
    tin = REAL(times)[i];
    tout = REAL(times)[i+1];
    repcount = 0;
    do
	{  /* error control */
 	    if (istate == -2) {
	      for (j = 0; j < lrtol; j++) Rtol[j] *= 10.0;
	      for (j = 0; j < latol; j++) Atol[j] *= 10.0;
	      warning("Excessive precision requested.  `rtol' and `atol' have been scaled upwards by the factor %g\n",10.0);
	      istate = 3;
	    }

      if (solver == 1) {
	      F77_CALL(dlsoda) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, Rtol, Atol, &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt, out, ipar); 
      } else if (solver == 2) {
        F77_CALL(dlsode) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, Rtol, Atol, &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt, out, ipar); 
      } else if (solver == 3) {
        F77_CALL(dlsodes) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, Rtol, Atol, &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jacvec, &jt, out, ipar); 
      } else if (solver == 4) {
        F77_CALL(dlsodar) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, Rtol, Atol,  &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt, root, &nroot, jroot, 
         out, ipar); 
      } else if (solver == 5) {
 	      F77_CALL(dvode) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, Rtol, Atol, &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt, out, ipar);
      }
    
	    if (istate == -1)  {
        warning("an excessive amount of work (> maxsteps ) was done, but integration was not successful - increase maxsteps");
      } else if (istate == 3 && solver == 4){
	      istate = -20;  repcount = 50;
      } else if (istate == -2)  {
	      warning("Excessive precision requested.  scale up `rtol' and `atol' e.g by the factor %g\n",10.0);
	    } else if (istate == -4)  {
        warning("repeated error test failures on a step, but integration was successful - singularity ?");
      } else if (istate == -5)  {
        warning("repeated convergence test failures on a step, but integration was successful - inaccurate Jacobian matrix?");
      } else if (istate == -6)  {
        warning("Error term became zero for some i: pure relative error control (ATOL(i)=0.0) for a variable which is now vanished");
      }

	    repcount ++;
	} while (tin < tout && istate >= 0 && repcount < maxit); 
	
  if (istate == -3)  {
    error("illegal input detected before taking any integration steps - see written message");
	  unprotect_all();
	}  else	{
	  REAL(yout)[(i+1)*(ntot+1)] = tin;
	  for (j = 0; j < neq; j++)
	    REAL(yout)[(i+1)*(ntot + 1) + j + 1] = xytmp[j];
	  if (isOut == 1)  {
      derivs (&neq, &tin, xytmp, dy, out, ipar) ;
	    for (j = 0; j < nout; j++)
	      REAL(yout)[(i+1)*(ntot + 1) + j + neq + 1] = out[j];
    }
	}
	  
/*                    ####  an error occurred   ####                          */    
   if (istate < 0 || tin < tout) {
	  if (istate != -20) warning("Returning early. Results are accurate, as far as they go\n");

	 /* need to redimension yout here, and add the attribute "istate" for */
	 /* the most recent value of `istate' from lsodx */
	  PROTECT(yout2 = allocMatrix(REALSXP,ntot+1,(i+2)));incr_N_Protect();
    for (k = 0; k < i+2; k++)
	    for (j = 0; j < ntot+1; j++)
	     REAL(yout2)[k*(ntot+1) + j] = REAL(yout)[k*(ntot+1) + j];
	    break;
    }
  }     /* end main time loop */

/*                   ####   returning output   ####                           */    
  if (istate == -20 && nroot > 0)  {
    isroot = 1   ;
    PROTECT(IROOT = allocVector(INTSXP, nroot));incr_N_Protect();
    for (k = 0;k<nroot;k++) INTEGER(IROOT)[k] = jroot[k];
  } else isroot = 0;

  PROTECT(ISTATE = allocVector(INTSXP, 23));incr_N_Protect();
  for (k = 0;k<22;k++) INTEGER(ISTATE)[k+1] = iwork[k];
        
  PROTECT(RWORK = allocVector(REALSXP, 5));incr_N_Protect();
  for (k = 0;k<5;k++) REAL(RWORK)[k] = rwork[k+10];

  INTEGER(ISTATE)[0] = istate;  
  if (istate == -20) INTEGER(ISTATE)[0] = 3; 	  
  if (istate > 0) {
    setAttrib(yout, install("istate"), ISTATE);
    setAttrib(yout, install("rstate"), RWORK);
    if (isroot==1) setAttrib(yout, install("iroot"), IROOT);
  }
  else  {
    setAttrib(yout2, install("istate"), ISTATE);
    setAttrib(yout2, install("rstate"), RWORK);
    if (isroot==1) setAttrib(yout2, install("iroot"), IROOT);
  }

/*                       ####   termination   ####                            */    
  unprotect_all();
  if (istate > 0)
    return(yout);
  else
    return(yout2);
}

