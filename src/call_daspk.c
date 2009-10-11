/* Patterned on code from package odesolve */
#include <time.h>
#include <string.h>
#include "deSolve.h"

/* definition of the call to the FORTRAN function ddaspk - in file ddaspk.f*/
void F77_NAME(ddaspk)(void (*)(double *, double *, double *, double*,
                               double *, int*, double *, int*),
		     int *, double *, double *, double *, double *, 
		     int *,double *, double *,  int *,  double *,  int *, 
		     int *, int *, double *, int *,
		     void (*)/*(double *, double *, double *, double *, double *, double *, int *)*/,
		     void (*)(int *, double *, double *, double *, double *, double *, 
                  double *, double *, double *, int *, double *, double *, 
                        int *, double *, int *));   

static void forc_daspk (double *t, double *y, double *yprime, double *cj,
                       double *delta, int *ires, double *yout, int *iout)
{
  updatedeforc(t);
  res_fun(t, y, yprime, cj, delta, ires, yout, iout);
}

static void daspk_psol (int *neq, double *t, double *y, double *yprime, 
                        double *savr, double *wk, double *cj, double* wght,
                        double *wp, int *iwp, double *b, double *eplin, 
                        int *ierr, double *RPAR, int *IPAR)
{
/* not yet implemented */
}

/* interface between FORTRAN function call and R function  */

static void daspk_res (double *t, double *y, double *yprime, double *cj, 
                       double *delta, int *ires, double *yout, int *iout)
{                             
  int i;
  SEXP R_fcall, ans;

  REAL(Time)[0] = *t;
  for (i = 0; i < n_eq; i++)
    {
      REAL(Y)[i] = y[i];
      REAL (YPRIME)[i] = yprime[i];
    }
  PROTECT(R_fcall = lang4(daspk_res_func,Time, Y, YPRIME));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, daspk_envir));                  incr_N_Protect();

  for (i = 0; i < n_eq; i++)  	delta[i] = REAL(ans)[i];

  my_unprotect(2);
}

/* deriv output function with rearrangement of state variables and rate of change */

static void daspk_out (int *nout, double *t, double *y, 
                       double *yprime, double *yout)
{
  int i;
  SEXP R_fcall, ans;

  REAL(Time)[0] = *t;
  for (i = 0; i < n_eq; i++)  
    {
      REAL(Y)[i] = y[i];
      REAL (YPRIME)[i] = yprime[i];      
    }
     
  PROTECT(R_fcall = lang4(daspk_res_func,Time, Y, YPRIME));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, daspk_envir));                  incr_N_Protect();

  for (i = 0; i < *nout; i++) yout[i] = REAL(ans)[i + n_eq];

  my_unprotect(2);
}      

/* interface between FORTRAN call to jacobian and R function */

static void daspk_jac (double *t, double *y, double *yprime, 
                       double *pd,  double *cj, double *RPAR, int *IPAR)
{
  int i;
  SEXP R_fcall, ans;

  REAL(Rin)[0] = *t;  
  REAL(Rin)[1] = *cj;  

  for (i = 0; i < n_eq; i++)
    {
      REAL(Y)[i] = y[i];
      REAL (YPRIME)[i] = yprime[i];      
    }
  PROTECT(R_fcall = lang4(daspk_jac_func, Rin, Y, YPRIME));  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, daspk_envir));                 incr_N_Protect();
  for (i = 0; i < n_eq * nrowpd; i++)  pd[i] = REAL(ans)[i];

  my_unprotect(2);
}

/* give name to data types */

typedef void jac_func(double *, double *, double *, double *, double *,
                      double *, int *);
typedef void psol_func(int *, double *, double *, double *, double *,
                      double *, double *, double *, double *, int*, double *,
                      double *, int*, double *, int*);
typedef void kryljac_func(double *, int *, int *, double *, double *,
                          double *, double *, double *,
           double *, double *, double *, double *, int*, int*, double *, int*);

/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_daspk(SEXP y, SEXP yprime, SEXP times, SEXP res, SEXP parms, 
		SEXP rtol, SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc, 
		SEXP psolfunc, SEXP verbose, SEXP info, SEXP iWork, SEXP rWork,  
    SEXP nOut, SEXP maxIt, SEXP bu, SEXP bd, SEXP nRowpd, SEXP Rpar,
    SEXP Ipar, SEXP flist)
{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP   yout, yout2=NULL, ISTATE, RWORK;
  int    i, j, k, nt, ny, repcount, latol, lrtol, lrw, liw, isDll;
  int    maxit, isForcing;
  double *xytmp,  *xdytmp, *rwork, tin, tout, *Atol, *Rtol;
  double *delta=NULL, cj;
  int    *Info,  ninfo, idid, *iwork, mflag, ires;

  res_func  *Resfun;
  jac_func  *jac=NULL;
  psol_func *psol=NULL;
  kryljac_func *kryljac=NULL;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    

  init_N_Protect();

  ny   = LENGTH(y);  
  n_eq = ny;                          /* n_eq is a global variable */
  nt = LENGTH(times);  
  mflag = INTEGER(verbose)[0];        
  nout  = INTEGER(nOut)[0];
  ntot  = n_eq+nout;

  ninfo=LENGTH(info);
  ml = INTEGER(bd)[0]; 
  mu = INTEGER(bu)[0]; 
  nrowpd = INTEGER(nRowpd)[0];  
  maxit = INTEGER(maxIt)[0];

/* function is a dll ?*/
  if (inherits(res, "NativeSymbol")) {
   isDll = 1;
  } else {
   isDll = 0;
  }

/* initialise output ... */
  if (isDll==1)  /* function is a dll */
  {
   lrpar = nout + LENGTH(Rpar);       /* length of rpar */
   lipar = 3    + LENGTH(Ipar);       /* length of ipar */

  } else                              /* function is not a dll */
  {
   lipar = 3;
   lrpar = nout;
  }
   out   = (double *) R_alloc(lrpar, sizeof(double));
   ipar  = (int *)    R_alloc(lipar, sizeof(int));

   if (isDll ==1)
   {
    ipar[0] = nout;          /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;
    /* other elements of ipar are set in R-function lsodx via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar);j++) ipar[j+3] = INTEGER(Ipar)[j];

    /* first nout elements of rpar reserved for output variables
      other elements are set in R-function lsodx via argument *rpar* */
    for (j = 0; j < nout;        j++) out[j] = 0.;
    for (j = 0; j < LENGTH(Rpar);j++) out[nout+j] = REAL(Rpar)[j];
   }

/*  output always done here in C-code (<-> lsode, vode)... */

  /* copies of all variables that will be changed in the FORTRAN subroutine */
  Info  = (int *) R_alloc(ninfo,sizeof(int));
   for (j = 0; j < ninfo; j++) Info[j] = INTEGER(info)[j];  
  
  xytmp = (double *) R_alloc(n_eq, sizeof(double));
   for (j = 0; j < n_eq; j++) xytmp[j] = REAL(y)[j];

  xdytmp = (double *) R_alloc(n_eq, sizeof(double));
   for (j = 0; j < n_eq; j++) xdytmp[j] = REAL(yprime)[j];

  latol = LENGTH(atol);
  Atol  = (double *) R_alloc((int) latol, sizeof(double));
    for (j = 0; j < latol; j++) Atol[j] = REAL(atol)[j];

  lrtol = LENGTH(rtol);
  Rtol  = (double *) R_alloc((int) lrtol, sizeof(double));
    for (j = 0; j < lrtol; j++) Rtol[j] = REAL(rtol)[j];
  
  liw = LENGTH(iWork);
  iwork = (int *) R_alloc(liw, sizeof(int));   
    for (j = 0; j < liw; j++) iwork[j] = INTEGER(iWork)[j];  

  lrw = LENGTH(rWork);
  rwork = (double *) R_alloc(lrw, sizeof(double));
    for (j = 0; j < lrw; j++) rwork[j] = REAL(rWork)[j];

  /* initialise global variables... */
  
  PROTECT(Time = NEW_NUMERIC(1));                    incr_N_Protect();
  PROTECT(Rin  = NEW_NUMERIC(2));                    incr_N_Protect();
  PROTECT(Y = allocVector(REALSXP,(n_eq)));          incr_N_Protect();
  PROTECT(YPRIME = allocVector(REALSXP,(n_eq)));     incr_N_Protect();
  PROTECT(yout = allocMatrix(REALSXP,ntot+1,nt));    incr_N_Protect();

  /**************************************************************************/
  /****** Initialization of Parameters and Forcings (DLL functions)    ******/
  /**************************************************************************/
  initParms(initfunc, parms);

  isForcing = initForcings(flist);

 /* pointers to functions res, psol and jac, passed to the FORTRAN subroutine */

  if (isDll == 1)  {       /* DLL address passed to FORTRAN */
      Resfun = (res_func *) R_ExternalPtrAddr(res);
      /* no need to communicate with R - but output variables set here */            
      delta = (double *) R_alloc(n_eq, sizeof(double));
      for (j = 0; j < n_eq; j++) delta[j] = 0.;

      if(isForcing==1) {
        res_fun = (res_func *) R_ExternalPtrAddr(res);
        Resfun = (res_func *) forc_daspk;
      }

    } else {
      /* interface function between FORTRAN and R passed to FORTRAN */
      Resfun = (res_func *) daspk_res;
      /* needed to communicate with R */      
      daspk_res_func = res; 
      daspk_envir = rho; 

    }
  if (!isNull(jacfunc))
    {
      if (inherits(jacfunc,"NativeSymbol"))
     	{
     	if (Info[11] ==0) {        /*ordinary jac*/
	      jac = (jac_func *) R_ExternalPtrAddr(jacfunc);
	      } else {                /*krylov*/
	      kryljac = (kryljac_func *) R_ExternalPtrAddr(jacfunc);
	      }
	    }
      else  {
	    daspk_jac_func = jacfunc;
	    jac = daspk_jac;
	    }
    }
  if (!isNull(psolfunc))
    {
      if (inherits(psolfunc,"NativeSymbol"))
     	{
	    psol = (psol_func *) R_ExternalPtrAddr(psolfunc);
	    }
      else  {
	    daspk_psol_func = psolfunc;
	    psol = daspk_psol;
	    }
    }

/*                      #### initial time step ####                           */    
  idid = 1;
  REAL(yout)[0] = REAL(times)[0];
  for (j = 0; j < n_eq; j++)
    {
      REAL(yout)[j+1] = REAL(y)[j];
    }
  if (nout>0)
    {
     tin = REAL(times)[0];

	   if (isDll == 1) Resfun (&tin, xytmp, xdytmp, &cj, delta, &ires, out, ipar) ;
	   else daspk_out(&nout,&tin,xytmp,xdytmp,out);
	      for (j = 0; j < nout; j++)
	       REAL(yout)[j + n_eq + 1] = out[j]; 
               }
               
/*                     ####   main time loop   ####                           */    
               
  for (i = 0; i < nt-1; i++)
  {
      tin = REAL(times)[i];
      tout = REAL(times)[i+1];

     repcount = 0;
     do  /* iterations in case maxsteps>500*/
	   {
     	if (Info[11] ==0) {        /*ordinary jac*/
	       F77_CALL(ddaspk) (Resfun, &ny, &tin, xytmp, xdytmp, &tout,
			   Info, Rtol, Atol, &idid, 
			   rwork, &lrw, iwork, &liw, out, ipar, jac, psol);

	      } else {                /*krylov*/
      	 F77_CALL(ddaspk) (Resfun, &ny, &tin, xytmp, xdytmp, &tout,
			   Info, Rtol, Atol, &idid, 
			   rwork, &lrw, iwork, &liw, out, ipar, kryljac, psol);
        }
	  repcount ++;
	  if (idid == -1) 
      {Info[0]=1;
       } else     if (idid == -2)   {
	      warning("Excessive precision requested.  scale up `rtol' and `atol' e.g. by the factor %g\n",10.0);
       Info[0]=1;          
	      repcount=maxit+2;
	    }   else    if (idid == -3)   {
       warning("Error term became zero for some i: pure relative error control (ATOL(i)=0.0) for a variable which is now vanished");
       repcount=maxit+2;
      }   else    if (idid == -5)   {
	      warning("jacfun routine failed with the Krylov method"); 
        repcount = maxit+2;     
      }   else    if (idid == -6)   {
       warning("repeated error test failures on a step - singularity ?");
        repcount = maxit+2;     
      }  else    if (idid == -7)    {
       warning("repeated convergence test failures on a step - inaccurate Jacobian or preconditioner?");
       repcount = maxit+2; 
      }  else    if (idid == -8)    {
       warning("matrix of partial derivatives is singular with direct method-some equations redundant");
       repcount = maxit+2; 
      }  else    if (idid == -9)    {
       warning("repeated convergence test failures and error test failures ?");
       repcount = maxit+2; 
      }  else    if (idid == -10)   {
       warning("repeated convergence test failures on a step, because ires was -1");
       repcount = maxit+2; 
      }  else    if (idid == -11)   {
       warning("unrecoverable error from inside noninear solver, ires=-2 ");
       repcount = maxit+2; 
      }  else    if (idid == -12)   {
       warning("failed to compute initial y and yprime vectors");
       repcount = maxit+2; 
      }  else    if (idid == -13)   {
       warning("unrecoverable error inside the PSOL routine");
       repcount = maxit+2; 
      }  else    if (idid == -14)   {
       warning("Krylov linear system solver failed to converge");
       repcount = maxit+2; 
      }  else    if (idid == -33)   {
       warning("fatal error");
       repcount = maxit+2; 
      }

	} while (tin < tout && repcount < maxit);

 	  REAL(yout)[(i+1)*(ntot+1)] = tin;
	  for (j = 0; j < n_eq; j++)
	    REAL(yout)[(i+1)*(ntot + 1) + j + 1] = xytmp[j];

	  if (nout>0) {
	    if (isDll == 1) Resfun (&tin, xytmp, xdytmp, &cj, delta, &ires, out, ipar) ;
 	    else daspk_out(&nout,&tin,xytmp,xdytmp,out);
      for (j = 0; j < nout; j++)
	       REAL(yout)[(i+1)*(ntot + 1) + j + n_eq + 1] = out[j]; 
               }
               
/*                    ####  an error occurred   ####                          */                     
    if (repcount > maxit || tin < tout) {
	   warning("Returning early from daspk  Results are accurate, as far as they go\n");

	/* redimension yout */
	PROTECT(yout2 = allocMatrix(REALSXP,ntot+1,(i+2)));incr_N_Protect();
	for (k = 0; k < i+2; k++)
	  for (j = 0; j < ntot+1; j++)
	    REAL(yout2)[k*(ntot+1) + j] = REAL(yout)[k*(ntot+1) + j];
	break;
      }
  }    /* end main time loop */

/*                   ####   returning output   ####                           */    

  PROTECT(ISTATE = allocVector(INTSXP, 23));incr_N_Protect();
  for (k = 0;k<21;k++) INTEGER(ISTATE)[k+1] = iwork[k];

  PROTECT(RWORK = allocVector(REALSXP, 3));incr_N_Protect();
  for (k = 1;k<3;k++) REAL(RWORK)[k] = rwork[k+1];
  REAL(RWORK)[0] = rwork[6];

  INTEGER(ISTATE)[0] = idid;  
  if (idid > 0)
    {
      setAttrib(yout, install("istate"), ISTATE);
      setAttrib(yout, install("rstate"), RWORK);    
    }
  else
    {
      setAttrib(yout2, install("istate"), ISTATE);
      setAttrib(yout2, install("rstate"), RWORK);   
    }
    
/*                       ####   termination   ####                            */       
  unprotect_all();
  if (idid > 0)
    return(yout);
  else
    return(yout2);
}
