/* Patterned on code call_lsoda.c from package odesolve */
#include <time.h>
#include <string.h>
#include "deSolve.h"   
                           
/* definition of the call to the fortran function dvode - in file dvode.f*/
void F77_NAME(dvode)(void (*)(int *, double *, double *, double *,
                              double *, int *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *,int *,int *, int *,
		     void (*)(int *, double *, double *, int *,
			            int *, double *, int *, double*, int*),
		     int *, double *, int *);

/* interface between fortran function call and R function 
   Fortran code calls vode_derivs(N, t, y, ydot, yout, iout) 
   R code called as vode_deriv_func(time, y) and returns ydot 
   Note: passing of parameter values and "..." is done in R-function vode*/

static void vode_derivs (int *neq, double *t, double *y, 
                         double *ydot, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;     

                              REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(vode_deriv_func,Time,Y)) ;incr_N_Protect();
  PROTECT(ans = eval(R_fcall, vode_envir))         ;incr_N_Protect();

  for (i = 0; i < *neq; i++)	ydot[i] = REAL(VECTOR_ELT(ans,0))[i];

  my_unprotect(2);      
}

/* interface between fortran call to jacobian and R function */

static void vode_jac (int *neq, double *t, double *y, int *ml,
		    int *mu, double *pd, int *nrowpd, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;

                              REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(vode_jac_func,Time,Y));  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, vode_envir));        incr_N_Protect();

  for (i = 0; i < *neq * *nrowpd; i++)  pd[i ] = REAL(ans)[i ];

  my_unprotect(2);
}

/* give name to data types */
typedef void deriv_func(int *, double *, double *,double *,double *, int *);
typedef void jac_func(int *, double *, double *, int *,
		                  int *, double *, int *, double *, int *);
typedef void init_func(void (*)(int *, double *));

/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_dvode(SEXP y, SEXP times, SEXP func, SEXP parms, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc, 
		SEXP verbose, SEXP iTask, SEXP rWork, SEXP iWork, SEXP jT, SEXP nOut, 
    SEXP lIw, SEXP lRw, SEXP Rpar, SEXP Ipar)
    
{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP   yout, yout2=NULL, ISTATE, RWORK;

  int    i, j, k, nt, latol, lrtol, lrw, liw, isOut, maxit;
  double *xytmp, *rwork, tin, tout, *Atol, *Rtol, *out, *dy=NULL, ss;
  int    neq, itol, itask, istate, iopt, *iwork, jt, mflag, nout, 
         lrpar, lipar, ntot, is, *ipar, isDll;

  deriv_func *derivs;
  jac_func   *jac=NULL;
  init_func  *initializer;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    

  init_N_Protect();

  jt = INTEGER(jT)[0];        
  neq = LENGTH(y);
  nt = LENGTH(times);

  maxit = 1;  /* iterations not allowed... */
  mflag = INTEGER(verbose)[0];
  
  nout  = INTEGER(nOut)[0];
  
/* The output:
    out and ipar are used to pass output variables (number set by nout)
    followed by other input (e.g. forcing functions) provided 
    by R-arguments rpar, ipar
    ipar[0]: number of output variables, ipar[1]: length of rpar, 
    ipar[2]: length of ipar */

  isOut = 0;
  if (inherits(func, "NativeSymbol"))  /* function is a dll */
  {
   isDll = 1;
   if (nout > 0) isOut = 1; 
   ntot  = neq + nout;          /* length of yout */
   lrpar = nout + LENGTH(Rpar); /* length of rpar; LENGTH(Rpar) is always >0 */
   lipar = 3 + LENGTH(Ipar);    /* length of ipar */

  } else                              /* function is not a dll */
  {
   isDll = 0;
   isOut = 0;
   ntot = neq;
   lipar = 1;
   lrpar = 1; 
  }
 
   out   = (double *) R_alloc(lrpar, sizeof(double));
   ipar  = (int *)    R_alloc(lipar, sizeof(int));

   if (isDll ==1)
   {
    ipar[0] = nout;             /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;
    /* other elements of ipar are set in R-function lsodx via argument *ipar* */    
    for (j = 0; j < LENGTH(Ipar);j++) ipar[j+3] = INTEGER(Ipar)[j];
    
    /* first nout elements of rpar reserved for output variables 
      other elements are set in R-function lsodx via argument *rpar* */
    for (j = 0; j < nout; j++) out[j] = 0.;  
    for (j = 0; j < LENGTH(Rpar);j++) out[nout+j] = REAL(Rpar)[j];
   }
   
  /* copies of all variables that will be changed in the FORTRAN subroutine */
 
  xytmp = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < neq; j++) xytmp[j] = REAL(y)[j];

  latol = LENGTH(atol);
  Atol = (double *) R_alloc((int) latol, sizeof(double));
    for (j = 0; j < latol; j++) Atol[j] = REAL(atol)[j];

  lrtol = LENGTH(rtol);
  Rtol = (double *) R_alloc((int) lrtol, sizeof(double));
    for (j = 0; j < lrtol; j++) Rtol[j] = REAL(rtol)[j];

  liw = INTEGER (lIw)[0];
  iwork = (int *) R_alloc(liw, sizeof(int));   
    for (j = 0; j < 30; j++) iwork[j] = INTEGER(iWork)[j];  

  lrw = INTEGER (lRw)[0];
  rwork = (double *) R_alloc(lrw, sizeof(double));
    for (j = 0; j < 20; j++) rwork[j] = REAL(rWork)[j];

  /* initialise global variables... */
            
  PROTECT(Time = NEW_NUMERIC(1))                  ;incr_N_Protect(); 
  PROTECT(Y = allocVector(REALSXP,( neq)))        ;incr_N_Protect();        
  PROTECT(yout = allocMatrix(REALSXP,ntot+1,nt));  incr_N_Protect();
  PROTECT(de_gparms = parms);                      incr_N_Protect();  

 /* The initialisation routine */
  if (!isNull(initfunc))
    	{
	     initializer = (init_func *) R_ExternalPtrAddr(initfunc);
	     initializer(Initdeparms); 	}

 /* pointers to functions derivs and jac, passed to the FORTRAN subroutine */

  if (isDll==1) 
    {/* DLL address passed to fortran */
      derivs = (deriv_func *) R_ExternalPtrAddr(func);
      /* no need to communicate with R - but output variables set here */      
      if (isOut) {dy = (double *) R_alloc(neq, sizeof(double));
                  for (j = 0; j < neq; j++) dy[j] = 0.; }

    } else {  
      /* interface function between fortran and R passed to Fortran*/     
      derivs = (deriv_func *) vode_derivs;  
      /* needed to communicate with R */
      vode_deriv_func = func; 
      vode_envir = rho;       
  
    }
    
   if (!isNull(jacfunc))
    {
      if (inherits(jacfunc,"NativeSymbol"))
     	{
	    jac = (jac_func *) R_ExternalPtrAddr(jacfunc);
	    } else {
	    vode_jac_func = jacfunc;
	    jac = vode_jac;
	    }
    }
/* tolerance specifications */
  if (latol == 1 && lrtol == 1 ) itol = 1;
  if (latol  > 1 && lrtol == 1 ) itol = 2;
  if (latol == 1 && lrtol  > 1 ) itol = 3;
  if (latol  > 1 && lrtol  > 1 ) itol = 4;

  itask = INTEGER(iTask)[0]; 
  istate = 1;

  iopt = 0;
  ss = 0.;
  is = 0;
  for (i = 5; i < 8 ; i++) ss = ss+rwork[i];
  for (i = 5; i < 10; i++) is = is+iwork[i];
  if (ss >0 || is > 0) iopt = 1;  /* non-standard input */

/*                      #### initial time step ####                           */    

  REAL(yout)[0] = REAL(times)[0];
  for (j = 0; j < neq; j++)
    {
      REAL(yout)[j+1] = REAL(y)[j];
    }
	  if (isOut == 1) { /* function in DLL and output */
        tin = REAL(times)[0];
        derivs (&neq, &tin, xytmp, dy, out, ipar) ;
	      for (j = 0; j < nout; j++)
	       REAL(yout)[j + neq + 1] = out[j]; 
               }

/*                     ####   main time loop   ####                           */    

  for (i = 0; i < nt-1; i++)
  {
    tin = REAL(times)[i];
    tout = REAL(times)[i+1];
      
 	  F77_CALL(dvode) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, Rtol, Atol, &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt, out, ipar);
	  if (istate == -1) {
      warning("an excessive amount of work (> mxstep ) was done, but integration was not successful - increase maxsteps ?");
      }
	  else if (istate == -2)  {
	      warning("Excessive precision requested.  scale up `rtol' and `atol' e.g by the factor %g\n",10.0);
	    }
	  else if (istate == -4)  {
       warning("repeated error test failures on a step, but integration was successful - singularity ?");
      }
    else if (istate == -5)  {
       warning("repeated convergence test failures on a step, but integration was successful - inaccurate Jacobian matrix?");
      }
    else if (istate == -6)  {
       warning("Error term became zero for some i: pure relative error control (ATOL(i)=0.0) for a variable which is now vanished");
      }      
    
    if (istate == -3) 	{
	      error("illegal input detected before taking any integration steps - see written message"); 
    	  unprotect_all();
  	}
      else
	  {
 	   REAL(yout)[(i+1)*(ntot+1)] = tin;
	   for (j = 0; j < neq; j++)
	    REAL(yout)[(i+1)*(ntot + 1) + j + 1] = xytmp[j];
   
	   if (isOut == 1) 
     {
        derivs (&neq, &tin, xytmp, dy, out, ipar) ;
	      for (j = 0; j < nout; j++)
	       REAL(yout)[(i+1)*(ntot + 1) + j + neq + 1] = out[j]; 
     }
    } 

/*                    ####  an error occurred   ####                          */      
  if (istate < 0 || tin < tout) {
	  warning("Returning early from dvode  Results are accurate, as far as they go\n");

	 /* redimension yout */
	 PROTECT(yout2 = allocMatrix(REALSXP,ntot+1,(i+2)));incr_N_Protect();
	 for (k = 0; k < i+2; k++)
	   for (j = 0; j < ntot+1; j++)
	     REAL(yout2)[k*(ntot+1) + j] = REAL(yout)[k*(ntot+1) + j];
	     break;
    }
  }  /* end main time loop */

/*                   ####   returning output   ####                           */    
       
  PROTECT(ISTATE = allocVector(INTSXP, 23));incr_N_Protect();
  for (k = 0;k<22;k++) INTEGER(ISTATE)[k+1] = iwork[k];

  PROTECT(RWORK = allocVector(REALSXP, 4));incr_N_Protect();
  for (k = 0;k<4;k++) REAL(RWORK)[k] = rwork[k+10];

  INTEGER(ISTATE)[0] = istate;  
  if (istate > 0)
    {
      setAttrib(yout, install("istate"), ISTATE);
      setAttrib(yout, install("rstate"), RWORK);    }
  else
    {
      setAttrib(yout2, install("istate"), ISTATE);
      setAttrib(yout2, install("rstate"), RWORK);    }

/*                       ####   termination   ####                            */         
  unprotect_all();
  if (istate > 0)
    return(yout);
  else
    return(yout2);
}

