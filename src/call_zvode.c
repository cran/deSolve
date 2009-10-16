/* complex number vode */
#include <time.h>
#include <string.h>
#include "deSolve.h"   
#include "zvode.h"   

SEXP cvode_deriv_func;
SEXP cvode_jac_func;
SEXP vode_envir;
                           
/* definition of the call to the FORTRAN function dvode - in file zvode.f*/
void F77_NAME(zvode)(void (*)(int *, double *, Rcomplex *, Rcomplex *,
                              Rcomplex *, int *),
		     int *, Rcomplex *, double *, double *,
		     int *, double *,  double *, int *, int *,
		     int *, Rcomplex *, int*, double *, int *,int *, int *,
		     void (*)(int *, double *, Rcomplex *, int *,
			            int *, Rcomplex *, int *, Rcomplex*, int*),
		     int *, Rcomplex *, int *);

/* interface between FORTRAN function call and R function
   Fortran code calls cvode_derivs(N, t, y, ydot, yout, iout) 
   R code called as cvode_deriv_func(time, y) and returns ydot 
   Note: passing of parameter values and "..." is done in R-function zvode*/

static void cvode_derivs (int *neq, double *t, Rcomplex *y, 
                         Rcomplex *ydot, Rcomplex *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;     

                              REAL(Time)[0]   = *t;
  for (i = 0; i < *neq; i++)  COMPLEX (cY)[i] = y[i];

  PROTECT(R_fcall = lang3(cvode_deriv_func,Time,cY)) ;incr_N_Protect();
  PROTECT(ans = eval(R_fcall, vode_envir))           ;incr_N_Protect();

  for (i = 0; i < *neq; i++)	ydot[i] = COMPLEX(VECTOR_ELT(ans,0))[i];

  my_unprotect(2);      
}

/* interface between FORTRAN call to jacobian and R function */

static void cvode_jac (int *neq, double *t, Rcomplex *y, int *ml,
		    int *mu, Rcomplex *pd, int *nrowpd, Rcomplex *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;

                              REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  COMPLEX(cY)[i] = y[i];

  PROTECT(R_fcall = lang3(cvode_jac_func,Time,cY));  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, vode_envir));        incr_N_Protect();

  for (i = 0; i < *neq * *nrowpd; i++)  pd[i ] = COMPLEX(ans)[i ];

  my_unprotect(2);
}

/* wrapper above the derivate function that first estimates the
values of the forcing functions */

static void forc_zvode (int *neq, double *t, Rcomplex *y,
                         Rcomplex *ydot, Rcomplex *yout, int *iout)
{
  updatedeforc(t);
  cderfun(neq, t, y, ydot, yout, iout);
}


/* give name to data types */
typedef void cjac_func(int *, double *, Rcomplex *, int *,
		                  int *, Rcomplex *, int *, Rcomplex *, int *);

/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_zvode(SEXP y, SEXP times, SEXP func, SEXP parms, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc, 
		SEXP verbose, SEXP iTask, SEXP rWork, SEXP iWork, SEXP jT, SEXP nOut, 
    SEXP lZw, SEXP lRw, SEXP lIw, SEXP Rpar, SEXP Ipar, SEXP flist)
    
{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP   yout, yout2=NULL, ISTATE, RWORK;

  int    i, j, k, nt, latol, lrtol, lrw, liw, lzw;
  double*rwork, tin, tout, *Atol, *Rtol, ss;
  int    neq, itol, itask, istate, iopt, *iwork, jt, mflag, nout, 
         is, isDll, isForcing;
  Rcomplex  *xytmp, *dy=NULL, *zwork;
  
  cderiv_func *derivs;
  cjac_func   *jac=NULL;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    

  init_N_Protect();

  jt = INTEGER(jT)[0];        
  neq = LENGTH(y);
  nt = LENGTH(times);

  mflag = INTEGER(verbose)[0];
  nout  = INTEGER(nOut)[0];
  
/* The output:
    zout and ipar are used to pass output variables (number set by nout)
    followed by other input (e.g. forcing functions) provided 
    by R-arguments rpar, ipar
    ipar[0]: number of output variables, ipar[1]: length of rpar, 
    ipar[2]: length of ipar */

/* is function a dll ?*/
  if (inherits(func, "NativeSymbol")) {
   isDll = 1;
  } else {
   isDll = 0;
  }

/* initialise output for Complex variables ... */
  initOutC(isDll, neq, nOut, Rpar, Ipar);

/* copies of all variables that will be changed in the FORTRAN subroutine */
 
  xytmp = (Rcomplex *) R_alloc(neq, sizeof(Rcomplex));
    for (j = 0; j < neq; j++) xytmp[j] = COMPLEX(y)[j];

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

  lzw = INTEGER (lZw)[0];
  zwork = (Rcomplex *) R_alloc(lzw, sizeof(Rcomplex));

  /* initialise global R-variables... */
  
  PROTECT(Time = NEW_NUMERIC(1))                  ;incr_N_Protect(); 
  PROTECT(cY = allocVector(CPLXSXP , neq) )       ;incr_N_Protect();        
  PROTECT(yout = allocMatrix(CPLXSXP,ntot+1,nt))  ;incr_N_Protect();
  
  /**************************************************************************/
  /****** Initialization of Parameters and Forcings (DLL functions)    ******/
  /**************************************************************************/
  initParms(initfunc, parms);
  isForcing = initForcings(flist);

/* pointers to functions derivs and jac, passed to the FORTRAN subroutine */

  if (isDll==1) 
    {/* DLL address passed to FORTRAN */
      derivs = (cderiv_func *) R_ExternalPtrAddr(func);
      /* no need to communicate with R - but output variables set here */      
      if (isOut) {dy = (Rcomplex *) R_alloc(neq, sizeof(Rcomplex));
                  // for (j = 0; j < neq; j++) dy[j] =  i0; 
                  }
	  /* here overruling derivs if forcing */
      if (isForcing) {
        cderfun = (cderiv_func *) R_ExternalPtrAddr(func);
        derivs = (cderiv_func *) forc_zvode;
      }

    } else {  
      /* interface function between FORTRAN and R passed to FORTRAN*/
      derivs = (cderiv_func *) cvode_derivs;  
      /* needed to communicate with R */
      cvode_deriv_func = func; 
      vode_envir = rho;       
  
    }
    
   if (!isNull(jacfunc))
    {
      if (isDll == 1)
     	{
	    jac = (cjac_func *) R_ExternalPtrAddr(jacfunc);
	    } else {
	    cvode_jac_func = jacfunc;
	    jac = cvode_jac;
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

/*  COMPLEX(yout)[0] = COMPLEX(times)[0];*/
  for (j = 0; j < neq; j++)
    {
      COMPLEX(yout)[j+1] = COMPLEX(y)[j];
    }      /* function in DLL and output */

	  if (isOut == 1) {
        tin = REAL(times)[0];
        derivs (&neq, &tin, xytmp, dy, zout, ipar) ;
	      for (j = 0; j < nout; j++)
	       COMPLEX(yout)[j + neq + 1] = zout[j]; 
               }  

/*                     ####   main time loop   ####                           */    

  for (i = 0; i < nt-1; i++)
  {
    tin = REAL(times)[i];
    tout = REAL(times)[i+1];
      
 	  F77_CALL(zvode) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, Rtol, Atol, &itask, &istate, &iopt, zwork, &lzw, rwork,
			   &lrw, iwork, &liw, jac, &jt, zout, ipar);
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
 	/*   REAL(yout)[(i+1)*(ntot+1)] = tin;*/
	   for (j = 0; j < neq; j++)
	    COMPLEX(yout)[(i+1)*(ntot + 1) + j + 1] = xytmp[j];
   
	  if (isOut == 1) 
     {
        derivs (&neq, &tin, xytmp, dy, zout, ipar) ;
	      for (j = 0; j < nout; j++)
	       COMPLEX(yout)[(i+1)*(ntot + 1) + j + neq + 1] = zout[j]; 
      }
    } 

/*                    ####  an error occurred   ####                          */      
  if (istate < 0 || tin < tout) {
	  warning("Returning early from dvode  Results are accurate, as far as they go\n");

	 /* redimension yout */
	 PROTECT(yout2 = allocMatrix(CPLXSXP,ntot+1,(i+2)));incr_N_Protect();

	 for (k = 0; k < i+2; k++)
	   for (j = 0; j < ntot+1; j++)
	     COMPLEX(yout2)[k*(ntot+1) + j] = COMPLEX(yout)[k*(ntot+1) + j];
	     break;
    }
  }  /* end main time loop */

/*                   ####   returning output   ####                           */    
       
  PROTECT(ISTATE = allocVector(INTSXP, 23));incr_N_Protect();
  for (k = 0;k<22;k++) INTEGER(ISTATE)[k+1] = iwork[k];

  PROTECT(RWORK = allocVector(REALSXP, 4));incr_N_Protect();
  for (k = 0;k<4;k++) REAL(RWORK)[k] = rwork[k+10]; 
  /*  */

  INTEGER(ISTATE)[0] = istate;  
  if (istate > 0)
    {
      setAttrib(yout, install("istate"), ISTATE);
       setAttrib(yout, install("rstate"), RWORK);       
  } else {
      setAttrib(yout2, install("istate"), ISTATE);
      setAttrib(yout2, install("rstate"), RWORK);    
  }
/*                       ####   termination   ####                            */         
  unprotect_all();
  if (istate > 0)
    return(yout);
  else
    return(yout2);
}



