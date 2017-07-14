/* complex number vode */
#include <time.h>
#include <string.h>
#include "deSolve.h"   
#include "zvode.h"
#include "externalptr.h"
   
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Ordinary differential equation solver for complex state variables, zvode.
   
   The C-wrappers that provide the interface between FORTRAN codes and R-code 
   are: C_zderiv_func: interface with R-code "func", passes derivatives  
        C_zjac_func  : interface with R-code "jacfunc", passes jacobian (except lsodes)
  
   C_zderiv_func_forc provides the interface between the function specified in
   a DLL and the integrator, in case there are forcing functions.
   
   Events and roots are not implemented for zvode
  
   changes since 1.4
   karline: version 1.5: added forcing functions in DLL
            improving names
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

SEXP R_zderiv_func;
SEXP R_zjac_func;
SEXP R_vode_envir;
                           
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
   R code called as R_zderiv_func(time, y) and returns ydot 
   Note: passing of parameter values and "..." is done in R-function zvode*/

static void C_zderiv_func (int *neq, double *t, Rcomplex *y, 
                         Rcomplex *ydot, Rcomplex *yout, int *iout)
{
  int i;
  SEXP R_fcall, Time, ans;     

  for (i = 0; i < *neq; i++)  COMPLEX(cY)[i] = y[i];

  PROTECT(Time = ScalarReal(*t));                  incr_N_Protect();
  PROTECT(R_fcall = lang3(R_zderiv_func,Time,cY)) ;incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_vode_envir))           ;incr_N_Protect();

  for (i = 0; i < *neq; i++)	ydot[i] = COMPLEX(VECTOR_ELT(ans,0))[i];

  my_unprotect(3);      
}

/* interface between FORTRAN call to jacobian and R function */

static void C_zjac_func (int *neq, double *t, Rcomplex *y, int *ml,
		    int *mu, Rcomplex *pd, int *nrowpd, Rcomplex *yout, int *iout)
{
  int i;
  SEXP R_fcall, Time, ans;

  for (i = 0; i < *neq; i++)  COMPLEX(cY)[i] = y[i];

  PROTECT(Time = ScalarReal(*t));                 incr_N_Protect();
  PROTECT(R_fcall = lang3(R_zjac_func,Time,cY));  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_vode_envir));     incr_N_Protect();

  for (i = 0; i < *neq * *nrowpd; i++)  pd[i ] = COMPLEX(ans)[i ];

  my_unprotect(3);
}

/* wrapper above the derivate function that first estimates the
values of the forcing functions */

static void C_zderiv_func_forc (int *neq, double *t, Rcomplex *y,
                         Rcomplex *ydot, Rcomplex *yout, int *iout)
{
  updatedeforc(t);
  DLL_cderiv_func(neq, t, y, ydot, yout, iout);
}


/* give name to data types */
typedef void C_zjac_func_type(int *, double *, Rcomplex *, int *,
		                  int *, Rcomplex *, int *, Rcomplex *, int *);

/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_zvode(SEXP y, SEXP times, SEXP derivfunc, SEXP parms, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc, 
		SEXP iTask, SEXP rWork, SEXP iWork, SEXP jT, SEXP nOut,
    SEXP lZw, SEXP lRw, SEXP lIw, SEXP Rpar, SEXP Ipar, SEXP flist)
    
{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

  int    i, j, k, nt, latol, lrtol, lrw, liw, lzw;
  double tin, tout, *Atol, *Rtol, ss;
  int    neq, itol, itask, istate, iopt, jt, //mflag, 
         is, isDll, isForcing;
  Rcomplex  *xytmp, *dy = NULL, *zwork;
  int    *iwork, it, ntot, nout;   
  double *rwork;  
  C_zderiv_func_type *zderiv_func;
  C_zjac_func_type   *zjac_func = NULL;

/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

  lock_solver(); /* prevent nested call of solvers that have global variables */

/*                      #### initialisation ####                              */    

  //init_N_Protect();
  long int old_N_Protect = save_N_Protected();  

  jt = INTEGER(jT)[0];        
  neq = LENGTH(y);
  nt = LENGTH(times);

  nout  = INTEGER(nOut)[0];
  
/* The output:
    zout and ipar are used to pass output variables (number set by nout)
    followed by other input (e.g. forcing functions) provided 
    by R-arguments rpar, ipar
    ipar[0]: number of output variables, ipar[1]: length of rpar, 
    ipar[2]: length of ipar */

/* is function a dll ?*/
  if (inherits(derivfunc, "NativeSymbol")) {
    isDll = 1;
  } else {
    isDll = 0;
  }

/* initialise output for Complex variables ... */
  initOutComplex(isDll, &nout, &ntot, neq, nOut, Rpar, Ipar);

/* copies of all variables that will be changed in the FORTRAN subroutine */
 
  xytmp = (Rcomplex *) R_alloc(neq, sizeof(Rcomplex));
  for (j = 0; j < neq; j++) xytmp[j] = COMPLEX(y)[j];

  latol = LENGTH(atol);
  Atol = (double *) R_alloc((int) latol, sizeof(double));
  for (j = 0; j < latol; j++) Atol[j] = REAL(atol)[j];

  lrtol = LENGTH(rtol);
  Rtol = (double *) R_alloc((int) lrtol, sizeof(double));
  for (j = 0; j < lrtol; j++) Rtol[j] = REAL(rtol)[j];

  liw = INTEGER(lIw)[0];
  iwork = (int *) R_alloc(liw, sizeof(int));   
  for (j = 0; j < 30; j++) iwork[j] = INTEGER(iWork)[j];  

  lrw = INTEGER(lRw)[0];
  rwork = (double *) R_alloc(lrw, sizeof(double));
  for (j = 0; j < 20; j++) rwork[j] = REAL(rWork)[j];

  /* global variable */
  //timesteps = (double *) R_alloc(2, sizeof(double));
  for (j=0; j<2; j++) timesteps[j] = 0.;

  lzw = INTEGER(lZw)[0];
  zwork = (Rcomplex *) R_alloc(lzw, sizeof(Rcomplex));

  /* initialise global R-variables... */
  
  PROTECT(cY = allocVector(CPLXSXP , neq) )       ;incr_N_Protect();        
  PROTECT(YOUT = allocMatrix(CPLXSXP,ntot+1,nt))  ;incr_N_Protect();
  
  /**************************************************************************/
  /****** Initialization of Parameters and Forcings (DLL functions)    ******/
  /**************************************************************************/
  initParms(initfunc, parms);
  isForcing = initForcings(flist);

/* pointers to functions zderiv_func and zjac_func, passed to the FORTRAN subroutine */

  if (isDll == 1) { /* DLL address passed to FORTRAN */
    zderiv_func = (C_zderiv_func_type *) R_ExternalPtrAddrFn_(derivfunc);
    /* no need to communicate with R - but output variables set here */      
    if (isOut) {
      dy = (Rcomplex *) R_alloc(neq, sizeof(Rcomplex));
      /* for (j = 0; j < neq; j++) dy[j] =  i0; */
    }
	  /* here overruling zderiv_func if forcing */
    if (isForcing) {
      DLL_cderiv_func = (C_zderiv_func_type *) R_ExternalPtrAddrFn_(derivfunc);
      zderiv_func = (C_zderiv_func_type *) C_zderiv_func_forc;
    }
  } else {  
    /* interface function between FORTRAN and R passed to FORTRAN*/
    zderiv_func = (C_zderiv_func_type *) C_zderiv_func;  
    /* needed to communicate with R */
    R_zderiv_func = derivfunc; 
    R_vode_envir = rho;       
  }
  
  if (!isNull(jacfunc)) {
    if (isDll == 1) {
	    zjac_func = (C_zjac_func_type *) R_ExternalPtrAddrFn_(jacfunc);
    } else {
	    R_zjac_func = jacfunc;
	    zjac_func = C_zjac_func;
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

/*  COMPLEX(YOUT)[0] = COMPLEX(times)[0];*/
  for (j = 0; j < neq; j++) {
    COMPLEX(YOUT)[j+1] = COMPLEX(y)[j];
  }      /* function in DLL and output */

  if (isOut == 1) {
    tin = REAL(times)[0];
    zderiv_func (&neq, &tin, xytmp, dy, zout, ipar) ;
    for (j = 0; j < nout; j++)
      COMPLEX(YOUT)[j + neq + 1] = zout[j]; 
  }  
/*                     ####   main time loop   ####                           */    
  for (it = 0; it < nt-1; it++) {
    tin = REAL(times)[it];
    tout = REAL(times)[it+1];
      
 	  F77_CALL(zvode) (zderiv_func, &neq, xytmp, &tin, &tout,
			   &itol, Rtol, Atol, &itask, &istate, &iopt, zwork, &lzw, rwork,
			   &lrw, iwork, &liw, zjac_func, &jt, zout, ipar);
	  
    /* in case size of timesteps is called for */
    timesteps [0] = rwork[10];
    timesteps [1] = rwork[11];
    
    if (istate == -1) {
      warning("an excessive amount of work (> mxstep ) was done, but integration was not successful - increase maxsteps ?");
    } else if (istate == -2)  {
	    warning("Excessive precision requested.  scale up `rtol' and `atol' e.g by the factor %g\n",10.0);
	  } else if (istate == -4)  {
      warning("repeated error test failures on a step, but integration was successful - singularity ?");
    } else if (istate == -5)  {
      warning("repeated convergence test failures on a step, but integration was successful - inaccurate Jacobian matrix?");
    } else if (istate == -6)  {
      warning("Error term became zero for some i: pure relative error control (ATOL(i)=0.0) for a variable which is now vanished");
    }      
    
    if (istate == -3) {
	    error("illegal input detected before taking any integration steps - see written message"); 
      unprotect_all();
  	} else {
    	/*   REAL(YOUT)[(it+1)*(ntot+1)] = tin;*/
      for (j = 0; j < neq; j++)
	    COMPLEX(YOUT)[(it+1)*(ntot + 1) + j + 1] = xytmp[j];
   
	    if (isOut == 1) {
        zderiv_func (&neq, &tin, xytmp, dy, zout, ipar) ;
	      for (j = 0; j < nout; j++)
        COMPLEX(YOUT)[(it+1)*(ntot + 1) + j + neq + 1] = zout[j]; 
      }
    } 

/*                    ####  an error occurred   ####                          */      
    if (istate < 0 || tin < tout) {
	    warning("Returning early from dvode  Results are accurate, as far as they go\n");

    	/* redimension YOUT */
	    PROTECT(YOUT2 = allocMatrix(CPLXSXP,ntot+1,(it+2)));incr_N_Protect();

  	  for (k = 0; k < it+2; k++)
  	    for (j = 0; j < ntot+1; j++)
  	      COMPLEX(YOUT2)[k*(ntot+1) + j] = COMPLEX(YOUT)[k*(ntot+1) + j];
      break;
    }
  }  /* end main time loop */

/*                   ####   returning output   ####                           */    
  terminate(istate, iwork, 23, 0, rwork, 4, 10);      
  
  unlock_solver();
  //unprotect_all();
  restore_N_Protected(old_N_Protect);
  
  if (istate > 0)
    return(YOUT);
  else
    return(YOUT2);
}



