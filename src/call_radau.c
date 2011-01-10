#include <time.h>
#include <string.h>
#include "deSolve.h"

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   RADAU: Implicit runge-Kutta of order 5
   due to Hairer and Wanner, with stepsize control and dense output
   
   The C-wrappers that provide the interface between FORTRAN codes and R-code 
   are: C_deriv_func_rad: interface with R-code "func", passes derivatives  
        C_deriv_out_rad : interface with R-code "func", passes derivatives + output variables  
  
   C_deriv_func_forc_rad provides the interface between the function specified in
   a DLL and the integrator, in case there are forcing functions.
   
   karline soetaert
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
  int maxt, it, nout, isDll, ntot;   
  C_deriv_func_type *deriv_func;
  double *xdytmp, *ytmp, *tt;
  
/* definition of the calls to the FORTRAN subroutines in file radau.f */

void F77_NAME(radau5)( int *,
         void (*)(int *, double *, double *, double *,
                              double *, int *),           // func
		     double *, double *, double *, double *,
		     double *, double *, int *,  
 	       void (*)(int *, double *, double *, int *, int *, 
                    double *, int *, double *, int *),    // jac
		     int *, int *, int *,
 	       void (*)(int *, double *, int *,
                              double *, int *),           // mas
		     int *, int *, int *,
         void (*)(int *, double *, double *, double *, double *,  
			            int *, int *, double *, int *, int *, double *),   // soloutrad
		     int *, double *, int *, int *, int*, double *, int*, int*);

/* continuous output formula */
void F77_NAME (contr5) (int *, double *, double *, int *, double *);

/* wrapper above the derivate function in a dll that first estimates the
values of the forcing functions */

static void C_deriv_func_forc_rad (int *neq, double *t, double *y,
                         double *ydot, double *yout, int *iout)
{
  updatedeforc(t);
  DLL_deriv_func(neq, t, y, ydot, yout, iout);
}

/* interface between FORTRAN function call and R function
   Fortran code calls C_deriv_func_rad(N, t, y, ydot, yout, iout) 
   R code called as R_deriv_func(time, y) and returns ydot */

static void C_deriv_func_rad (int *neq, double *t, double *y,
                          double *ydot, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;

                              REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(R_deriv_func,Time,Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));           incr_N_Protect();

  for (i = 0; i < *neq; i++)   ydot[i] = REAL(ans)[i];

  my_unprotect(2);
}

static void C_mas_func_rad (int *neq, double *am, int *lmas,
                             double *yout, int *iout)
{
  int i;
  SEXP NEQ, LM, R_fcall, ans;

  PROTECT(NEQ = NEW_INTEGER(1));                  incr_N_Protect();
  PROTECT(LM = NEW_INTEGER(1));                   incr_N_Protect();

                              INTEGER(NEQ)[0] = *neq;
                              INTEGER(LM) [0] = *lmas;
  PROTECT(R_fcall = lang3(R_mas_func,NEQ,LM));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));           incr_N_Protect();

  for (i = 0; i <*lmas * *neq; i++)   am[i] = REAL(ans)[i];

  my_unprotect(4);
}

/* deriv output function - for ordinary output variables */

static void C_deriv_out_rad (int *nOut, double *t, double *y, 
                       double *ydot, double *yout)
{
  int i;
  SEXP R_fcall, ans;
  
  REAL(Time)[0] = *t;
  for (i = 0; i < n_eq; i++)  
      REAL(Y)[i] = y[i];
     
  PROTECT(R_fcall = lang3(R_deriv_func,Time, Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));            incr_N_Protect();

  for (i = 0; i < *nOut; i++) yout[i] = REAL(ans)[i + n_eq];

  my_unprotect(2);                                  
}      

/* save output in R-variables */

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

/* function called by Fortran to check for output */

static void C_soloutrad(int * nr, double * told, double *t, double * y, 
  double * con, int *lrc, int * neq, double * rpar, int * ipar, int * irtrn, double *xout)
{
  timesteps[0] = *told-*t;
  timesteps[1] = *told-*t;

  while (*told <= tt[it] && tt[it] < *t) {
    F77_CALL(contr5) (neq, &tt[it], con, lrc, ytmp);
    saveOut(tt[it], ytmp);	     
    it++;
    if ( it >= maxt) break;
  }
}

/* interface to jacobian function */
static void C_jac_func_rad(int *neq, double *t, double *y, int *ml,
		    int *mu, double *pd, int *nrowpd, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;

                             REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++) REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(R_jac_func,Time,Y));    incr_N_Protect();
  PROTECT(ans = eval(R_fcall, R_envir));          incr_N_Protect();

  for (i = 0; i < *neq * *nrowpd; i++)  pd[i] = REAL(ans)[i];

  my_unprotect(2);
}

/* give name to data types */
typedef void C_solout_type (int *, double *, double *, double *,
  double *, int *, int *, double *, int *, int *, double *) ;

typedef void C_mas_type (int *, double *, int *, double *, int *);

// to be changed...
typedef void C_jac_func_type_rad(int *, double *, double *, int *, 
                     int *, double *, int*, double *, int *);

/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_radau(SEXP y, SEXP times, SEXP derivfunc, SEXP masfunc, SEXP jacfunc,
    SEXP parms, SEXP rtol, SEXP atol, 
    SEXP Nrjac, SEXP Nrmas,
		SEXP rho, SEXP initfunc, SEXP verbose, SEXP rWork, SEXP iWork, 
    SEXP nOut, SEXP lRw, SEXP lIw, 
    SEXP Rpar, SEXP Ipar, SEXP Hini, SEXP flist, SEXP Type)

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

  int  j, nt, latol, lrtol, lrw, liw, type, 
       ijac, mljac, mujac, imas, mlmas, mumas;
  int  isForcing;
  double *xytmp, tin, tout, *Atol, *Rtol, hini=0;
  int itol, iout, idid, mflag;
  int    *iwork;   
  double *rwork;
  

  /* pointers to functions passed to FORTRAN */
  C_solout_type         *solout = NULL;
  C_jac_func_type_rad   *jac_func = NULL;
  C_mas_type            *mas_func = NULL;
    
/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/
  lock_solver(); /* prevent nested call of solvers that have global variables */
/*                      #### initialisation ####                              */    
  //init_N_Protect();
  long int old_N_Protect = save_N_Protected();  

  n_eq = LENGTH(y);             /* number of equations */ 
  nt   = LENGTH(times);         /* number of output times */
  maxt = nt; 
  
  tt = (double *) R_alloc(nt, sizeof(double));
  for (j = 0; j < nt; j++) tt[j] = REAL(times)[j];
  
  mflag = INTEGER(verbose)[0];
  type  = INTEGER(Type)[0];     /* 1 = dopri 8, 2 = dopri5 */

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
  initglobals (nt, ntot);
  //timesteps = (double *) R_alloc(2, sizeof(double));
  for (j=0; j<2; j++) timesteps[j] = 0.;
  
  /* Initialization of Parameters, Forcings (DLL) */
  initParms (initfunc, parms);
  isForcing = initForcings(flist);
  
  if (nout > 0 ) {
     xdytmp= (double *) R_alloc(n_eq, sizeof(double));
     for (j = 0; j < n_eq; j++) xdytmp[j] = 0.; 
  }
  
 /* pointers to functions deriv_func, jac_func, passed to FORTRAN */
  if (isDll)  { /* DLL address passed to FORTRAN */
      deriv_func = (C_deriv_func_type *) R_ExternalPtrAddr(derivfunc);  
	  
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
      R_envir = rho;
  }

  if (!isNull(jacfunc))   {
      if (isDll)
	      jac_func = (C_jac_func_type_rad *) R_ExternalPtrAddr(jacfunc);
	    else  {
	      R_jac_func = jacfunc;
	      jac_func= C_jac_func_rad;
	    }
    }
  if (!isNull(masfunc))   {
	      R_mas_func = masfunc;
	      mas_func= C_mas_func_rad;
    }


 	solout = C_soloutrad;              
 
  iout = 2;                           /* solout called after each step OR 1???*/
  idid = 0;

/*                   ####      integration     ####                           */    
  it   = 0;
  tin  = REAL(times)[0];
  tout = REAL(times)[nt-1];
  saveOut (tin, xytmp);               /* save initial condition */ 
  
    F77_CALL(radau5) ( &n_eq, deriv_func, &tin, xytmp, &tout, &hini, 
		     Rtol, Atol, &itol, jac_func, &ijac, &mljac, &mujac, 
         mas_func, &imas, &mlmas, &mumas, solout, &iout,
		     rwork, &lrw, iwork, &liw, out, ipar, &idid);
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
  if (idid < 0 ) {
    it = it-1;
    returnearly (1, it, ntot);
  }  

/*                   ####   returning output   ####                           */    
  rwork[0] = hini;
  rwork[1] = tin ; 
  terminate(idid,iwork,7,13,rwork,5,0);       
  
/*                   ####     termination      ####                           */    
  unlock_solver();
  restore_N_Protected(old_N_Protect);                           
  //unprotect_all();
  if (idid > 0)
    return(YOUT);
  else
    return(YOUT2);
}
 
