/* Define some global variables and functions that operate on some of them */

#include <R.h>
#include <Rdefines.h>

#ifndef  R_INTERNALS_H_
#include <Rinternals.h>
#endif

#include <R_ext/Applic.h>
#include <R_ext/Rdynload.h>
#include "deSolve.h"

/*==================================================
some functions for keeping track of how many SEXPs
are PROTECTed, and UNPROTECTing them in the
case of a FORTRAN stop.
==================================================*/
 
long int N_Protected = 0; /* initialize this with zero at the first time */

int solver_locked = 0;   /* prevent nested calls of odepack solvers */

void init_N_Protect(void) { N_Protected = 0; }

void incr_N_Protect(void) { N_Protected++; }

void unprotect_all(void) { UNPROTECT((int) N_Protected); }

long int save_N_Protected(void) {
  int saved_N = N_Protected;
  init_N_Protect();
  return saved_N; 
}

void restore_N_Protected(long int n) {
  unprotect_all();
  N_Protected = n; 
}

void my_unprotect(int n) {
    UNPROTECT(n);
    N_Protected -= n;
}

void lock_solver(void) {
  if (solver_locked) {
    /* important: unlock for the next call *after* error */
    solver_locked = 0; 
    error("The used combination of solvers cannot be nested.\n");
  }
  solver_locked = 1;
}

void unlock_solver(void) {
  solver_locked = 0;
  timesteps[0] = 0;
  timesteps[1] = 0;
}


/* Globals :*/
SEXP R_deriv_func;
SEXP R_jac_func;
SEXP R_jac_vec;
SEXP R_root_func;
SEXP R_event_func;

SEXP R_envir;
SEXP odesolve_gparms;

SEXP R_res_func;
SEXP R_daejac_func;
SEXP R_psol_func;

SEXP R_mas_func;

SEXP de_gparms;

/*======================================================
SEXP initialisation functions
=======================================================*/

void initglobals(int nt, int ntot) {
/*  PROTECT(Time = NEW_NUMERIC(1));                  incr_N_Protect(); */
  PROTECT(Y = allocVector(REALSXP,(n_eq)));        incr_N_Protect();
  PROTECT(YOUT = allocMatrix(REALSXP,ntot+1,nt));  incr_N_Protect();
}

void initdaeglobals(int nt, int ntot) {
/*  PROTECT(Time = NEW_NUMERIC(1));                    incr_N_Protect(); */
  PROTECT(Rin  = NEW_NUMERIC(2));                    incr_N_Protect();
  PROTECT(Y = allocVector(REALSXP,n_eq));            incr_N_Protect();
  PROTECT(YPRIME = allocVector(REALSXP,n_eq));       incr_N_Protect();
  PROTECT(YOUT = allocMatrix(REALSXP,ntot+1,nt));    incr_N_Protect();
}

/*======================================================
Parameter initialisation functions
note: forcing initialisation function is in forcings.c
=======================================================*/

void initParms(SEXP Initfunc, SEXP Parms) {
  if (Initfunc == NA_STRING) return;
  if (inherits(Initfunc, "NativeSymbol")) {
    init_func_type *initializer;
    PROTECT(de_gparms = Parms); incr_N_Protect();
    initializer = (init_func_type *) R_ExternalPtrAddr(Initfunc);
    initializer(Initdeparms);
  }
}


void Initdeparms(int *N, double *parms) {
  int i, Nparms;
  Nparms = LENGTH(de_gparms);
  if ((*N) != Nparms) {
    warning("Number of parameters passed to solver, %i; number in DLL, %i\n",
      Nparms, *N);
    PROBLEM "Confusion over the length of parms"
    ERROR;
  } else {
    for (i = 0; i < *N; i++) parms[i] = REAL(de_gparms)[i];
  }
}
  
SEXP get_deSolve_gparms(void) {
  return de_gparms;
}

/*=========================================================================== 
  C-equivalent of R-function timestep: gets the past and new time step
  =========================================================================== */
SEXP getTimestep() {
  SEXP value;
  PROTECT(value = NEW_NUMERIC(2));
  if (timesteps == NULL) {         /* integration not yet started... */
    for (int i = 0; i < 2; i++) 
      NUMERIC_POINTER(value)[i] = 0.0;
  } else
    for (int i = 0; i < 2; i++) 
      NUMERIC_POINTER(value)[i] = timesteps[i];
  UNPROTECT(1);
  return(value);
}
  
/*============================ ======================
 Termination 
===================================================*/

/* an error occurred - save output in YOUT2 */
void returnearly (int Print, int it, int ntot) {
  int j, k;
  if (Print) 
    warning("Returning early. Results are accurate, as far as they go\n");
  PROTECT(YOUT2 = allocMatrix(REALSXP,ntot+1,(it+2))); incr_N_Protect();
  for (k = 0; k < it+2; k++)
    for (j = 0; j < ntot+1; j++)
      REAL(YOUT2)[k*(ntot+1) + j] = REAL(YOUT)[k*(ntot+1) + j];
}   

/* add ISTATE and RSTATE */
void terminate(int istate, int * iwork, int ilen, int ioffset, 
  double * rwork, int rlen, int roffset) {

  int k;
  
  PROTECT(ISTATE = allocVector(INTSXP, ilen)); incr_N_Protect();
  for (k = 0; k < ilen-1; k++) INTEGER(ISTATE)[k+1] = iwork[k +ioffset];
  INTEGER(ISTATE)[0] = istate;  
        
  PROTECT(RWORK = allocVector(REALSXP, rlen)); incr_N_Protect();
  for (k = 0; k < rlen; k++) REAL(RWORK)[k] = rwork[k+roffset];
  if (istate > 0) {
    setAttrib(YOUT, install("istate"), ISTATE);
    setAttrib(YOUT, install("rstate"), RWORK);
  }
  else  {
    setAttrib(YOUT2, install("istate"), ISTATE);
    setAttrib(YOUT2, install("rstate"), RWORK);
  }
  /* timestep = 0 - for use in getTimestep */
  timesteps[0] = 0;
  timesteps[1] = 0;
}

/*==================================================
 extracting elements from a list
===================================================*/

SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  for (i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

/*==================================================
 output initialisation function

 out and ipar are used to pass output variables
 (number set by nout) followed by other input
 by R-arguments rpar, ipar
 ipar[0]: number of output variables,
 ipar[1]: length of rpar,
 ipar[2]: length of ipar
===================================================*/

/* Initialise output - output variables calculated in R-code ... */

void initOutR(int isDll, int *nout, int *ntot, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar) {

  int j, lrpar, lipar;
  *nout = INTEGER(nOut)[0];       /* number of output variables */
  if (isDll) {                   /* function is a dll */
    if (*nout > 0) isOut = 1;
    *ntot  = neq + *nout;          /* length of yout */
    lrpar = *nout + LENGTH(Rpar); /* length of rpar; LENGTH(Rpar) is always >0 */
    lipar = 3 + LENGTH(Ipar);    /* length of ipar */
  } else {                       /* function is not a dll */
    isOut = 0;
    *ntot = neq;
    lipar = 1;
    lrpar = 1;
  }
  out  = (double*) R_alloc(lrpar, sizeof(double));
  ipar = (int*)    R_alloc(lipar, sizeof(int));
  if (isDll ==1) {
    ipar[0] = *nout;              /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;

    /* other elements of ipar are set in R-function lsodx via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar); j++) ipar[j+3] = INTEGER(Ipar)[j];

    /* first nout elements of rpar reserved for output variables
       other elements are set in R-function lsodx via argument *rpar* */
    for (j = 0; j < *nout; j++)        out[j] = 0.;
    for (j = 0; j < LENGTH(Rpar); j++) out[*nout+j] = REAL(Rpar)[j];
   }
}

/* Initialise output - output variables calculated in C-code ... */

void initOutC(int isDll, int *nout, int *ntot, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar) {
  int j, lrpar, lipar;
  /* initialise output when a dae ... */   
  /*  output always done here in C-code (<-> lsode, vode)... */

  *nout  = INTEGER(nOut)[0];
  *ntot  = n_eq+*nout;
  
  if (isDll == 1) {                /* function is a dll */
    lrpar = *nout + LENGTH(Rpar);   /* length of rpar */
    lipar = 3    + LENGTH(Ipar);   /* length of ipar */
  } else {                         /* function is not a dll */
    lipar = 3;
    lrpar = *nout;
  }
  out   = (double*) R_alloc(lrpar, sizeof(double));
  ipar  = (int*)    R_alloc(lipar, sizeof(int));

  if (isDll == 1) {
    ipar[0] = *nout;                /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;

    /* other elements of ipar are set in R-function lsodx via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar); j++) ipar[j+3] = INTEGER(Ipar)[j];

    /* first nout elements of rpar reserved for output variables
       other elements are set in R-function lsodx via argument *rpar* */
    for (j = 0; j < *nout;         j++) out[j] = 0.;
    for (j = 0; j < LENGTH(Rpar); j++) out[*nout+j] = REAL(Rpar)[j];
   }
}
/*==================================================
 1-D, 2-D and 3-D sparsity structure
================================================== */
void sparsity1D (SEXP Type, int* iwork, int neq, int liw) {

    int nspec, nx, ij, i, j, k, l;
    nspec = INTEGER(Type)[1]; /* number of components*/
    nx    = INTEGER(Type)[2]; /* dimension x*/

    ij    = 31 + neq;
    iwork[30] = 1;
    k = 1;
    for( i = 0; i < nspec; i++) {
      for( j = 0; j < nx; j++) {
        if (ij > liw-3-nspec)  error ("not enough memory allocated in iwork - increase liw %i ",liw);
        iwork[ij++] = k;
        if (j < nx-1) iwork[ij++] = k+1 ;
        if (j > 0)    iwork[ij++] = k-1 ;

        for(l = 0; l < nspec; l++)
          if (l != i) iwork[ij++] = l*nx+j+1;

        iwork[30+k] = ij-30-neq;
        k = k+1;
      }
    }
    iwork[ij] = 0;
}

/*==================================================*/

void sparsity2D (SEXP Type, int* iwork, int neq, int liw) {

    int nspec, nx, ny, bndx, bndy, Nt, ij, isp, i, j, k, l, m;

    nspec = INTEGER(Type)[1]; /* number components*/
    nx    = INTEGER(Type)[2]; /* dimension x*/
    ny    = INTEGER(Type)[3]; /* dimension y*/
    bndx  = INTEGER(Type)[4]; /* cyclic boundary x*/
    bndy  = INTEGER(Type)[5]; /* cyclic boundary y*/
    Nt    = nx*ny;
    ij    = 31 + neq;
    iwork[30] = 1;
    m = 1;
    for( i = 0; i < nspec; i++) {
      isp = i*Nt;
      for( j = 0; j < nx; j++) {
        for( k = 0; k < ny; k++) {
          if (ij > liw-8-nspec)  
            error("not enough memory allocated in iwork - increase liw %i ",liw);
          iwork[ij++] = m;
          if (k < ny-1)     iwork[ij++] = m+1;
          if (j < nx-1)     iwork[ij++] = m+ny;
          if (j > 0)        iwork[ij++] = m-ny;
          if (k > 0)        iwork[ij++] = m-1;
          if (bndx == 1) {
            if (j == 0)     iwork[ij++] = isp+(nx-1)*ny+k+1;
            if (j == nx-1)  iwork[ij++] = isp+k+1;
          }
          if (bndy == 1) {
            if (k == 0)    iwork[ij++] = isp+(j+1)*ny;
            if (k == ny-1) iwork[ij++] = isp + j*ny +1;
          }
          for(l = 0; l < nspec; l++)
            if (l != i)    iwork[ij++] = l*Nt+j*ny+k+1;

          iwork[30+m] = ij-30-neq;
          m = m+1;
        }
      }
    }
}
void interact (int *ij, int nnz, int *iwork, int is, int ival) {
  int i, isave;

  isave = 1;
/* check if not yet present for current state */
	for (i = is; i < *ij; i++) 
	  if (iwork[i] == ival) {
	     isave = 0;
	     break;
	  }

/* save */
	if (isave == 1) {
    if (*ij > nnz) 
       error ("not enough memory allocated in iwork - increase liw %i ", nnz);
     iwork[(*ij)++] = ival;
  }
}

/*==================================================*/
/* an element in C-array A(I,J,K), i=0,dim(1)-1 etc... is positioned at 
   j*dim(2)*dim(3) + k*dim(3) + l + 1 in FORTRAN VECTOR! 
   includes check on validity

   dimens and boundary are reversed ... 
*/

void sparsity3D (SEXP Type, int* iwork, int neq, int liw) {
    int nspec, nx, ny, nz, bndx, bndy, bndz, Nt, ij, is, isp, i, j, k, l, m, ll;

    nspec = INTEGER(Type)[1]; 
    nx    = INTEGER(Type)[2]; 
    ny    = INTEGER(Type)[3]; 
    nz    = INTEGER(Type)[4]; 
    bndx  = INTEGER(Type)[5];
    bndy  = INTEGER(Type)[6]; 
    bndz  = INTEGER(Type)[7]; 

    Nt    = nx*ny*nz;
    ij    = 31+neq;
    iwork[30] = 1;
    m = 1;
    for( i = 0; i < nspec; i++) {
      isp = i*Nt;
      for( j = 0; j < nx; j++) {
        for( k = 0; k < ny; k++) {
          for( ll = 0; ll < nz; ll++) {
            is = ij;
            
            if (ij > liw-6-nspec)  
              error ("not enough memory allocated in iwork - increase liw %i ", liw);
                            interact (&ij, liw, iwork, is, m);
            if (ll < nz-1)  
              interact (&ij, liw, iwork, is, m+1);
            else if (bndz == 1)
              interact (&ij, liw, iwork, is, isp + j*ny*nz + k*nz + 1);              
            
            if (k  < ny-1) 
              interact (&ij, liw, iwork, is, m+nz);
            else  if (bndy == 1)
              interact (&ij, liw, iwork, is, isp + j*ny*nz + ll + 1);
              
            if (j  < nx-1)  
              interact (&ij, liw, iwork, is, m+ny*nz);
            else if (bndx == 1)  
              interact (&ij, liw, iwork, is, isp + k*nz + ll + 1);

            if (j > 0)      
              interact (&ij, liw, iwork, is, m-ny*nz);
            else if (bndx == 1)
              interact (&ij, liw, iwork, is, isp+(nx-1)*ny*nz+k*nz+ll+1);
                                         
            if (k > 0)  
              interact (&ij, liw, iwork, is, m-nz);
            else  if (bndy == 1)
              interact (&ij, liw, iwork, is, isp + j*ny*nz+(ny-1)*nz+ll+1);              
            
            if (ll > 0) 
              interact (&ij, liw, iwork, is, m-1);
            else if (bndz == 1)  
              interact (&ij, liw, iwork, is, isp + j*ny*nz+k*nz+nz);

            for(l = 0; l < nspec; l++)
              if (l != i) interact (&ij, liw, iwork, is, l*Nt+j*ny*nz+k*nz+ll+1);
            iwork[30+m] = ij-30-neq;
            m = m+1;
          }
        }
      }
    }
}

