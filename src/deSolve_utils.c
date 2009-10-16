/* Define some global variables and functions that operate on some of them */

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#include "deSolve.h"

/*==================================================
some functions for keeping track of how many SEXPs
are PROTECTed, and UNPROTECTing them in the
case of a FORTRAN stop.
==================================================*/
 
long int N_Protected;

void init_N_Protect(void) { N_Protected = 0; }

void incr_N_Protect(void) { N_Protected++; }

void unprotect_all(void) { UNPROTECT((int) N_Protected); }

void my_unprotect(int n)
{
    UNPROTECT(n);
    N_Protected -= n;
}

/* Globals :*/
SEXP odesolve_deriv_func;
SEXP odesolve_jac_func;
SEXP odesolve_jac_vec;
SEXP odesolve_root_func;
SEXP odesolve_envir;
SEXP odesolve_gparms;

SEXP daspk_res_func;
SEXP daspk_jac_func;
SEXP daspk_psol_func;
SEXP daspk_envir;

SEXP de_gparms;

/*======================================================
Parameter initialisation functions
note: forcing initialisation function is in forcings.c
=======================================================*/

void initParms(SEXP Initfunc, SEXP Parms) {
  // ks: added this to prevent entering this if initfunc does not exist
  if (Initfunc == NA_STRING) return;
  if (inherits(Initfunc, "NativeSymbol"))  {
    init_func *initializer;

    PROTECT(de_gparms = Parms);     incr_N_Protect();
    initializer = (init_func *) R_ExternalPtrAddr(Initfunc);
    initializer(Initdeparms);
  }

}


void Initdeparms(int *N, double *parms)
{
  int i, Nparms;

  Nparms = LENGTH(de_gparms);
  if ((*N) != Nparms)
    {
      warning("Number of parameters passed to solver, %i; number in DLL, %i\n",Nparms, *N);
      PROBLEM "Confusion over the length of parms"

      ERROR;
    } 
  else
    {
      for (i = 0; i < *N; i++) parms[i] = REAL(de_gparms)[i];
    }
}
  
SEXP get_deSolve_gparms(void)
{
  return de_gparms;
}

/*==================================================
 extracting elements from a list
===================================================*/

SEXP getListElement(SEXP list, const char *str)
{
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

void initOut(int isDll, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar) {

  int j;
  nout   = INTEGER(nOut)[0];    /* number of output variables */
  if (isDll)                    /* function is a dll */
  {
   if (nout > 0) isOut = 1;
   ntot  = neq + nout;          /* length of yout */
   lrpar = nout + LENGTH(Rpar); /* length of rpar; LENGTH(Rpar) is always >0 */
   lipar = 3 + LENGTH(Ipar);    /* length of ipar */

  } else                              /* function is not a dll */
  {
   isOut = 0;
   ntot = neq;
   lipar = 1;
   lrpar = 1;
  }

   out   = (double *) R_alloc(lrpar, sizeof(double));
   ipar  = (int *)    R_alloc(lipar, sizeof(int));

   if (isDll ==1)
   {
    ipar[0] = nout;              /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;
    /* other elements of ipar are set in R-function lsodx via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar);j++) ipar[j+3] = INTEGER(Ipar)[j];

    /* first nout elements of rpar reserved for output variables
      other elements are set in R-function lsodx via argument *rpar* */
    for (j = 0; j < nout; j++)        out[j] = 0.;
    for (j = 0; j < LENGTH(Rpar);j++) out[nout+j] = REAL(Rpar)[j];
   }

}

/*==================================================
 1-D, 2-D and 3-D sparsity structure
================================================== */
void sparsity1D (SEXP Type, int* iwork, int neq, int liw) {

    int nspec, nx, ij, i, j, k, l;
    nspec = INTEGER(Type)[1]; /* number of components*/
    nx    = INTEGER(Type)[2]; /* dimension x*/

    ij     = 31+neq;
    iwork[30] = 1;
    k = 1;
    for( i = 0; i<nspec; i++) {
      for( j = 0; j<nx; j++) {
        if (ij > liw-4)  error ("not enough memory allocated in iwork - increase liw %i ",liw);
                    iwork[ij++] = k;
        if (j<nx-1) iwork[ij++] = k+1 ;
        if (j>0)    iwork[ij++] = k-1 ;

        for(l = 0; l<nspec;l++)
          if (l != i) iwork[ij++] = l*nx+j+1;

        iwork[30+k] = ij-30-neq;
        k=k+1;
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
    ij     = 31+neq;
    iwork[30] = 1;
    m = 1;
    for( i = 0; i<nspec; i++) {
       isp = i*Nt;
       for( j = 0; j<nx; j++) {
         for( k = 0; k<ny; k++) {
           if (ij > liw-4)  error ("not enough memory allocated in iwork - increase liw %i ",liw);
                              iwork[ij++] = m;
           if (k<ny-1)        iwork[ij++] = m+1;

           if (j<nx-1)        iwork[ij++] = m+ny;
           if (j >0)          iwork[ij++] = m-ny;
           if (k >0)          iwork[ij++] = m-1;
           if (bndx == 1) {
               if (j == 0)      iwork[ij++] = isp+(nx-1)*ny+k+1;
               if (j == nx-1)   iwork[ij++] = isp+k+1;
           }
           if (bndy == 1) {
               if (k == 0)      iwork[ij++] = isp+(j+1)*ny;
               if (k == ny-1)   iwork[ij++] = isp + j*ny +1;
           }
           for(l = 0; l<nspec;l++)
               if (l != i)      iwork[ij++] = l*Nt+j*ny+k+1;

           iwork[30+m] = ij-30-neq;
           m = m+1;
           }
        }
    }
}

/*==================================================*/

void sparsity3D (SEXP Type, int* iwork, int neq, int liw) {

    int nspec, nx, ny, nz,  Nt, ij, isp, i, j, k, l, m, ll;

    nspec = INTEGER(Type)[1]; /* number components*/
    nx    = INTEGER(Type)[2]; /* dimension x*/
    ny    = INTEGER(Type)[3]; /* dimension y*/
    nz    = INTEGER(Type)[4]; /* dimension y*/
/*     bndx  = INTEGER(Type)[5];
       bndy  = INTEGER(Type)[6];  cyclic boundary NOT yet implemented*/
    Nt    = nx*ny*nz;
    ij     = 31+neq;
    iwork[30] = 1;
    m = 1;
    for( i = 0; i<nspec; i++) {
      isp = i*Nt;
      for( j = 0; j<nx; j++) {
        for( k = 0; k<ny; k++) {
           for( ll = 0; ll<nz; ll++) {
              if (ij > liw-4)  error ("not enough memory allocated in iwork - increase liw %i ",liw);
                                 iwork[ij++] = m;
              if (ll<nz-1)       iwork[ij++] = m+1;
              if (k<ny-1)        iwork[ij++] = m+nz;
              if (j<nx-1)        iwork[ij++] = m+ny*nz;

              if (j >0)          iwork[ij++] = m-ny*nz;
              if (k >0)          iwork[ij++] = m-nz;
              if (ll >0)         iwork[ij++] = m-1;
              for(l = 0; l<nspec;l++)
                if (l != i)       iwork[ij++] = l*Nt+j*ny*nz+k*nz+ll+1;

              iwork[30+m] = ij-30-neq;
              m = m+1;
            }
         }
      }
    }
}
