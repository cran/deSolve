#include <R.h>
#include <Rdefines.h>
/* global variables */

typedef void C_zderiv_func_type (int *, double *, Rcomplex *,Rcomplex *, 
  Rcomplex *, int *);
C_zderiv_func_type *DLL_cderiv_func;  

SEXP cY;

/* livermore solver globals */
extern SEXP cvode_deriv_func;
extern SEXP cvode_jac_func;
extern SEXP vode_envir;

Rcomplex *zout;

void initOutComplex(int isDll, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar) {

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

   zout   = (Rcomplex *) R_alloc(lrpar, sizeof(Rcomplex));
   ipar  = (int *)    R_alloc(lipar, sizeof(int));

   if (isDll ==1)
   {
    ipar[0] = nout;              /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;
    /* other elements of ipar are set in R-function lsodx via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar);j++) ipar[j+3] = INTEGER(Ipar)[j];

    /* first nout elements of rpar reserved for output variables
      other elements are set in R-function zvode via argument *rpar*  */
//     for (j = 0; j < nout; j++)        zout[j] = 0+0i;                
    for (j = 0; j < LENGTH(Rpar);j++) zout[nout+j] = COMPLEX(Rpar)[j];
   }

}

