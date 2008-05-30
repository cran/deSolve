/* Define some global variables and functions that operate on some of them */
/* Patterned on code odesolve_utils.c from package odesolve */
#include <R.h>
#include <Rdefines.h>

/* some functions for keeping track of how many SEXPs 
 * 	are PROTECTed, and UNPROTECTing them in the case of a fortran stop.
 */
 
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

SEXP vode_deriv_func;
SEXP vode_jac_func;
SEXP vode_envir;
SEXP de_gparms;

void Initdeparms(int *N, double *parms)
{
  int i, Nparms;

  Nparms = LENGTH(de_gparms);
  if ((*N) != Nparms)
    {
      PROBLEM "Confusion over the length of parms"
      ERROR;
    } 
  else
    {
      for (i = 0; i < *N; i++) parms[i] = REAL(de_gparms)[i];
    }
}
  
