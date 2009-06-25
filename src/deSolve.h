#include <R.h>
#include <Rdefines.h>
/* global variables */
SEXP Time, Y, YPRIME , Rin;
extern SEXP de_gparms;

/* vode globals */
extern SEXP vode_deriv_func;
extern SEXP vode_jac_func;
extern SEXP vode_envir;

/* lsoda globals */
extern SEXP odesolve_deriv_func;
extern SEXP odesolve_jac_func;
extern SEXP odesolve_jac_vec;
extern SEXP odesolve_root_func;
extern SEXP odesolve_envir;

/* daspk globals */
extern SEXP daspk_res_func;
extern SEXP daspk_jac_func;
extern SEXP daspk_psol_func;
extern SEXP daspk_envir;

/* utilities */
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);

/* declarations for initideparms;*/
void Initdeparms(int *, double *);

/* use in daspk */
long int n_eq;
long int mu;
long int ml;
long int nrowpd;
