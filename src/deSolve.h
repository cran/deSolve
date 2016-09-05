#include <R.h>
#include <Rdefines.h>

#ifndef  R_INTERNALS_H_
#include <Rinternals.h>
#endif



/*============================================================================
  global R variables 
============================================================================*/

double *timesteps; /* see also: R_init_deSolve.c */

SEXP YOUT, YOUT2, ISTATE, RWORK, IROOT;    /* returned to R */
SEXP Y, YPRIME , Rin;

int     n_eq; 


/* use in daspk */
long int nrowpd;

/* output in DLL globals */
int  isOut, *ipar;
double *out;

/* forcings  */
long int nforc;  /* the number of forcings */
double *tvec;
double *fvec;
int    *ivec;
int    fmethod;

int    *findex;
double *intpol;
int    *maxindex;

double *forcings;

/* events */
double tEvent;
int iEvent, nEvent, typeevent, rootevent, Rootsave;
double *troot, *valroot;
int *nrroot, *termroot;

double *timeevent, *valueevent;
int *svarevent, *methodevent;

/* time delays */
int interpolMethod;  /* for time-delays : 1 = hermite; 2=dense */

/*============================================================================
 type definitions for C functions
============================================================================*/
typedef void C_deriv_func_type(int*, double*, double*, double*, double*, int*);
C_deriv_func_type* DLL_deriv_func;

typedef void C_res_func_type(double*, double*, double*, double*, double*,
                             int*, double*, int*);
C_res_func_type* DLL_res_func;


/* this is for use in compiled code */
typedef void init_func_type (void (*)(int*, double*));

/*============================================================================
  solver R- global functions 
============================================================================*/
extern SEXP R_deriv_func;
extern SEXP R_jac_func;
extern SEXP R_jac_vec;
extern SEXP R_root_func;
extern SEXP R_event_func;
extern SEXP R_envir;

/* DAE globals */
extern SEXP R_res_func;
extern SEXP R_daejac_func;
extern SEXP R_psol_func;
extern SEXP R_mas_func;

extern SEXP de_gparms;
SEXP getListElement(SEXP list, const char* str);

SEXP getTimestep();

/*============================================================================ 
  C- utilities, functions 
============================================================================*/
void init_N_Protect(void);
void incr_N_Protect(void);
long int save_N_Protected(void);
void restore_N_Protected(long int);
void unprotect_all(void);
void my_unprotect(int);

void lock_solver(void);
void unlock_solver(void);

void returnearly (int, int, int);
void terminate(int, int*, int, int, double *, int, int);

/* declarations for initialisations */
void initParms(SEXP Initfunc, SEXP Parms);
void Initdeparms(int*, double*);
void Initdeforc(int*, double*);
void initOutR(int isDll, int *nout, int *ntot, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar);
void initOutC(int isDll, int *nout, int *ntot, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar);

/* sparsity of Jacobian */
void sparsity1D(SEXP Type, int* iwork, int neq, int liw);
void sparsity2D(SEXP Type, int* iwork, int neq, int liw);
void sparsity3D(SEXP Type, int* iwork, int neq, int liw);
void sparsity2Dmap(SEXP Type, int* iwork, int neq, int liw);  /* testing, since version 1.10.4*/
void sparsity3Dmap(SEXP Type, int* iwork, int neq, int liw);  /* testing, since version 1.10.4*/
void interactmap (int *ij, int nnz, int *iwork, int *ipres, int ival);
void initglobals(int, int);
void initdaeglobals(int, int);

/* the forcings and event functions */
void updatedeforc(double*);
int initForcings(SEXP list);
int initEvents(SEXP list, SEXP, int);
void updateevent(double*, double*, int*);


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                         DECLARATIONS for time lags
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*==========================================
  R-functions
==========================================*/

SEXP getPastValue   (SEXP T, SEXP nr);
SEXP getPastGradient(SEXP T, SEXP nr);

/*==========================================
  C- utilities, functions
==========================================*/
/* Hermitian interpolation */
double Hermite (double t0, double t1, double y0, double y1, double dy0,
                double dy1, double t);

double dHermite(double t0, double t1, double y0, double y1, double dy0,
                double dy1, double t);

int initLags(SEXP elag, int solver, int nroot);

/* history vectors  */
void inithist(int max, int maxlags, int solver, int nroot);

void updatehistini(double t, double *y, double *dY, double *rwork, int *iwork);
void updatehist(double t, double *y, double *dy, double *rwork, int *iwork);

int nexthist(int i);
double interpolate(int i, int k, double t0, double t1, double t, 
  double *Yh, int nq); 


/*==========================================
  Global variables for history arrays
==========================================*/

int indexhist, indexlag, endreached, starthist;
double *histvar, *histdvar, *histtime, *histhh, *histsave;
int    *histord;
int    histsize, offset;
int    initialisehist, lyh, lhh, lo;

