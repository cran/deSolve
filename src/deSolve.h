#ifndef R_R_H
#  include <R.h>
#endif

#ifndef R_DEFINES_H
#  include <Rdefines.h>
#endif

#ifndef R_INTERNALS_H_
#  include <Rinternals.h>
#endif


/*============================================================================
  global R variables 
============================================================================*/

#ifndef EXTERN
# define EXTERN extern
#endif

EXTERN double *timesteps; /* see also: R_init_deSolve.c */

EXTERN SEXP YOUT, YOUT2, ISTATE, RWORK, IROOT;    /* returned to R */
EXTERN SEXP Y, YPRIME , Rin;

EXTERN int     n_eq; 


/* use in daspk */
EXTERN long int nrowpd;

/* output in DLL globals */
EXTERN int  isOut, *ipar;
EXTERN double *out;

/* forcings  */
EXTERN long int nforc;  /* the number of forcings */
EXTERN double *tvec;
EXTERN double *fvec;
EXTERN int    *ivec;
EXTERN int    fmethod;

EXTERN int    *findex;
EXTERN double *intpol;
EXTERN int    *maxindex;

EXTERN double *forcings;

/* events */
EXTERN double tEvent;
EXTERN int iEvent, nEvent, typeevent, rootevent, Rootsave;
EXTERN double *troot, *valroot;
EXTERN int *nrroot, *termroot;

EXTERN double *timeevent, *valueevent;
EXTERN int *svarevent, *methodevent;

/* time delays */
EXTERN int interpolMethod;  /* for time-delays : 1 = hermite; 2=dense */

/*============================================================================
 type definitions for C functions
============================================================================*/
typedef void C_deriv_func_type(int*, double*, double*, double*, double*, int*);
EXTERN C_deriv_func_type* DLL_deriv_func;

typedef void C_res_func_type(double*, double*, double*, double*, double*,
                             int*, double*, int*);
EXTERN C_res_func_type* DLL_res_func;


/* this is for use in compiled code */
typedef void init_func_type (void (*)(int*, double*));

/*============================================================================
  solver R- global functions 
============================================================================*/
EXTERN SEXP R_deriv_func;
EXTERN SEXP R_jac_func;
EXTERN SEXP R_jac_vec;
EXTERN SEXP R_root_func;
EXTERN SEXP R_event_func;
EXTERN SEXP R_envir;

/* DAE globals */
EXTERN SEXP R_res_func;
EXTERN SEXP R_daejac_func;
EXTERN SEXP R_psol_func;
EXTERN SEXP R_mas_func;

EXTERN SEXP de_gparms;
SEXP getListElement(SEXP list, const char* str);

SEXP getTimestep();

/*============================================================================ 
  C- utilities, functions 
============================================================================*/

void lock_solver(void);
void unlock_solver(void);

void returnearly (int, int, int);
void terminate(int, int*, int, int, double *, int, int);

/* declarations for initialisations */
// void initParms(SEXP Initfunc, SEXP Parms);
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
//void initglobals(int, int);
//void initdaeglobals(int, int);

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

EXTERN int indexhist, indexlag, endreached, starthist;
EXTERN double *histvar, *histdvar, *histtime, *histhh, *histsave;
EXTERN int    *histord;
EXTERN int    histsize, offset;
EXTERN int    initialisehist, lyh, lhh, lo;

#undef EXTERN
