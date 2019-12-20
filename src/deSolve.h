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

extern double *timesteps; /* see also: R_init_deSolve.c */

extern SEXP YOUT, YOUT2, ISTATE, RWORK, IROOT;    /* returned to R */
extern SEXP Y, YPRIME , Rin;

extern int     n_eq; 


/* use in daspk */
extern long int nrowpd;

/* output in DLL globals */
extern int  isOut, *ipar;
extern double *out;

/* forcings  */
extern long int nforc;  /* the number of forcings */
extern double *tvec;
extern double *fvec;
extern int    *ivec;
extern int    fmethod;

extern int    *findex;
extern double *intpol;
extern int    *maxindex;

extern double *forcings;

/* events */
extern double tEvent;
extern int iEvent, nEvent, typeevent, rootevent, Rootsave;
extern double *troot, *valroot;
extern int *nrroot, *termroot;

extern double *timeevent, *valueevent;
extern int *svarevent, *methodevent;

/* time delays */
extern int interpolMethod;  /* for time-delays : 1 = hermite; 2=dense */

/*============================================================================
 type definitions for C functions
============================================================================*/
typedef void C_deriv_func_type(int*, double*, double*, double*, double*, int*);
extern C_deriv_func_type* DLL_deriv_func;

typedef void C_res_func_type(double*, double*, double*, double*, double*,
                             int*, double*, int*);
extern C_res_func_type* DLL_res_func;


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

extern int indexhist, indexlag, endreached, starthist;
extern double *histvar, *histdvar, *histtime, *histhh, *histsave;
extern int    *histord;
extern int    histsize, offset;
extern int    initialisehist, lyh, lhh, lo;

