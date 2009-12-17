#include <R.h>
#include <Rdefines.h>

/*============================================================================
  global R variables 
============================================================================*/
SEXP YOUT, YOUT2, ISTATE, RWORK, IROOT;    /* returned to R */
SEXP Time, Y, YPRIME , Rin;

int    it, n_eq; 
int    *iwork;   
double *rwork;

/* use in daspk */
long int mu;
long int ml;
long int nrowpd;

/* output in DLL globals */
int nout, ntot, isOut, lrpar, lipar, *ipar;
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
int iEvent, nEvent, typeevent, rootevent;

double *timeevent, *valueevent;
int *svarevent, *methodevent;

/*============================================================================
 type definitions for C functions
============================================================================*/
typedef void C_deriv_func_type(int*, double*, double*, double*, double*, int*);
C_deriv_func_type* DLL_deriv_func;

typedef void C_res_func_type(double*, double*, double*, double*, double*,
                             int*, double*, int*);
C_res_func_type* DLL_res_func;


/* this is in compiled code */
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

extern SEXP de_gparms;
SEXP getListElement(SEXP list, const char* str);

/*============================================================================ 
  C- utilities, functions 
============================================================================*/
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);
void returnearly (int);
void terminate(int, int, int, int, int);

/* declarations for initialisations */
void initParms(SEXP Initfunc, SEXP Parms);
void Initdeparms(int*, double*);
void Initdeforc(int*, double*);
void initOut(int isDll, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar);
void initOutdae(int isDll, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar);

/* sparsity of Jacobian */
void sparsity1D(SEXP Type, int* iwork, int neq, int liw);
void sparsity2D(SEXP Type, int* iwork, int neq, int liw);
void sparsity3D(SEXP Type, int* iwork, int neq, int liw);

void initglobals(int);
void initdaeglobals(int);

/* the forcings and event functions */
void updatedeforc(double*);
int initForcings(SEXP list);
int initEvents(SEXP list, SEXP);
void updateevent(double*, double*, int*);
