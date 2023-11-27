#define USE_FC_LEN_T
#include <R.h>

#ifdef FC_LEN_T
# include <stddef.h> // for size_t if needed
# define FCLEN ,FC_LEN_T msg_len
#else
# define FCLEN
#endif


void F77_SUB(rprintf)(const char* msg FCLEN) {
   Rprintf("%s", msg);
   Rprintf("\n");
}

void F77_SUB(rprintfid)(const char* msg, int *i, double *d FCLEN) {
   Rprintf(msg, *i, *d);
   Rprintf("\n");
}

void F77_SUB(rprintfdi)(const char* msg, double *d, int *i FCLEN) {
   Rprintf(msg, *d, *i);
   Rprintf("\n");
}

void F77_SUB(rprintfdid)(const char* msg, double *d1, int *i, double *d2 FCLEN) {
   Rprintf(msg, *d1, *i, *d2);
   Rprintf("\n");
}

void F77_SUB(rprintfd1)(const char* msg, double *d FCLEN) {
   Rprintf(msg, *d);
   Rprintf("\n");
}

void F77_SUB(rprintfd2)(const char* msg, double *d1, double *d2 FCLEN) {
   Rprintf(msg, *d1, *d2);
   Rprintf("\n");
}

void F77_SUB(rprintfi1)(const char* msg, int *i FCLEN) {
   Rprintf(msg, *i);
   Rprintf("\n");
}

void F77_SUB(rprintfi2)(const char* msg, int *i1, int *i2 FCLEN) {
   Rprintf(msg, *i1, *i2);
   Rprintf("\n");
}

void F77_SUB(rprintfi3)(const char* msg, int *i1, int *i2, int* i3 FCLEN) {
   Rprintf(msg, *i1, *i2, *i3);
   Rprintf("\n");
}
