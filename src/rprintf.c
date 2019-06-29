#include <R.h>

void F77_SUB(rprintf)(const char* msg) {
   Rprintf(msg);
   Rprintf("\n");
}

// may be redundant
void F77_SUB(rprintf2)(const char* msg) {
   Rprintf(msg);
   Rprintf("\n");
}

void F77_SUB(rprintfid)(const char* msg, int *i, double *d) {
   Rprintf(msg, *i, *d);
   Rprintf("\n");
}

void F77_SUB(rprintfdi)(const char* msg, double *d, int *i) {
   Rprintf(msg, *d, *i);
   Rprintf("\n");
}

void F77_SUB(rprintfdid)(const char* msg, double *d1, int *i, double *d2) {
   Rprintf(msg, *d1, *i, *d2);
   Rprintf("\n");
}

void F77_SUB(rprintfd1)(const char* msg, double *d) {
   Rprintf(msg, *d);
   Rprintf("\n");
}

void F77_SUB(rprintfd2)(const char* msg, double *d1, double *d2) {
   Rprintf(msg, *d1, *d2);
   Rprintf("\n");
}

void F77_SUB(rprintfi1)(const char* msg, int *i) {
   Rprintf(msg, *i);
   Rprintf("\n");
}

void F77_SUB(rprintfi2)(const char* msg, int *i1, int *i2) {
   Rprintf(msg, *i1, *i2);
   Rprintf("\n");
}

void F77_SUB(rprintfi3)(const char* msg, int *i1, int *i2, int* i3) {
   Rprintf(msg, *i1, *i2, *i3);
   Rprintf("\n");
}
