/*----------------------------------------------------------------
 The chemical model example of daspk but with the
 production rate a forcing function rather than
 a parameter...
----------------------------------------------------------------*/
#include <R.h>

/* -------- ChemicalDAE.c -> ChemicalDAE.dll ------
c compile in R with: system("g77 -shared -o ChemicalDAE.dll ChemicalDAE.c")
c or with system("R CMD SHLIB ChemicalDAE.c") */


/* A trick to address the parameters and forcings by name */
static double parms[3];
static double forc[1];

#define K parms[0]
#define ka parms[1]
#define r parms[2]
#define prod forc[0]

/*----------------------------------------------------------------
 Initialiser for parameters
----------------------------------------------------------------*/
void initparms(void (* daspkparms)(int *, double *)) {

  int N=3;
  daspkparms(&N, parms);

}

/*----------------------------------------------------------------
c Initialiser for forcings
----------------------------------------------------------------*/
void initforcs(void (* daspkforcs)(int *, double *))  {

  int N=1;
  daspkforcs(&N, forc);

}

/*----------------------------------------------------------------
 Derivatives
----------------------------------------------------------------*/
void chemres (double *t, double *y, double *ydot, double *cj, double *delta,
              int *ires, double *out, int *ip)  {

    double ra, rb;
    if (ip[0] <2) error("nout should be at least 2");

    ra  = ka* y[2];            /* forward rate */
    rb  = ka/K *y[0] * y[1];   /* backward rate */

    /* residuals of rates of changes */
    delta[2] = -ydot[2]  -  ra + rb + prod;
    delta[0] = -ydot[0]  +  ra - rb;
    delta[1] = -ydot[1]  +  ra - rb - r*y[1];
    out[0]   = y[0] + y[1] + y[2];
    out[1]   = prod;

}
