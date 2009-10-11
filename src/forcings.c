/* deals with forcing functions that are passed via arguments in the call
to the integration routines */

#include "deSolve.h"

int    finit = 0;

/*         -----     Check for presence of forcing functions     -----        */


int initForcings(SEXP flist) {

    SEXP Tvec, Fvec, Ivec, initforc;
    int i, j, isForcing = 0;
    init_func  *initforcings;


    initforc = getListElement(flist,"ModelForc");
    if (!isNull(initforc))
    	{
    	 Tvec = getListElement(flist,"tmat");
       Fvec = getListElement(flist,"fmat");
       Ivec = getListElement(flist,"imat");
       nforc =LENGTH(Ivec)-2; /* nforc, fvec, ivec =globals */

       i = LENGTH(Fvec);
       fvec = (double *) R_alloc((int) i, sizeof(double));
       for (j = 0; j < i; j++) fvec[j] = REAL(Fvec)[j];

       tvec = (double *) R_alloc((int) i, sizeof(double));
       for (j = 0; j < i; j++) tvec[j] = REAL(Tvec)[j];

       i = LENGTH (Ivec)-1; /* last element: the interpolation method...*/
       ivec = (int *) R_alloc(i, sizeof(int));
       for (j = 0; j < i; j++) ivec[j] = INTEGER(Ivec)[j];

       fmethod =INTEGER(Ivec)[i];

	     initforcings = (init_func *) R_ExternalPtrAddr(initforc);
	     initforcings(Initdeforc);

       isForcing = 1;
       }
       return(isForcing);
}

/*         -----     INITIALISATION  called from compiled code   -----
   1. Check the length of forcing functions in solver call and code in DLL
   2. Initialise the forcing function vectors
   3. set pointer to DLL; FORTRAN common block or C globals /
*/


/* THESE ARE CLEANER VERSIONS ; the other versions are in isnt/removed.txt*/

void Initdeforc(int *N, double *forc)
{
  int i, ii;

  if ((*N) != nforc) {
    warning("Number of forcings passed to solver, %i; number in DLL, %i\n",nforc, *N);

    PROBLEM "Confusion over the length of forc"
    ERROR;
  }

/* for each forcing function: index to current position of data,
 current value, interpolation factor, current forcing time, next forcing time,..
*/
   finit = 1;
   findex   = (int    *) R_alloc(nforc, sizeof(int));
   intpol   = (double *) R_alloc(nforc, sizeof(double));
   maxindex = (int    *) R_alloc(nforc, sizeof(int));

/* Input is in three vectors:
   tvec, fvec: time and value;
   ivec : index to each forcing in tvec and fvec
*/
   for (i = 0; i<nforc; i++) {
     ii = ivec[i]-1;
     findex[i] = ii;
     maxindex[i] = ivec[i+1]-2;
     if (fmethod == 1) {
       intpol[i] = (fvec[ii+1]-fvec[ii])/(tvec[ii+1]-tvec[ii]);
     } else  intpol[i] = 0;
     forc[i] = fvec[ii];
   }
   forcings = forc;      /* set pointer to C globals or FORTRAN common block */
}

void updatedeforc(double *time)
{
  int i, ii,  zerograd;

/* check if initialised? */
   if (finit == 0)
     error ("error in forcing function: not initialised");

   for (i=0; i<nforc; i++) {
     ii = findex[i];

     zerograd=0;
     while (*time > tvec[ii+1]){
         if (ii+2 > maxindex[i]) {   /* this probably redundant...*/
           zerograd=1;
           break;
         }
         ii = ii+1;
       }
    while (*time < tvec[ii]){       /* test here for ii < 1 ?...*/
         ii = ii-1;
    }
    if (ii != findex[i]) {
       findex[i] = ii;
       if ((zerograd == 0) & (fmethod == 1)) {  /* fmethod 1=linear */
         intpol[i] = (fvec[ii+1]-fvec[ii])/(tvec[ii+1]-tvec[ii]); }
       else {
         intpol[i] = 0;
       }
     }

     forcings[i]=fvec[ii]+intpol[i]*(*time-tvec[ii]);
   }
}

