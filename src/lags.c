#include <time.h>
#include <string.h>
#include "deSolve.h"

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   time lags and delay-differential equations; from deSolve version 1.7
   
   For delay-differential equations, a history of past values, past derivatives 
   and past times, is kept (time-lags). 
   
   They are in ring-vectors "histvar", "histdvar" and "histtime" respectively. 
   These vectors are initialised at the start of the integration ("inithist")
   and then updated with new values every accepted timestep ("updatehist").
   When the end of the history vectors is reached, new values are stored at the 
   start (it is a ringbuffer); 
   function "nexthist" finds the next position in this ringbuffer.
   
   The history buffers can be interrogated in the R-code, via R-functions
   "lagvalue(t,nr)" and "lagderiv(t,nr)", where nr can be one index or a vector
   containing the nr of the variable whose lag has to be computed at time t.
   
   These R-functions call C-functions "getLagValue" and "getLagDeriv" which
   first find the interval in the history vectors in which the lagged value is 
   to be found ("findHistInt"), and then use hermite interpolation to interpolate
   to the requested time (functions "Hermite" and "dHermite" for values and 
   derivatives). 
   
   findHistInt finds interval by bisectioning; only marginally
   more/less efficient than straightforward findHistInt2...
   
   
   to do:    
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*=========================================================================== 
  Hermitian interpolation of y to x
  =========================================================================== */

double Hermite (double t0, double t1, double y0, double y1, double dy0,
                double dy1, double t) {
  double tt0, tt1, tt12, tt02, hh, res;
  
  tt0 = t-t0;
  tt1 = t-t1;
  tt12 = tt1*tt1;
  tt02 = tt0*tt0;
  hh   = t1-t0;
  if (hh)
    res=( dy0* tt0* tt12 + dy1* tt1* tt02
		   + ( y0* (2.0* tt0 + hh)* tt12 
			   	-y1* (2.0* tt1 - hh)* tt02 )/hh) / (hh * hh);
  else
    res=y0;
  return(res);
}

/*=========================================================================== 
  Hermitian interpolation of dy to x
  =========================================================================== */

double dHermite (double t0, double t1, double y0, double y1, double dy0,
                double dy1, double t) {
  double tt0, tt1, tt12, tt02, hh, res;
  
	tt0 = t-t0;
  tt1 = t-t1;
  tt12 = tt1*tt1;
  tt02 = tt0*tt0;
  hh   = t1-t0;

  if (hh)
    res=( dy0 * (tt12+2.0* tt0* tt1) + dy1 * ( tt02+2.0* tt0* tt1)
			 + ( y0 *2.0* tt1*(2.0* tt0+ hh +  tt1)
				  -y1 *2.0* tt0*(2.0* tt1- hh +  tt0))/ hh ) / ( hh* hh) ;
  else
    res= dy0;
  return(res);
}

/*=========================================================================== 
  initialise history arrays + indices at start integration
  =========================================================================== */

void inithist(int max, int maxlags) {

  histsize = max;
  initialisehist = 1;
  indexhist  = -1; // indexhist+1 = next time in circular buffer.
  starthist  = 0;  // start time in circular buffer.
  endreached = 0;  // if end of buffer reached and new values added at start
  
  histvar  = (double *) R_alloc (n_eq * histsize, sizeof(double));
  histdvar = (double *) R_alloc (n_eq * histsize, sizeof(double));
  histtime = (double *) R_alloc (histsize, sizeof(double));

  /*lagindex = (int *) R_alloc ( maxlags, sizeof(int));  // not yet used
   for (j = 0; j < maxlags; j++) lagindex[j] = 0; */
}

/*=========================================================================== 
  given the maximum size of the history arrays; finds the next index
  =========================================================================== */

int nexthist(int i) {
  if (i < histsize-1)
    return(i+1);
  else {
    endreached = 1;
    return(0);
  }  
}

/*=========================================================================== 
  update history arrays each time step
  =========================================================================== */

void updatehist(double t, double *y, double *dy) {
  int j, ii;

  indexhist = nexthist(indexhist);
  
  ii = indexhist * n_eq;    /* offset */
  
  for (j = 0; j < n_eq; j++) {
    histvar [ii  + j ] = y[j];
    histdvar[ii  + j ] = dy[j];
  }
  histtime [indexhist] = t;

  if (endreached == 1) {      /* starthist stays 0 until end reached... */
   starthist = nexthist(starthist);
  }
}

/* ============================================================================
   not yet used - to be used at start of each timestep, to keep up with the 
   previous position of the time-lagged variable 
void initlag() {
  indexlag = 0;
}
*/

/*=========================================================================== 
  find a past value (val==1) or a past derivative (val == 0) 
  =========================================================================== */

double past(int i, int interval, double t, int val)

	/* Interrogates the history ringbuffers. Not very efficient...*/

{ int j, jn;
  double t0, t1, y0, y1, dy0, dy1, res;

  /* error checking */
  if ( i > n_eq)
    error("illegal input in lagvalue - var nr too high, %i", i+1);

  if ( interval == indexhist) {    
    if (val == 1)
      res = histvar [interval * n_eq  + i ];
    else 
      res = histdvar [interval * n_eq  + i ];   
  }  else {
   
    j  = interval;
    jn = nexthist(j);

    t0  = histtime[j];
    t1  = histtime[jn];
    y0  = histvar [j * n_eq  + i ];
    y1  = histvar [jn * n_eq  + i ];
    dy0 = histdvar [j * n_eq  + i ];
    dy1 = histdvar [jn * n_eq  + i ];
    if (val == 1)
      res = Hermite (t0, t1, y0, y1, dy0, dy1, t);
    else
      res = dHermite (t0, t1, y0, y1, dy0, dy1, t);
  }
  return(res);
}

/*=========================================================================== 
  Find interval in history ring buffers, corresponding to "t"
  two alternatives; only findHistInt used
  =========================================================================== */
int findHistInt2 (double t) {
  int j, jn;
  
  if ( t >= histtime[indexhist]) 
    return(indexhist);
	if ( t < histtime[starthist])
	  error("illegal input in lagvalue - lag, %g, too large, at time = %g\n",
      t, histtime[starthist]);
 
   /* find embracing time starting from beginning  */
    j  = starthist;
    jn = nexthist(j);

    while (histtime[jn]<t) {
      j = jn;
      jn = nexthist(j);
    }
    return(j);  
}
/* alternative: bisectioning... */

int findHistInt (double t) {
  int ilo, ihi, imid, ii, n;
  
  if ( t >= histtime[indexhist]) 
    return(indexhist);
	if ( t < histtime[starthist])
	  error("illegal input in lagvalue - lag, %g, too large, at time = %g\n",
      t, histtime[starthist]);

  if (endreached == 0) {  /* still filling buffer; not yet wrapped */
    ilo = 0;
    ihi = indexhist;
    for(;;) {
   	  imid = (ilo + ihi) / 2;
	    if (imid == ilo) return ilo;
      if (t >= histtime[imid])
	      ilo = imid;
    	else
	      ihi = imid;
     }
  }
  n = histsize -1;
  ilo = 0;
  ihi = n;
  for(;;) {
   	imid = (ilo + ihi) / 2;
	
    ii = imid + starthist;
    if (ii > n) ii = ii - n - 1;

	  if (imid == ilo) return ii;

    if (t >= histtime[ii])
	    ilo = imid;
	  else
	    ihi = imid;
  }
}

/*=========================================================================== 
  C-equivalent of R-function lagvalue
  =========================================================================== */
SEXP getLagValue(SEXP T, SEXP nr)
{
	SEXP value;
	int i, ilen, interval;
	double t;

  ilen = LENGTH(nr);
	if (initialisehist == 0)
    error("pastvalue can only be called from `func` or 'res' when triggered by appropriate integrator.");
	if (!isNumeric(T)) error("‘t’ should be numeric");

  t = *NUMERIC_POINTER(T);
  interval = findHistInt (t);

  if ((ilen ==1) && (INTEGER(nr)[0] == 0)) {
  	PROTECT(value=NEW_NUMERIC(n_eq));
  	for(i=0; i<n_eq; i++) {
	  	NUMERIC_POINTER(value)[i] = past(i, interval, t, 1);
  }
  } else {
	  PROTECT(value=NEW_NUMERIC(ilen));
  	for(i=0; i<ilen; i++) {
		NUMERIC_POINTER(value)[i] = past(INTEGER(nr)[i]-1, interval, t, 1);
	  }
	}
	UNPROTECT(1);
	return(value);
}

/*===========================================================================
  C-equivalent of R-function lagderiv 
  =========================================================================== */
SEXP getLagDeriv(SEXP T, SEXP nr)
{
	SEXP value;
	int i, ilen, interval;
	double t;

  ilen = LENGTH(nr);
	if (initialisehist == 0)
    error("pastgradient can only be called from `func` or 'res' when triggered by appropriate integrator.");
	if (!isNumeric(T)) error("‘t’ should be numeric");

  t = *NUMERIC_POINTER(T);
  interval = findHistInt (t);

  if ((ilen ==1) && (INTEGER(nr)[0] == 0)) {
  	PROTECT(value=NEW_NUMERIC(n_eq));
  	for(i=0; i<n_eq; i++) {
	  	NUMERIC_POINTER(value)[i] = past(i, interval, t, 0);
	  }
  } else {
	  PROTECT(value=NEW_NUMERIC(ilen));
  	for(i=0; i<ilen; i++) {
	  	NUMERIC_POINTER(value)[i] = past(INTEGER(nr)[i]-1, interval, t, 0);
	  }
	}
  UNPROTECT(1);
	return(value);
}

/* ============================================================================
  Interrogate the lag settings as in an R-list   
   ==========================================================================*/

int initLags(SEXP elag) {

  SEXP Mxhist, Islag;       
  int mxhist, islag;
    
  Islag = getListElement(elag, "islag");
  islag = INTEGER(Islag)[0];

  if (islag ==1) { 
   Mxhist = getListElement(elag, "mxhist");
   mxhist = INTEGER(Mxhist)[0];
   inithist(mxhist,1);
  } else mxhist = 0;

  return(islag);
}





