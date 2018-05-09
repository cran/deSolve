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
   to be found ("findHistInt"), and then either use hermite interpolation 
   to the requested time (functions "Hermite" and "dHermite" for values and 
   derivatives), or use the Nordsieck history array. 
   
   Note: findHistInt finds interval by bisectioning; only marginally
   more/less efficient than straightforward findHistInt2...
   
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/*=========================================================================== 
  Higher-order interpolation of y to x, based on Nordsieck history array 
  if interpolMethod ==2
  =========================================================================== */

/* definition of call to FORTRAN function INTERPOLY, as derived from dintdy   */

void F77_NAME(interpoly)(double *, int *, int *, double *, int *, double *, 
   int *, double *, double *);
   
double interpolate(int i, int k, double t0, double hh, double t, 
  double *Yh, int nq) {
  
  double  res;
 
  if (nq > 12)
    error("illegal nq in interpolate, %i, at time %g", nq, t);
  if (k > nq)
    error("illegal k %i, nq in interpolate, %i, at time %g", k, nq, t);
                
  if (i > n_eq || i <1)
    error("illegal i %i, n_eq %i, at time %g", i, n_eq, t);

  F77_CALL(interpoly) (&t, &k, &i, Yh, &n_eq, &res, &nq, &t0, &hh); 
  return(res);
}  

/* continuous output formula for radau                                        */
void F77_NAME (contr5alone) (int *, int *, double *, double *, int *, double *,
                             double *, int *);
void F77_NAME (getconra) (double *);

/*=========================================================================== 
  Hermitian interpolation of y to x (interpolMethod==1)
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

void inithist(int max, int maxlags, int solver, int nroot) {
  int maxord;  
  
  histsize = max;
  initialisehist = 1;
  indexhist  = -1; /* indexhist+1 = next time in circular buffer.  */
  starthist  = 0;  /* start time in circular buffer.               */
  endreached = 0;  /* if end of buffer reached and new values added at start  */
  
  /* interpolMethod = Hermite */
  if (interpolMethod == 1) {
    offset   = n_eq; /* size needed for saving one time-step in histvar*/

  /* interpolMethod = HigherOrder, Livermore solvers */
  } else if (interpolMethod == 2) {
    if (solver == 0) 
      error("illegal input in lags - cannot combine interpol=2 with chosen solver");
    maxord  = 12;   /* 5(bdf) or 12 (adams) */
    lyh     = 20;   /* position of history array in rwork (C-index) */
    lhh     = 11;   /* position of h in rwork (C-index)
                       Note: for lsodx this is NEXT time step! */
    lo      = 13;   /* position of method order in iwork (C-index) */
    if (solver == 5) {  /* different for vode!  uses current time step*/  
      lhh = 10;        
      lo = 13;             
    }
    if (solver == 4 || solver == 6 || solver == 7)  /* lsodar or lsoder */
      lyh = 20+3*nroot;

    offset  = n_eq*(maxord+1);       
    histord = (int *) R_alloc (histsize, sizeof(int));
    histhh  = (double *) R_alloc (histsize, sizeof(double));

  /* interpolMethod = 3; HigherOrder, radau */
  } else {
    offset  = n_eq * 4 + 2;
    histsave = (double *) R_alloc (2, sizeof(double));
  }

  histtime = (double *) R_alloc (histsize, sizeof(double));
  histvar  = (double *) R_alloc (offset * histsize, sizeof(double));
  histdvar = (double *) R_alloc (n_eq * histsize, sizeof(double));
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

/* first time: just store y, (dy) and t */
void updatehistini(double t, double *y, double *dY, double *rwork, int *iwork){
  int intpol;

  intpol = interpolMethod;
  interpolMethod = 1; 
  updatehist(t, y, dY, rwork, iwork);
  interpolMethod = intpol; 
  if (interpolMethod == 2){
    histord[0] = 0;
    histhh[0] = timesteps[0];    
  }  
}

void updatehist(double t, double *y, double *dY, double *rwork, int *iwork) {
  int j, ii;
  double ss[2];
  
  indexhist = nexthist(indexhist);
  ii = indexhist * offset;     

  /* interpolMethod = Hermite */
  if (interpolMethod == 1) {
    for (j = 0; j < n_eq; j++)  
      histvar [ii  + j ] = y[j];

  /* higherOrder, livermores */
  } else if (interpolMethod == 2) {
    histord[indexhist] = iwork[lo];    

    for (j = 0; j < offset; j++)
      histvar[ii + j] = rwork[lyh + j];
    histhh [indexhist] = rwork[lhh];   

  /* higherOrder, radau */
  }  else if (interpolMethod == 3) {
    for (j = 0; j < 4 * n_eq; j++)
      histvar[ii + j] = rwork[j];
    F77_CALL(getconra) (ss);
    for (j = 0; j < 2; j++)
      histvar[ii + 4*n_eq + j] = ss[j];
  }

  ii = indexhist * n_eq;     
 
  for (j = 0; j < n_eq; j++)
      histdvar[ii + j] = dY[j];

  histtime [indexhist] = t;

  if (endreached == 1)       /* starthist stays 0 until end reached... */
    starthist = nexthist(starthist);
}

/*=========================================================================== 
  find a past value (val=1) or a past derivative (val = 2)
  =========================================================================== */

double past(int i, int interval, double t, int val)

  /* finds past values (val=1) or past derivatives (val=2)*/

{ int j, jn, nq, ip;
  double t0, t1, y0, y1, dy0, dy1, res, hh;
  double *Yh;

  /* error checking */
  if ( i >= n_eq)
    error("illegal input in lagvalue - var nr too high, %i", i+1);
  
  /* equal to current value... */   
  if ( interval == indexhist && t == histtime[interval]) {   
    if (val == 1)
      res = histvar [interval * offset  + i ];
    else 
      res = histdvar [interval * offset  + i ];   
  
  /* within last interval - for now: just extrapolate last value */
  } else if ( interval == indexhist && interpolMethod == 1) {
    if (val == 1) {
      t0  = histtime[interval];
      y0  = histvar [interval * offset  + i ];
      dy0 = histdvar [interval * n_eq  + i ];
      res = y0 + dy0*(t-t0);
    }
    else 
      res = histdvar [interval * n_eq  + i ];

  /* Hermite interpolation */
  }  else if (interpolMethod == 1) {
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
  
  /* dense interpolation - livermore solvers */
  } else if (interpolMethod == 2) {
    j  = interval;
    jn = nexthist(j);

    t0  = histtime[j];
    t1  = histtime[jn];
    nq  = histord [j];
    if (nq == 0) {
      y0  = histvar [j  * offset  + i ];
      y1  = histvar [jn * offset  + i ];
      dy0 = histdvar [j  * n_eq  + i ];
      dy1 = histdvar [jn * n_eq  + i ];
      if (val == 1)
        res = Hermite (t0, t1, y0, y1, dy0, dy1, t);
      else
        res = dHermite (t0, t1, y0, y1, dy0, dy1, t);
    } else { 
      Yh  = &histvar [j * offset];
      hh = histhh[j];
      res = interpolate(i+1, val-1, t0, hh, t, Yh, nq); 
    }  
  /* dense interpolation - radau - gets all values (i not used) */
  } else {
 //   if (val == 2)
 //     error("radau interpol = 2 does not work for lagderiv");
    j  = interval;
    Yh  = &histvar [j * offset];
    histsave  = &histvar [j * offset + 4*n_eq];
    ip = i+1;
    F77_CALL(contr5alone) (&ip, &n_eq, &t, Yh, &offset, histsave, &res, &val);
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
      t, histtime[indexhist]);
 
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
      t, histtime[indexhist]);

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
    error("pastvalue can only be called from 'func' or 'res' when triggered by appropriate integrator.");
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
    error("pastgradient can only be called from 'func' or 'res' when triggered by appropriate integrator.");
  if (!isNumeric(T)) error("'t' should be numeric");

  t = *NUMERIC_POINTER(T);
  interval = findHistInt (t);

  if ((ilen ==1) && (INTEGER(nr)[0] == 0)) {
    PROTECT(value=NEW_NUMERIC(n_eq));
    for(i=0; i<n_eq; i++) {
      NUMERIC_POINTER(value)[i] = past(i, interval, t, 2);
    }
  } else {
    PROTECT(value=NEW_NUMERIC(ilen));
    for(i=0; i<ilen; i++) {                                              
      NUMERIC_POINTER(value)[i] = past(INTEGER(nr)[i]-1, interval, t, 2);
    }
  }
  UNPROTECT(1);
  return(value);
}


/* ============================================================================
  Interrogate the lag settings as in an R-list   
   ==========================================================================*/

int initLags(SEXP elag, int solver, int nroot) {

  SEXP Mxhist, Islag, Interpol ;       
  int mxhist, islag;
    
  Islag = getListElement(elag, "islag");
  islag = INTEGER(Islag)[0];
    
  if (islag == 1) {
   Mxhist = getListElement(elag, "mxhist");
   mxhist = INTEGER(Mxhist)[0];
   Interpol = getListElement(elag, "interpol");
   interpolMethod = INTEGER(Interpol)[0];
   if (interpolMethod < 1) interpolMethod = 1;
   if ((interpolMethod == 2) && (solver == 10)) interpolMethod = 3; /* radau */
//   if((solver == 7 || solver == 3) && interpolMethod == 2)
//     error("cannot combine lags in lsodes, with interpol=2");
   inithist(mxhist, 1, solver, nroot);
  } else {
    mxhist = 0;
    interpolMethod = 1;
  }  
  return(islag);
}

/*===========================================================================
  lagderiv and lagderiv versions for use in compiled (.so/.dll) models
  see interface in R_init_deSolve.c
  tested with C , but not yet with Fortran
  thpe 2013-03-21
  =========================================================================== */

void lagvalue(double T, int *nr, int N, double *ytau) {
  int i, interval;

  if (initialisehist == 0)
    error("pastvalue can only be called from 'func' or 'res' when triggered by appropriate integrator.");

  interval = findHistInt(T);
  for(i = 0; i < N; i++)  ytau[i] = past(nr[i], interval, T, 1);
}

void lagderiv(double T, int *nr, int N, double *ytau) {
  int i, interval;

  if (initialisehist == 0)
    error("pastvalue can only be called from 'func' or 'res' when triggered by appropriate integrator.");

  interval = findHistInt(T);

  for(i = 0; i < N; i++)  ytau[i] = past(nr[i], interval, T, 2);
}

