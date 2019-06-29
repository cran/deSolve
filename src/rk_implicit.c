/*==========================================================================*/
/* Implicit RK Solver with fixed step size                                  */
/*==========================================================================*/

#include "rk_util.h"
void F77_NAME(dgefa)(double*, int*, int*, int*, int*);
void F77_NAME(dgesl)(double*, int*, int*, int*, double*, int*);
/*
void lu_solve(double, int, int, double);
void kfunc(int, int, double, double, double, double, double, double,
      double, SEXP, SEXP, SEXP, double, double, double, int, int, int);
void dkfunc(int, int, double, double, double, double, double, double, double,
      SEXP, SEXP, SEXP, double, double, double, double, int, int, int, double);
*/

/* lower upper decomposition - no error checking */
void lu_solve(double *alfa, int n, int *index, double *bet) {
  int info;

  F77_CALL(dgefa)(alfa, &n, &n, index, &info);
	if (info != 0)
    error("error during factorisation of matrix (dgefa), singular matrix");
  F77_CALL(dgesl)(alfa, &n, &n, index, bet, &info);
	if (info != 0)
    error("error during backsubstitution");
}

/* function that returns -k + dt*derivs(t+c[i]*dt, y+sum(a[i,)*k
   this is the function whose roots should be found in the implicit method */

void kfunc(int stage, int neq, double t, double dt,
   double *FF, double *Fj, double *A, double *cc, double *y0 ,
   SEXP Func, SEXP Parms, SEXP Rho, double *tmp, double *tmp2,
   double *out, int *ipar, int isDll, int isForcing){

   int i, j, k;
   /******  Prepare Coefficients from Butcher table ******/
   for (j = 0; j < stage; j++) {
     for (i = 0; i < neq; i++) Fj[i] = 0.;
     for (k =0; k < stage; k++) { /* implicit part */
       for(i = 0; i < neq; i++)
         Fj[i] = Fj[i] + A[j + stage * k] * FF[i + neq * k] * dt;
     }
     for (int i = 0; i < neq; i++) {
       tmp[i] = Fj[i] + y0[i];
     }
     /******  Compute Derivatives ******/
     /* pass option to avoid unnecessary copying in derivs note:tmp2 rather than FF */
     derivs(Func, t + dt * cc[j], tmp, Parms, Rho, tmp2, out, j, neq,
              ipar, isDll, isForcing);
   }
   for (i = 0; i< neq*stage;i++)
     tmp[i] = FF[i] - tmp2[i];       /* tmp should be = 0 at root */
}

/* function that returns the Jacobian of kfunc; df[i,j] should contain:
   dkfunc_i/dFFj CHECK */
void dkfunc(int stage, int neq, double t, double dt,
   double *FF, double *Fj, double *A, double *cc, double *y0,
   SEXP Func, SEXP Parms, SEXP Rho, double *tmp, double *tmp2, double *tmp3,
   double *out, int *ipar, int isDll, int isForcing, double *df){

   int i, j, nroot;
   double d1, d2;

   nroot = neq*stage;

   /* function reference value in tmp2 */
   kfunc(stage, neq, t, dt, FF, Fj, A, cc, y0, Func, Parms, Rho,
         tmp2, tmp3, out, ipar, isDll, isForcing);

   for (i = 0; i < nroot; i++) {
     d1 = FF[i];                      /* copy */
     d2 = fmax(1e-8, FF[i] * 1e-8);     /* perturb */
     FF[i] = FF[i] + d2;
     kfunc(stage, neq, t, dt, FF, Fj, A, cc, y0, Func, Parms, Rho,
        tmp, tmp3, out, ipar, isDll, isForcing);
     for (j = 0; j < nroot; j++)
       df[nroot * i + j] = (tmp[j] - tmp2[j])/d2;   //df[j,i] j,i=1:nroot
     FF[i] = d1;                      /* restore */
   }
}

/* ks: check if tmp3 necessary ... */
void rk_implicit( double * alfa,  /* neq*stage * neq*stage */
       int *index,                /* neq*stage */
       /* integers */
       int fsal, int neq, int stage,
       int isDll, int isForcing, int verbose,
       int nknots, int interpolate, int maxsteps, int nt,
       /* int pointers */
       int* _iknots, int* _it, int* _it_ext, int* _it_tot,
       int* istate,  int* ipar,
       /* double */
        double t, double tmax, double hini,
       /* double pointers */
       double* _dt,
       /* arrays */
       double* tt, double* y0, double* y1, double* dy1,
       double* f, double* y, double* Fj,
       double* tmp, double* tmp2, double* tmp3,
       double* FF, double* rr, double* A, double* out,
       double* bb1, double* cc,
       double* yknots, double* yout,
       /* SEXPs */
       SEXP Func, SEXP Parms, SEXP Rho
  )
{
  int i = 0, one = 1;
  int iknots = *_iknots, it = *_it, it_ext = *_it_ext, it_tot = *_it_tot;
  double t_ext;
  double dt = *_dt;
  int iter, maxit = 100;
  double errf, errx;
  int nroot = neq * stage;

  /*------------------------------------------------------------------------*/
  /* Main Loop                                                              */
  /*------------------------------------------------------------------------*/
  do {
    /* select time step (possibly irregular) */
    if (hini > 0.0)
      dt = fmin(hini, tmax - t); /* adjust dt for step-by-step-mode */
    else
      dt = tt[it] - tt[it-1];

    timesteps[0] = timesteps[1];
    timesteps[1] = dt;

    /* Newton-Raphson steps */
    for (iter = 0; iter < maxit; iter++) {
      /* function value and Jacobian*/
      kfunc(stage, neq, t, dt, FF, Fj, A, cc, y0, Func, Parms, Rho,
        tmp, tmp2, out, ipar, isDll, isForcing);
      it_tot++; /* count total number of time steps */
      errf = 0.;
      for ( i = 0; i < nroot; i++) errf = errf + fabs(tmp[i]);
      if (errf < 1e-8) break;
      dkfunc(stage, neq, t, dt, FF, Fj, A, cc, y0, Func, Parms, Rho,
        tmp, tmp2, tmp3, out, ipar, isDll, isForcing, alfa);
      it_tot = it_tot + nroot + 1;
      lu_solve (alfa, nroot, index, tmp);
      errx = 0;
      for (i = 0; i < nroot; i++) {
        errx = errx + fabs(tmp[i]);
        FF[i] = FF[i] - tmp[i];
      }
      //  Rprintf("iter %i errf %g errx %g\n",iter, errf, errx);
      if (errx < 1e-8) break;
    }

    /*====================================================================*/
    /* Estimation of new values                                           */
    /*====================================================================*/

    /* use BLAS with reduced error checking */
    blas_matprod1(FF, neq, stage, bb1, stage, one, dy1);

    for (i = 0; i < neq; i++) {
      y1[i] = y0[i] +  dt * dy1[i];
    }

    /*====================================================================*/
    /*      Interpolation and Data Storage                                */
    /*====================================================================*/
    if (interpolate) {
      /*------------------------------------------------------------------*/
      /* Neville-Aitken-Interpolation                                     */
      /* the fixed step integrators have no dense output                  */
      /*------------------------------------------------------------------*/
      /* (1) collect number "nknots" of knots in advanve */
      yknots[iknots] = t + dt;   /* time in first column */
      for (i = 0; i < neq; i++) yknots[iknots + nknots * (1 + i)] = y1[i];
      if (iknots < (nknots - 1)) {
        iknots++;
      } else {
       /* (2) do polynomial interpolation */
       t_ext = tt[it_ext];
       while (t_ext <= t + dt) {
        neville(yknots, &yknots[nknots], t_ext, tmp, nknots, neq);
        /* (3) store outputs */
        if (it_ext < nt) {
          yout[it_ext] = t_ext;
          for (i = 0; i < neq; i++)
            yout[it_ext + nt * (1 + i)] = tmp[i];
        }
        if(it_ext < nt-1) t_ext = tt[++it_ext]; else break;
       }
       shiftBuffer(yknots, nknots, neq + 1);
      }
    } else {
      /*--------------------------------------------------------------------*/
      /* No interpolation mode for step to step integration                 */
      /*         results are stored after the call                          */
      /*--------------------------------------------------------------------*/
    }
    /*--------------------------------------------------------------------*/
    /* next time step                                                     */
    /*--------------------------------------------------------------------*/
    t = t + dt;
    it++;
    for (i = 0; i < neq; i++) y0[i] = y1[i];
    if (it_ext > nt) {
      Rprintf("error in RK solver rk_implicit.c: output buffer overflow\n");
      break;
    }
    if (it_tot > maxsteps) {
      istate[0] = -1;
      warning("Number of time steps %i exceeded maxsteps at t = %g\n", it, t);
      break;
    }
    /* tolerance to avoid rounding errors */
  } while (t < (tmax - 100.0 * DBL_EPSILON * dt)); /* end of rk main loop */

  /* return reference values */
  *_iknots = iknots; *_it = it; *_it_ext = it_ext; *_it_tot = it_tot;
}

