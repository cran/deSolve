/*==========================================================================*/
/* Runge-Kutta Solvers, (C) Th. Petzoldt, License: GPL >=2                  */
/* General RK Solver for methods with adaptive step size                    */
/* -- main loop == core function --                                         */
/*==========================================================================*/

#include "rk_util.h"

void rk_fixed(
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
       double* f, double* y, double* Fj, double* tmp,
       double* FF, double* rr, double* A, double* out,
       double* bb1, double* cc,
       double* yknots, double* yout,
       /* SEXPs */
       SEXP Func, SEXP Parms, SEXP Rho
  ) {

  int i = 0, j = 0, one = 1;
  int iknots = *_iknots, it = *_it, it_ext = *_it_ext, it_tot = *_it_tot;
  double t_ext;
  double dt = *_dt;

  /*------------------------------------------------------------------------*/
  /* Main Loop                                                              */
  /*------------------------------------------------------------------------*/
  //Rprintf("1: dt, hini = %g , %g\n", dt, hini);
  do {
    /* select time step (possibly irregular) */
    if (fabs(hini) < (DBL_EPSILON * 100.0))
      dt = tt[it] - tt[it-1];
    else
      dt = fmin(fabs(hini), fabs(tmax - t)) * sign(hini);
    //Rprintf("dt, hini = %g , %g\n", dt, hini);
    timesteps[0] = timesteps[1];
    timesteps[1] = dt;

    /******  Prepare Coefficients from Butcher table ******/
    /* NOTE: the fixed-step solver needs coefficients as vector, not matrix!  */
    for (j = 0; j < stage; j++) {
      if (j == 0)
        for(i = 0; i < neq; i++) Fj[i] = 0;
      else
        for(i = 0; i < neq; i++)
          Fj[i] = A[j] * FF[i + neq * (j - 1)] * dt;
      for (int i = 0; i < neq; i++) {
        tmp[i] = Fj[i] + y0[i];
      }
      /******  Compute Derivatives ******/
      derivs(Func, t + dt * cc[j], tmp, Parms, Rho, FF, out, j, neq,
        ipar, isDll, isForcing);
    }

    /*====================================================================*/
    /* Estimation of new values                                           */
    /*====================================================================*/

    /* use BLAS with reduced error checking */
    blas_matprod1(FF, neq, stage, bb1, stage, one, dy1);

    it_tot++; /* count total number of time steps */
    for (i = 0; i < neq; i++) {
      y1[i] = y0[i] +  dt * dy1[i];
    }

    /*====================================================================*/
    /*      Interpolation and Data Storage                                */
    /*====================================================================*/
    if (interpolate) {
      /*------------------------------------------------------------------*/
      /* "Neville-Aitken-Interpolation";                                  */
      /* the fixed step integrators have no dense output                  */
      /*------------------------------------------------------------------*/
      /* (1) collect number "nknots" of knots in advance */
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
      /* No interpolation mode(for step to step integration);               */
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
      Rprintf("error in RK solver rk_fixed.c: output buffer overflow\n");
      break;
    }
    if (it_tot > maxsteps) {
      istate[0] = -1;
      warning("Number of time steps %i exceeded maxsteps at t = %g\n", it, t);
      break;
    }
    /* tolerance to avoid rounding errors */
  } while (fabs(t - tmax) > 100.0 * DBL_EPSILON); /* end of rk main loop */

  /* return reference values */
  *_iknots = iknots; *_it = it; *_it_ext = it_ext; *_it_tot = it_tot;
}

