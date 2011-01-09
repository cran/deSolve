/*==========================================================================*/
/* Runge-Kutta Solvers, (C) Th. Petzoldt, License: GPL >= 2                 */
/* General RK Solver for methods with adaptive step size                    */
/* -- main loop == core function --                                         */
/* Parts inspired by Press et al., 2002, 2007;                              */
/*   see vignette for full references                                       */
/*==========================================================================*/

#include "rk_util.h"

void rk_auto(
       /* integers */
       int fsal, int neq, int stage,
       int isDll, int isForcing, int verbose,
       int nknots, int interpolate, int densetype, int maxsteps, int nt,
       /* int pointers */
       int* _iknots, int* _it, int* _it_ext, int* _it_tot, int* _it_rej,
       int* istate,  int* ipar,
       /* double */
       double t, double tmax, double hmin, double hmax, 
       double alpha, double beta,
       /* double pointers */
       double* _dt, double* _errold,
       /* arrays */
       double* tt, double* y0, double* y1, double* y2, double* dy1, double* dy2,
       double* f, double* y, double* Fj, double* tmp,
       double* FF, double* rr, double* A, double* out, 
       double* bb1, double* bb2, double* cc, double* dd, 
       double* atol, double* rtol, double* yknots, double* yout,
       /* SEXPs */
       SEXP Func, SEXP Parms, SEXP Rho
  ) 
{

  int i = 0, j = 0, j1 = 0, k = 0, accept = FALSE, nreject = *_it_rej, one = 1; 
  int iknots = *_iknots, it = *_it, it_ext = *_it_ext, it_tot = *_it_tot;
  double err, dtnew, t_ext;
  double dt = *_dt, errold = *_errold;

  /* todo: make this user adjustable */
  static const double minscale = 0.2, maxscale = 10.0, safe = 0.9;

  /*------------------------------------------------------------------------*/
  /* Main Loop                                                              */
  /*------------------------------------------------------------------------*/
  do {
    if (accept) timesteps[0] = timesteps[1];
    timesteps[1] = dt;

    /*  save former results of last step if the method allows this
       (first same as last)                                             */
    /* Karline: improve by saving "accepted" FF, use this when rejected */
    if (fsal && accept){
      j1 = 1;
      for (i = 0; i < neq; i++) FF[i] = FF[i + neq * (stage - 1)];
    } else {
      j1 = 0;
    }
    /******  Prepare Coefficients from Butcher table ******/
    for (j = j1; j < stage; j++) {
      for(i = 0; i < neq; i++) Fj[i] = 0;
        k = 0;
        while(k < j) {
          for(i = 0; i < neq; i++)
            Fj[i] = Fj[i] + A[j + stage * k] * FF[i + neq * k] * dt;
          k++;
        }
        for (int i = 0; i < neq; i++) {
          tmp[i] = Fj[i] + y0[i];
        }
        /******  Compute Derivatives ******/
        /* pass option to avoid unnecessary copying in derivs */
        derivs(Func, t + dt * cc[j], tmp, Parms, Rho, FF, out, j, neq, 
               ipar, isDll, isForcing);
    }

    /*====================================================================*/
    /* Estimation of new values                                           */
    /*====================================================================*/

    /* use BLAS wrapper with reduced error checking */
    blas_matprod1(FF, neq, stage, bb1, stage, one, dy1);
    blas_matprod1(FF, neq, stage, bb2, stage, one, dy2);

    it_tot++; /* count total number of time steps */
    for (i = 0; i < neq; i++) {
      y1[i] = y0[i] + dt * dy1[i];
      y2[i] = y0[i] + dt * dy2[i];
    }

    /*====================================================================*/
    /*      stepsize adjustment                                           */
    /*====================================================================*/
    
    err = maxerr(y0, y1, y2, atol, rtol, neq);
    dtnew = dt;
    if (err == 0) {  /* use max scale if all tolerances are zero */
      dtnew  = fmin(dt * 10, hmax);
      errold = fmax(err, 1e-4); /* 1e-4 taken from Press et al. */
      accept = TRUE;
    } else if (err < 1.0) {
      /* increase step size only if last one was accepted */
      if (accept) 
        dtnew = fmin(hmax, dt * 
          fmin(safe * pow(err, -alpha) * pow(errold, beta), maxscale));
      errold = fmax(err, 1e-4); /* 1e-4 taken from Press et al. */
      accept = TRUE;
    } else if (err > 1.0) {
      nreject++;    /* count total number of rejected steps */
      accept = FALSE;
      dtnew = dt * fmax(safe * pow(err, -alpha), minscale);
    }

    if (dtnew < hmin) {
      accept = TRUE;
      if (verbose) Rprintf("warning, h < Hmin\n");
      istate[0] = -2;
      dtnew = hmin;
    }
    /*====================================================================*/
    /*      Interpolation and Data Storage                                */
    /*====================================================================*/
    if (accept) {
      if (interpolate) {
      /*--------------------------------------------------------------------*/
      /* case A1) "dense output type 1": built-in polynomial interpolation          */
      /* available for certain rk formulae, e.g. for rk45dp7                */
      /*--------------------------------------------------------------------*/
      if (densetype == 1) {
        denspar(FF, y0, y2, dt, dd, neq, stage, rr);
        t_ext = tt[it_ext];
        while (t_ext <= t + dt) {
          densout(rr, t, t_ext, dt, tmp, neq);
          /* store outputs */
          if (it_ext < nt) {
            yout[it_ext] = t_ext;
            for (i = 0; i < neq; i++)
              yout[it_ext + nt * (1 + i)] = tmp[i];
          }
          if(it_ext < nt-1) t_ext = tt[++it_ext]; else break;
        }

        /*--------------------------------------------------------------------*/
        /* case A2) dense output type 2: the Cash-Karp method                */
        /*--------------------------------------------------------------------*/
      } else if (densetype == 2)  {   /* dense output method 2 = Cash-Karp */
        derivs(Func, t + dt, y2, Parms, Rho, dy2, out, 0, neq, 
               ipar, isDll, isForcing);
        
        t_ext = tt[it_ext];

        while (t_ext <= t + dt) {
          densoutck(t, t_ext, dt, y0, FF, dy2, tmp, neq);
          /* store outputs */
          if (it_ext < nt) {
            yout[it_ext] = t_ext;
            for (i = 0; i < neq; i++)
              yout[it_ext + nt * (1 + i)] = tmp[i];
          }
          if(it_ext < nt-1) t_ext = tt[++it_ext]; else break;
       }
       /* fsal trick for Cash-Karp */
       for (i = 0; i < neq; i++) FF[i + neq * (stage - 1)] = dy2[i] ;

        /*--------------------------------------------------------------------*/
        /* case B) "Neville-Aitken-Interpolation" for integrators             */
        /* without dense output                                               */
        /*--------------------------------------------------------------------*/
      } else {
          /* (1) collect number "nknots" of knots in advance */
          yknots[iknots] = t + dt;   /* time is first column */
          for (i = 0; i < neq; i++) yknots[iknots + nknots * (1 + i)] = y2[i];
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
        }
      } else {
        /*--------------------------------------------------------------------*/
        /* Case C) no interpolation at all (for step to step integration);    */
        /*         results are stored after the call                          */
        /*--------------------------------------------------------------------*/
      }
      /*--------------------------------------------------------------------*/
      /* next time step                                                     */
      /*--------------------------------------------------------------------*/
      t = t + dt;
      it++;
      for (i=0; i < neq; i++) y0[i] = y2[i];
    } /* else rejected time step */
    dt = fmin(dtnew, tmax - t);
    if (it_ext > nt) {
      Rprintf("error in RK solver rk_auto.c: output buffer overflow\n");
      break;
    }
    if (it_tot > maxsteps) {
      if (verbose) Rprintf("Max. number of steps exceeded\n");
      istate[0] = -1;
      break;
    }
    /* tolerance to avoid rounding errors */
  } while (t < (tmax - 100.0 * DBL_EPSILON * dt)); /* end of rk main loop */

  /* return reference values */
  *_iknots = iknots; *_it = it; *_it_ext = it_ext; *_it_rej = nreject;
  *_it_tot = it_tot; *_dt = dtnew; *_errold = errold;
}
