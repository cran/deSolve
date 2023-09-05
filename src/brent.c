/* brent's rootfinding method, based on R_Zeroin_2, itself based on
  NETLIB c/brent.shar */

/*************************************************************************
 *			    C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,f,info,Tol,Maxit)
 *	double ax;			Root will be seeked for within
 *	double bx;			a range [ax,bx]
 *	double (f)(double x, void *info); Name of the function whose zero
 *					will be seeked for
 *	double *rw; int *iw;	Additional real and integer vector
 *	double tol;			Acceptable tolerance for the root
 *	int maxit;			Max. iterations
 *
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition

 ************************************************************************
 */

#include <float.h>
#include <math.h>
#include <R_ext/Applic.h>
#include <R_ext/BLAS.h>

#define EPSILON DBL_EPSILON

double brent(		    	/* An estimate of the root */
    double ax,				/* Left border | of the range	*/
    double bx,				/* Right border| the root is seeked*/
    double fa, double fb,		/* f(a), f(b) */
    double f (double x, double *rw, int *iw),	/* Function under investigation	*/
    double *rw,
    int *iw,
    double tol,			  /* Acceptable tolerance		*/
    int maxit)				/* Max # of iterations */
{
    double a,b,c, fc;

    a = ax;  b = bx;
    c = a;   fc = fa;
    maxit = maxit + 1;

    /* First test if we have a root at an endpoint */
    if(fa == 0.0)  return a;
    if(fb == 0.0)  return b;

    /* Main iteration loop	*/
    while(maxit--)  {
 	    double prev_step = b-a;
	    double tol_act;			/* Actual tolerance		*/
    	double p;			      /* Interpolation step in the form p/q; */
	    double q;
	    double new_step;		/* Step at this iteration	*/

	    if( fabs(fc) < fabs(fb) ){				/* Swap data for b to be the	*/
	      a = b;  b = c;  c = a;  	      /* best approximation		*/
	      fa=fb;  fb=fc;  fc=fa;
    	}
	    tol_act = 2*EPSILON*fabs(b) + tol/2;
	    new_step = (c-b)/2;

    	if( fabs(new_step) <= tol_act || fb == (double)0 )  return b;

	    /* Decide if the interpolation can be tried	*/
	    if( fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb) ) {
	      register double t1,cb,t2;
	      cb = c-b;
	      if( a == c ) {    /* linear interpolation*/
      		t1 = fb/fa;
		      p = cb*t1;
		      q = 1.0 - t1;
	      }  else {			    /* Quadric inverse interpolation*/
	 	      q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
		      p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
		      q = (q-1.0) * (t1-1.0) * (t2-1.0);
	      }
	      if( p > (double)0 )
		       q = -q;
	      else
		       p = -p;

	      if( p < (0.75*cb*q-fabs(tol_act*q)/2)
		        && p < fabs(prev_step*q/2) )
		      new_step = p/q;
	    }

	    if( fabs(new_step) < tol_act) {	/* Adjust step to be not less than tol*/
	      if( new_step > (double)0 )
		      new_step = tol_act;
	      else
		      new_step = -tol_act;
	    }
	    a = b;	fa = fb;			/* Save the previous approx. */
    	b += new_step;	fb = f (b, rw, iw);
	    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
	      c = a;  fc = fa;  /* Adjust c to have a sign opposite to that of b */
	    }

    }
    /* failed! */
    return b;
}
