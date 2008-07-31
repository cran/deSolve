\name{vode}
\alias{vode}
\title{General solver for ordinary differential equations (ODE)}
\description{Solves the initial value problem for stiff or nonstiff systems of 
    ordinary differential equations (ODE) in the form:
    \deqn{dy/dt = f(t,y)} 
    \cr The \R function \code{vode} provides an interface to the
    Fortran ODE solver of the same name, written by Peter N. Brown, Alan C. Hindmarsh and George D. Byrne.  
    
    The system of ODE's is written as an \R function or be defined in
    compiled code that has been dynamically loaded. 
  
    In contrast to \code{\link{lsoda}}, the user has to specify whether or not the problem is stiff
    and choose the appropriate solution method.
    
  \code{vode} is very similar to \code{\link{lsode}}, but uses a variable-coefficient method 
  rather than the fixed-step-interpolate methods in \code{\link{lsode}}.
  In addition, in vode it is possible to choose whether or not a copy
  of the Jacobian is saved for reuse in the corrector iteration algorithm; 
  In \code{lsode}, a copy is not kept. }
\usage{vode(y, times, func, parms, rtol=1e-6, atol=1e-8,  
  jacfunc=NULL, jactype="fullint", mf=NULL, verbose=FALSE,   
  tcrit=NULL, hmin=0, hmax=NULL, hini=0, ynames=TRUE, maxord=NULL, 
  bandup=NULL, banddown=NULL, maxsteps=5000, dllname=NULL, 
  initfunc=dllname, initpar=parms, rpar=NULL, 
  ipar=NULL, nout=0, outnames=NULL, ...)}
\arguments{
  \item{y }{the initial (state) values for the ODE system. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{times }{time sequence for which output is wanted; the first value of \code{times} must be the initial time; if only one step is to be taken; set \code{times} = NULL}
  \item{func }{either an \R-function that computes the values of the
      derivatives in the ODE system (the \emph{model definition}) at time
      t, or a character string giving the name of a compiled function in a dynamically loaded
      shared library.
  
      If \code{func}  is an \R-function, it must be defined as:
      \code{yprime = func(t, y, parms,...)}.  \code{t} is the current time point
      in the integration, \code{y} is the current estimate of the variables
      in the ODE system.  If the initial values \code{y} has a names
      attribute, the names will be available inside \code{func}.  \code{parms} is
      a vector or list of parameters; ... (optional) are any other arguments passed to the function. 
  
      The return value of \code{func} should be a list, whose first element is a
      vector containing the derivatives of \code{y} with respect to
      \code{time}, and whose next elements are global values that are required at
      each point in \code{times}.
          If \code{func} is a string, then \code{dllname} must give the name
      of the shared library (without extension) which must be loaded
      before \code{vode()} is called. See package vignette for more details.
      }
  \item{parms }{vector or list of parameters used in \code{func} or \code{jacfunc}.}
  \item{rtol }{relative error tolerance, either a scalar or an array as
      long as \code{y}.  See details. }
  \item{atol }{absolute error tolerance, either a scalar or an array as
      long as \code{y}.  See details.}
  \item{jacfunc }{if not \code{NULL}, an \R function that computes
      the jacobian of the system of differential equations
      dydot(i)/dy(j), or a string giving the name of a function or
      subroutine in \file{dllname} that computes the jacobian (see Details
      below for more about this option).  In some circumstances, supplying
      \code{jacfunc} can speed up the computations, if the system is stiff.  The \R calling sequence for
      \code{jacfunc} is identical to that of \code{func}.  
      \cr If the jacobian is a full matrix, \code{jacfunc} should return a matrix dydot/dy, where the ith
      row contains the derivative of \eqn{dy_i/dt} with respect to \eqn{y_j},
      or a vector containing the matrix elements by columns (the way \R and Fortran store matrices).
      \cr If the jacobian is banded, \code{jacfunc} should return a matrix containing only the 
      nonzero bands of the jacobian, rotated row-wise. See first example of lsode.}
  \item{jactype }{the structure of the jacobian, one of "fullint", "fullusr", "bandusr" or "bandint" - either full or banded and estimated internally or by user; overruled if \code{mf} is not NULL}
  \item{mf }{the "method flag" passed to function vode - overrules \code{jactype} - provides more options than \code{jactype} - see details}
  \item{verbose }{if TRUE: full output to the screen, e.g. will output the settings of vectors *istate* and *rstate* - see details}
  \item{tcrit }{if not \code{NULL}, then \code{vode} cannot integrate past \code{tcrit}. The Fortran routine \code{dvode} overshoots its targets (times points in the vector \code{times}), and interpolates values
      for the desired time points.  If there is a time beyond which integration should not proceed (perhaps because of a singularity),
      that should be provided in \code{tcrit}.}
  \item{hmin }{an optional minimum value of the integration stepsize. In special situations this parameter may speed up computations with the cost of precision. Don't use hmin if you don't know why!}
  \item{hmax }{an optional maximum value of the integration stepsize. If not specified, hmax is set to the largest difference in \code{times}, to avoid that the simulation possibly ignores short-term events. If 0, no maximal size is specified}
  \item{hini }{initial step size to be attempted; if 0, initial step size is determined by the solver}
  \item{ynames }{if FALSE: names of state variables are not passed to function \code{func} ; this may speed up the simulation especially for multi-D models}
  \item{maxord }{the maximum order to be allowed. NULL uses the default, i.e. order 12 if implicit Adams method (meth=1), order 5 if BDF method (meth=2). Reduce maxord to save storage space}
  \item{bandup }{number of non-zero bands above the diagonal, in case the Jacobian is banded}
  \item{banddown }{number of non-zero bands below the diagonal, in case the Jacobian is banded}
  \item{maxsteps }{maximal number of steps during one call to the solver}
  \item{dllname }{a string giving the name of the shared library (without extension) that contains all the compiled function or subroutine definitions refered to in \code{func} and \code{jacfunc}.  See package vignette}
  \item{initfunc }{if not NULL, the name of the initialisation function (which initialises values of parameters), as provided in \file{dllname}. See package vignette}
  \item{initpar }{only when \file{dllname} is specified and an initialisation function \code{initfunc} is in the dll: the parameters passed to the initialiser, to initialise the common blocks (fortran) or global variables (C, C++)}
  \item{rpar }{only when \file{dllname} is specified: a vector with double precision values passed to the dll-functions whose names are specified by \code{func} and \code{jacfunc}}
  \item{ipar }{only when \file{dllname} is specified: a vector with integer values passed to the dll-functions whose names are specified by \code{func} and \code{jacfunc}}
  \item{nout }{only used if \code{dllname} is specified and the model is defined in compiled code: the number of output variables calculated in the compiled function \code{func}, present in the shared library. Note:
      it is not automatically checked whether this is indeed the number of output variables calculed in the dll - you have to perform this check in the code - See package vignette}
  \item{outnames }{only used if \file{dllname} is specified and \code{nout} > 0: the names of output variables calculated in the compiled function \code{func}, present in the shared library}
  \item{... }{additional arguments passed to \code{func} and \code{jacfunc} allowing this to be a generic function}
}
\value{
  A matrix with up to as many rows as elements in times and as many columns as elements in \code{y} plus the number of "global" values returned 
      in the next elements of the return from \code{func}, plus an additional column (the first) for the time value. 
      There will be one row for each element in \code{times} unless the Fortran routine `vode' returns with an unrecoverable error. 
      If \code{y} has a names attribute, it will be used to label the columns of the output value. 
      
      The output will have the attributes \code{istate}, and \code{rstate}, two vectors with several useful elements. 
      See details.
      The first element of istate returns the conditions under which the last call to lsoda returned. Normal is \code{istate[1] = 2}.
      If \code{verbose} = TRUE, the settings of istate and rstate will be written to the screen

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# The famous Lorenz equations: chaos in the earth's atmosphere
# Lorenz 1963. J. Atmos. Sci. 20, 130-141.

chaos<-function(t,state,parameters)
  {
  with(as.list(c(state)),{

    dx     <- -8/3*x+y*z
    dy     <- -10*(y-z)
    dz     <- -x*y+28*y-z

    list(c(dx,dy,dz))            })

 }  # end of model

state <-c(x=1, y=1, z=1)
times <-seq(0,100,0.01)
out   <-as.data.frame(vode(state,times,chaos,0))

plot(out$x,out$y,type="l",main="Lorenz butterfly")
}
\references{\item P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, 1989. 
    VODE: A Variable Coefficient ODE Solver," SIAM J. Sci. Stat. Comput., 10, pp. 1038-1051.  
    \cr Also, LLNL Report UCRL-98412, June 1988.
    \item G. D. Byrne and A. C. Hindmarsh, 1975. A Polyalgorithm for the Numerical Solution of Ordinary Differential Equations.
      ACM Trans. Math. Software, 1, pp. 71-96.
    \item A. C. Hindmarsh and G. D. Byrne, 1977. EPISODE: An Effective Package for the Integration of Systems of Ordinary Differential
   Equations. LLNL Report UCID-30112, Rev. 1.
    \item G. D. Byrne and A. C. Hindmarsh, 1976. EPISODEB: An Experimental Package for the Integration of Systems of Ordinary Differential
   Equations with Banded Jacobians. LLNL Report UCID-30132, April 1976.
    \item A. C. Hindmarsh, 1983. ODEPACK, a Systematized Collection of ODE Solvers. in Scientific Computing, R. S. Stepleman et al., eds.,
    North-Holland, Amsterdam, pp. 55-64.
   \item K. R. Jackson and R. Sacks-Davis, 1980. An Alternative Implementation of Variable Step-Size Multistep Formulas for Stiff ODEs.
    ACM Trans. Math. Software, 6, pp. 295-318.
    
      Netlib: \url{http://www.netlib.org}}
\details{Before using the integrator \code{vode}, the user has to decide whether or not the problem is stiff. 
  
  If the problem is nonstiff, use method flag \code{mf} = 10, which selects a nonstiff (Adams) method, no Jacobian used.  
  \cr If the problem is stiff, there are four standard choices 
     which can be specified with \code{jactype} or \code{mf}.
  
  The options for \bold{jactype} are
      \item jac = "fullint" : a full jacobian, calculated internally by vode, corresponds to \code{mf}=22
      \item jac = "fullusr" : a full jacobian, specified by user function \code{jacfunc}, corresponds to \code{mf}=21 
      \item jac = "bandusr" : a banded jacobian, specified by user function \code{jacfunc}; the size of the bands specified by \code{bandup} and \code{banddown}, corresponds to \code{mf}=24
      \item jac = "bandint" : a banded jacobian, calculated by vode; the size of the bands specified by \code{bandup} and \code{banddown}, corresponds to \code{mf}=25
  
  More options are available when specifying \bold{mf} directly. 
     \cr The legal values of \code{mf} are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25, -11, -12, -14, -15, -21, -22, -24, -25.
     \cr \code{mf} is a signed two-digit integer, \code{mf} = JSV*(10*METH + MITER),where
  \item JSV = SIGN(mf) indicates the Jacobian-saving strategy: JSV =  1 means a copy of the Jacobian is saved for reuse
  in the corrector iteration algorithm. JSV = -1 means a copy of the Jacobian is not saved.
  \item METH indicates the basic linear multistep method: METH = 1 means the implicit Adams method. METH = 2 means the method based on backward
  differentiation formulas (BDF-s).
  \item MITER indicates the corrector iteration method: MITER = 0 means functional iteration (no Jacobian matrix is involved).
   MITER = 1 means chord iteration with a user-supplied full (NEQ by NEQ) Jacobian.
   MITER = 2 means chord iteration with an internally generated (difference quotient) full Jacobian
   (using NEQ extra calls to \code{func} per df/dy value).
  MITER = 3 means chord iteration with an internally generated diagonal Jacobian approximation
   (using 1 extra call to \code{func} per df/dy evaluation).
  MITER = 4 means chord iteration with a user-supplied banded Jacobian.
  MITER = 5 means chord iteration with an internally generated banded Jacobian (using ML+MU+1 extra calls to \code{func} per df/dy evaluation).
  
  If MITER = 1 or 4, the user must supply a subroutine \code{jacfunc}.
  
  The example for integrator \code{\link{lsode}} demonstrates how to specify both a banded and full jacobian.
  
  The input parameters \code{rtol}, and \code{atol} determine the \bold{error
  control} performed by the solver. 
  If the request for precision exceeds the capabilities of the machine, 
  vode will return an error code. See \code{\link{lsoda}} for details.
  
  \bold{Models} may be defined in compiled C or Fortran code, as well as in an R-function. See package vignette for details.
  
  
   The output will have the \bold{attributes *istate*, and *rstate*}, two vectors with several useful elements. 
      \cr if \code{verbose} = TRUE, the settings of istate and rstate will be written to the screen. 
      
      the following elements of \bold{istate} are meaningful:    
      \item el 1 : returns the conditions under which the last call to vode returned. 
      2  if DVODE was successful, -1 if excess work done, -2 means excess accuracy requested. (Tolerances too small), 
  -3 means illegal input detected. (See printed message.), -4 means repeated error test failures. (Check all input),
  -5 means repeated convergence failures. (Perhaps bad Jacobian supplied or wrong choice of MF or tolerances.),
  -6 means error weight became zero during problem. (Solution component i vanished, and atol or atol(i) = 0.)
      \item el 12 : The number of steps taken for the problem so far.
      \item el 13 : The number of function evaluations for the problem so far.",
      \item el 14 : The number of Jacobian evaluations so far.,
      \item el 15 : The method order last used (successfully).,
      \item el 16 : The order to be attempted on the next step.,
      \item el 17 : if el 1 =-4,-5: the largest component in the error vector,
      \item el 20 : The number of matrix LU decompositions so far.,
      \item el 21 : The number of nonlinear (Newton) iterations so far.,
      \item el 22 : The number of convergence failures of the solver so far ,
      \item el 23 : The number of error test failures of the integrator so far.)
  
      \bold{rstate} contains the following:
      \item 1: The step size in t last used (successfully).
      \item 2: The step size to be attempted on the next step.
      \item 3: The current value of the independent variable which the solver has actually reached, i.e. the
  current internal mesh point in t.  
      \item 4: A tolerance scale factor, greater than 1.0, computed when a request for too much accuracy was detected.
      
      For more information, see the comments in the original code dvode.f}
\seealso{\code{\link{ode}}, \code{\link{lsoda}}, \code{\link{lsode}}, \code{\link{lsodes}}, 
    \code{\link{lsodar}}, \code{\link{daspk}}, \code{\link{rk}}.}
\keyword{math}

