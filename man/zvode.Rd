\name{zvode}
\alias{zvode}                                                        

\title{Solver for Ordinary Differential Equations (ODE) for COMPLEX variables}

\description{
  Solves the initial value problem for stiff or nonstiff systems of 
  ordinary differential equations (ODE) in the form:

  \deqn{dy/dt = f(t,y)}
  
  where \eqn{dy} and \eqn{y} are complex variables. 
    
  The \R function \code{zvode} provides an interface to the FORTRAN ODE
  solver of the same name, written by Peter N. Brown, Alan C. Hindmarsh
  and George D. Byrne.
    
}

\usage{zvode(y, times, func, parms, rtol = 1e-6, atol = 1e-6,  
  jacfunc = NULL, jactype = "fullint", mf = NULL, verbose = FALSE,   
  tcrit = NULL, hmin = 0, hmax = NULL, hini = 0, ynames = TRUE,
  maxord = NULL, bandup = NULL, banddown = NULL, maxsteps = 5000,
  dllname = NULL, initfunc = dllname, initpar = parms, rpar = NULL,
  ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, ...)
}
\arguments{
  \item{y }{the initial (state) values for the ODE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix. \emph{y has to be complex}
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time; if only one step is
    to be taken; set \code{times = NULL}.
  }
  \item{func }{either an \R-function that computes the values of the
    derivatives in the ODE system (the \emph{model definition}) at time
    \code{t}, or a character string giving the name of a compiled function in a
    dynamically loaded shared library.
  
    If \code{func} is an \R-function, it must be defined as:
     \code{func <- function(t, y, parms, ...)}.  \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system.  If the initial values \code{y} has a
    \code{names} attribute, the names will be available inside \code{func}.
    \code{parms} is a vector or list of parameters; ... (optional) are
    any other arguments passed to the function.
  
    The return value of \code{func} should be a list, whose first
    element is a vector containing the derivatives of \code{y} with
    respect to \code{time}, and whose next elements are global values
    that are required at each point in \code{times}.  The derivatives
    must be specified in the \bold{same order} as the state variables \code{y}.
    They should be \emph{complex numbers}.
    
    If \code{func} is
    a string, then \code{dllname} must give the name of the shared
    library (without extension) which must be loaded before
    \code{zvode()} is called. See package vignette \code{"compiledCode"}
    for more details.
  }
  \item{parms }{vector or list of parameters used in \code{func} or
    \code{jacfunc}.
  }
  \item{rtol }{relative error tolerance, either a scalar or an array as
    long as \code{y}.  See details.
  }
  \item{atol }{absolute error tolerance, either a scalar or an array as
    long as \code{y}.  See details.
  }
  \item{jacfunc }{if not \code{NULL}, an \R function that computes the
    Jacobian of the system of differential equations
    \eqn{\partial\dot{y}_i/\partial y_j}{dydot(i)/dy(j)}, or
    a string giving the name of a function or subroutine in
    \file{dllname} that computes the Jacobian (see vignette
    \code{"compiledCode"} for more about this option).

    In some circumstances, supplying
    \code{jacfunc} can speed up the computations, if the system is
    stiff.  The \R calling sequence for \code{jacfunc} is identical to
    that of \code{func}.

    If the Jacobian is a full matrix, \code{jacfunc} should return a
    matrix \eqn{\dot{dy}/dy}{dydot/dy}, where the ith row contains the derivative of
    \eqn{dy_i/dt} with respect to \eqn{y_j}, or a vector containing the
    matrix elements by columns (the way \R and FORTRAN store matrices).
    Its elements should be \emph{complex numbers}.

    If the Jacobian is banded, \code{jacfunc} should return a matrix
    containing only the nonzero bands of the Jacobian, rotated
    row-wise. See first example of \code{lsode}.
  }
  \item{jactype }{the structure of the Jacobian, one of
    \code{"fullint"}, \code{"fullusr"}, \code{"bandusr"} or
    \code{"bandint"} - either full or banded and estimated internally or
    by user; overruled if \code{mf} is not \code{NULL}.
  }
  \item{mf }{the "method flag" passed to function \code{zvode} - overrules
    \code{jactype} - provides more options than \code{jactype} - see
    details.
  }
  \item{verbose }{if TRUE: full output to the screen, e.g. will
    print the \code{diagnostiscs} of the integration - see details.
  }
  \item{tcrit }{if not \code{NULL}, then \code{zvode} cannot integrate
    past \code{tcrit}. The FORTRAN routine \code{dvode} overshoots its
    targets (times points in the vector \code{times}), and interpolates
    values for the desired time points.  If there is a time beyond which
    integration should not proceed (perhaps because of a singularity),
    that should be provided in \code{tcrit}.
  }
  \item{hmin }{an optional minimum value of the integration stepsize. In
    special situations this parameter may speed up computations with the
    cost of precision. Don't use hmin if you don't know why!
  }
  \item{hmax }{an optional maximum value of the integration stepsize. If
    not specified, hmax is set to the largest difference in
    \code{times}, to avoid that the simulation possibly ignores
    short-term events. If 0, no maximal size is specified.
  }
  \item{hini }{initial step size to be attempted; if 0, the initial step
    size is determined by the solver.
  }
  \item{ynames }{logical; if \code{FALSE}: names of state variables are not
    passed to function \code{func} ; this may speed up the simulation
    especially for multi-D models.
  }
  \item{maxord }{the maximum order to be allowed. \code{NULL} uses the default,
    i.e. order 12 if implicit Adams method (\code{meth = 1}), order 5 if BDF
    method (\code{meth = 2}). Reduce maxord to save storage space.
  }
  \item{bandup }{number of non-zero bands above the diagonal, in case
    the Jacobian is banded.
  }
  \item{banddown }{number of non-zero bands below the diagonal, in case
    the Jacobian is banded.
  }
  \item{maxsteps }{maximal number of steps per output interval taken by the
    solver.
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions refered to in \code{func} and \code{jacfunc}.
    See package vignette \code{"compiledCode"}.
  }
  \item{initfunc }{if not \code{NULL}, the name of the initialisation function
    (which initialises values of parameters), as provided in
    \file{dllname}. See package vignette \code{"compiledCode"}.
  }
  \item{initpar }{only when \file{dllname} is specified and an
    initialisation function \code{initfunc} is in the dll: the
    parameters passed to the initialiser, to initialise the common
    blocks (FORTRAN) or global variables (C, C++).
  }
  \item{rpar }{only when \file{dllname} is specified: a vector with
    double precision values passed to the DLL-functions whose names are
    specified by \code{func} and \code{jacfunc}.
  }
  \item{ipar }{only when \file{dllname} is specified: a vector with
    integer values passed to the dll-functions whose names are specified
    by \code{func} and \code{jacfunc}.
  }
  \item{nout }{only used if \code{dllname} is specified and the model is
    defined in compiled code: the number of output variables calculated
    in the compiled function \code{func}, present in the shared
    library. Note: it is not automatically checked whether this is
    indeed the number of output variables calculated in the DLL - you have
    to perform this check in the code - See package vignette \code{"compiledCode"}.
  }
  \item{outnames }{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{func}, present in the shared library.
    These names will be used to label the output matrix.
  }
  \item{forcings }{only used if \file{dllname} is specified: a list with
    the forcing function data sets, each present as a two-columned matrix,
    with (time,value); interpolation outside the interval
    [min(\code{times}), max(\code{times})] is done by taking the value at
    the closest data extreme.

    See \link{forcings} or package vignette \code{"compiledCode"}.
  }
  \item{initforc }{if not \code{NULL}, the name of the forcing function
    initialisation function, as provided in
    \file{dllname}. It MUST be present if \code{forcings} has been given a
    value.
    See \link{forcings} or package vignette \code{"compiledCode"}.
  }
  \item{fcontrol }{A list of control parameters for the forcing functions.
    \link{forcings} or package vignette \code{"compiledCode"}
  }
  \item{... }{additional arguments passed to \code{func} and
    \code{jacfunc} allowing this to be a generic function.
  }
}
\value{
  A matrix of class \code{deSolve} with up to as many rows as elements
  in \code{times} and as many columns as elements in \code{y} plus the 
  number of "global" values returned in the next elements of the return 
  from \code{func},
  plus and additional column for the time value.  There will be a row
  for each element in \code{times} unless the FORTRAN routine `zvode'
  returns with an unrecoverable error. If \code{y} has a names
  attribute, it will be used to label the columns of the output value.
}
\author{Karline Soetaert <karline.soetaert@nioz.nl>}
\examples{

## =======================================================================
## Example 1 - very simple example 
## df/dt = 1i*f, where 1i is the imaginary unit
## The initial value is f(0) = 1 = 1+0i
## =======================================================================

ZODE <- function(Time, f, Pars) {
  df <-  1i*f
  return(list(df))
}

pars    <- NULL
yini    <- c(f = 1+0i)
times   <- seq(0, 2*pi, length = 100)
out     <- zvode(func = ZODE, y = yini, parms = pars, times = times,
  atol = 1e-10, rtol = 1e-10)

# The analytical solution to this ODE is the exp-function:
# f(t) = exp(1i*t)
#      = cos(t)+1i*sin(t)  (due to Euler's equation)

analytical.solution  <- exp(1i * times) 

## compare numerical and analytical solution
tail(cbind(out[,2], analytical.solution))


## =======================================================================
## Example 2 - example in "zvode.f",  
## df/dt = 1i*f        (same as above ODE)
## dg/dt = -1i*g*g*f   (an additional ODE depending on f)
##
## Initial values are
## g(0) = 1/2.1 and
## z(0) = 1  
## =======================================================================

ZODE2<-function(Time,State,Pars) {
  with(as.list(State), {
    df <- 1i * f
    dg <- -1i * g*g * f
    return(list(c(df, dg)))
  })
}

yini    <- c(f = 1 + 0i, g = 1/2.1 + 0i)
times   <- seq(0, 2*pi, length = 100)
out     <- zvode(func = ZODE2, y = yini, parms = NULL, times = times,
  atol = 1e-10, rtol = 1e-10)


## The analytical solution is
## f(t) = exp(1i*t)   (same as above)
## g(t) = 1/(f(t) + 1.1)

analytical <- cbind(f = exp(1i * times), g = 1/(exp(1i * times) + 1.1))

## compare numerical solution and the two analytical ones:
tail(cbind(out[,2], analytical[,1]))

}
\references{
  
  P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, 1989.  VODE: A Variable
  Coefficient ODE Solver, SIAM J. Sci. Stat. Comput., 10, pp. 1038-1051.
  \cr Also, LLNL Report UCRL-98412, June 1988.
  \doi{10.1137/0910062}

  G. D. Byrne and A. C. Hindmarsh, 1975. A Polyalgorithm for the
  Numerical Solution of Ordinary Differential Equations.  ACM
  Trans. Math. Software, 1, pp. 71-96.
  \doi{10.1145/355626.355636}

  A. C. Hindmarsh and G. D. Byrne, 1977. EPISODE: An Effective Package
  for the Integration of Systems of Ordinary Differential
  Equations. LLNL Report UCID-30112, Rev. 1.

  G. D. Byrne and A. C. Hindmarsh, 1976. EPISODEB: An Experimental
  Package for the Integration of Systems of Ordinary Differential
  Equations with Banded Jacobians. LLNL Report UCID-30132, April 1976.

  A. C. Hindmarsh, 1983. ODEPACK, a Systematized Collection of ODE
  Solvers. in Scientific Computing, R. S. Stepleman et al., eds.,
  North-Holland, Amsterdam, pp. 55-64.
  
  K. R. Jackson and R. Sacks-Davis, 1980. An Alternative Implementation
  of Variable Step-Size Multistep Formulas for Stiff ODEs.  ACM
  Trans. Math. Software, 6, pp. 295-318.
  \doi{10.1145/355900.355903}
  
  Netlib: \url{https://netlib.org}
}
\details{
    see \code{\link{vode}}, the double precision version, for details.
    }
\note{
From version 1.10.4, the default of atol was changed from 1e-8 to 1e-6,
    to be consistent with the other solvers.

The following text is adapted from the zvode.f source code:

When using \code{zvode} for a stiff system, it should only be used for
the case in which the function f is analytic, that is, when each f(i)
is an analytic function of each y(j).  Analyticity means that the
partial derivative df(i)/dy(j) is a unique complex number, and this
fact is critical in the way \code{zvode} solves the dense or banded linear
systems that arise in the stiff case.  For a complex stiff ODE system
in which f is not analytic, \code{zvode} is likely to have convergence
failures, and for this problem one should instead use \code{ode} on the
equivalent real system (in the real and imaginary parts of y).
}

\seealso{
\code{\link{vode}} for the double precision version
}

\keyword{math}

