\name{lsodar}
\alias{lsodar}
\title{Solver for Ordinary Differential Equations (ODE),
  Switching Automatically Between Stiff and Non-stiff Methods and With
  Root Finding
}
\description{Solving initial value problems for stiff or non-stiff
  systems of first-order ordinary differential equations (ODEs) and
  including root-finding.

  The \R function \code{lsodar} provides an interface to the FORTRAN ODE
  solver of the same name, written by Alan C. Hindmarsh and Linda
  R. Petzold.
  
  The system of ODE's is written as an \R function or be defined in
  compiled code that has been dynamically loaded. - see description of
  \code{\link{lsoda}} for details.

  \code{lsodar} differs from \code{lsode} in two respects.
  \itemize{
    \item It switches automatically between stiff and nonstiff methods
      (similar as lsoda).
    \item It finds the root of at least one of a set of constraint
      functions g(i) of the independent and dependent variables.
  }
}
\usage{lsodar(y, times, func, parms, rtol = 1e-6, atol = 1e-6, 
  jacfunc = NULL, jactype = "fullint", rootfunc = NULL,
  verbose = FALSE, nroot = 0, tcrit = NULL, hmin = 0,
  hmax = NULL, hini = 0, ynames = TRUE, maxordn = 12,
  maxords = 5, bandup = NULL, banddown = NULL, maxsteps = 5000,
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, forcings=NULL,
  initforc = NULL, fcontrol=NULL, ...)
}
\arguments{
  \item{y }{the initial (state) values for the ODE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{times }{times at which explicit estimates for \code{y} are
    desired.  The first value in \code{times} must be the initial time.
  }
  \item{func }{either an \R-function that computes the values of the
    derivatives in the ODE system (the \emph{model definition}) at time
    t, or a character string giving the name of a compiled function in a
    dynamically loaded shared library.

    If \code{func} is an \R-function, it must be defined as:
     \code{func <- function(t, y, parms,...)}.  \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system.  If the initial values \code{y} has a
    \code{names} attribute, the names will be available inside \code{func}.
    \code{parms} is a vector or list of parameters; ... (optional) are
    any other arguments passed to the function.

    The return value of \code{func} should be a list, whose first
    element is a vector containing the derivatives of \code{y} with
    respect to \code{time}, and whose next elements are global values
    that are required at each point in \code{times}.  The derivatives
    should be specified in the same order as the state variables \code{y}.

    If \code{func} is a string, then \code{dllname} must give the name
    of the shared library (without extension) which must be loaded
    before \code{lsodar()} is called. See package vignette
    \code{"compiledCode"} for more
    details.
  }
  \item{parms }{vector or list of parameters used in \code{func} or
    \code{jacfunc}.
  }
  \item{rtol }{relative error tolerance, either a scalar or an array as
    long as \code{y}. See details.
  }
  \item{atol }{absolute error tolerance, either a scalar or an array as
    long as \code{y}. See details.
  }
  \item{jacfunc }{if not \code{NULL}, an \R function, that computes the
    Jacobian of the system of differential equations dydot(i)/dy(j), or
    a string giving the name of a function or subroutine in
    \file{dllname} that computes the Jacobian (see vignette
    \code{"compiledCode"} for more about this option).

    In some circumstances, supplying
    \code{jacfunc} can speed up the computations, if the system is
    stiff.  The \R calling sequence for \code{jacfunc} is identical to
    that of \code{func}.

    If the Jacobian is a full matrix, \code{jacfunc} should return a
    matrix dydot/dy, where the ith row contains the derivative of
    \eqn{dy_i/dt} with respect to \eqn{y_j}, or a vector containing the
    matrix elements by columns (the way \R and FORTRAN store matrices).
    
    If the Jacobian is banded, \code{jacfunc} should return a matrix
    containing only the nonzero bands of the Jacobian, rotated
    row-wise. See first example of \code{lsode}.
  }
  \item{jactype }{the structure of the Jacobian, one of
    \code{"fullint"}, \code{"fullusr"}, \code{"bandusr"} or
    \code{"bandint"} - either full or banded and estimated internally or
    by user.
  }
  \item{rootfunc }{if not \code{NULL}, an \R function that computes the
    function whose root has to be estimated or a string giving the name
    of a function or subroutine in \file{dllname} that computes the root
    function.  The \R calling sequence for \code{rootfunc} is identical
    to that of \code{func}.  \code{rootfunc} should return a vector with
    the function values whose root is sought.
  }
  \item{verbose }{a logical value that, when \code{TRUE},  will
    print the \code{diagnostiscs} of the integration - see details.
  }
  \item{nroot }{only used if \file{dllname}  is specified: the number of
    constraint functions whose roots are desired during the integration;
    if \code{rootfunc} is an R-function, the solver estimates the number
    of roots.
  }
  \item{tcrit }{if not \code{NULL}, then \code{lsodar} cannot integrate
    past \code{tcrit}. The FORTRAN routine \code{lsodar} overshoots its
    targets (times points in the vector \code{times}), and interpolates
    values for the desired time points.  If there is a time beyond which
    integration should not proceed (perhaps because of a singularity),
    that should be provided in \code{tcrit}.
  }
  \item{hmin }{an optional minimum value of the integration stepsize. In
    special situations this parameter may speed up computations with the
    cost of precision. Don't use \code{hmin} if you don't know why!
  }
  \item{hmax }{an optional maximum value of the integration stepsize. If
    not specified, \code{hmax} is set to the largest difference in
    \code{times}, to avoid that the simulation possibly ignores
    short-term events. If 0, no maximal size is specified.
  }
  \item{hini }{initial step size to be attempted; if 0, the initial step
    size is determined by the solver.
  }

  \item{ynames }{logical, if \code{FALSE}: names of state variables are not
    passed to function \code{func}; this may speed up the simulation
    especially for large models.
  }
  \item{maxordn }{the maximum order to be allowed in case the method is
    non-stiff. Should be <= 12. Reduce \code{maxord} to save storage space.
  }
  \item{maxords }{the maximum order to be allowed in case the method is
    stiff. Should be <= 5. Reduce maxord to save storage space.
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
    subroutine definitions refered to in \code{func} and
    \code{jacfunc}. See package vignette \code{"compiledCode"}.
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
    double precision values passed to the dll-functions whose names are
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
    indeed the number of output variables calculed in the dll - you have
    to perform this check in the code - See package vignette
    \code{"compiledCode"}.
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
    See \link{forcings} or vignette \code{compiledCode}.
  }
  \item{... }{additional arguments passed to \code{func} and
    \code{jacfunc} allowing this to be a generic function.
  }
}
\value{
  A matrix of class \code{deSolve} with up to as many rows as elements
  in \code{times} and as many columns as elements in \code{y} plus the number of "global"
  values returned in the next elements of the return from \code{func},
  plus and additional column for the time value.  There will be a row
  for each element in \code{times} unless the FORTRAN routine `lsoda'
  returns with an unrecoverable error. If \code{y} has a names
  attribute, it will be used to label the columns of the output value.

  If a root has been found, the output will have the attribute
  \code{iroot}, an integer indicating which root has been found.
  
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
## =======================================================================
## Example 1:
##   from lsodar source code
## =======================================================================

Fun <- function (t, y, parms)
{
  ydot <- vector(len = 3)
  ydot[1] <- -.04*y[1] + 1.e4*y[2]*y[3]
  ydot[3] <- 3.e7*y[2]*y[2]
  ydot[2] <- -ydot[1]-ydot[3]
  return(list(ydot,ytot = sum(y)))
}

rootFun <- function (t, y, parms)
{
  yroot <- vector(len = 2)
  yroot[1] <- y[1] - 1.e-4
  yroot[2] <- y[3] - 1.e-2
  return(yroot)
}

y     <- c(1, 0, 0)
times <- c(0, 0.4*10^(0:8))
Out   <- NULL
ny    <- length(y)

out   <- lsodar(y = y, times = times, fun = Fun, rootfun = rootFun,
                rtol = 1e-4, atol = c(1e-6, 1e-10, 1e-6), parms = NULL)
print(paste("root is found for eqn", which(attributes(out)$iroot == 1)))
print(out[nrow(out),])

diagnostics(out)
  
## =======================================================================
## Example 2:
##   using lsodar to estimate steady-state conditions
## =======================================================================

## Bacteria (Bac) are growing on a substrate (Sub)
model <- function(t, state, pars)
{
  with (as.list(c(state, pars)), {
    ##       substrate uptake             death  respiration
    dBact =  gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact
    dSub  = -gmax    *Sub/(Sub+ks)*Bact + dB*Bact            + input

    return(list(c(dBact,dSub)))
  })
}

## root is the condition where sum of |rates of change|
## is very small

rootfun <- function (t, state, pars)
{
  dstate <- unlist(model(t, state, pars)) # rate of change vector
  return(sum(abs(dstate)) - 1e-10)
}

pars <- list(Bini = 0.1, Sini = 100, gmax = 0.5, eff = 0.5,
             ks = 0.5, rB = 0.01, dB = 0.01, input = 0.1)

tout    <- c(0, 1e10)
state   <- c(Bact = pars$Bini, Sub = pars$Sini)
out     <- lsodar(state, tout, model, pars, rootfun = rootfun)
print(out)
}
\references{

  Alan C. Hindmarsh, ODEPACK, A Systematized Collection of ODE Solvers,
  in Scientific Computing, R. S. Stepleman et al. (Eds.), North-Holland,
  Amsterdam, 1983, pp. 55-64.
    
  Linda R. Petzold, Automatic Selection of Methods for Solving Stiff and
  Nonstiff Systems of Ordinary Differential Equations, Siam
  J. Sci. Stat. Comput. 4 (1983), pp. 136-148.

 Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined Output
 Points for Solutions of ODEs, Sandia Report SAND80-0180, February 1980.
    
 Netlib: \url{http://www.netlib.org}
 
}
\details{
  The work is done by the FORTRAN subroutine \code{lsodar}, whose
  documentation should be consulted for details (it is included as
  comments in the source file \file{src/opkdmain.f}).  The
  implementation is based on the November, 2003 version of lsodar, from
  Netlib.
    
  \code{lsodar} switches automatically between stiff and nonstiff
  methods (similar as \code{lsoda}).  This means that the user does not
  have to determine whether the problem is stiff or not, and the solver
  will automatically choose the appropriate method.  It always starts
  with the nonstiff method.
  
  It finds the root of at least one of a set of constraint functions
  g(i) of the independent and dependent variables.  It then returns the
  solution at the root if that occurs sooner than the specified stop
  condition, and otherwise returns the solution according the specified
  stop condition.
  
  The form of the \bold{Jacobian} can be specified by \code{jactype}
  which can take the following values:
  \describe{
    \item{jactype = "fullint":}{a full Jacobian, calculated internally
      by lsodar, the default,
    }
    \item{jactype = "fullusr":}{a full Jacobian, specified by user
      function \code{jacfunc},
    }
    \item{jactype = "bandusr":}{a banded Jacobian, specified by user
      function \code{jacfunc}; the size of the bands specified by
      \code{bandup} and \code{banddown},
    }
    \item{jactype = "bandint":}{banded Jacobian, calculated by lsodar;
      the size of the bands specified by \code{bandup} and
      \code{banddown}.
    }
  }
  If \code{jactype} = "fullusr" or "bandusr" then the user must supply a
  subroutine \code{jacfunc}.

  The input parameters \code{rtol}, and \code{atol} determine the
  \bold{error control} performed by the solver. See \code{\link{lsoda}}
  for details.
  
  The output will have the attribute \bold{iroot}, if a root was found
  \bold{iroot} is a vector, its length equal to the number of constraint
  functions it will have a value of 1 for the constraint function whose
  root that has been found and 0 otherwise.

  The diagnostics of the integration can be printed to screen
  by calling \code{\link{diagnostics}}. If \code{verbose} = \code{TRUE},
  the diagnostics will written to the screen at the end of the integration.

  See vignette("deSolve") for an explanation of each element in the vectors
  containing the diagnostic properties and how to directly access them.

  \bold{Models} may be defined in compiled C or FORTRAN code, as well as
  in an R-function. See package vignette \code{"compiledCode"} for details.

  More information about models defined in compiled code is in the package
  vignette ("compiledCode"); information about linking forcing functions
  to compiled code is in \link{forcings}.

  Examples in both C and FORTRAN are in the \file{dynload} subdirectory
  of the \code{deSolve} package directory.
}
\seealso{
  \itemize{
    \item \code{\link{rk}}, \item \code{\link{rk4}} and \code{\link{euler}} for
      Runge-Kutta integrators.
    \item \code{\link{lsoda}}, \code{\link{lsode}},
      \code{\link{lsodes}}, \code{\link{vode}},
      \code{\link{daspk}} for other solvers of the Livermore family,
    \item \code{\link{ode}} for a general interface to most of the ODE solvers,
    \item \code{\link{ode.band}} for solving models with a banded
      Jacobian,
    \item \code{\link{ode.1D}} for integrating 1-D models,
    \item \code{\link{ode.2D}} for integrating 2-D models,
    \item \code{\link{ode.3D}} for integrating 3-D models,
  }

  \code{\link{diagnostics}} to print diagnostic messages.
}
\keyword{math}

