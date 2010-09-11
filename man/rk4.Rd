\name{rk4}
\alias{rk4}
\alias{euler}
\title{Solve System of ODE (Ordinary Differential Equation)s by
  Euler's Method or Classical Runge-Kutta 4th Order Integration.
}
\description{Solving initial value problems for systems of first-order
  ordinary differential equations (ODEs) using Euler's method or the
  classical Runge-Kutta 4th order integration.
}
\usage{
euler(y, times, func, parms, verbose = FALSE, ynames = TRUE,
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL, ipar = NULL, nout = 0, outnames = NULL,
  forcings = NULL, initforc = NULL, fcontrol = NULL, ...)
rk4(y, times, func, parms, verbose = FALSE, ynames = TRUE,
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL, ipar = NULL, nout = 0, outnames = NULL,
  forcings = NULL, initforc = NULL, fcontrol = NULL, ...)
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
     \code{func <- function(t, y, parms,...)}.  \code{t} is the current
    time point in the integration, \code{y} is the current estimate of
    the variables in the ODE system.  If the initial values \code{y} has
    a \code{names} attribute, the names will be available inside \code{func}.
    \code{parms} is a vector or list of parameters; ... (optional) are
    any other arguments passed to the function.

    The return value of \code{func} should be a list, whose first
    element is a vector containing the derivatives of \code{y} with
    respect to \code{time}, and whose next elements are global values
    that are required at each point in \code{times}.  The derivatives
    should be specified in the same order as the state variables \code{y}.
    
    If \code{func} is a string, then \code{dllname} must give the name
    of the shared library (without extension) which must be loaded
    before \code{rk4} is called. See package vignette
    \code{"compiledCode"} for more details.
  }
  \item{parms }{vector or list of parameters used in \code{func}.
  }
  \item{verbose }{a logical value that, when \code{TRUE}, triggers more
    verbose output from the ODE solver.
  }
  \item{ynames }{if \code{FALSE}: names of state variables are not passed
    to function \code{func} ; this may speed up the simulation especially
    for large models.
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions refered to in \code{func}.
    See package vignette \code{"compiledCode"}.
  }
  \item{initfunc }{if not \code{NULL}, the name of the initialisation function
    (which initialises values of parameters), as provided in
    \file{dllname}. See package vignette \code{"compiledCode"},
  }
  \item{initpar }{only when \file{dllname} is specified and an
    initialisation function \code{initfunc} is in the DLL: the
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
    to perform this check in the code. See package vignette
    \code{"compiledCode"}.
  }
  \item{outnames }{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{func}, present in the shared library.
  }
  \item{forcings }{only used if \file{dllname} is specified: a list with
    the forcing function data sets, each present as a two-columned matrix,
    with (time, value); interpolation outside the interval
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
   \item{... }{additional arguments passed to \code{func} allowing this
    to be a generic function.
  }
}
\author{Thomas Petzoldt \email{thomas.petzoldt@tu-dresden.de}}
\details{
  \code{rk4} and \code{euler} are special versions of the two fixed step
  solvers with less overhead and less functionality (e.g. no interpolation
  and no events) compared to the generic Runge-Kutta codes called by 
  \code{\link{rk}}.
   
  If you need different internal and external time steps or want to use events, 
  please use:
  \code{rk(y, times, func, parms, method = "rk4")} or
  \code{rk(y, times, func, parms, method = "euler")}.

  See help pages of \code{\link{rk}} and \code{\link{rkMethod}}
  for details.
}
\note{
  For most practical cases, solvers with flexible timestep
  (e.g. \code{rk(method = "ode45")} and especially solvers of the
  Livermore family (ODEPACK, e.g. \code{\link{lsoda}}) are superior.
}

\value{
  A matrix of class \code{deSolve} with up to as many rows as elements
  in \code{times} and as many columns as elements in \code{y} plus the
  number of "global" values returned in the next elements of the return
  from \code{func}, plus and additional column for the time value.
  There will be a row for each element in \code{times} unless the
  integration routine returns with an unrecoverable error. If \code{y}
  has a names attribute, it will be used to label the columns of the
  output value.
}

\examples{
## =======================================================================
## Example: Analytical and numerical solutions of logistic growth
## =======================================================================

## the derivative of the logistic
logist <- function(t, x, parms) {
  with(as.list(parms), {
    dx <- r * x[1] * (1 - x[1]/K)
    list(dx)
  })
}

time  <- 0:100
N0    <- 0.1; r <- 0.5; K <- 100
parms <- c(r = r, K = K)
x <- c(N = N0)

## analytical solution
plot(time, K/(1 + (K/N0-1) * exp(-r*time)), ylim = c(0, 120),
  type = "l", col = "red", lwd = 2)

## reasonable numerical solution
time <- seq(0, 100, 2)
out <- as.data.frame(rk4(x, time, logist, parms))
points(out$time, out$N, pch = 16, col = "blue", cex = 0.5)

## same time step, systematic under-estimation
time <- seq(0, 100, 2)
out <- as.data.frame(euler(x, time, logist, parms))
points(out$time, out$N, pch = 1)

## unstable result
time <- seq(0, 100, 4)
out <- as.data.frame(euler(x, time, logist, parms))
points(out$time, out$N, pch = 8, cex = 0.5)

## method with automatic time step
out <- as.data.frame(lsoda(x, time, logist, parms))
points(out$time, out$N, pch = 1, col = "green")

legend("bottomright",
  c("analytical","rk4, h=2", "euler, h=2",
    "euler, h=4", "lsoda"),
  lty = c(1, NA, NA, NA, NA), lwd = c(2, 1, 1, 1, 1),
  pch = c(NA, 16, 1, 8, 1),
  col = c("red", "blue", "black", "black", "green"))
}
\seealso{
  \itemize{
    \item \code{\link{rkMethod}} for a list of available Runge-Kutta
      parameter sets,
    \item \code{\link{rk}} for the more general Runge-Code,
    \item \code{\link{lsoda}}, \code{\link{lsode}},
      \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}},
      \code{\link{daspk}} for solvers of the Livermore family,
    \item \code{\link{ode}} for a general interface to most of the ODE solvers,
    \item \code{\link{ode.band}} for solving models with a banded
      Jacobian,
    \item \code{\link{ode.1D}} for integrating 1-D models,
    \item \code{\link{ode.2D}} for integrating 2-D models,
    \item \code{\link{ode.3D}} for integrating 3-D models,
    \item \code{\link{dede}} for integrating models with delay
      differential equations,
  }
  \code{\link{diagnostics}} to print diagnostic messages.
}
\keyword{math}

