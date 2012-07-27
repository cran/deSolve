\name{DLLfunc}
\alias{DLLfunc}
\title{Evaluates a Derivative Function Represented in a DLL}
\description{Calls a function, defined in a compiled language as a DLL}
\usage{DLLfunc(func, times, y, parms, dllname,
  initfunc = dllname, rpar = NULL, ipar = NULL, nout = 0,
  outnames = NULL, forcings = NULL, initforc = NULL,
  fcontrol = NULL)

}
\arguments{
  \item{func }{the name of the function in the dynamically loaded
    shared library,
  }
  \item{times }{first value = the time at which the function needs to be
    evaluated,
  }
  \item{y }{the values of the dependent variables for which the function
    needs to be evaluated,
  }
  \item{parms }{the parameters that are passed to the initialiser function,
  }
  \item{dllname }{a string giving the name of the shared library (without
    extension) that contains the compiled function or subroutine definitions
    referred to in \code{func},
  }
  \item{initfunc }{if not \code{NULL}, the name of the initialisation function
    (which initialises values of parameters), as provided in \file{dllname}.
    See details.
  }
  \item{rpar }{a vector with double precision values passed to the
    DLL-function \code{func} and \code{jacfunc} present in the DLL, via
    argument rpar,
  }
  \item{ipar }{a vector with integer values passed to the dll-function
    \code{func} and \code{jacfunc} present in the DLL, via function argument
    ipar,
  }
  \item{nout }{the number of output variables.
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

    See package vignette \code{"compiledCode"}.
  }
  \item{initforc }{if not \code{NULL}, the name of the forcing function
    initialisation function, as provided in
    \file{dllname}. It MUST be present if \code{forcings} has been given a
    value.
    See package vignette \code{"compiledCode"}.
  }
  \item{fcontrol }{A list of control parameters for the forcing functions.
    See package vignette \code{"compiledCode"}.
  }
}
\value{
  a list containing:
  \item{dy }{the rate of change estimated by the function,
  }
  \item{var }{the ordinary output variables of the function.
  }
}
\details{
  This function is meant to help developing FORTRAN or C models that are
  to be used to solve ordinary differential equations (ODE) in packages
  \code{deSolve} and/or \code{rootSolve}.
}
\author{Karline Soetaert <karline.soetaert@nioz.nl>}
\examples{
## ==========================================================================
## ex. 1
## ccl4model
## ==========================================================================
## Parameter values and initial conditions
## see example(ccl4model) for a more comprehensive implementation

Parms <- c(0.182, 4.0, 4.0, 0.08, 0.04, 0.74, 0.05, 0.15, 0.32,
        16.17, 281.48, 13.3, 16.17, 5.487, 153.8, 0.04321671,
        0.4027255, 1000, 0.02, 1.0, 3.8)

yini <- c(AI = 21, AAM = 0, AT = 0, AF = 0, AL = 0, CLT = 0,  AM = 0)

## the rate of change
DLLfunc(y = yini, dllname = "deSolve", func = "derivsccl4",
        initfunc = "initccl4", parms = Parms, times = 1,
        nout = 3, outnames = c("DOSE", "MASS", "CP")  )

## ==========================================================================
## ex. 2
## SCOC model, in fortran  - to see the FORTRAN code:
## ==========================================================================

## Forcing function "data"
Flux  <- matrix(ncol = 2, byrow = TRUE, data = c(1, 0.654, 2, 0.167))
parms <- c(k = 0.01)
Yini  <- 60

DLLfunc(y=Yini, times=1, func = "scocder",
    parms = parms, dllname = "deSolve",
    initforc = "scocforc",  forcings = Flux,
    initfunc = "scocpar", nout = 2,
    outnames = c("Mineralisation","Depo"))
## correct value = dy = flux - k * y = 0.654 - 0.01 * 60

DLLfunc(y = Yini, times = 2, func = "scocder",
    parms = parms, dllname = "deSolve",
    initforc = "scocforc",  forcings = Flux,
    initfunc = "scocpar", nout = 2,
    outnames = c("Mineralisation", "Depo"))
}
\keyword{utilities}

\seealso{
  \code{\link{ode}} for a general interface to most of the ODE solvers
}
