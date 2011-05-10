\name{DLLres}
\alias{DLLres}
\title{Evaluates a Residual Derivative Function Represented in a DLL
}
\description{
  Calls a residual function, \eqn{F(t,y,y')} of a DAE system
  (differential algebraic equations) defined in a compiled language as a
  DLL.
  
  To be used for testing the implementation of DAE problems in compiled code
}
\usage{DLLres(res, times, y, dy, parms, dllname,
  initfunc = dllname, rpar = NULL, ipar = NULL, nout = 0,
  outnames = NULL, forcings = NULL, initforc = NULL,
  fcontrol = NULL)
}
\arguments{
  \item{res }{the name of the function in the dynamically loaded
    shared library,
  }
  \item{times }{first value = the time at which the function needs to be
    evaluated,
  }
  \item{y }{the values of the dependent variables for which the function needs
    to be evaluated,
  }
  \item{dy }{the derivative of the values of the dependent variables for
    which the function needs to be evaluated,
  }
  \item{parms }{the parameters that are passed to the initialiser function,
  }
  \item{dllname }{a string giving the name of the shared library (without
    extension) that contains the compiled function or subroutine definitions
    referred to in \code{func},
  }
  \item{initfunc }{if not NULL, the name of the initialisation function
    (which initialises values of parameters), as provided in \file{dllname}.
    See details,
  }
  \item{rpar }{a vector with double precision values passed to the
    DLL-function \code{func} and \code{jacfunc} present in the DLL, via
    argument \code{rpar},
  }
  \item{ipar }{a vector with integer values passed to the DLL-function
    \code{func} and \code{jacfunc} present in the DLL, via function argument
    \code{ipar},
  }
  \item{nout }{the number of output variables.
  }
  \item{outnames }{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{func}, present in the shared library.
  }
  \item{forcings }{only used if \file{dllname} is specified: a list with
    the forcing function data sets, each present as a two-columned matrix,
    with (time,value); interpolation outside the interval
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
  \item{res }{the residual of derivative estimated by the function}
  \item{var }{the ordinary output variables of the function}
}
\details{
This function is meant to help developing FORTRAN or C models that are to be
used to solve differential algebraic equations (DAE) in
package \code{deSolve}.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\keyword{utilities}

\examples{
## =========================================================================
## Residuals from the daspk chemical model, production a forcing function
## =========================================================================
## Parameter values and initial conditions
## see example(daspk) for a more comprehensive implementation

pars <- c(K = 1, ka = 1e6, r = 1)

## Initial conc; D is in equilibrium with A,B
y <- c(A = 2, B = 3, D = 2 * 3/pars["K"])

## Initial rate of change
dy <- c(dA = 0, dB = 0, dD = 0)

## production increases with time
prod <- matrix(ncol = 2,
          data = c(seq(0, 100, by = 10), seq(0.1, 0.5, len = 11)))

DLLres(y = y, dy = dy, times = 5, res = "chemres",
          dllname = "deSolve", initfunc = "initparms",
          initforc = "initforcs", parms = pars, forcings = prod,
          nout = 2, outnames = c("CONC", "Prod"))

}
\seealso{
 \link{daspk} to solve DAE problems
}
