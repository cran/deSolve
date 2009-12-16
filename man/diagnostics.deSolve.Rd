\name{diagnostics.deSolve}
\alias{diagnostics.deSolve}
\title{Print Diagnostic Characteristics of ODE and DAE Solvers}
\description{
  Prints several diagnostics of the simulation to the screen, e.g.
  number of steps taken, the last step size, ...
}
\usage{
\method{diagnostics}{deSolve}(obj, Full = FALSE, ...)
}
\arguments{
  \item{obj}{is the output matrix as produced by one of the
    integration routines.
  }
  \item{Full}{when \code{TRUE} then all messages will be printed, including
    the ones that are not relevant for the solver. If \code{FALSE}, then only
    the relevant messages will be printed.
  }
  \item{...}{optional arguments allowing to extend \code{diagnostics}
    as a generic function.
  }
}
\value{
  The integer and real vector with diagnostic values; for function \code{lsodar}
  also the root information.
  
  See tables 2 and 3 in vignette("deSolve") for what these vectors contain.
}
\details{
  When the integration output is saved as a \code{data.frame}, then the required
  attributes are lost and method \code{diagnostics} will not work anymore.
}
\examples{
## The famous Lorenz equations: chaos in the earth's atmosphere
## Lorenz 1963. J. Atmos. Sci. 20, 130-141.

chaos <- function(t, state, parameters) {
  with(as.list(c(state)),{

    dx     <- -8/3*x+y*z
    dy     <- -10*(y-z)
    dz     <- -x*y+28*y-z

    list(c(dx, dy, dz))
  })
}

state <- c(x = 1, y = 1, z = 1)
times <- seq(0, 50, 0.01)
out   <- vode(state, times, chaos, 0)
pairs(out, pch=".")
diagnostics(out)
}
\keyword{ utilities }