\name{diagnostics}
\alias{diagnostics}
\alias{diagnostics.default}
\title{Print Diagnostic Characteristics of Solvers}
\description{
  Prints several diagnostics of the simulation to the screen, e.g.
  number of steps taken, the last step size, ...
}
\usage{
diagnostics(obj, ...)
\method{diagnostics}{default}(obj, ...)
}
\arguments{
  \item{obj}{is an output data structure produced by one of the
    solver routines.
  }
  \item{...}{optional arguments allowing to extend \code{diagnostics}
    as a generic function.
  }
}
\details{
  Detailed information obout the success of a simulation is printed,
  if a \code{diagnostics} function exists for a specific solver routine.
  A warning is printed, if no class-specific diagnostics exists.

  Please consult the class-specific help page for details.
}
\seealso{
  \code{\link{diagnostics.deSolve}} for diagnostics of differential
  equaton solvers.
%% enable this when bvpSolve is on CRAN
%  \code{\link[bvpSolve:diagnostics]{diagnostics.bvpSolve}} for
%  diagnostics of boundary value problem solvers.
}

\keyword{ utilities }