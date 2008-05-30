\name{rkInternal}
\alias{rkAuto}
\alias{rkFixed}

\title{Internal Functions for the Runge-Kutta Family of ODE Solvers.}

\description{
 These are the internal workhorses of \code{\link{rk}}, which are normally not
 intended to be called by the end user.

 Though internal, both functions are visible in the user's workspace to allow for
 easy inspection and fine-tuning of the algorithms by experienced users.
}

\usage{
rkFixed(y, times, func, parms, tcrit = NULL, verbose = FALSE, hini = 0,
  method = rkMethod("rk4", ...), ...)

rkAuto(y, times, func, parms, rtol = 1e-06, atol = 1e-06,
  tcrit = NULL, verbose = FALSE, hmin = 0, hmax = NULL, hini = 0,
  method = rkMethod("rk45dp7", ...), maxsteps = 5000, ...)
}

\details{
 Arguments and value are described in the manual page of the top-level function
 \code{\link{rk}}. Note that  \code{rkFixed} is used for the methods
 with fixed time step and \code{rkAuto} for the embedded methods with automatic step
 size determination. Both functions return raw values of the integration
 without external interpolation and without the "global" return values of the model.
}


\keyword{ internal }