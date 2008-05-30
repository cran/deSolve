\name{rk4}
\alias{rk4}
\alias{euler}
\title{Solve System of ODE (Ordinary Differential Equation)s by
    Euler's Method or Classical Runge-Kutta 4th Order Integration.}
\description{Solving initial value problems for 
    systems of first-order ordinary differential equations
    (ODEs) using the classical Runge-Kutta 4th order integration.
    }
\usage{
rk4(y, times, func, parms, hini=min(diff(times)), ...)
euler(y, times, func, parms, hini=min(diff(times)), ...)
}
\arguments{
  \item{y }{the initial values for the ODE system.  If \code{y} has a
      name attribute, the names will be used to label the output matrix.}
  \item{times }{times at which explicit estimates for \code{y} are
      desired.  The first value in \code{times} must be the initial time.}
  \item{func }{an \R-function that computes the values of the
      derivatives in the ode system (the \emph{model defininition}) at time
      t.
      
      The \R-function  \code{func} must be defined as:
      \code{yprime = func(t, y, parms,...)}.  \code{t} is the current time point
      in the integration, \code{y} is the current estimate of the variables
      in the ode system, and \code{parms} is a vector or list of parameters.
      ... (optional) are any other arguments passed to the function.
        
      The return value of func should be a list, whose first element is a
      vector containing the derivatives of \code{y} with respect to
      \code{time}, and whose second element is a vector of global values that are required at
      each point in \code{times}.}
  \item{parms }{vector or list holding the parameters used in \code{func}.}
  \item{hini}{length of the (fixed) internal time step.}
  \item{... }{additional arguments, allowing this to be a generic function}
}
\author{Thomas Petzoldt \email{thomas.petzoldt@tu-dresden.de}}

\details{See help pages of \code{\link{rk}} and \code{\link{rkMethod}} 
  for details.
}

\note{\code{rk4} and \code{euler} are primarily provided for didactic reasons. 
     For most practical cases, solvers with flexible timestep (e.g. "rk45dp7")
     and especially solvers of the Livermore family  (e.g. \code{\link{lsoda}}) are superior.
}

\examples{
## numeric solution of logistic growth
logist <- function(t, x, parms) {
  with(as.list(parms), {
    dx <- r * x[1] * (1 - x[1]/K)
    list(dx)
  })
}

time  <- 0:100
N0    <- 0.1; r <- 0.5; K <- 100
parms <- c(r = r, K = K)
x<- c(N = N0)

## analytical solution
plot(time, K/(1+(K/N0-1)*exp(-r*time)), ylim=c(0, 120), type="l", col="red", lwd=2)

## good solution
out <- as.data.frame(rk4(x, time, logist, parms, hini=2))
points(out$time, out$N, pch=16, col="blue", cex=0.5)

## same time step, systematic under-estimation
out <- as.data.frame(euler(x, time, logist, parms, hini=2))
points(out$time, out$N, pch=1)

## unstable result
out <- as.data.frame(euler(x, time, logist, parms, hini=4))
points(out$time, out$N, pch=8, cex=0.5)

## method with automatic time step
out <- as.data.frame(lsoda(x, time, logist, parms))
points(out$time, out$N, pch=1, col="green")
}


\seealso{
  \code{\link{rk}}, 
  \code{\link{ode}}, 
  \code{\link{rkMethod}}, 
  \code{\link{lsoda}},\code{\link{lsode}}, \code{\link{lsodes}}, 
    \code{\link{lsodar}}, \code{\link{daspk}}}
\keyword{math}

