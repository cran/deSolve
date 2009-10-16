\name{plot.deSolve}
\alias{plot.deSolve}
\alias{hist.deSolve}
\title{
 Plot and Histogram Method for deSolve Objects
}
\description{
Plot the output of numeric integration routines.
}
\usage{
\method{plot}{deSolve}(x, which = 1:(ncol(x)-1), ask = NULL, ...)
\method{hist}{deSolve}(x, which = 1:(ncol(x)-1), ask = NULL, ...)
}

\arguments{
  \item{x }{an object of class \code{deSolve}, as returned by the
    integrators, and to be plotted.}

  \item{which }{the name(s) or the index to the variables that should be
    plotted.  Default = all variables. }

  \item{ask }{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, if \code{NULL} the user is only asked if more than one
    page of plots is necessary and the current graphics device is set
    interactive, see \code{\link{par}(ask=.)} and
    \code{\link{dev.interactive}}.}
  
  \item{\dots}{additional graphics arguments passed to
    \code{\link{plot.default}} or \code{\link{hist}}}
}

\details{
  The number of panels per page is automatically determined up to 3 x 3
  (\code{par(mfrow=c(3, 3))}). This default can be overwritten by
  specifying user-defined settings for \code{mfrow} or \code{mfcol}.

  Other graphical parameters can be passed as well. Parameters
  \code{xlab} and \code{ylab} are vectorized, so it is possible to
  assign specific axis labels to individual plots.
}

\seealso{
 \code{\link{print.deSolve}}, \code{\link{ode}},  \code{\link{deSolve}}  
}
\examples{
## A Predator-Prey model with 4 species in matrix formulation
LVmatrix <- function(t, n, parms) {
  with(parms, {
    dn <- r * n + n * (A \%*\% n)
    return(list(c(dn)))
  })
}
parms <- list(
  r = c(r1 = 0.1, r2 = 0.1, r3 = -0.1, r4 = -0.1),
  A = matrix(c(0.0, 0.0, -0.2, 0.01,      # prey 1
               0.0, 0.0, 0.02, -0.1,      # prey 2
               0.2, 0.02, 0.0, 0.0,       # predator 1; prefers prey 1
               0.01, 0.1, 0.0, 0.0),      # predator 2; prefers prey 2
               nrow = 4, ncol = 4, byrow=TRUE)
)
times <- seq(from = 0, to = 500, by = 0.1)
y     <- c(prey1 = 1, prey2 = 1, pred1 = 2, pred2 = 2)

out <- ode(y, times, LVmatrix, parms)

## Basic line plot
plot(out, type = "l")

## User-specified axis labels
plot(out, type = "l", ylab = c("Prey 1", "Prey 2", "Pred 1", "Pred 2"),
  xlab = "Time (d)", main = "Time Series") 

hist(out, col="darkblue", breaks = 50)
}

\keyword{ hplot }

