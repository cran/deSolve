\name{plot.deSolve}
\alias{plot.deSolve}
\alias{hist.deSolve}
\alias{image.deSolve}
\title{
 Plot and Histogram Method for deSolve Objects
}
\description{
Plot the output of numeric integration routines.
}
\usage{
\method{plot}{deSolve}(x, which = 1:(ncol(x)-1), ask = NULL, ...)
\method{hist}{deSolve}(x, which = 1:(ncol(x)-1), ask = NULL, ...)
\method{image}{deSolve}(x, which = NULL, ask = NULL, 
  add.contour = FALSE, grid = NULL, ...)  
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
  
  \item{add.contour }{if \code{TRUE}, will add contours to the image plot. }

  \item{grid }{only for \code{image} plots: the 1-D grid as a vector (output 
    generated with \code{ode.1D}, or the x- and y-grid, as a \code{list} for 
    output generated with \code{ode.2D}. }

  \item{\dots}{additional graphics arguments passed to
    \code{\link{plot.default}}, \code{\link{image}} or \code{\link{hist}}}
}

\details{
  The number of panels per page is automatically determined up to 3 x 3
  (\code{par(mfrow=c(3, 3))}). This default can be overwritten by
  specifying user-defined settings for \code{mfrow} or \code{mfcol}.

  Other graphical parameters can be passed as well. Parameters
  \code{xlab} and \code{ylab} are vectorized, so it is possible to
  assign specific axis labels to individual plots.
  
  Image plots will only work for 1-D and 2-D variables, as solved with 
  \code{\link{ode.1D}} and \code{\link{ode.2D}}. In the first case, an
  image with \code{times} as x- and \code{the grid} as y-axis will be created.
  in second case, an x-y plot will be created, for all times. Unless 
  \code{ask=FALSE}, the user will be asked to confirm page changes. 
}

\seealso{
 \code{\link{print.deSolve}}, \code{\link{ode}},  \code{\link{deSolve}}  
}
\examples{
## =======================================================================
## A Predator-Prey model with 4 species in matrix formulation
## =======================================================================

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

## =======================================================================
## The Aphid model from Soetaert and Herman, 2009.
## A practical guide to ecological modelling.
## Using R as a simulation platform. Springer.
## =======================================================================

## 1-D diffusion model

## ================
## Model equations
## ================
Aphid <- function(t, APHIDS, parameters)
{
  deltax  <- c (0.5, rep(1, numboxes-1), 0.5)
  Flux    <- -D*diff(c(0, APHIDS, 0))/deltax
  dAPHIDS <- -diff(Flux)/delx + APHIDS*r

  list(dAPHIDS)   # the output
}
  
## ==================
## Model application
## ==================

## the model parameters:
D         <- 0.3    # m2/day  diffusion rate
r         <- 0.01   # /day    net growth rate
delx      <- 1      # m       thickness of boxes
numboxes  <- 60 

## distance of boxes on plant, m, 1 m intervals
Distance  <- seq(from = 0.5, by = delx, length.out = numboxes)

## Initial conditions, ind/m2
## aphids present only on two central boxes
APHIDS        <- rep(0, times = numboxes)
APHIDS[30:31] <- 1
state         <- c(APHIDS = APHIDS)      # initialise state variables 
                  
## RUNNING the model:
times <- seq(0, 200, by = 1)   # output wanted at these time intervals
out   <- ode.1D(state, times, Aphid, parms = 0, nspec = 1)

image(out, grid = Distance, main="Aphid model", ylab="distance, m")
}

\keyword{ hplot }

