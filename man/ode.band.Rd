\name{ode.band}
\alias{ode.band}
\title{Solver for Ordinary Differential Equations; Assumes a Banded
  Jacobian
}
\description{
  Solves a system of ordinary differential equations.
  
  Assumes a banded Jacobian matrix, but does not rearrange the state
  variables (in contrast to ode.1D).  Suitable for 1-D models that
  include transport only between adjacent layers and that model only one
  species.
}

\usage{ode.band(y, times, func, parms, nspec = NULL, dimens = NULL, 
  bandup = nspec, banddown = nspec, method = "lsode", names = NULL, 
  ...)
}
\arguments{
  \item{y }{the initial (state) values for the ODE system, a vector. If
    \code{y} has a name attribute, the names will be used to label the
    output matrix.
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time.
  }
  \item{func }{either an \R-function that computes the values of the
    derivatives in the ODE system (the model definition) at time
    \code{t}, or a character string giving the name of a compiled
    function in a dynamically loaded shared library.

    If \code{func} is an \R-function, it must be defined as:
    \code{func <- function(t, y, parms, ...)}.  \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system.  If the initial values \code{y} has a
    \code{names} attribute, the names will be available inside \code{func}.
    \code{parms} is a vector or list of parameters; \code{...} (optional) are
    any other arguments passed to the function.

    The return value of \code{func} should be a list, whose first
    element is a vector containing the derivatives of \code{y} with
    respect to \code{time}, and whose next elements are global values
    that are required at each point in \code{times}.The derivatives
    should be specified in the same order as the state variables \code{y}.

  }
  \item{parms }{parameters passed to \code{func}.
  }
  \item{nspec }{the number of *species* (components) in the model.
  }
  \item{dimens}{the number of \bold{boxes} in the model. If \code{NULL}, then
    \code{nspec} should be specified.
  }
  \item{bandup }{the number of nonzero bands above the Jacobian
    diagonal.
  }
  \item{banddown }{the number of nonzero bands below the Jacobian
    diagonal.
  }
  \item{method }{the integrator to use, one of \code{"vode"},
    \code{"lsode"}, \code{"lsoda"}, \code{"lsodar"}, \code{"radau"}.
  }
  \item{names }{the names of the components; used for plotting.
  }
  \item{... }{additional arguments passed to the integrator.}
}
\value{
  A matrix with up to as many rows as elements in \code{times} and as
  many columns as elements in \code{y} plus the number of "global"
  values returned in the second element of the return from \code{func},
  plus an additional column (the first) for the time value.  There will
  be one row for each element in \code{times} unless the integrator
  returns with an unrecoverable error.  If \code{y} has a names
  attribute, it will be used to label the columns of the output value.
  
  The output will have the attributes \code{istate} and \code{rstate},
  two vectors with several elements. See the help for the selected
  integrator for details.  the first element of istate returns the
  conditions under which the last call to the integrator returned. Normal is
  \code{istate = 2}.  If \code{verbose = TRUE}, the settings of
  \code{istate} and \code{rstate} will be written to the screen.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{

## =======================================================================
## The Aphid model from Soetaert and Herman, 2009.
## A practical guide to ecological modelling.
## Using R as a simulation platform. Springer.
## =======================================================================

## 1-D diffusion model

## ================
## Model equations
## ================
Aphid <- function(t, APHIDS, parameters) {
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
out   <- ode.band(state, times, Aphid, parms = 0, 
                  nspec = 1, names = "Aphid")

## ================
## Plotting output
## ================
image(out, grid = Distance, method = "filled.contour", 
      xlab = "time, days", ylab = "Distance on plant, m",
      main = "Aphid density on a row of plants")

matplot.1D(out, grid = Distance, type = "l", 
   subset = time \%in\% seq(0, 200, by = 10))

# add an observed dataset to 1-D plot (make sure to use correct name):
data <- cbind(dist  = c(0,10, 20,  30,  40, 50, 60), 
              Aphid = c(0,0.1,0.25,0.5,0.25,0.1,0))

matplot.1D(out, grid = Distance, type = "l", 
   subset = time \%in\% seq(0, 200, by = 10), 
   obs = data, obspar = list(pch = 18, cex = 2, col="red"))
\dontrun{
plot.1D(out, grid = Distance, type = "l")

}
}
\details{
  This is the method of choice for single-species 1-D reactive transport
  models.
  
  For multi-species 1-D models, this method can only be used if the
  state variables are arranged per box, per species (e.g. A[1], B[1],
  A[2], B[2], A[3], B[3], ... for species A, B).  By default, the
  \bold{model} function will have the species arranged as A[1], A[2],
  A[3], ... B[1], B[2], B[3], ...  in this case, use \code{ode.1D}.
  
  See the selected integrator for the additional options.
}
\seealso{
  \itemize{
    \item \code{\link{ode}} for a general interface to most of the ODE solvers,
    \item  \code{\link{ode.1D}} for integrating 1-D models
    \item  \code{\link{ode.2D}} for integrating 2-D models
    \item  \code{\link{ode.3D}} for integrating 3-D models
    \item  \code{\link{lsode}}, \code{\link{lsoda}}, \code{\link{lsodar}},
    \code{\link{vode}} for the integration options.
  }
  \code{\link{diagnostics}} to print diagnostic messages.
}
\keyword{math}

