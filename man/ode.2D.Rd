\name{ode.2D}
\alias{ode.2D}

\title{Solver for 2-Dimensional Ordinary Differential Equations}

\description{
  Solves a system of ordinary differential equations resulting from
  2-Dimensional reaction-transport models that include transport only
  between adjacent layers.
}

\usage{ode.2D(y, times, func, parms, nspec = NULL, dimens, ...)}
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
    \code{yprime = func(t, y, parms,...)}.  \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system.  If the initial values \code{y} has a
    names attribute, the names will be available inside \code{func}.
    \code{parms} is a vector or list of parameters; ... (optional) are
    any other arguments passed to the function.
  
    The return value of \code{func} should be a list, whose first
    element is a vector containing the derivatives of \code{y} with
    respect to \code{time}, and whose next elements are global values
    that are required at each point in \code{times}.
  }
  \item{parms }{parameters passed to \code{func}.}
  \item{nspec }{the number of \bold{species} (components) in the model.}
  \item{dimens}{2-valued vector with the number of \bold{boxes} in two
    dimensions in the model. 
  }
  \item{... }{additional arguments passed to \code{lsodes}.}
}
\value{
  
  A matrix with up to as many rows as elements in times and as many
  columns as elements in \code{y} plus the number of "global" values
  returned in the second element of the return from \code{func}, plus an
  additional column (the first) for the time value.  There will be one
  row for each element in \code{times} unless the integrator returns
  with an unrecoverable error.  If \code{y} has a names attribute, it
  will be used to label the columns of the output value.
      
  The output will have the attributes \code{istate}, and \code{rstate},
  two vectors with several useful elements.  The first element of istate
  returns the conditions under which the last call to the integrator
  returned. Normal is \code{istate = 2}.  If \code{verbose = TRUE}, the
  settings of istate and rstate will be written to the screen. See the
  help for the selected integrator for details.
  
}
\note{
  It is advisable though not mandatory to specify \bold{both}
  \code{nspec} and \code{dimens}. In this case, the solver can check
  whether the input makes sense (as \code{nspec * dimens[1] * dimens[2]
  == length(y)}).

  Do \bold{not} use this method for problems that are not 2D!
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
## ============================================================
## A Lotka-Volterra predator-prey model with predator and prey
## dispersing in 2 dimensions
## ============================================================

## ==================
## Model definitions
## ==================

lvmod2D <- function (time, state, pars, N, Da, dx)
{
  NN <- N*N
  Prey <- matrix(nr = N,nc = N,state[1:NN])
  Pred <- matrix(nr = N,nc = N,state[(NN+1):(2*NN)])

  with (as.list(pars),
  {
    ## Biology
    dPrey   <- rGrow* Prey *(1- Prey/K) - rIng* Prey *Pred
    dPred   <- rIng* Prey *Pred*assEff -rMort* Pred

    zero <- rep(0,N)

    ## 1. Fluxes in x-direction; zero fluxes near boundaries
    FluxPrey <- -Da * rbind(zero,(Prey[2:N,]-Prey[1:(N-1),]), zero)/dx
    FluxPred <- -Da * rbind(zero,(Pred[2:N,]-Pred[1:(N-1),]), zero)/dx

    ## Add flux gradient to rate of change
    dPrey    <- dPrey - (FluxPrey[2:(N+1),]-FluxPrey[1:N,])/dx
    dPred    <- dPred - (FluxPred[2:(N+1),]-FluxPred[1:N,])/dx

    ## 2. Fluxes in y-direction; zero fluxes near boundaries
    FluxPrey <- -Da * cbind(zero,(Prey[,2:N]-Prey[,1:(N-1)]), zero)/dx
    FluxPred <- -Da * cbind(zero,(Pred[,2:N]-Pred[,1:(N-1)]), zero)/dx

    ## Add flux gradient to rate of change
    dPrey    <- dPrey - (FluxPrey[,2:(N+1)]-FluxPrey[,1:N])/dx
    dPred    <- dPred - (FluxPred[,2:(N+1)]-FluxPred[,1:N])/dx

    return (list(c(as.vector(dPrey), as.vector(dPred))))
 })
}


## ===================
## Model applications
## ===================


pars    <- c(rIng   = 0.2,    # /day, rate of ingestion
             rGrow  = 1.0,    # /day, growth rate of prey
             rMort  = 0.2 ,   # /day, mortality rate of predator
             assEff = 0.5,    # -, assimilation efficiency
             K      = 5  )    # mmol/m3, carrying capacity

R  <- 20                      # total length of surface, m
N  <- 50                      # number of boxes in one direction
dx <- R/N                     # thickness of each layer
Da <- 0.05                    # m2/d, dispersion coefficient

NN <- N*N                     # total number of boxes

## initial conditions
yini    <- rep(0, 2*N*N)
cc      <- c((NN/2):(NN/2+1)+N/2, (NN/2):(NN/2+1)-N/2)
yini[cc] <- yini[NN+cc] <- 1

## solve model (5000 state variables...
times   <- seq(0, 50, by = 1)
out <- ode.2D(y = yini, times = times, func = lvmod2D, parms = pars,
              dimens = c(N, N), N = N, dx = dx, Da = Da, lrw = 5000000)

## plot results
Col <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

for (i in seq(1, length(times), by = 1))
   image(matrix(nr = N, nc = N, out[i, 2:(NN+1)]),
   col = Col(100), xlab = "x", ylab = "y", zlim = range(out[,2:(NN+1)]))

}
\details{
  This is the method of choice for 2-dimensional models, that are only
  subjected to transport between adjacent layers.

  Based on the dimension of the problem, the method first calculates the
  sparsity pattern of the Jacobian, under the assumption that transport
  is onely occurring between adjacent layers. Then \code{lsodes} is
  called to solve the problem.
  
  As \code{lsodes} is used to integrate, it will probably be necessary
  to specify the length of the real work array, \code{lrw}.

  Although a reasonable guess of \code{lrw} is made, it is likely that
  this will be too low.  In this case, \code{ode.2D} will return with an
  error message telling the size of the work array actually needed. In
  the second try then, set \code{lrw} equal to this number.

  See \code{lsodes} for the additional options.
}
\seealso{
  \itemize{
    \item \code{\link{ode}}
    \item \code{\link{ode.band}} for integrating models with a banded Jacobian
    \item  \code{\link{ode.1D}} for integrating 1-D models
    \item  \code{\link{lsodes}} for the integration options.
  }
}
\keyword{math}

