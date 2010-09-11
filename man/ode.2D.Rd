\name{ode.2D}
\alias{ode.2D}

\title{Solver for 2-Dimensional Ordinary Differential Equations}

\description{
  Solves a system of ordinary differential equations resulting from
  2-Dimensional partial differential equations
  that have been converted to ODEs by numerical differencing.
}

\usage{ode.2D(y, times, func, parms, nspec = NULL, dimens,
  method= c("lsodes", "euler", "rk4", "ode23", "ode45", "adams"),
  names = NULL, cyclicBnd = NULL, ...)}
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
    \code{func <- function(t, y, parms, ...)}. \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system.  If the initial values \code{y} has a
    \code{names} attribute, the names will be available inside \code{func}.
    \code{parms} is a vector or list of parameters; \code{...} (optional) are
    any other arguments passed to the function.
  
    The return value of \code{func} should be a list, whose first
    element is a vector containing the derivatives of \code{y} with
    respect to \code{time}, and whose next elements are global values
    that are required at each point in \code{times}. The derivatives
    should be specified in the same order as the state variables \code{y}.

  }
  \item{parms }{parameters passed to \code{func}.}
  \item{nspec }{the number of \bold{species} (components) in the model.}
  \item{dimens}{2-valued vector with the number of \bold{boxes} in two
    dimensions in the model. 
  }
  \item{cyclicBnd }{if not \code{NULL} then a number or a 2-valued vector
    with the dimensions where a cyclic boundary is used - \code{1}: x-dimension,
    \code{2}: y-dimension; see details.
  }
  \item{names }{the names of the components; used for plotting.
  }
  \item{method }{the integrator. Use \code{"lsodes"} if the model is very stiff;
     \code{"impAdams"} may be best suited for mildly stiff problems; 
     \code{"euler", "rk4", "ode23", "ode45", "adams"} are most
     efficient for non-stiff problems. Also allowed is to pass an integrator
     \code{function}. Use one of the other Runge-Kutta methods via 
     \code{rkMethod}. For instance, \code{method = rkMethod("ode45ck")} will
     trigger the Cash-Karp method of order 4(5).
     
     If  \code{"lsodes"} is used, then also the size of the work array should
     be specified (\code{lrw}) (see \link{lsodes}).
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
## =======================================================================
## A Lotka-Volterra predator-prey model with predator and prey
## dispersing in 2 dimensions
## =======================================================================

## ==================
## Model definitions
## ==================

lvmod2D <- function (time, state, pars, N, Da, dx) {
  NN <- N*N
  Prey <- matrix(nr = N,nc = N,state[1:NN])
  Pred <- matrix(nr = N,nc = N,state[(NN+1):(2*NN)])

  with (as.list(pars), {
    ## Biology
    dPrey <- rGrow * Prey * (1- Prey/K) - rIng  * Prey * Pred
    dPred <- rIng  * Prey * Pred*assEff - rMort * Pred

    zero <- rep(0, N)

    ## 1. Fluxes in x-direction; zero fluxes near boundaries
    FluxPrey <- -Da * rbind(zero,(Prey[2:N,] - Prey[1:(N-1),]), zero)/dx
    FluxPred <- -Da * rbind(zero,(Pred[2:N,] - Pred[1:(N-1),]), zero)/dx

    ## Add flux gradient to rate of change
    dPrey    <- dPrey - (FluxPrey[2:(N+1),] - FluxPrey[1:N,])/dx
    dPred    <- dPred - (FluxPred[2:(N+1),] - FluxPred[1:N,])/dx

    ## 2. Fluxes in y-direction; zero fluxes near boundaries
    FluxPrey <- -Da * cbind(zero,(Prey[,2:N] - Prey[,1:(N-1)]), zero)/dx
    FluxPred <- -Da * cbind(zero,(Pred[,2:N] - Pred[,1:(N-1)]), zero)/dx

    ## Add flux gradient to rate of change
    dPrey    <- dPrey - (FluxPrey[,2:(N+1)] - FluxPrey[,1:N])/dx
    dPred    <- dPred - (FluxPred[,2:(N+1)] - FluxPred[,1:N])/dx

    return(list(c(as.vector(dPrey), as.vector(dPred))))
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

## solve model (5000 state variables...  use Cash-Karp Runge-Kutta method
times   <- seq(0, 50, by = 1)
out <- ode.2D(y = yini, times = times, func = lvmod2D, parms = pars,
              dimens = c(N, N), names = c("Prey", "Pred"),
              N = N, dx = dx, Da = Da, , method = rkMethod("rk45ck"))

## plot results
Col <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                          "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

\dontrun{
for (i in seq(1, length(times), by = 1))
   image(matrix(nr = N, nc = N, out[i, 2:(NN+1)]),
   col = Col(100), xlab = , zlim = range(out[,2:(NN+1)]))

## similar:
image(out, xlab = "x", ylab = "y")
}

## =======================================================================
## An example with a cyclic boundary condition.
## Diffusion in 2-D; extra flux on 2 boundaries,
## cyclic boundary in y
## =======================================================================


diffusion2D <- function(t, Y, par) {
  y    <- matrix(nrow = nx, ncol = ny, data = Y)  # vector to 2-D matrix
  dY   <- -r * y       # consumption
  BNDx <- rep(1, nx)   # boundary concentration
  BNDy <- rep(1, ny)   # boundary concentration

  ## diffusion in X-direction; boundaries=imposed concentration
  Flux <- -Dx * rbind(y[1,]-BNDy,(y[2:nx,]-y[1:(nx-1),]),BNDy-y[nx,])/dx
  dY   <- dY - (Flux[2:(nx+1),]-Flux[1:nx,])/dx

  ## diffusion in Y-direction
  Flux <- -Dy * cbind(y[,1]-BNDx,(y[,2:ny]-y[,1:(ny-1)]),BNDx-y[,ny])/dy
  dY   <- dY - (Flux[,2:(ny+1)]-Flux[,1:ny])/dy

  ## extra flux on two sides
  dY[,1] <- dY[,1] + 10
  dY[1,] <- dY[1,] + 10

  ## and exchange between sides on y-direction
  dY[,ny] <- dY[,ny] + (y[,1] - y[,ny]) * 10
  return(list(as.vector(dY)))
}

## parameters
dy <- dx <- 1  # grid size
Dy <- Dx <- 1  # diffusion coeff, X- and Y-direction
r  <- 0.05     # consumption rate

nx <- 50
ny <- 100
y  <- matrix(nrow = nx, ncol = ny, 1)

## model most efficiently solved with lsodes - need to specify lrw

print(system.time(
  ST3 <- ode.2D(y, times = 1:100, func = diffusion2D, parms = NULL,
                dimens = c(nx, ny), verbose = TRUE,
                lrw = 400000, atol = 1e-10, rtol = 1e-10, cyclicBnd = 2)
))

\dontrun{
zlim <- range(ST3[,-1])
for (i in 2:nrow(ST3)) {
  y <- matrix(nr=nx,nc=ny,data=ST3[i,-1])
  filled.contour(y,zlim=zlim,main=i)
}
}

}
\details{
  This is the method of choice for 2-dimensional models, that are only
  subjected to transport between adjacent layers.

  Based on the dimension of the problem, and if \code{lsodes} is used as 
  the integrator, the method first calculates the
  sparsity pattern of the Jacobian, under the assumption that transport
  is only occurring between adjacent layers. Then \code{lsodes} is
  called to solve the problem.
  
  If \code{lsodes} is used to integrate, it will probably be necessary
  to specify the length of the real work array, \code{lrw}.

  Although a reasonable guess of \code{lrw} is made, it is likely that
  this will be too low.  In this case, \code{ode.2D} will return with an
  error message telling the size of the work array actually needed. In
  the second try then, set \code{lrw} equal to this number.

  For instance, if you get the error:
\preformatted{   
DLSODES- RWORK length is insufficient to proceed.                               
  Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)                    
  In above message,  I1 =     27627   I2 =     25932 
}      
  set \code{lrw} equal to 27627 or a higher value

  In some cases, a cyclic boundary condition exists. This is when the first
  boxes in x-or y-direction interact with the last boxes. In this case, there
  will be extra non-zero fringes in the Jacobian which need to be taken
  into account. The occurrence of cyclic boundaries can be
  toggled on by specifying argument \code{cyclicBnd}. For innstance,
  \code{cyclicBnd = 1} indicates that a cyclic boundary is required only for
  the x-direction, whereas \code{cyclicBnd = c(1,2)} imposes a cyclic boundary
  for both x- and y-direction. The default is no cyclic boundaries.

  See \code{lsodes} for the additional options.
  
  If the model is not stiff, then it is more efficient to use one of the 
  explicit integration routines
}
\seealso{
  \itemize{
    \item \code{\link{ode}} for a general interface to most of the ODE solvers,
    \item \code{\link{ode.band}} for integrating models with a banded Jacobian
    \item  \code{\link{ode.1D}} for integrating 1-D models
    \item  \code{\link{ode.3D}} for integrating 3-D models
    \item  \code{\link{lsodes}} for the integration options.
  }
  \code{\link{diagnostics}} to print diagnostic messages.
}
\keyword{math}

