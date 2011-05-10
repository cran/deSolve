\name{ode.3D}
\alias{ode.3D}

\title{Solver for 3-Dimensional Ordinary Differential Equations}

\description{
  Solves a system of ordinary differential equations resulting from
  3-Dimensional partial differential equations
  that have been converted to ODEs by numerical differencing.
}

\usage{ode.3D(y, times, func, parms, nspec = NULL, dimens, 
  method = c("lsodes", "euler", "rk4", "ode23", "ode45", "adams", "iteration"),
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
    \code{func <- function(t, y, parms, ...)}.  \code{t} is the current time
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
  \item{dimens}{3-valued vector with the number of \bold{boxes} in three
    dimensions in the model. 
  }
  \item{names }{the names of the components; used for plotting.
  }
  \item{cyclicBnd }{if not \code{NULL} then a number or a 3-valued vector
    with the dimensions where a cyclic boundary is used - \code{1}: x-dimension,
    \code{2}: y-dimension; \code{3}: z-dimension.
  }
  \item{method }{the integrator. Use \code{"lsodes"} if the model is very stiff;
     "impAdams" may be best suited for mildly stiff problems; 
     \code{"euler", "rk4", "ode23", "ode45", "adams"} are most
     efficient for non-stiff problems. Also allowed is to pass an integrator
     \code{function}. Use one of the other Runge-Kutta methods via 
     \code{rkMethod}. For instance, \code{method = rkMethod("ode45ck")} will
     trigger the Cash-Karp method of order 4(5).
     
    Method \code{"iteration"} is special in that here the function \code{func} should
  return the new value of the state variables rather than the rate of change.
  This can be used for individual based models, for difference equations,
  or in those cases where the integration is performed within \code{func})

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
  whether the input makes sense (as \code{nspec*dimens[1]*dimens[2]*dimens[3]
  == length(y)}).

  Do \bold{not} use this method for problems that are not 3D!
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
## =======================================================================
## Diffusion in 3-D; imposed boundary conditions
## =======================================================================
diffusion3D <- function(t, Y, par) {

  ## function to bind two matrices to an array
  mbind <- function (Mat1, Array, Mat2, along = 1) {
    dimens <- dim(Array) + c(0, 0, 2)
    if (along == 3)
      array(dim = dimens, data = c(Mat1, Array, Mat2))
    else if (along == 1)
      aperm(array(dim = dimens,
        data=c(Mat1, aperm(Array, c(3, 2, 1)), Mat2)), c(3, 2, 1))
    else if (along == 2)
      aperm(array(dim = dimens,
        data = c(Mat1, aperm(Array, c(1, 3, 2)), Mat2)), c(1, 3, 2))
  }

  yy  <- array(dim=c(n, n, n), data = Y)      # vector to 3-D array
  dY  <- -r*yy                                # consumption
  BND <- matrix(nrow = n, ncol = n, data = 1) # boundary concentration

  ## diffusion in x-direction
  ## new array including boundary concentrations in X-direction
  BNDx <- mbind(BND, yy, BND, along = 1)
  ## diffusive Flux
  Flux <- -Dx * (BNDx[2:(n+2),,] - BNDx[1:(n+1),,])/dx
  ## rate of change = - flux gradient
  dY[] <- dY[] - (Flux[2:(n+1),,] - Flux[1:n,,])/dx

  ## diffusion in y-direction
  BNDy <- mbind(BND, yy, BND, along = 2)
  Flux <- -Dy * (BNDy[,2:(n+2),] - BNDy[,1:(n+1),])/dy
  dY[] <- dY[] - (Flux[,2:(n+1),] - Flux[,1:n,])/dy

  ## diffusion in z-direction
  BNDz <- mbind(BND, yy, BND, along = 3)
  Flux <- -Dz * (BNDz[,,2:(n+2)] - BNDz[,,1:(n+1)])/dz
  dY[] <- dY[] - (Flux[,,2:(n+1)] - Flux[,,1:n])/dz

  return(list(as.vector(dY)))
}

## parameters
dy <- dx <- dz <-1   # grid size
Dy <- Dx <- Dz <-1   # diffusion coeff, X- and Y-direction
r  <- 0.025          # consumption rate

n  <- 10
y  <- array(dim=c(n,n,n),data=10.)

## use lsodes, the default (for n>20, Runge-Kutta more efficient)
print(system.time(
  RES <- ode.3D(y, func = diffusion3D, parms = NULL, dimens = c(n, n, n),
                times = 1:20, lrw = 120000, atol = 1e-10,
                rtol = 1e-10, verbose = TRUE)
))

y <- array(dim = c(n, n, n), data = RES[nrow(RES), -1])
filled.contour(y[, , n/2], color.palette = terrain.colors)

summary(RES)

\dontrun{
for (i in 2:nrow(RES)) {
  y <- array(dim=c(n,n,n),data=RES[i,-1])
  filled.contour(y[,,n/2],main=i,color.palette=terrain.colors)
}
}
}
\details{
  This is the method of choice for 3-dimensional models, that are only
  subjected to transport between adjacent layers.

  Based on the dimension of the problem, the method first calculates the
  sparsity pattern of the Jacobian, under the assumption that transport
  is only occurring between adjacent layers. Then \code{lsodes} is
  called to solve the problem.
  
  As \code{lsodes} is used to integrate, it will probably be necessary
  to specify the length of the real work array, \code{lrw}.

  Although a reasonable guess of \code{lrw} is made, it is likely that
  this will be too low.  
  
  In this case, \code{ode.2D} will return with an
  error message telling the size of the work array actually needed. In
  the second try then, set \code{lrw} equal to this number.

  For instance, if you get the error:
\preformatted{   
DLSODES- RWORK length is insufficient to proceed.                               
  Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)                    
  In above message,  I1 =     27627   I2 =     25932 
}      
  set \code{lrw} equal to 27627 or a higher value.
    
  See \code{lsodes} for the additional options.
}
\seealso{
  \itemize{
    \item \code{\link{ode}} for a general interface to most of the ODE solvers,
    \item \code{\link{ode.band}} for integrating models with a banded Jacobian
    \item  \code{\link{ode.1D}} for integrating 1-D models
    \item  \code{\link{ode.2D}} for integrating 2-D models
    \item  \code{\link{lsodes}} for the integration options.
  }
  \code{\link{diagnostics}} to print diagnostic messages.
}
\keyword{math}

