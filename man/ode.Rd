\name{ode}
\alias{ode}
\alias{print.deSolve}

\title{General Solver for Ordinary Differential Equations}

\description{Solves a system of ordinary differential equations; a wrapper
around the implemented ODE solvers}

\usage{ode(y, times, func, parms, 
method = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk",
           "euler", "rk4", "ode23", "ode45"), ...)

\method{print}{deSolve}(x, \dots)
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
    derivatives in the ODE system (the model definition) at time t, or a
    character string giving the name of a compiled function in a
    dynamically loaded shared library.

    If \code{func} is an \R-function, it must be defined as:
     \code{func <- function(t, y, parms,...)}.  \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system.  If the initial values \code{y} has a
    \code{names} attribute, the names will be available inside \code{func}.
    \code{parms} is a vector or list of parameters; ... (optional) are
    any other arguments passed to the function.
  
    The return value of \code{func} should be a list, whose first
    element is a vector containing the derivatives of \code{y} with
    respect to \code{time}, and whose next elements are global values
    that are required at each point in \code{times}.The derivatives
    should be specified in the same order as the state variables \code{y}.

    If \code{func} is a string, then \code{dllname} must give the name
    of the shared library (without extension) which must be loaded
    before \code{ode} is called. See package vignette
    \code{"compiledCode"} for more details.
  }
  \item{parms }{parameters passed to \code{func}.}
  \item{method }{the integrator to use, either a string (\code{"lsoda"},
    \code{"lsode"}, \code{"lsodes"}, \code{"lsodar"}, \code{"vode"},
    \code{"daspk"}, \code{"euler"}, \code{"rk4"}, \code{"ode23"} or
    \code{"ode45"}) or a function that performs integration, or a list
    of class \code{rkMethod}.
    The default integrator used is \link{lsoda}
  }
   \item{x }{an object of class \code{deSolve}, as returned by the
    integrators, and to be printed.}
  \item{... }{additional arguments passed to the integrator or to the
    methods.}
}

\value{
  A matrix of class \code{deSolve} with up to as many rows as elements in
  \code{times} and as many
  columns as elements in \code{y} plus the number of "global" values
  returned in the second element of the return from \code{func}, plus an
  additional column (the first) for the time value.  There will be one
  row for each element in \code{times} unless the integrator returns
  with an unrecoverable error.  If \code{y} has a names attribute, it
  will be used to label the columns of the output value.

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}

\details{
  This is simply a wrapper around the various ode solvers.
  
  See package vignette for information about specifying the model in
  compiled code.
  
  See the selected integrator for the additional options.
  
  The default integrator used is \link{lsoda}
}

\seealso{
  \itemize{
    \item \code{\link{plot.deSolve}} for plotting the outputs,
    \item \code{\link{ode.band}} for solving models with a banded
      Jacobian,
    \item \code{\link{ode.1D}} for integrating 1-D models,
    \item \code{\link{ode.2D}} for integrating 2-D models,
    \item \code{\link{ode.3D}} for integrating 3-D models,
    \item \code{\link{aquaphy}}, \code{\link{ccl4model}}, where
      \code{ode} is used,
    \item \code{\link{lsoda}}, \code{\link{lsode}},
      \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}},
      \code{\link{daspk}},
    \item  \code{\link{rk}}, \code{\link{rkMethod}}
  }
  \code{\link{diagnostics}} to print diagnostic messages.
}
\keyword{math}
\examples{

## =======================================================================
## Example1: Predator-Prey Lotka-Volterra model
## =======================================================================

LVmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Ingestion    <- rIng  * Prey*Predator
    GrowthPrey   <- rGrow * Prey*(1-Prey/K)
    MortPredator <- rMort * Predator

    dPrey        <- GrowthPrey - Ingestion
    dPredator    <- Ingestion*assEff -MortPredator

    return(list(c(dPrey, dPredator)))
  })
}

pars    <- c(rIng   = 0.2,    # /day, rate of ingestion
             rGrow  = 1.0,    # /day, growth rate of prey
             rMort  = 0.2 ,   # /day, mortality rate of predator
             assEff = 0.5,    # -, assimilation efficiency
             K      = 10)     # mmol/m3, carrying capacity

yini    <- c(Prey = 1, Predator = 2)
times   <- seq(0, 200, by = 1)
out     <- as.data.frame(ode(func = LVmod, y = yini,
                         parms = pars, times = times))

matplot(out$time,out[,2:3],type = "l",xlab = "time",ylab = "Conc",
        main = "Lotka-Volterra",lwd = 2)
legend("topright", c("prey", "predator"), col = 1:2, lty = 1:2)

## =======================================================================
## Example2: Substrate-Producer-Consumer Lotka-Volterra model
## =======================================================================

## Note:
## Function sigimp passed as an argument (input) to model
##   (see also lsoda and rk examples)

SPCmod <- function(t, x, parms, input)  {
  with(as.list(c(parms, x)), {
    import <- input(t)
    dS <- import - b*S*P + g*C    # substrate
    dP <- c*S*P  - d*C*P          # producer
    dC <- e*P*C  - f*C            # consumer
    res <- c(dS, dP, dC)
    list(res)
  })
}

## The parameters 
parms  <- c(b = 0.001, c = 0.1, d = 0.1, e = 0.1, f = 0.1, g = 0.0)

## vector of timesteps
times  <- seq(0, 200, length = 101)

## external signal with rectangle impulse
signal <- as.data.frame(list(times = times,
                            import = rep(0, length(times))))

signal$import[signal$times >= 10 & signal$times <= 11] <- 0.2

sigimp <- approxfun(signal$times, signal$import, rule = 2)


## Start values for steady state
xstart <- c(S = 1, P = 1, C = 1)

## Solve model
out <- ode(y = xstart,times = times,
           func = SPCmod, parms, input = sigimp)

## Default plot method
plot(out, type="l")

## User specified plotting
mf <- par(mfrow = c(1, 2))
matplot(out[,1], out[,2:4], type = "l", xlab = "time", ylab = "state")
legend("topright", col = 1:3, lty = 1:3, legend = c("S", "P", "C"))
plot(out[,"P"], out[,"C"], type = "l", lwd = 2, xlab = "producer",
  ylab = "consumer")
par(mfrow=mf)

}

