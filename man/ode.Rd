\name{ode}
\alias{ode}
\alias{print.deSolve}
\alias{summary.deSolve}
\title{General Solver for Ordinary Differential Equations}

\description{Solves a system of ordinary differential equations; a wrapper
around the implemented ODE solvers}

\usage{ode(y, times, func, parms, 
method = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk",
           "euler", "rk4", "ode23", "ode45", "radau", 
           "bdf", "bdf_d", "adams", "impAdams", "impAdams_d", "iteration"), ...)

\method{print}{deSolve}(x, \dots)
\method{summary}{deSolve}(object, select = NULL, which = select, 
                 subset = NULL, \dots)
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
    that are required at each point in \code{times}. The derivatives
    should be specified in the same order as the state variables \code{y}.

    If \code{func} is a string, then \code{dllname} must give the name
    of the shared library (without extension) which must be loaded
    before \code{ode} is called. See package vignette
    \code{"compiledCode"} for more details.
  }
  \item{parms }{parameters passed to \code{func}.}
  \item{method }{the integrator to use, either a \bold{function} that performs 
    integration, or a \bold{list} of class \code{rkMethod}, or a \bold{string} 
    (\code{"lsoda"},
    \code{"lsode"}, \code{"lsodes"},\code{"lsodar"},\code{"vode"},
    \code{"daspk"}, \code{"euler"}, \code{"rk4"},   \code{"ode23"},
    \code{"ode45"}, \code{"radau"}, \code{"bdf"},   \code{"bdf_d"}, \code{"adams"}, 
    \code{"impAdams"} or \code{"impAdams_d"}  ,"iteration").
    Options "bdf", "bdf_d", "adams", "impAdams" or "impAdams_d" are the backward
    differentiation formula, the BDF with diagonal representation of the Jacobian,
    the (explicit) Adams and the implicit Adams method, and the implicit Adams
    method with diagonal representation of the Jacobian respectively (see details).
    The default integrator used is \link{lsoda}.

    Method \code{"iteration"} is special in that here the function \code{func} should
    return the new value of the state variables rather than the rate of change.
    This can be used for individual based models, for difference equations,
    or in those cases where the integration is performed within \code{func}).
    See last example.
  }
  \item{x }{an object of class \code{deSolve}, as returned by the
    integrators, and to be printed or to be subsetted.
  }
  \item{object }{an object of class \code{deSolve}, as returned by the
    integrators, and whose summary is to be calculated. In contrast to R's default,
    this returns a data.frame. It returns one summary column for a multi-dimensional variable.
  }
  \item{which }{the name(s) or the index to the variables whose summary should be
    estimated.  Default = all variables. }

  \item{select }{which variable/columns to be selected.
  }
  \item{subset }{logical expression indicating elements or rows to keep when 
    calculating a \code{summary}: missing values are taken as \code{FALSE}
  }
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
  
  The default integrator used is \code{\link{lsoda}}.
  
  The option \code{method = "bdf"} provdes a handle to the backward 
  differentiation formula (it is equal to using \code{method = "lsode"}). 
  It is best suited to solve stiff (systems of) equations.
  
  The option \code{method = "bdf_d"} selects the backward 
  differentiation formula that uses Jacobi-Newton iteration (neglecting the 
  off-diagonal elements of the Jacobian (it is equal to using 
  \code{method = "lsode", mf = 23}). 
  It is best suited to solve stiff (systems of) equations.

  \code{method = "adams"} triggers the Adams method that uses functional
  iteration (no Jacobian used); 
  (equal to \code{method = "lsode", mf = 10}. It is often the best 
  choice for solving non-stiff (systems of) equations. Note: when functional
  iteration is used, the method is often said to be explicit, although it is  
  in fact implicit.
  
  \code{method = "impAdams"} selects the implicit Adams method that uses Newton-
  Raphson iteration (equal to \code{method = "lsode", mf = 12}.

  \code{method = "impAdams_d"} selects the implicit Adams method that uses Jacobi-
  Newton iteration, i.e. neglecting all off-diagonal elements (equal to 
  \code{method = "lsode", mf = 13}.

  For very stiff systems, \code{method = "daspk"} may outperform 
  \code{method = "bdf"}.
  
}

\seealso{
  \itemize{
    \item \code{\link{plot.deSolve}} for plotting the outputs,
    \item \code{\link{dede}} general solver for delay differential equations
    \item \code{\link{ode.band}} for solving models with a banded
      Jacobian,
    \item \code{\link{ode.1D}} for integrating 1-D models,
    \item \code{\link{ode.2D}} for integrating 2-D models,
    \item \code{\link{ode.3D}} for integrating 3-D models,
    \item \code{\link{aquaphy}}, \code{\link{ccl4model}}, where
      \code{ode} is used,
    \item \code{\link{lsoda}}, \code{\link{lsode}},
      \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}},
      \code{\link{daspk}}, \code{\link{radau}},
    \item  \code{\link{rk}}, \code{\link{rkMethod}},
  }
  \code{\link{diagnostics}} to print diagnostic messages.
}
\keyword{math}
\examples{

## =======================================================================
## Example1: Predator-Prey Lotka-Volterra model (with logistic prey)
## =======================================================================

LVmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Ingestion    <- rIng  * Prey * Predator
    GrowthPrey   <- rGrow * Prey * (1 - Prey/K)
    MortPredator <- rMort * Predator

    dPrey        <- GrowthPrey - Ingestion
    dPredator    <- Ingestion * assEff - MortPredator

    return(list(c(dPrey, dPredator)))
  })
}

pars  <- c(rIng   = 0.2,    # /day, rate of ingestion
           rGrow  = 1.0,    # /day, growth rate of prey
           rMort  = 0.2 ,   # /day, mortality rate of predator
           assEff = 0.5,    # -, assimilation efficiency
           K      = 10)     # mmol/m3, carrying capacity

yini  <- c(Prey = 1, Predator = 2)
times <- seq(0, 200, by = 1)
out   <- ode(yini, times, LVmod, pars)
summary(out)

## Default plot method
plot(out)

## User specified plotting
matplot(out[ , 1], out[ , 2:3], type = "l", xlab = "time", ylab = "Conc",
        main = "Lotka-Volterra", lwd = 2)
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
parms <- c(b = 0.001, c = 0.1, d = 0.1, e = 0.1, f = 0.1, g = 0.0)

## vector of timesteps
times <- seq(0, 200, length = 101)

## external signal with rectangle impulse
signal <- data.frame(times = times,
                     import = rep(0, length(times)))

signal$import[signal$times >= 10 & signal$times <= 11] <- 0.2

sigimp <- approxfun(signal$times, signal$import, rule = 2)

## Start values for steady state
xstart <- c(S = 1, P = 1, C = 1)

## Solve model
out <- ode(y = xstart, times = times,
           func = SPCmod, parms = parms, input = sigimp)

## Default plot method
plot(out)

## User specified plotting
mf <- par(mfrow = c(1, 2))
matplot(out[,1], out[,2:4], type = "l", xlab = "time", ylab = "state")
legend("topright", col = 1:3, lty = 1:3, legend = c("S", "P", "C"))
plot(out[,"P"], out[,"C"], type = "l", lwd = 2, xlab = "producer",
  ylab = "consumer")
par(mfrow = mf)

## =======================================================================
## Example3: Discrete time model - using method = "iteration"
##           The host-parasitoid model from Soetaert and Herman, 2009, 
##           Springer - p. 284.
## =======================================================================

Parasite <- function(t, y, ks) {
  P <- y[1]
  H <- y[2]
  f    <- A * P / (ks + H)
  Pnew <- H * (1 - exp(-f))
  Hnew <- H * exp(rH * (1 - H) - f)
  
  list (c(Pnew, Hnew))
}
rH <- 2.82 # rate of increase
A  <- 100  # attack rate
ks <- 15   # half-saturation density

out <- ode(func = Parasite, y = c(P = 0.5, H = 0.5), times = 0:50, parms = ks,
           method = "iteration")
            
out2<- ode(func = Parasite, y = c(P = 0.5, H = 0.5), times = 0:50, parms = 25,
           method = "iteration")

out3<- ode(func = Parasite, y = c(P = 0.5, H = 0.5), times = 0:50, parms = 35,
           method = "iteration")

## Plot all 3 scenarios in one figure
plot(out, out2, out3, lty = 1, lwd = 2)

## Same like "out", but *output* every two steps
## hini = 1 ensures that the same *internal* timestep of 1 is used
outb <- ode(func = Parasite, y = c(P = 0.5, H = 0.5),
            times = seq(0, 50, 2), hini = 1, parms = ks,
            method = "iteration")
plot(out, outb, type = c("l", "p"))
}

