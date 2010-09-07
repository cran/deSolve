\name{rk}
\alias{rk}

\title{Explicit One-Step Solvers for Ordinary Differential Equations
(ODE)}

\description{Solving initial value problems for non-stiff systems of
  first-order ordinary differential equations (ODEs).

  The \R function \code{rk} is a top-level function that provides
  interfaces to a collection of common explicit one-step solvers of the
  Runge-Kutta family with fixed or variable time steps.

  The system of ODE's is written as an \R function (which may, of
  course, use \code{\link{.C}}, \code{\link{.Fortran}},
  \code{\link{.Call}}, etc., to call foreign code) or be defined in
  compiled code that has been dynamically loaded.  A vector of
  parameters is passed to the ODEs, so the solver may be used as part of
  a modeling package for ODEs, or for parameter estimation using any
  appropriate modeling tool for non-linear models in \R such as
  \code{\link{optim}}, \code{\link{nls}}, \code{\link{nlm}} or
  \code{\link[nlme]{nlme}}
}
\usage{
rk(y, times, func, parms, rtol = 1e-6, atol = 1e-6,
  verbose = FALSE, tcrit = NULL, hmin = 0, hmax = NULL,
  hini = hmax, ynames = TRUE, method = rkMethod("rk45dp7", ... ),
  maxsteps = 5000, dllname = NULL, initfunc = dllname,
  initpar = parms, rpar = NULL, ipar = NULL,
  nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, events = NULL, ...)
}
\arguments{
  \item{y }{the initial (state) values for the ODE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{times }{times at which explicit estimates for \code{y} are
    desired.  The first value in \code{times} must be the initial time.
  }
  \item{func }{either an \R-function that computes the values of the
    derivatives in the ODE system (the \emph{model definition}) at time
    t, or a character string giving the name of a compiled function in a
    dynamically loaded shared library.

    If \code{func} is an \R-function, it must be defined as:
     \code{func <- function(t, y, parms,...)}.  \code{t} is the current
    time point in the integration, \code{y} is the current estimate of
    the variables in the ODE system.  If the initial values \code{y} has
    a \code{names} attribute, the names will be available inside \code{func}.
    \code{parms} is a vector or list of parameters; ... (optional) are
    any other arguments passed to the function.

    The return value of \code{func} should be a list, whose first
    element is a vector containing the derivatives of \code{y} with
    respect to \code{time}, and whose next elements are global values
    that are required at each point in \code{times}.  The derivatives
    should be specified in the same order as the state variables \code{y}.

    If \code{func} is
    a string, then \code{dllname} must give the name of the shared
    library (without extension) which must be loaded before
    \code{rk} is called. See package vignette \code{"compiledCode"}
    for more details.
    
  }
  \item{parms }{vector or list of parameters used in \code{func}.
  }
  \item{rtol }{relative error tolerance, either a scalar or an array as
    long as \code{y}. Only applicable to methods with variable time
    step, see details.
  }
  \item{atol }{absolute error tolerance, either a scalar or an array as
    long as \code{y}. Only applicable to methods with variable time
    step, see details.
  }
  \item{tcrit }{if not \code{NULL}, then \code{rk} cannot integrate past
    \code{tcrit}. This parameter is for compatibility with other solvers.
  }
  \item{verbose }{a logical value that, when TRUE, triggers more
    verbose output from the ODE solver.
  }
  \item{hmin }{an optional minimum value of the integration stepsize. In
    special situations this parameter may speed up computations with the
    cost of precision. Don't use \code{hmin} if you don't know why!
  }
  \item{hmax }{an optional maximum value of the integration stepsize. If
    not specified, \code{hmax} is set to the maximum of \code{hini} and
    the largest difference in \code{times}, to avoid that the simulation
    possibly ignores short-term events.  If 0, no maximal size is
    specified. Note that \code{hmin} and \code{hmax} are ignored by
    fixed step methods like \code{"rk4"} or \code{"euler"}.
  }
  \item{hini }{initial step size to be attempted; if 0, the initial step
    size is determined automatically by solvers with flexible time step.
    Setting \code{hini = 0} for fixed step methods forces setting of
    internal time steps identically to external time steps provided by
    \code{times}.
  }
  \item{ynames }{if \code{FALSE}: names of state variables are not passed
    to function \code{func} ; this may speed up the simulation especially
    for large models.
  }
  
  \item{method }{the integrator to use. This can either be a string
    constant naming one of the pre-defined methods or a call to function
    \code{\link{rkMethod}} specifying a user-defined method.  The most
    common methods are the fixed-step methods \code{"euler"}, second and
    fourth-order Runge Kutta (\code{"rk2"}, \code{"rk4"}), or the
    variable step methods Bogacki-Shampine \code{"rk23bs"},
    Runge-Kutta-Fehlberg \code{"rk34f"}, the fifth-order Cash-Karp
    method \code{"rk45ck"} or the fifth-order Dormand-Prince method with
    seven stages \code{"rk45dp7"}.
    As a suggestion, one may use \code{"rk23"} (alias \code{"ode23"}) for
    simple problems and \code{"rk45dp7"} (alias \code{"ode45"}) for
    rough problems.
  }
  \item{maxsteps }{average maximal number of steps per output interval
    taken by the solver. This argument is defined such as to ensure
    compatibility with the Livermore-solvers. \code{rk} only accepts the maximal
    number of steps for the entire integration, and this is calculated 
    as \code{length(times) * maxsteps}.
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions refered to in \code{func} and
    \code{jacfunc}. See package vignette \code{"compiledCode"}.
  }
  \item{initfunc }{if not \code{NULL}, the name of the initialisation function
    (which initialises values of parameters), as provided in
    \file{dllname}. See package vignette \code{"compiledCode"}.
  }
  \item{initpar }{only when \file{dllname} is specified and an
    initialisation function \code{initfunc} is in the dll: the
    parameters passed to the initialiser, to initialise the common
    blocks (FORTRAN) or global variables (C, C++).
  }
  \item{rpar }{only when \file{dllname} is specified: a vector with
    double precision values passed to the dll-functions whose names are
    specified by \code{func} and \code{jacfunc}.
  }
  \item{ipar }{only when \file{dllname} is specified: a vector with
    integer values passed to the dll-functions whose names are specified
    by \code{func} and \code{jacfunc}.
  }
  \item{nout }{only used if \code{dllname} is specified and the model is
    defined in compiled code: the number of output variables calculated
    in the compiled function \code{func}, present in the shared
    library. Note: it is not automatically checked whether this is
    indeed the number of output variables calculated in the dll - you have
    to perform this check in the code. See package vignette
    \code{"compiledCode"}.
  }
  \item{outnames }{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{func}, present in the shared library.
  }
  \item{forcings }{only used if \file{dllname} is specified: a list with
    the forcing function data sets, each present as a two-columned matrix,
    with (time,value); interpolation outside the interval
    [min(\code{times}), max(\code{times})] is done by taking the value at
    the closest data extreme.

    See \link{forcings} or package vignette \code{"compiledCode"}.
  }
  \item{initforc }{if not \code{NULL}, the name of the forcing function
    initialisation function, as provided in
    \file{dllname}. It MUST be present if \code{forcings} has been given a
    value.
    See \link{forcings} or package vignette \code{"compiledCode"}.
  }
  \item{fcontrol }{A list of control parameters for the forcing functions.
    See \link{forcings} or vignette \code{compiledCode}.
  }
  \item{events }{A list that specifies events, i.e. when the value of a 
    state variable is suddenly changed. See \link{events} for more information.
    Not also that if events are specified, then polynomial interpolation
    is switched off and integration takes place from one external time step 
    to the next, with an internal step size less than or equal the difference
    of two adjacent points of \code{times}.
  }
  \item{... }{additional arguments passed to \code{func} allowing this
    to be a generic function.
  }
}

\details{
  Function \code{rk} is a generalized implementation that can be used to
  evaluate different solvers of the Runge-Kutta family of explicit ODE
  solvers. A pre-defined set of common method parameters is in function
  \code{\link{rkMethod}} which also allows to supply user-defined
  Butcher tables.

  The input parameters \code{rtol}, and \code{atol} determine the error
  control performed by the solver.  The solver will control the vector
  of estimated local errors in \bold{y}, according to an inequality of
  the form max-norm of ( \bold{e}/\bold{ewt} ) \eqn{\leq}{ <= } 1, where
  \bold{ewt} is a vector of positive error weights.  The values of
  \code{rtol} and \code{atol} should all be non-negative.  The form of
  \bold{ewt} is:

  \deqn{\mathbf{rtol} \times \mathrm{abs}(\mathbf{y}) +
  \mathbf{atol}}{\bold{rtol} * abs(\bold{y}) + \bold{atol}}

  where multiplication of two vectors is element-by-element.

  \bold{Models} can be defined in \R as a user-supplied
  \bold{R-function}, that must be called as: \code{yprime = func(t, y,
  parms)}.  \code{t} is the current time point in the integration,
  \code{y} is the current estimate of the variables in the ODE system.

  The return value of \code{func} should be a list, whose first element
  is a vector containing the derivatives of \code{y} with respect to
  time, and whose second element contains output variables that are
  required at each point in time. Examples are given below.  
}
\value{
  A matrix of class \code{deSolve} with up to as many rows as elements
  in \code{times} and as many columns as elements in \code{y} plus the
  number of "global" values returned in the next elements of the return
  from \code{func}, plus and additional column for the time value.
  There will be a row for each element in \code{times} unless the
  integration routine returns with an unrecoverable error. If \code{y}
  has a names attribute, it will be used to label the columns of the
  output value.
}
\note{  
  Arguments \code{rpar} and \code{ipar} are provided for compatibility
  with \code{lsoda}.

  Starting with version 1.8 implicit Runge-Kutta methods are also
  supported by this general \code{rk} interface, however their
  implementation is still experimental.  Instead of this you may
  consider \code{\link{radau}} for a specific full implementation of an
  implicit Runge-Kutta method.
}
\references{
  Butcher, J. C. (1987) The numerical analysis of ordinary differential
  equations, Runge-Kutta and general linear methods, Wiley, Chichester
  and New York.

  Engeln-Muellges, G. and Reutter, F. (1996) Numerik Algorithmen:
  Entscheidungshilfe zur Auswahl und Nutzung. VDI Verlag, Duesseldorf.

  Hindmarsh, Alan C. (1983) ODEPACK, A Systematized Collection of ODE
  Solvers; in p.55--64 of Stepleman, R.W. et al.[ed.] (1983)
  \emph{Scientific Computing}, North-Holland, Amsterdam.

  Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery,
  B. P. (2007) Numerical Recipes in C. Cambridge University Press.
}

\author{Thomas Petzoldt \email{thomas.petzoldt@tu-dresden.de}}

\seealso{  
  For most practical cases, solvers of the Livermore family (i.e. the
  ODEPACK solvers, see below) are superior. Some of them are also
  suitable for stiff ODEs, differential algebraic equations (DAEs), or
  partial differential equations (PDEs).

  \itemize{
    \item \code{\link{rkMethod}} for a list of available Runge-Kutta
      parameter sets,
    \item \code{\link{rk4}} and \code{\link{euler}} for special
      versions without interpolation (and less overhead),
    \item \code{\link{lsoda}}, \code{\link{lsode}},
      \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}},
      \code{\link{daspk}} for solvers of the Livermore family,
    \item \code{\link{ode}} for a general interface to most of the ODE solvers,
    \item \code{\link{ode.band}} for solving models with a banded
      Jacobian,
    \item \code{\link{ode.1D}} for integrating 1-D models,
    \item \code{\link{ode.2D}} for integrating 2-D models,
    \item \code{\link{ode.3D}} for integrating 3-D models,
  }
  \code{\link{diagnostics}} to print diagnostic messages.
}


\examples{
## =======================================================================
## Example: Resource-producer-consumer Lotka-Volterra model
## =======================================================================

## Note:
## parameters are a list, names accessible via "with" function
## (see also ode and lsoda examples)

SPCmod <- function(t, x, parms) {
  S <- x[1] # substrate
  P <- x[2] # producer
  C <- x[3] # consumer

  with(parms, {
    import <- approx(signal$times, signal$import, t)$y
    dS <- import - b * S * P + g * C
    dP <- c * S * P  - d * C * P
    dC <- e * P * C  - f * C
    res <- c(dS, dP, dC)
    list(res)
  })
}

## vector of timesteps
times  <- seq(0, 100, length = 101)

## external signal with rectangle impulse
signal <- as.data.frame(list(times = times,
                            import = rep(0,length(times))))

signal$import[signal$times >= 10 & signal$times <= 11] <- 0.2

## Parameters for steady state conditions
parms <- list(b = 0.0, c = 0.1, d = 0.1, e = 0.1, f = 0.1, g = 0.0)

## Start values for steady state
y <- xstart <- c(S = 1, P = 1, C = 1)

system.time({
## Euler method
out1  <- as.data.frame(rk(xstart, times, SPCmod, parms,
                          hini = 0.1, method = "euler"))

## classical Runge-Kutta 4th order
out2 <- as.data.frame(rk(xstart, times, SPCmod, parms,
                         hini = 1, method = "rk4"))

## Dormand-Prince method of order 5(4)
out3 <- as.data.frame(rk(xstart, times, SPCmod, parms,
                         hmax = 1, method = "rk45dp7"))
})

mf <- par(mfrow = c(2,2))
plot (out1$time, out1$S, type = "l", ylab = "Substrate")
lines(out2$time, out2$S, col = "red", lty = "dotted", lwd = 2)
lines(out3$time, out3$S, col = "green", lty = "dotted")

plot (out1$time, out1$P, type = "l", ylab = "Producer")
lines(out2$time, out2$P, col = "red", lty = "dotted")
lines(out3$time, out3$P, col = "green", lty = "dotted")

plot (out1$time, out1$C, type = "l", ylab = "Consumer")
lines(out2$time, out2$C, col = "red", lty = "dotted", lwd = 2)
lines(out3$time, out3$C, col = "green", lty = "dotted")

plot (out1$P, out1$C, type = "l", xlab = "Producer", ylab = "Consumer")
lines(out2$P, out2$C, col = "red",   lty = "dotted", lwd = 2)
lines(out3$P, out3$C, col = "green", lty = "dotted")

legend("center", legend = c("euler", "rk4", "rk45dp7"),
  lty = c(1, 3, 3), lwd = c(1, 2, 1),
  col = c("black", "red", "green"))
par(mfrow = mf)
}

\keyword{ math }