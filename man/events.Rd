\name{events}
\alias{events}
\alias{roots}
\title{
  Implementing Events and Roots in Differential Equation Models.
}
\description{
  An \code{event} occurs when the value of a state variable is suddenly
  changed, e.g. because a value is added, subtracted, or multiplied. The
  integration routines cannot deal easily with such state variable
  changes. Typically these events occur only at specific times.  In
  \code{deSolve}, events can be imposed by means of an input data.frame,
  that specifies at which time and how a certain state variable is altered,
  or via an event function.

  Roots occur when a root function becomes zero. By default when a root is
  found, the simulation either stops (no event), or triggers an event.
}
\details{

  The \code{events} are specified by means of argument \code{events}
  passed to the integration routines.

  \code{events} should be a list that contains one of the following:

  \describe{
    \item{func: }{an R-function or the name of a function in compiled
      code that specifies the event,
    }
    \item{data: }{a data.frame that specifies the state variables, times, values and
      types of the events. Note that the event times must also
      be part of the integration output times, else the event will not take
      place. As from version 1.9.1, this is checked by the solver,
      and a warning message is produced if event times are missing in times;
      see also \code{\link{cleanEventTimes}} for utility functions
      to check and solve such issues.
    }
    \item{time: }{when events are specified by an event function: the times at
      which the events take place. Note that these event times must also
      be part of the integration output times exactly, else the event
      would not take place. As from version 1.9.1 this is checked by the solver,
      and an error message produced if event times are missing in times;
      see also \code{\link{cleanEventTimes}} for utility functions
      to check and solve such issues.
    }
    \item{root: }{when events are specified by a function and triggered
      by a root, this logical should be set equal to \code{TRUE}
    }
    \item{terminalroot: }{when events are triggered by a root, the default is
      that the simulation continues after the event is executed. In
      \code{terminalroot}, we can specify which roots should terminate the
      simulation.
    }
    \item{maxroot: }{when \code{root = TRUE}, the maximal number of
      times at with a root is found and that are kept; defaults to 100.  If
      the number of roots > \code{maxroot}, then only the first
      \code{maxroot} will be outputted.
    }
    \item{ties: }{if events, as specified by a data.frame are "ordered",
       set to "ordered", the default is "notordered". This will save
       some computational time.
     }
   }

   In case the events are specified by means of an \R \bold{function}
   (argument \code{events$func}),
   it must be defined as: \code{function(t, y, parms, ...)}.
   \code{t} is the current time point in the integration,
   \code{y} is the current estimate of the variables in the ODE system.
   If the initial values \code{y} has a \code{names} attribute, the
   names will be available inside \code{events$func}.  \code{parms} is a
   vector or list of parameters; \code{...} (optional) are any other
   arguments passed to the function via the call to the integration method.
   The event function should return the y-values (some of which modified),
   as a \emph{vector}.

   If \code{events$func} is a string, this indicates that the events are
   specified by a \code{function} in compiled code. This function has as
   arguments, the number of state variables, the time, and the state
   variable vector. See package vignette "compiledCode" for more details.

   In case events are specified by an \R-function,
   this requires either: input of the \emph{time} of the events, a vector as
   defined in \code{events$time}
   OR the specification of a \emph{root} function.  In the
   latter case, the model must be solved with an integration routine
   with root-finding capability

   The root function itself should be  specified with argument \code{rootfunc}.
   In this case, the integrator is informed that the simulation it to be
   continued after a root is found by
   setting \code{events$root} equal to \code{TRUE}.


   If the events are specified by a \bold{data frame}
   (argument \code{events$data}), this should
   contain the following columns (and in that order):
  \describe{
    \item{var: }{the state variable \emph{name} or \emph{number} that is affected
      by the event}
    \item{time: }{the time at which the event is to take place; the solvers will
      check if the time is embraced by the simulation time}
    \item{value: }{the value, magnitude of the event}
    \item{method: }{which event is to take place; should be one of ("replace",
      "add", "multiply"); also allowed is to specify the number
      (1 = replace, 2 = add, 3 = multiply)
    }
  }

  For instance, the following line

  \code{"v1"  10  2  "add"}

  will cause the value 2 to be added to a state variable, called \code{"v1"} at
  \code{time = 10}.


   From deSolve version 1.9.1 the following routines have \bold{root-finding} capability:
   \link{lsoda}, \link{lsode}, \link{lsodes}, and \link{radau}.
   For the first 3 integration methods, the
   root finding algorithm is based on the algorithm in solver LSODAR, and is
   implemented in FORTRAN. For radau, the root solving algorithm is written
   in C-code, and it works slightly different. Thus, some problems involving roots
   may be more efficiently solved with either lsoda, lsode, or lsodes, while other
   problems are more efficiently solved with radau.

   If a root function is defined, but not an event function, then by default the
   solver will stop at a root. If this is not desirable, e.g. because we want
   to record the position of many roots, then a dummy "event" function can be
   defined which returns the values of the state variables - unaltered.

   If roots and events are combined, and roots are found, then the output will have attribute
   \code{troot} which will contain the \code{times} at which a root was found (and
   the event trigerred).
   There will be at most \code{events$maxroot} such values. The default is 100.

   See two last examples; also see example of \code{\link{ccl4model}}.
}
\author{
  Karline Soetaert
}
\seealso{
  \link{forcings}, for how to implement forcing functions.

  \link{lsodar}, for more examples of roots
}
\examples{
## =============================================================================
## 1. EVENTS in a data.frame
## =============================================================================

## derivative function: derivatives set to 0
derivs <- function(t, var, parms) {
  list(dvar = rep(0, 2))
}

yini <- c(v1 = 1, v2 = 2)
times <- seq(0, 10, by = 0.1)

eventdat <- data.frame(var = c("v1", "v2", "v2", "v1"),
                       time = c(1, 1, 5, 9) ,
                       value = c(1, 2, 3, 4),
                       method = c("add", "mult", "rep", "add"))
eventdat

out <- vode(func = derivs, y = yini, times = times, parms = NULL,
            events = list(data = eventdat))
plot(out)

##
eventdat <- data.frame(var = c(rep("v1", 10), rep("v2", 10)),
                       time = c(1:10, 1:10),
                       value = runif(20),
                       method = rep("add", 20))
eventdat

out <- ode(func = derivs, y = yini, times = times, parms = NULL,
           events = list(data = eventdat))

plot(out)

## =============================================================================
## 2. EVENTS in a function
## =============================================================================

## derivative function: rate of change v1 = 0, v2 reduced at first-order rate
derivs <- function(t, var, parms) {
   list(c(0, -0.5 * var[2]))
}


# events: add 1 to v1, multiply v2 with random number
eventfun <- function(t, y, parms){
  with (as.list(y),{
    v1 <- v1 + 1
    v2 <- 5 * runif(1)
    return(c(v1, v2))
  })
}

yini <- c(v1 = 1, v2 = 2)
times <- seq(0, 10, by = 0.1)

out <- ode(func = derivs, y = yini, times = times, parms = NULL,
           events = list(func = eventfun, time = c(1:9, 2.2, 2.4)) )
plot(out, type = "l")

## =============================================================================
## 3. EVENTS triggered by a root function
## =============================================================================

## derivative: simple first-order decay
derivs <- function(t, y, pars) {
  return(list(-0.1 * y))
}

## event triggered if state variable = 0.5
rootfun <- function (t, y, pars) {
  return(y - 0.5)
}

## sets state variable = 1
eventfun <- function(t, y, pars) {
  return(y = 1)
}

yini <- 2
times <- seq(0, 100, 0.1)

## uses ode to solve; root = TRUE specifies that the event is
## triggered by a root.
out <- ode(times = times, y = yini, func = derivs, parms = NULL,
           events = list(func = eventfun, root = TRUE),
           rootfun = rootfun)

plot(out, type = "l")

## time of the root:
troot <- attributes(out)$troot
points(troot, rep(0.5, length(troot)))


## =============================================================================
## 4. More ROOT examples: Rotation function
## =============================================================================
Rotate <- function(t, x, p )
  list(c( x[2],
         -x[1]  ))

## Root = when second state variable = 0
rootfun <- function(t, x, p) x[2]

## "event" returns state variables unchanged
eventfun <- function(t, x, p) x
times <- seq(from = 0, to = 15, by = 0.1)

## 1. No event: stops at first root
out1 <- ode(func = Rotate, y = c(5, 5), parms = 0,
           times = times, rootfun = rootfun)
tail(out1)

## 2. Continues till end of times and records the roots
out <- ode(func = Rotate, y = c(5, 5), parms = 0,
           times = times, rootfun = rootfun,
           events = list(func = eventfun, root = TRUE) )

plot(out)
troot <- attributes(out)$troot  # time of roots
points(troot,rep(0, length (troot)))

## Multiple roots:  either one of the state variables = 0
root2 <- function(t, x, p) x

out2 <- ode(func = Rotate, y = c(5, 5), parms = 0,
           times = times, rootfun = root2,
           events = list(func = eventfun, root = TRUE) )

plot(out2, which = 2)
troot <- attributes(out2)$troot
indroot <- attributes(out2)$indroot  # which root was found
points(troot, rep(0, length (troot)), col = indroot, pch = 18, cex = 2)

## Multiple roots and stop at first time root 1.
out3 <- ode(func = Rotate, y = c(5, 5), parms = 0,
      times = times, rootfun = root2,
      events = list(func = eventfun, root = TRUE, terminalroot = 1))


## =============================================================================
## 5. Stop at 5th root - only works with radau.
## =============================================================================
Rotate <- function(t, x, p )
  list(c( x[2],
         -x[1],
         0  ))

## Root = when second state variable = 0
root3  <- function(t, x, p)  c(x[2], x[3] - 5)
event3 <- function (t, x, p) c(x[1:2], x[3]+1)
times <- seq(0, 15, 0.1)
out3 <- ode(func = Rotate, y = c(x1 = 5, x2 = 5, nroot = 0),
      parms = 0, method = "radau",
      times = times, rootfun = root3,
      events = list(func = event3, root = TRUE, terminalroot = 2))
plot(out3)
attributes(out3)[c("troot", "nroot", "indroot")]

## =============================================================================
## 6 Event in R-code, model function in compiled code - based on vode example
## =============================================================================

times <- 1:365
Flux <- cbind(times, sin(pi*times/365)^2) # forcing function

# run without events
out <- ode(y = c(C = 1), times, func = "scocder", parms = c(k=0.01),
   dllname = "deSolve", initforc = "scocforc", forcings = Flux,
   initfunc = "scocpar", nout = 2, outnames = c("Mineralisation", "Depo"))

# Event halves the concentration
EventMin <- function(t, y , p) y/2

out2 <- ode(y = c(C = 1), times, func = "scocder", parms = c(k=0.01),
   dllname = "deSolve", initforc = "scocforc", forcings = Flux,
   initfunc = "scocpar", nout = 2, outnames = c("Mineralisation", "Depo"),
   events = list (func = EventMin, time = c(50.1, 200, 210.5)))

plot(out, out2)


}

\keyword{utilities}