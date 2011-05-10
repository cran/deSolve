\name{events}
\alias{events}
\title{
  Implementing Events in Differential Equation Models.
}
\description{
  An \code{event} occurs when the value of a state variable is suddenly
  changed, e.g. because a value is added, subtracted, or multiplied. The
  integration routines cannot deal easily with such state variable
  changes. Typically these events occur only at specific times.  In
  \code{deSolve}, events can be imposed by means of an input data.frame,
  that specifies at which time and how a certain state variable is altered,
  or via an event function.
}
\details{
  
  The \code{events} are specified by means of argument \code{events}
  passed to the integration routines.
   
  \code{events} should be a list that contains one of the following:

  \enumerate{
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
  
   In case the events are specified by means of an \R \emph{**FUNCTION**}
   (argument \code{events$func}),
   it must be defined as: \code{function(t, y, parms, ...)}.
   \code{t} is the current time point in the integration,
   \code{y} is the current estimate of the variables in the ODE system.
   If the initial values \code{y} has a \code{names} attribute, the
   names will be available inside \code{events$func}.  \code{parms} is a
   vector or list of parameters; \code{...} (optional) are any other
   arguments passed to the function.  The function should return the
   modified y-values, as a \emph{vector}.

   If \code{events$func} is a string, this indicates that the events are
   specified by a \code{function} in compiled code. This function has as
   arguments, the number of state variables, the time, and the state
   variable vector. See package vignette "compiledCode" for more details.

   Also, if events are specified by an \R-function,
   this requires either input of the \emph{time} of the events, a vector in
   \code{events$time} OR the specification of a \emph{root} function.  In the
   latter case, the model must be solved with an integration routine
   with root-finding capability
   (note that from version 1.9.1 the following routines have root-finding capability:
   lsoda, lsode, lsodes, radau).
   The root function itself should be  specified with argument \code{rootfunc}.
   In this case, the integrator is informed that the simulation it to be
   continued after a root is found by
   setting \code{events$root} equal to \code{TRUE}.
  
   If roots are found, then the output will have attribute
   \code{troot} which will contain the \code{times} at which a root was found.
   There will be at most \code{events$maxroot} such values. The default is 100.
  
   See two last examples; also see example of \code{\link{ccl4model}}.
  
   If the evenst are specified by a **\code{data.frame}**
   (argument \code{events$data}), this should
   contain the following columns (and in that order):
  \enumerate{
    \item{var }{the state variable \emph{name} or \emph{number} that is affected
      by the event}
    \item{time }{the time at which the event is to take place; the solvers will
      check if the time is embraced by the simulation time}
    \item{value }{the value, magnitude of the event}
    \item{method }{which event is to take place; should be one of ("replace", 
      "add", "multiply"); also allowed is to specify the number 
      (1=replace, 2=add, 3=multiply)
    }
  }

  For instance, the following line 
  
  "v1"  10  2  "add"
  
  will cause the value 2 to be added to a state variable, called "v1" at 
  time = 10.
  
}
\author{
  Karline Soetaert,
}
\seealso{
  \link{forcings}, for how to implement forcing functions.
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
           events = list(func = eventfun, time = 1:9) )
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
}

\keyword{utilities}