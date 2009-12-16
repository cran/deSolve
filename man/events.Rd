\name{events}
\alias{events}
\title{
  Implementing events in differential equation models.
}
\description{
  An \code{event} occurs when the value of a state variable is suddenly changed,
  e.g. because new substrate is added, or part is removed. The integration 
  routines cannot deal easily with such state variable changes. Typically these
  events occur only at specific times.   
  In \code{deSolve}, events can be imposed by means of an input file that 
  specifies at which time a certain state variable is altered, or via
  an event function. 
}
\details{
  The \code{events} are specified by means of argument \code{events} passed to
   the integration routines. 
   
   \code{events} should be a list that contains one
   of the following:
  \enumerate{
  \item{func: }{an R-function or the name of a function in compiled code that 
    specifies the event,}
  \item{data: }{a data.frame that specifies the variables, times, values and
    types of the events,}
  \item{time: }{when events are specified by a function: the times at which 
    the events take place,}
  \item{root: }{when events are specified by a function and triggered by a root,
    this logical should be set equal to \code{TRUE},}
  \item{ties: }{if events, as specified by a data.frame are "ordered", set to 
    "ordered", the default is "notordered". This will save some computational 
    time.}
  }
  
  If specified by an \code{R-function} (argument \code{events$func}), 
    this requires either input of the time of the events, a vector in 
    \code{events$time} OR the specification of a root function. 
    In the latter case, the model MUST be solved with integration routine 
    \link{lsodar}
  
  the \R-function, must be defined as:
     \code{function(t, y, parms,...)}.  \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system.  If the initial values \code{y} has a
    \code{names} attribute, the names will be available inside \code{events$func}.
    \code{parms} is a vector or list of parameters; ... (optional) are
    any other arguments passed to the function.  The function should return
    the modified y-values, as a \emph{vector}. 
  
  If \code{events$func} is a string, this indicates that the events are 
  specified by a \code{function in compiled code}. This function 
  has as arguments, the number of state variables, the time, and the 
  state variable vector. See package vignette
  "compiledCode" for more details.

  In addition, either the \code{time} at which the events take place should be 
  specified as a vector (\code{event$time}). When the model is solved with
  \link{lsodar}, an event can also be triggered by a root function, as 
  specified with lsodar argument \code{rootfunc}. In this case, \code{lsodar}
  is informed that the simulation it to be continued after a root is found
  by setting \code{events$root} equal to \code{TRUE}. 
  
  See two last examples; also see example of \link{ccl4model}.
  
  If specified by a \code{data.frame} (argument \code{events$data}), this should 
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
## EVENTS in a data.frame
## =============================================================================

# model: rate of change set to 0
eventmod <- function(t,var,parms) {
  list(dvar = rep(0,2))
}

yini <- c(v1 = 1, v2 = 2)
times <- seq(0,10, by=0.1)

#
eventdat <- data.frame(var = c("v1", "v2", "v2", "v1"), time = c(1,1,5,9) ,
  value = c(1,2,3,4), method =c("add", "mult","rep","add"))

eventdat
  
out <- vode(func=eventmod, y=yini, times=times, parms=NULL, 
  events=list(data=eventdat))
plot(out,type="l")

#
eventdat <- data.frame(var = c(rep("v1",10),rep( "v2", 10)), 
  time = c(1:10,1:10), value = runif(20), method =rep("add", 20))

eventdat
out <- ode(func=eventmod, y=yini, times=times, parms=NULL, 
  events=list(data=eventdat))

plot(out,type="l")

## =============================================================================
## EVENTS in a function
## =============================================================================

# model: rate of change v1 = 0, v2 consumed at first-order rate
eventmod <- function(t,var,parms) {
   list(c(0,-0.5*var[2]))
}


# events: add 1 to v1,multiply v2 with 0.8
eventfun <- function(t,y,parms){
  with (as.list(y),{
    v1 <- v1+1
    v2 <- 5*runif(1)
    return(c(v1,v2))
  })
}

yini <- c(v1 = 1, v2 = 2)
times <- seq(0,10, by=0.1)

out <- ode(func=eventmod, y=yini, times=times, parms=NULL, 
  events=list(func=eventfun, time=1:9) )
plot(out,type="l")

## =============================================================================
## EVENTS triggered by a root function
## =============================================================================

# derivative: simple first-order decay
func <- function(t, y, pars) {
  return(list(-0.1*y))
}

# event triggered if state variable =0.5
rootfun <- function (t, y, pars) {
  return(y-0.5)
}

# sets state vaiable = 1                                                  
eventfun <- function(t, y, pars) {
  return(y=1)
}

yini <- 2

times <- seq(0,100,0.1)

# uses lsodar to solve; root =TRUE specifies that the event is triggered by 
# a root.
out <- lsodar(times=times, y=yini, func = func, parms=NULL,
  events=list(func = eventfun, root = TRUE),  
  rootfun = rootfun)

plot(out,type="l")

}

\keyword{utilities}