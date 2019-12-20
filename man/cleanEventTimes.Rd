\name{cleanEventTimes}
\alias{cleanEventTimes}
\alias{nearestEvent}

\title{
Find Nearest Event for Each Time Step and Clean Time Steps to Avoid Doubles
}
\description{
These functions can be used for checking time steps and events used by 
ode solver functions. They are normally called internally within the solvers. 
}
\usage{
nearestEvent(times, eventtimes)
cleanEventTimes(times, eventtimes, eps = .Machine$double.eps * 10)
}
\arguments{
  \item{times}{the vector of output times,}
  \item{eventtimes}{a vector with the event times,}
  \item{eps}{relative tolerance value below which two numbers are assumed 
    to be numerically equal.}
}
\details{
  In floating point arithmetics, problems can occur if values have to be
  compared for 'equality' but are only close to each other and not exactly
  the same.
  
  The utility functions can be used to add all \code{eventtimes} to
  the output \code{times} vector, but without including times that are
  very close to an event.

  This means that all values of \code{eventtimes} are contained
  but only the subset of \code{times} that have no close neighbors in
  \code{eventtimes}.

  These checks are normally performed internally by the integration solvers.
}
\value{
    \code{nearestEvent} returns a vector with the closest events for
      each time step and
  
    \code{cleanEventTimes} returns a vector with the output times
      without all those that are 'very close' to an event.
}

\author{
  Thomas Petzoldt
}

\seealso{
  \code{\link{events}}
}

\examples{

events <- sort(c(0, 2, 3, 4 + 1e-10, 5, 7 - 1e-10,
                 7 + 6e-15, 7.5, 9, 24.9999, 25, 80, 1001, 1e300))
times  <- sort(c(0, 1:7, 4.5, 6.75, 7.5, 9.2, 9.0001, 25, 879, 1e3, 1e300+5))

nearest <- nearestEvent(times, events)
data.frame(times=times, nearest = nearest)

## typical usage: include all events in times after removing values that
## are numerically close together, events have priority 
times
unique_times <- cleanEventTimes(times, events)
newtimes <- sort(c(unique_times, events))
newtimes
}
\keyword{ misc }
