\name{deSolve-internal}
\alias{timestep}
\title{Internal deSolve Functions}

\description{
  Internal deSolve functions, these are not to be called by the user.
}

\usage{
timestep(prev = TRUE)
} 

\arguments{
  \item{prev }{if \code{TRUE} will return the timestep previously used;
    when \code{FALSE} will return the time step to be currently tried by 
    the integrator.
  }
} 

\details{
  Function \code{timestep} is intended to return the current or next
  timestep of the integration.  It works only under specific
  circumstances and should not be used by the end user.

  Instead of this, please see the example below for a pure
  \R solution.  
}

\seealso{
  \code{\link{diagnostics}} for information about the time steps used,\cr
  \code{\link{lagvalue}} and \code{\link{lagderiv}} that can be used for DDEs. 
}

\examples{
###################################################
### This example shows how to retrieve information 
### about the used time steps.
###################################################

## a function closure  ('lexical scoping')
modelClosure <- function(t0) {
  t.old <- t.act <- t0
  function(t, y, parms) {
    t.old  <<- t.act
    t.act  <<- t
    cat(t, "\t", t - t.old, "\n")
    with (as.list(c(y, parms)), {
      dP <- a * P      -  b * P * K
      dK <- b * P * K  -  c * K
      list(c(dP, dK))
    })
  }
}

model <- modelClosure(0) # initialization

parms <- c(a = 0.1, b = 0.1, c = 0.1)

y <- c(P = 1, K = 2)

out <- ode(y = y, func = model, times = c(0, 2),
 parms = parms, method = "lsoda")

ls() # prove that t.old and t.new are local within 'model' 
}

\keyword{ internal }