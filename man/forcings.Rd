\name{forcings}
\alias{forcings}
\title{
  Passing forcing functions to models written in R or compiled code.
}
\description{
  Models may be defined in compiled C or FORTRAN code, as well as
  in an R-function.
                                   
  If a model uses "external variables" (or "forcing variables"), their value
  at each time point needs to be estimated by interpolation of a data series.
  
  If the models are defined in R-code, then forcing function interpolation 
  is most efficiently done using \R's function \code{approxfun}. 
  
  If the models are defined in compiled C or FORTRAN code, it is possible to
  use \code{deSolve}s forcing function update algorithm. This is the
  compiled-code equivalent of \code{approxfun} or \code{approx}, 
}

\details{

  If the model is defined in \emph{R code}, it is most efficient to:
  
  1. define  a function that performs the linear interpolation,
  using \R's \link[stats]{approxfun}.
  
  2. call this function within the model's derivative function, 
  to interpolate at the current timestep.
  
  See first example
  
  If the model is defined in \emph{compiled code}, the forcing function 
  data series is provided by means of argument \code{forcings}, 
  \code{initforc} is the name of the forcing function initialisation function, 
  as provided in \file{dllname}, while \code{fcontrol} is a list used to 
  finetune how the forcing update should be performed.

  The \bold{fcontrol} argument is a list that can supply any of the following
  components (conform the definitions in the \link[stats]{approxfun} function):
  \describe{
    \item{method }{specifies the interpolation method to be used.
      Choices are "linear" or "constant",}
    \item{rule }{an integer describing how interpolation is to take place
      outside the interval [min(times), max(times)].
      If \code{rule} is \code{1} then an error will be triggered and the
      calculation will stop if \code{times} extends the interval of the
      forcing function data set. If it is \code{2}, the *default*, the
      value at the closest data extreme is used, a warning will be printed if
      \code{verbose} is TRUE,

      Note that the default differs from the \code{approx} default}
    \item{f }{For method=\code{"constant"} a number between \code{0} and
      \code{1} inclusive, indicating a compromise between left- and
      right-continuous step functions. If \code{y0} and \code{y1} are the
      values to the left and right of the point then the value is
      \code{y0*(1-f)+y1*f} so that \code{f=0} is right-continuous and
      \code{f=1} is left-continuous,
      }
    \item{ties }{Handling of tied \code{times} values. Either a function
      with a single vector argument returning a single number result or the
      string "ordered".

      Note that the default is "ordered", hence the existence of ties will
      NOT be investigated; in the \code{C} code this will mean that -if ties
      exist, the first value will be used; if the dataset is not ordered,
      then nonsense will be produced.

      Alternative values for \code{ties} are \code{mean}, \code{min} etc
      }
   }
   The defaults are:

   \code{fcontrol=list(method="linear", rule = 2,  f = 0, ties = "ordered")}

   Note that only ONE specification is allowed, even if there is more than
   one forcing function data set.

   More information about models defined in compiled code is in the package
   vignette ("compiledCode").
}
\note{
  How to write compiled code is described in package vignette
  \code{"compiledCode"}, which should be referred to for details.

  This vignette also contains examples on how to pass forcing functions.
}
\author{
  Karline Soetaert,
  
  Thomas Petzoldt,

  R. Woodrow Setzer (Maintainer)
}
\seealso{
  \code{approx} or \link[stats]{approxfun}, the \R function
}
\examples{
## =============================================================================
# The sediment oxygen consumption example - R-code:
## =============================================================================

# Forcing function data
Flux <- matrix(ncol=2,byrow=TRUE,data=c(
  1, 0.654, 11, 0.167,   21, 0.060, 41, 0.070, 73,0.277, 83,0.186,
  93,0.140,103, 0.255,  113, 0.231,123, 0.309,133,1.127,143,1.923,
  153,1.091,163,1.001,  173, 1.691,183, 1.404,194,1.226,204,0.767,
  214, 0.893,224,0.737, 234,0.772,244, 0.726,254,0.624,264,0.439,
  274,0.168,284 ,0.280, 294,0.202,304, 0.193,315,0.286,325,0.599,
  335, 1.889,345, 0.996,355,0.681,365,1.135))

parms <- c(k=0.01)

times <- 1:365


# the model
sediment <- function( t, O2, k) 
  list (c(Depo(t) - k*O2), depo = Depo(t))

# the forcing functions
Depo <- approxfun(x=Flux[,1], y= Flux[,2], method="linear")

Out <- ode(times = times, func= sediment, y = c(O2=63), parms=parms)
  
# same forcing functions, now constant interpolation
Depo <- approxfun(x=Flux[,1], y= Flux[,2], method="constant", f=0.5)

Out2 <- ode(times = times, func= sediment, y = c(O2=63), parms=parms)

mf <- par(mfrow=c(2,1))
plot (Out, which = "depo", type="l", lwd=2, mfrow=NULL)
lines(Out2[,"time"], Out2[,"depo"],col="red",lwd=2)

plot (Out, which = "O2", type="l", lwd=2, mfrow=NULL)
lines(Out2[,"time"], Out2[,"O2"],col="red",lwd=2)

## =============================================================================
## SCOC is the same model, as implemented in FORTRAN
## =============================================================================

out<- SCOC(times,parms=parms,Flux=Flux)

plot(out$time,out$Depo,type="l",col="red")
lines(out$time,out$Mineralisation,col="blue")

# Constant interpolation of forcing function - left side of interval
fcontrol<-list(method="constant")

out2 <- SCOC(times,parms=parms,Flux=Flux, fcontrol=fcontrol)

plot(out2$time,out2$Depo,type="l",col="red")
lines(out2$time,out2$Mineralisation,col="blue")

\dontrun{
## =============================================================================
## show examples (see respective help pages for details)
## =============================================================================

example(aquaphy)

## show package vignette with tutorial about how to use compiled models
## + source code of the vignette
## + directory with C and FORTRAN sources
vignette("compiledCode")
edit(vignette("compiledCode"))
browseURL(paste(system.file(package = "deSolve"), "/doc", sep = ""))
}

}

\keyword{utilities}