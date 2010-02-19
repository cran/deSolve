\name{SCOC}
\alias{SCOC}
\title{A Sediment Model of Oxygen Consumption}
\description{A model that describes oxygen consumption in a marine sediment

  One state variable:
  \itemize{
    \item sedimentary organic carbon,
  }

  Organic carbon settles on the sediment surface (forcing function Flux)
  and decays at a constant rate.

  The equation is simple:
  
  \deqn{\frac{dC}{dt}=Flux-k C}

  This model is written in \code{FORTRAN}.
}
\usage{SCOC(times, y=NULL, parms, Flux, ...)}
\arguments{
  \item{times}{time sequence for which output is wanted; the first value
    of times must be the initial time,}
  \item{y}{the initial value of the state variable; if \code{NULL} it
    will be estimated based on \code{Flux} and \code{parms},}
  \item{parms }{the model parameter, \code{k},}
  \item{Flux }{a data set with the organic carbon deposition rates, }
  \item{...}{any other parameters passed to the integrator \code{ode}
    (which solves the model).}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{

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
out<- SCOC(times,parms=parms,Flux=Flux)

plot(out$time,out$Depo,type="l",col="red")
lines(out$time,out$Mineralisation,col="blue")

# Constant interpolation of forcing function - left side of interval
fcontrol<-list(method="constant")

out2 <- SCOC(times,parms=parms,Flux=Flux, fcontrol=fcontrol)

plot(out2$time,out2$Depo,type="l",col="red")
lines(out2$time,out2$Mineralisation,col="blue")


}
\references{
Soetaert, K. and P.M.J. Herman, 2009.
A Practical Guide to Ecological Modelling. Using \R as a Simulation Platform.
Springer, 372 pp.
}
\details{
  The model is implemented primarily to demonstrate the linking of
  FORTRAN with R-code.

  The source can be found in the \file{dynload} subdirectory of the package.
}
\seealso{
  \code{\link{ccl4model}}, the CCl4 inhalation model.

  \code{\link{aquaphy}}, the algal growth model.

}
\keyword{models}

