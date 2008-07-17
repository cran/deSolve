\name{ode}
\alias{ode}
\title{General solver for ordinary differential equations}
\description{Solves a system of ordinary differential equations.}
\usage{ode(y, times, func, parms, 
method=c("lsoda","lsode","lsodes","lsodar","vode","daspk", "euler", "rk4", 
         "ode23", "ode45"), ...)}
\arguments{
  \item{y }{the initial (state) values for the ODE system, a vector. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{times }{time sequence for which output is wanted; the first value of \code{times} must be the initial time}
  \item{func }{either an \R-function that computes the values of the derivatives in the ODE system (the model definition) at time t, or a character string giving the name of a compiled function in a dynamically loaded shared library. 

      If \code{func}  is an \R-function, it must be defined as:
      \code{yprime = func(t, y, parms,...)}.  \code{t} is the current time point
      in the integration, \code{y} is the current estimate of the variables
      in the ODE system.  If the initial values \code{y} has a names
      attribute, the names will be available inside \code{func}.  \code{parms} is
      a vector or list of parameters; ... (optional) are any other arguments passed to the function. 
  
      The return value of \code{func} should be a list, whose first element is a
      vector containing the derivatives of \code{y} with respect to
      \code{time}, and whose next elements are global values that are required at
      each point in \code{times}. }
  \item{parms }{parameters passed to \code{func}}
  \item{method }{the integrator to use, either a string ("lsoda","lsode","lsodes",
    "lsodar","vode", "daspk", "euler", "rk4", "ode23" or "ode45") or a function
     that performs integration, or a list of class \code{rkMethod}.}
  \item{... }{additional arguments passed to the integrator}
}
\value{
  A matrix with up to as many rows as elements in times and as many columns as elements in \code{y} plus the number of "global" values returned 
    in the second element of the return from \code{func}, plus an additional column (the first) for the time value. 
    There will be one row for each element in \code{times} unless the integrator returns with an unrecoverable error. 
    If \code{y} has a names attribute, it will be used to label the columns of the output value. 
    
    The output will have the attributes \code{istate}, and \code{rstate}, two vectors with several useful elements. 
    The first element of istate returns the conditions under which the last call to the integrator returned. Normal is istate = 2.
    If \code{verbose} = TRUE, the settings of istate and rstate will be written to the screen. See the help for the selected integrator for details.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\details{This is simply a wrapper around the various ode solvers.
\cr See package vignette for information about specifying the model in compiled code.
\cr See the selected integrator for the additional options}
\seealso{\item \code{\link{ode.band}} for solving models with a banded Jacobian
\item  \code{\link{ode.1D}} for integrating 1-D models 
\item \code{\link{ode.2D}} for integrating 2-D models
 	\item \code{\link{aquaphy}}, \code{\link{ccl4model}}, where \code{ode} is used

	\item \code{\link{lsoda}}, \code{\link{lsode}}, \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}}, \code{\link{daspk}},
  \item  \code{\link{rk}}, \code{\link{rkMethod}}}
\keyword{math}
\examples{

#########################################
## Example1: Pred-Prey Lotka-volterra model
#########################################

LVmod <- function(Time,State,Pars)
 {

   with(as.list(c(State,Pars)),

    {
    Ingestion    <- rIng  * Prey*Predator
    GrowthPrey   <- rGrow * Prey*(1-Prey/K)
    MortPredator <- rMort * Predator

    dPrey        <- GrowthPrey - Ingestion
    dPredator    <- Ingestion*assEff -MortPredator

    return(list(c( dPrey, dPredator)))

    })
 }


pars    <- c(rIng   =0.2,    # /day, rate of ingestion
             rGrow  =1.0,    # /day, growth rate of prey
             rMort  =0.2 ,   # /day, mortality rate of predator
             assEff =0.5,    # -, assimilation efficiency
             K      =10  )   # mmol/m3, carrying capacity

yini    <- c(Prey=1,Predator=2)
times   <- seq(0,200,by=1)
out     <- as.data.frame(lsoda(func= LVmod, y=yini,
                         parms=pars, times=times))

matplot(out$time,out[,2:3],type="l",xlab="time",ylab="Conc",
        main="Lotka-Volterra",lwd=2)
legend("topright",c("prey", "predator"),col=1:2, lty=1:2)

#########################################
## Example2: Resource-producer-consumer Lotka-volterra model
#########################################

  ## Note:
  ## 1. parameter and state variable names made
  ## accessible via "with" statement
  ## 2. function sigimp passed as an argument (input) to model
  ## (see also lsoda and rk examples)

  lvmodel <- function(t, x, parms, input)  {

      with(as.list(c(parms,x)),  {
      
      import <- input(t)
      dS <- import - b*S*P + g*K    #substrate
      dP <- c*S*P  - d*K*P          #producer
      dK <- e*P*K  - f*K            #consumer
      res<-c(dS, dP, dK)
      list(res)                   
                                  })
    }
  
  ## The parameters 
  parms  <- c(b=0.0, c=0.1, d=0.1, e=0.1, f=0.1, g=0.0)

  ## vector of timesteps
  times  <- seq(0, 100, length=101)
  
  ## external signal with rectangle impulse
  signal <- as.data.frame(list(times = times,
                              import = rep(0,length(times))))
  
  signal$import[signal$times >= 10 & signal$times <=11] <- 0.2
  
  sigimp <- approxfun(signal$times, signal$import, rule=2)
  
  
  ## Start values for steady state
  xstart <- c(S=1, P=1, K=1)
  
  ## Solve model
  out <- as.data.frame(ode(y=xstart,times= times, 
                       func=lvmodel, parms, input =sigimp))
  plot(out$P,out$K,type="l",lwd=2,xlab="producer",ylab="consumer")
}

