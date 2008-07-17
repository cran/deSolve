\name{ode.1D}
\alias{ode.1D}
\title{Solver for multicomponent 1-D ordinary differential equations}
\description{Solves a system of ordinary differential equations resulting from 
1-Dimensional  multi-component transport-reaction models that include transport only between adjacent layers.}
\usage{ode.1D(y, times, func, parms, nspec=NULL, dimens=NULL, 
       method="lsode", ...)}
\arguments{
  \item{y }{the initial (state) values for the ODE system, a vector. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{times }{time sequence for which output is wanted; the first value of \code{times} must be the initial time}
  \item{func }{either an \R-function that computes the values of the derivatives in the ODE system (the model definition) at time t, 
  or a character string giving the name of a compiled function in a dynamically loaded shared library.
  
      If \code{func}  is an \R-function, it must be defined as:
      \code{yprime = func(t, y, parms,...)}.  \code{t} is the current time point
      in the integration, \code{y} is the current estimate of the variables
      in the ODE system.  If the initial values \code{y} has a names
      attribute, the names will be available inside \code{func}.  \code{parms} is
      a vector or list of parameters; ... (optional) are any other arguments passed to the function. 
  
      The return value of \code{func} should be a list, whose first element is a
      vector containing the derivatives of \code{y} with respect to
      \code{time}, and whose next elements are global values that are required at
      each point in \code{times}.
   
       		If \code{func} is a character string then integrator \code{lsodes} will be used. See details}
  \item{parms }{parameters passed to \code{func}}
  \item{nspec }{the number of *species* (components) in the model. If NULL, then \code{dimens} should be specified}
  \item{dimens}{the number of *boxes* in the model. If NULL, then \code{nspec} should be specified}
  \item{method }{the integrator to use, one of "vode", "lsode", "lsoda", "lsodar", "lsodes" }
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
\note{It is advisable though not mandatory to specify BOTH \code{nspec} and \code{dimens}. In this case, the solver can check whether the input makes sense 
(i.e. if nspec*dimens = length(y))}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# example 1

  #=======================================================
  # a predator and its prey diffusing on a flat surface
  # in concentric circles
  # 1-D model with using cylindrical coordinates
  # Lotka-Volterra type biology
  #=======================================================

  #==================#
  # Model equations  #
  #==================#

  lvmod <- function (time, state, parms,N,rr,ri,dr,dri)

  {
    with (as.list(parms),{
    PREY <- state[1:N]
    PRED <- state[(N+1):(2*N)]
    
    # Fluxes due to diffusion 
    # at internal and external boundaries: zero gradient
    FluxPrey <- -Da * diff(c(PREY[1],PREY,PREY[N]))/dri   
    FluxPred <- -Da * diff(c(PRED[1],PRED,PRED[N]))/dri   

    # Biology: Lotka-Volterra model
    Ingestion     <- rIng * PREY*PRED
    GrowthPrey    <- rGrow* PREY*(1-PREY/cap)
    MortPredator  <- rMort* PRED

    # Rate of change = Flux gradient + Biology   
    dPREY    <- -diff(ri * FluxPrey)/rr/dr   +
                GrowthPrey - Ingestion
    dPRED    <- -diff(ri * FluxPred)/rr/dr   +
                Ingestion*assEff -MortPredator

    return (list(c(dPREY,dPRED)))
   }) 
  }
  
  #==================#
  # Model application#
  #==================#
  # model parameters: 

  R  <- 20                    # total radius of surface, m
  N  <- 100                   # 100 concentric circles
  dr <- R/N                   # thickness of each layer
  r  <- seq(dr/2,by=dr,len=N) # distance of center to mid-layer
  ri <- seq(0,by=dr,len=N+1)  # distance to layer interface
  dri<- dr                    # dispersion distances

  parms <- c( Da     =0.05,   # m2/d, dispersion coefficient 
              rIng   =0.2,    # /day, rate of ingestion
              rGrow  =1.0,    # /day, growth rate of prey
              rMort  =0.2 ,   # /day, mortality rate of pred
              assEff =0.5,    # -, assimilation efficiency
              cap    =10  )   # density, carrying capacity

  # Initial conditions: both present in central circle (box 1) only
  state <- rep(0,2*N)
  state[1] <- state[N+1] <- 10
                  
  # RUNNING the model:   #
  times  <-seq(0,200,by=1)   # output wanted at these time intervals           

  # the model is solved by the two implemented methods:
  # 1. Default: banded reformulation
  print(system.time(
  out    <- ode.1D(y=state,times=times,func=lvmod,parms=parms,nspec=2,
                    N=N,rr=r,ri=ri,dr=dr,dri=dri)  
                    ))

  # 2. Using sparse method
  print(system.time(
  out2   <- ode.1D(y=state,times=times,func=lvmod,parms=parms,nspec=2,
                    N=N,rr=r,ri=ri,dr=dr,dri=dri,method="lsodes")  
                    ))

  #==================#
  # Plotting output  #
  #==================#
  # the data in 'out' consist of: 1st col times, 2-N+1: the prey
  # N+2:2*N+1: predators

  PREY   <- out[,2:(N  +1)]

  filled.contour(x=times,y=r,PREY,color= topo.colors,
                 xlab="time, days", ylab= "Distance, m",
                 main="Prey density")
 
  # Example 2.
  #=======================================================
  # Biochemical Oxygen Demand (BOD) and oxygen (O2) dynamics
  # in a river
  #=======================================================
  
  #==================#
  # Model equations  #
  #==================#
  O2BOD <- function(t,state,pars)
  
  {
    BOD <- state[1:N]
    O2  <- state[(N+1):(2*N)]
  
  # BOD dynamics
    FluxBOD <-  v*c(BOD_0,BOD)  # fluxes due to water transport
    FluxO2  <-  v*c(O2_0,O2)
    
    BODrate <- r*BOD            # 1-st order consumption
  
  #rate of change = flux gradient - consumption  + reaeration (O2)
    dBOD         <- -diff(FluxBOD)/dx  - BODrate
    dO2          <- -diff(FluxO2)/dx   - BODrate + p*(O2sat-O2)
  
    return(list(c(dBOD=dBOD,dO2=dO2)))
  
  }    # END O2BOD
   
   
  #==================#
  # Model application#
  #==================#
  # parameters
  dx      <- 25        # grid size of 25 meters
  v       <- 1e3       # velocity, m/day
  x       <- seq(dx/2,5000,by=dx)  # m, distance from river
  N       <- length(x)
  r       <- 0.05      # /day, first-order decay of BOD
  p       <- 0.5       # /day, air-sea exchange rate 
  O2sat   <- 300       # mmol/m3 saturated oxygen conc
  O2_0    <- 200       # mmol/m3 riverine oxygen conc
  BOD_0   <- 1000      # mmol/m3 riverine BOD concentration
  
  # initial conditions:
  state <- c(rep(200,N),rep(200,N))
  times     <- seq(0,20,by=1)
  
  # running the model
  #  step 1  : model spinup
  out       <- ode.1D (y=state,times,O2BOD,parms=NULL,nspec=2)
  
  #==================#
  # Plotting output  #
  #==================#
  # select oxygen (first column of out:time, then BOD, then O2
  O2   <- out[,(N+2):(2*N+1)]
  color= topo.colors
  
  filled.contour(x=times,y=x,O2,color= color,nlevels=50,
                 xlab="time, days", ylab= "Distance from river, m",
                 main="Oxygen")}
\details{This is the method of choice for multi-species 1-dimensional models, that are only subjected to transport between adjacent layers.
    \cr More specifically, this method is to be used if the state variables are arranged per species: 
    \cr A[1],A[2],A[3],....B[1],B[2],B[3],.... (for species A, B)) 

Two methods are implemented.
\item The default method rearranges the state variables as A[1],B[1],...A[2],B[2],...A[3],B[3],.... This reformulation leads to a banded Jacobian with (upper and lower) half bandwidth = number of species. 
Then the selected integrator solves the banded problem.
\item The second method uses \code{lsodes}. Based on the dimension of the problem, the method first calculates the sparsity pattern of the Jacobian, under the assumption 
  that transport is only occurring between adjacent layers. Then \code{lsodes} is called to solve the problem.
  \cr As \code{lsodes} is used to integrate, it may be necessary to specify the length of the real work array, \code{lrw}.
  \cr Although a reasonable guess of \code{lrw} is made, it is possible that this will be too low. 
  In this case, \code{ode.1D} will return with an error message telling
  the size of the work array actually needed. In the second try then, set \code{lrw} equal to this number.

If the model is specified in compiled code (in a DLL), then option 2, based on \code{lsodes} is the only solution method.
 
For single-species 1-D models, use \code{\link{ode.band}}. 
  
  See the selected integrator for the additional options}
\seealso{\code{\link{ode}},
\item \code{\link{ode.band}} for solving models with a banded Jacobian
\item \code{\link{ode.2D}} for integrating 2-D models 
\item \code{\link{lsoda}}, \code{\link{lsode}}, \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}}, \code{\link{daspk}}.}
\keyword{math}

