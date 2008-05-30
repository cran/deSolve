\name{rk}
\alias{rk}


\title{Explicit One-Step Solvers for Ordinary Differential Equations (ODE)}
\description{
  Solving initial value problems for 
  non-stiff systems of first-order ordinary differential equations
  (ODEs). 
  
  The \R function \code{rk} is a top-level function that provides interfaces
  to a collection of common explicit one-step solvers of the
  Runge-Kutta family with fixed step or variable time steps.
  
  The system of ODE's is written as an \R function (which
  may, of course, use \code{\link{.C}}, \code{\link{.Fortran}},
  \code{\link{.Call}}, etc., to call foreign code). A vector of
  parameters is passed to the ODEs, so the solver may be used as part of
  a modeling package for ODEs, or for parameter estimation using any
  appropriate modeling tool for non-linear models in \R such as
  \code{\link{optim}}, \code{\link[nls]{nls}}, \code{\link{nlm}} or
  \code{\link[nlme]{nlme}}
    
}
\usage{
rk(y, times, func, parms, rtol = 1e-06, atol = 1e-06, 
  tcrit = NULL, verbose = FALSE, hmin = 0, hmax = NULL, hini = hmax, 
  method = rkMethod("rk45dp7", ...), maxsteps = 5000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y }{the initial (state) values for the ODE system. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{times }{times at which explicit estimates for \code{y} are desired.  The first value in \code{times} must be the initial time.}
  \item{func }{an \R-function that computes the values of the
      derivatives in the ODE system (the \emph{model definition}) at time t.
  
      The \R-function  \code{func} must be defined as:
      \code{yprime = func(t, y, parms,...)}.  \code{t} is the current time point
      in the integration, \code{y} is the current estimate of the variables
      in the ode system, and \code{parms} is a vector or list of parameters.
      ... (optional) are any other arguments passed to the function.
  
      The return value of \code{func} should be a list, whose first element is a
      vector containing the derivatives of \code{y} with respect to
      \code{time}, and whose next elements are global values that are required at
      each point in \code{times}.}
  \item{parms }{any parameters used in \code{func}}
  \item{rtol }{relative error tolerance, either a scalar or an array as
      long as \code{y}. Only applicable to methods with variable time step, see details.}
  \item{atol }{absolute error tolerance, either a scalar or an array as
      long as \code{y}. Only applicable to methods with variable time step, see details.}
  \item{tcrit }{if not \code{NULL}, then \code{rk} cannot integrate past 
      \code{tcrit}. The solver routines may overshoot their targets
      (times points in the vector \code{times}), and interpolates values
      for the desired time points.  If there is a time beyond which
      integration should not proceed (perhaps because of a singularity),
      that should be provided in \code{tcrit}. }
  \item{verbose }{a logical value that, when TRUE, triggers more
      verbose output from the ODE solver.  
      %% Will output the settings of vectors *istate* and *rstate* - see details
      }
  \item{hmin }{an optional minimum value of the integration
      stepsize. In special situations this parameter may speed up computations with
      the cost of precision. Don't use \code{hmin} if you don't know why!}
  \item{hmax }{an optional maximum value of the integration stepsize. If not specified, 
       \code{hmax} is set to the largest difference in \code{times}, to 
       avoid that the simulation possibly ignores short-term events. If 0, no maximal size is specified}
  \item{hini }{initial step size to be attempted; if 0, the initial step size is determined by the solver}

  \item{method}{the integrator to use. This can either be a string constant naming one of the
    pre-defined methods or a call to function \code{\link{rkMethod}} specifying a user-defined method.
    The most common methods are the fixed-step methods "euler", "rk2", "rk4" 
    or the variable step methods "rk23bs", "rk34f", "rk45f" or "rk45dp7".}
  \item{maxsteps }{maximal number of steps during one call to the solver}
  \item{... }{additional arguments passed to \code{func} allowing this to be a generic function}

}
\details{
  The Runge-Kutta solvers are primarily provided for didactic reasons. 
  For most practical cases, solvers of the Livermore family 
  (\code{\link{lsoda}}, \code{\link{lsode}}, \code{\link{lsodes}}, 
  \code{\link{lsodar}}, \code{\link{vode}}, \code{\link{daspk}})
  are superior because of higher efficiency and faster implementation 
  (FORTRAN and C). 
  In addition to this, some of the Livermore solvers are also suitable for stiff 
  ODEs, differential algebraic equations (DAEs), or partial differential equations 
  (PDEs).

  Function \code{rk} is a generalized implementation that can be used to evaluate
  different solvers of the Runge-Kutta family. A pre-defined set of common
  method parameters is in function \code{\link{rkMethod}} which also 
  allows to supply user-defined Butcher tables. 
    
  The input parameters \code{rtol}, and \code{atol} determine the error
  control performed by the solver.  The solver will control the vector
  of estimated local errors in \bold{y}, according to an
  inequality of the form max-norm of ( \bold{e}/\bold{ewt} )
  \eqn{\leq}{<=} 1, where \bold{ewt} is a vector of positive error
  weights.  The values of \code{rtol} and \code{atol} should all be
  non-negative.
  The form of \bold{ewt} is:
  
  \deqn{\mathbf{rtol} \times \mathrm{abs}(\mathbf{y}) + \mathbf{atol}}{\bold{rtol} * abs(\bold{y}) + \bold{atol}}
  
  where multiplication of two vectors is element-by-element.
     
  \bold{Models} can be defined in \R as a user-supplied \bold{R-function}, 
  that must be called as: yprime = func(t, y, parms). 
  t is the current time point in the integration, 
  y is the current estimate of the variables in the ODE system. 
  The return value of \code{func} should be a list, whose first element is a 
  vector containing the derivatives of y with respect to time, 
  and whose second element contains output variables that are required at each 
  point in time. 
  \cr An example is given below:
  
    \code{model<-function(t,Y,parameters)}\cr
    \code{\{}\cr
    \code{with (as.list(parameters),\{}\cr 
       \code{dy1 = -k1*Y[1] + k2*Y[2]*Y[3]}\cr
       \code{dy3 = k3*Y[2]*Y[2]}\cr
       \code{dy2 = -dy1 - dy3}\cr
      \code{list(c(dy1,dy2,dy3))}\cr         
      \code{\})}\cr
  \code{\}}\cr
 
}
\value{
  A matrix with up to as many rows as elements in \code{times} and as
  many columns as elements in \code{y} plus the number of "global"
  values returned in the next elements of the return from \code{func},
  plus and additional column for the time value.  There will be a row
  for each element in \code{times} unless the solver
  returns with an unrecoverable error.  If \code{y} has a names
  attribute, it will be used to label the columns of the output value.

  The output will have the attributes \code{istate}, and \code{rstate}, 
  two vectors with several useful elements. 
  See details.

  The output will have the \bold{attributes *istate*, and *rstate*}, 
  two vectors with several useful elements, whose interpretation is
  compatible with \code{\link{lsoda}}:
  
  \itemize{
    \item{el 1:}{0 for normal return,  -2 means excess accuracy requested.}
      (Tolerances too small),
    \item{el 12:}{The number of steps taken for the problem so far.}
    \item{el 13:}{The number of function evaluations for the problem so far.}
    \item{el 15:}{The order of the method.}
   }
}
\references{
  Press, W. H.,  Teukolsky, S. A., Vetterling, W. T. and
  Flannery, B. P. (2007) Numerical Recipes in C. Cambridge
  University Press.
  
  Butcher, J. C. (1987) The numerical analysis of ordinary differential equations,
  Runge-Kutta and general linear methods, Wiley, Chichester and New York.
  
  Engeln-Muellges, G. and Reutter, F. (1996) Numerik Algorithmen:
  Entscheidungshilfe zur Auswahl und Nutzung. VDI Verlag, Duesseldorf.
}
\author{Thomas Petzoldt \email{thomas.petzoldt@tu-dresden.de}}


\seealso{\code{\link{rkMethod}}, \code{\link{ode}}, 
         \code{\link{lsoda}}, 
         \code{\link{lsode}}, 
         \code{\link{lsodes}}, 
         \code{\link{lsodar}}, 
         \code{\link{vode}}, 
         \code{\link{daspk}}}
\examples{
  #########################################
  ## Example: Lotka-volterra model
  #########################################
  
  ## Note: 
  ## parameters are a list, names accessible via "with" statement
  ## (see also ode and lsoda examples)

  lvmodel <- function(t, x, parms) {
    S <- x[1] # substrate
    P <- x[2] # producer
    K <- x[3] # consumer
    
    with(parms,{
      import <- approx(signal$times, signal$import, t)$y
      dS <- import - b * S * P + g * K
      dP <- c * S * P  - d * K * P
      dK <- e * P * K  - f * K
      res<-c(dS, dP, dK)
      list(res)
    })
  }
  
  ## vector of timesteps
  times  <- seq(0, 100, length=101)
  
  ## external signal with rectangle impulse
  signal <- as.data.frame(list(times = times,
                              import = rep(0,length(times))))
  
  signal$import[signal$times >= 10 & signal$times <=11] <- 0.2
  
  ## Parameters for steady state conditions
  parms  <- list(b=0.0, c=0.1, d=0.1, e=0.1, f=0.1, g=0.0)
  
  ## Start values for steady state
  y<-xstart <- c(S=1, P=1, K=1)
  
  ## Euler method
  out1  <- as.data.frame(rk(xstart, times, lvmodel, parms, 
                            hini = 0.1, method="euler"))
  
  ## classical Runge-Kutta 4th order
  out2 <- as.data.frame(rk(xstart, times, lvmodel, parms, 
                           hini = 1, method="rk4"))
  
  ## Dormand-Prince method of order 5(4)
  out3 <- as.data.frame(rk(xstart, times, lvmodel, parms, 
                           hmax=1, method = "rk45dp7"))
  
  mf <- par(mfrow=c(2,2))
  plot (out1$time, out1$S, type="l",    ylab="Substrate")
  lines(out2$time, out2$S, col="red",   lty="dotted", lwd=2)
  lines(out3$time, out3$S, col="green", lty="dotted")
  
  plot (out1$time, out1$P, type="l",    ylab="Producer")
  lines(out2$time, out2$P, col="red",   lty="dotted")
  lines(out3$time, out3$P, col="green", lty="dotted")
  
  plot (out1$time, out1$K, type="l",    ylab="Consumer")
  lines(out2$time, out2$K, col="red",   lty="dotted", lwd=2)
  lines(out3$time, out3$K, col="green", lty="dotted")
  
  plot (out1$P, out1$K, type="l", xlab="Producer",ylab="Consumer")
  lines(out2$P, out2$K, col="red",   lty="dotted", lwd=2)
  lines(out3$P, out3$K, col="green", lty="dotted")
  legend("center",legend=c("euler","rk4","rk45dp7"),lty=c(1,3,3),
         lwd=c(1,2,1),col=c("black","red","green"))
  par(mfrow=mf)
           
}

\keyword{ math }