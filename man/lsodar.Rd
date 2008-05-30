\name{lsodar}
\alias{lsodar}
\title{General solver for ordinary differential equations (ODE), switching automatically between stiff and non-stiff methods and with root finding}
\description{Solving initial value problems for stiff or
    non-stiff systems of first-order ordinary differential equations
    (ODEs) and including root-finding. 
    \cr The \R function \code{lsodar} provides an interface to the
    Fortran ODE solver of the same name, written by Alan
    C. Hindmarsh and Linda R. Petzold.  
    
    The system of ODE's is written as an \R function or be defined in
    compiled code that has been dynamically loaded. - see description of \code{\link{lsoda}} for details. 
  
    \code{lsodar} differs from \code{lsode} in two respects.
    \item It switches automatically between stiff and nonstiff methods (similar as lsoda).
    \item It finds the root of at least one of a set of constraint
  functions g(i) of the independent and dependent variables.}
\usage{lsodar(y, times, func, parms, rtol=1e-6, atol=1e-6, 
  jacfunc=NULL, jactype="fullint", rootfunc=NULL, verbose=FALSE,   
  nroot=0, tcrit=NULL, hmin=0, hmax=NULL, hini=0, ynames=TRUE, 
  maxordn=12, maxords = 5, bandup=NULL, banddown=NULL, maxsteps=5000, 
  dllname=NULL, initfunc=dllname, initpar=parms, rpar=NULL, 
  ipar=NULL, nout=0, outnames=NULL, ...)}
\arguments{
  \item{y }{the initial (state) values for the ODE system. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{times }{times at which explicit estimates for \code{y} are desired.  The first value in \code{times} must be the initial time.}
  \item{func }{either an \R-function that computes the values of the
      derivatives in the ODE system (the \emph{model definition}) at time
      t, or a character string
      giving the name of a compiled function in a dynamically loaded
      shared library.
  
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
  
      If \code{func} is a string, then \code{dllname} must give the name
      of the shared library (without extension) which must be loaded
      before \code{lsodar()} is called. See package vignette for more details}
  \item{parms }{vector or list of parameters used in \code{func} or \code{jacfunc}.}
  \item{rtol }{relative error tolerance, either a scalar or an array as
      long as \code{y}.  See details. }
  \item{atol }{absolute error tolerance, either a scalar or an array as
      long as \code{y}.  See details.}
  \item{jacfunc }{if not \code{NULL}, an \R function, that computes
      the jacobian of the system of differential equations
      dydot(i)/dy(j), or a string giving the name of a function or
      subroutine in \file{dllname} that computes the jacobian (see Details
      below for more about this option).  In some circumstances, supplying
      \code{jacfunc} can speed up
      the computations, if the system is stiff.  The \R calling sequence for
      \code{jacfunc} is identical to that of \code{func}.  
      \cr If the jacobian is a full matrix, \code{jacfunc} should return a matrix dydot/dy, where the ith
      row contains the derivative of \eqn{dy_i/dt} with respect to \eqn{y_j},
      or a vector containing the matrix elements by columns (the way \R and Fortran store matrices).
      \cr If the jacobian is banded, \code{jacfunc} should return a matrix containing only the 
      nonzero bands of the jacobian, rotated row-wise. See first example of lsode.}
  \item{jactype }{the structure of the jacobian, one of "fullint", "fullusr", "bandusr" or "bandint" - either full or banded and estimated internally or by user}
  \item{rootfunc }{if not \code{NULL}, an \R function that computes
      the function whose root has to be estimated or a string giving the name of a function or
      subroutine in \file{dllname} that computes the root function.  The \R calling sequence for
      \code{rootfunc} is identical to that of \code{func}.  \code{rootfunc} should
      return a vector with the function values whose root is sought}
  \item{verbose }{a logical value that, when TRUE, triggers more
      verbose output from the ODE solver.  Will output the settings of vectors *istate* and *rstate* - see details}
  \item{nroot }{only used if \file{dllname}  is specified: the number of constraint functions whose roots are desired during the integration; if \code{rootfunc} is an R-function, the solver estimates the number of roots}
  \item{tcrit }{if not \code{NULL}, then \code{lsodar} cannot integrate past \code{tcrit}. The Fortran routine \code{lsodar} overshoots its targets
      (times points in the vector \code{times}), and interpolates values
      for the desired time points.  If there is a time beyond which
      integration should not proceed (perhaps because of a singularity),
      that should be provided in \code{tcrit}.}
  \item{hmin }{an optional minimum value of the integration
      stepsize. In special situations this parameter may speed up computations with
      the cost of precision. Don't use \code{hmin} if you don't know why!}
  \item{hmax }{an optional maximum value of the integration stepsize. If not specified, \code{hmax} is set to the largest difference in \code{times}, to avoid that the simulation possibly ignores short-term events. If 0, no maximal size is specified}
  \item{hini }{initial step size to be attempted; if 0, the initial step size is determined by the solver}
  \item{ynames }{if FALSE: names of state variables are not passed to function \code{func} ; this may speed up the simulation especially for multi-D models}
  \item{maxordn }{the maximum order to be allowed in case the method is non-stiff. Should be <=12. Reduce maxord to save storage space}
  \item{maxords }{the maximum order to be allowed in case the method is stiff. Should be <=5. Reduce maxord to save storage space}
  \item{bandup }{number of non-zero bands above the diagonal, in case the Jacobian is banded}
  \item{banddown }{number of non-zero bands below the diagonal, in case the Jacobian is banded}
  \item{maxsteps }{maximal number of steps during one call to the solver}
  \item{dllname }{a string giving the name of the shared library (without
      extension) that contains all the compiled function or subroutine
      definitions refered to in \code{func} and \code{jacfunc}. See package vignette}
  \item{initfunc }{if not NULL, the name of the initialisation function (which initialises values of parameters), as provided in \file{dllname}. See package vignette. }
  \item{initpar }{only when \file{dllname} is specified and an initialisation function \code{initfunc} is in the dll: the parameters passed to the initialiser, to initialise the common blocks (fortran) or global variables (C, C++)}
  \item{rpar }{only when \file{dllname} is specified: a vector with double precision values passed to the dll-functions whose names are specified by \code{func} and \code{jacfunc}}
  \item{ipar }{only when \file{dllname} is specified: a vector with integer values passed to the dll-functions whose names are specified by \code{func} and \code{jacfunc}}
  \item{nout }{only used if \code{dllname} is specified and the model is defined in compiled code: the number of output variables calculated in the compiled function \code{func}, present in the shared library. Note:
      it is not automatically checked whether this is indeed the number of output variables calculed in the dll - you have to perform this check in the code - See package vignette}
  \item{outnames }{only used if \file{dllname} is specified and \code{nout} > 0: the names of output variables calculated in the compiled function \code{func}, present in the shared library}
  \item{... }{additional arguments passed to \code{func} and \code{jacfunc} allowing this to be a generic function}
}
\value{
  A matrix with up to as many rows as elements in \code{times} and as
    many columns as elements in \code{y} plus the number of "global"
    values returned in the next elements of the return from \code{func},
    plus and additional column for the time value.  There will be a row
    for each element in \code{times} unless the Fortran routine `lsodar'
    returns with an unrecoverable error or has found a root, in which case the last row will contain the function value at the root.  
    \cr If \code{y} has a names attribute, it will be used to label the columns of the output value.
      
      The output will have the attributes \code{istate}, and \code{rstate}, two vectors with several useful elements. 
      See details.
      The first element of istate returns the conditions under which the last call to lsoda returned. Normal is \code{istate[1] = 2}.
      If \code{verbose} = TRUE, the settings of istate and rstate will be written to the screen
      \cr if a root has been found, the output will also have the attribute \code{iroot}, an integer indicating which root has been found.
  
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
#########################################
### example 1: from lsodar source code
#########################################
  
  Fun <- function (t,y,parms)
  {
   ydot <- vector(len=3)
   ydot[1] <- -.04*y[1] + 1.e4*y[2]*y[3]
   ydot[3] <- 3.e7*y[2]*y[2]
   ydot[2] <- -ydot[1]-ydot[3]
  
   return(list(ydot,ytot = sum(y)))
  }
  
  rootFun <- function (t,y,parms)
  {
   yroot <- vector(len=2)
   yroot[1] <- y[1] - 1.e-4
   yroot[2] <- y[3] - 1.e-2
   return(yroot)
  }
  
  y     <- c(1,0,0)
  times <- c(0,0.4*10^(0:8))
  Out   <- NULL
  ny    <- length(y)
  
  out   <- lsodar(y=y,times=times,fun=Fun,rootfun=rootFun,
         rtol=1e-4,atol=c(1e-6,1e-10,1e-6), parms=NULL)
  print(paste("root is found for eqn",which(attributes(out)$iroot==1)))
  print(out[nrow(out),])
  
#########################################
### example 2:
### using lsodar to estimate steady-state conditions
#########################################
  
  # Bacteria (Bac) are growing on a substrate (Sub)
  model <- function(t,state,pars)
  {
  with (as.list(c(state,pars)), {
  #       substrate uptake             death  respiration
  dBact = gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact
  dSub  =-gmax    *Sub/(Sub+ks)*Bact + dB*Bact          +input
  
  return(list(c(dBact,dSub)))
                                })
  }
  
  # root is the condition where sum of |rates of change|
  # is very small
  
  rootfun <- function (t,state,pars)
  {
  dstate <- unlist(model(t,state,pars)) #rate of change vector 
  return(sum(abs(dstate))-1e-10)
  }
  
  pars <- list(Bini=0.1,Sini=100,gmax =0.5,eff = 0.5,
                ks =0.5, rB =0.01, dB =0.01, input=0.1)
  
  tout    <- c(0,1e10)
  state   <- c(Bact=pars$Bini,Sub =pars$Sini)
  out     <- lsodar(state,tout,model,pars,rootfun=rootfun)
  print(out)}
\references{\item Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
    Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
    North-Holland, Amsterdam, 1983, pp. 55-64.
  \item Linda R. Petzold, Automatic Selection of Methods for Solving
      Stiff and Nonstiff Systems of Ordinary Differential Equations,
      Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148.
  \item Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined
      Output Points for Solutions of ODEs, Sandia Report SAND80-0180,
      February 1980.
    
    Netlib: \url{http://www.netlib.org}}
\details{The work is done by the Fortran subroutine \code{lsodar},
    whose documentation should be consulted for details (it is included as
    comments in the source file \file{src/opkdmain.f}).  The implementation is based on the
    November, 2003 version of lsodar, from Netlib.  
    
  \code{lsodar} switches automatically between stiff and nonstiff methods (similar as \code{lsoda}).
  This means that the user does not have to determine whether the
  problem is stiff or not, and the solver will automatically choose the
  appropriate method.  It always starts with the nonstiff method.  
  
  It finds the root of at least one of a set of constraint
  functions g(i) of the independent and dependent variables.
  It then returns the solution at the root if that occurs
  sooner than the specified stop condition, and otherwise returns
  the solution according the specified stop condition.
  
  The form of the \bold{jacobian} can be specified by \code{jactype} which can take the following values.
  
      \item jactype = "fullint" : a full jacobian, calculated internally by lsodar, the default
      \item jactype = "fullusr" : a full jacobian, specified by user function \code{jacfunc}
      \item jactype = "bandusr" : a banded jacobian, specified by user function \code{jacfunc}; the size of the bands specified by \code{bandup} and \code{banddown}
      \item jactype = "bandint" : a banded jacobian, calculated by lsodar; the size of the bands specified by \code{bandup} and \code{banddown}
  
  if \code{jactype}= "fullusr" or "bandusr" then the user must supply a subroutine \code{jacfunc}.
  
  The input parameters \code{rtol}, and \code{atol} determine the \bold{error
  control} performed by the solver. See \code{\link{lsoda}} for details.
  
  \bold{Models} may be defined in compiled C or Fortran code, as well as in an R-function.  See package vignette for details.
  Examples in Fortran are in the \file{dynload} subdirectory of
  the \code{deSolve} package directory.
  
  
  The output will have the \bold{attributes *istate*, *rstate*}, and if a root was found \bold{iroot}, three vectors with several useful elements. 
   \cr if \code{verbose} = TRUE, the settings of istate and rstate will be written to the screen. 
  
      the following elements of \bold{istate} are meaningful:    
  
      \item el 1 : returns the conditions under which the last call to lsodar returned.
      2  if lsodar was successful, 3 if lsodar was succesful and one or more roots were found - see \code{iroot}.
      -1 if excess work done, -2 means excess accuracy requested. (Tolerances too small),
  -3 means illegal input detected. (See printed message.), -4 means repeated error test failures. (Check all input),
  -5 means repeated convergence failures. (Perhaps bad Jacobian supplied or wrong choice of MF or tolerances.),
  -6 means error weight became zero during problem. (Solution component i vanished, and atol or atol(i) = 0.)
      \item el 12 : The number of steps taken for the problem so far.
      \item el 13 : The number of function evaluations for the problem so far.,
      \item el 14 : The number of Jacobian evaluations and LU decompositions so far.,
      \item el 15 : The method order last used (successfully).,
      \item el 16 : The order to be attempted on the next step.,
      \item el 17 : if el 1 =-4,-5: the largest component in the error vector,
      \item el 18 : The length of rwork actually required.,
      \item el 19 : The length of IUSER actually required.,
      \item el 20 : The method indicator for the last succesful step, 1=adams (nonstiff), 2= bdf (stiff),
      \item el 21 : The current method indicator to be attempted on th next step, 1=adams (nonstiff), 2= bdf (stiff),
  
      \bold{rstate} contains the following:
      \item 1: The step size in t last used (successfully).
      \item 2: The step size to be attempted on the next step.
      \item 3: A tolerance scale factor, greater than 1.0, computed when a request for too much accuracy was detected.
      \item 4: the value of t at the time of the last method switch, if any.
  
      \bold{iroot} is a vector, its length equal to the number of constraint functions; 
      it will have a value of 1 for the constraint function whose root that has been found and 0 otherwise.}
\seealso{\code{\code{\link{ode}}, \link{lsoda}}, \code{\link{lsode}}, \code{\link{lsodes}}, 
    \code{\link{vode}}, \code{\link{daspk}}, \code{\link{rk}}.}
\keyword{math}

