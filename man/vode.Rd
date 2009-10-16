\name{vode}
\alias{vode}

\title{Solver for Ordinary Differential Equations (ODE)}

\description{
  Solves the initial value problem for stiff or nonstiff systems of 
  ordinary differential equations (ODE) in the form:

  \deqn{dy/dt = f(t,y)}
    
  The \R function \code{vode} provides an interface to the FORTRAN ODE
  solver of the same name, written by Peter N. Brown, Alan C. Hindmarsh
  and George D. Byrne.
    
  The system of ODE's is written as an \R function or be defined in
  compiled code that has been dynamically loaded.
  
  In contrast to \code{\link{lsoda}}, the user has to specify whether or
  not the problem is stiff and choose the appropriate solution method.
    
  \code{vode} is very similar to \code{\link{lsode}}, but uses a
  variable-coefficient method rather than the fixed-step-interpolate
  methods in \code{\link{lsode}}.  In addition, in vode it is possible
  to choose whether or not a copy of the Jacobian is saved for reuse in
  the corrector iteration algorithm; In \code{lsode}, a copy is not
  kept.
}

\usage{vode(y, times, func, parms, rtol = 1e-6, atol = 1e-8,  
  jacfunc = NULL, jactype = "fullint", mf = NULL, verbose = FALSE,   
  tcrit = NULL, hmin = 0, hmax = NULL, hini = 0, ynames = TRUE,
  maxord = NULL, bandup = NULL, banddown = NULL, maxsteps = 5000,
  dllname = NULL, initfunc = dllname, initpar = parms, rpar = NULL,
  ipar = NULL, nout = 0, outnames = NULL, forcings=NULL,
  initforc = NULL, fcontrol=NULL, ...)
}
\arguments{
  \item{y }{the initial (state) values for the ODE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time; if only one step is
    to be taken; set \code{times = NULL}.
  }
  \item{func }{either an \R-function that computes the values of the
    derivatives in the ODE system (the \emph{model definition}) at time
    \code{t}, or a character string giving the name of a compiled function in a
    dynamically loaded shared library.
  
    If \code{func} is an \R-function, it must be defined as:
     \code{func <- function(t, y, parms,...)}.  \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system.  If the initial values \code{y} has a
    \code{names} attribute, the names will be available inside \code{func}.
    \code{parms} is a vector or list of parameters; ... (optional) are
    any other arguments passed to the function.
  
    The return value of \code{func} should be a list, whose first
    element is a vector containing the derivatives of \code{y} with
    respect to \code{time}, and whose next elements are global values
    that are required at each point in \code{times}.  The derivatives
    should be specified in the same order as the state variables \code{y}.
    
    If \code{func} is
    a string, then \code{dllname} must give the name of the shared
    library (without extension) which must be loaded before
    \code{vode()} is called. See package vignette \code{"compiledCode"}
    for more details.
  }
  \item{parms }{vector or list of parameters used in \code{func} or
    \code{jacfunc}.
  }
  \item{rtol }{relative error tolerance, either a scalar or an array as
    long as \code{y}.  See details.
  }
  \item{atol }{absolute error tolerance, either a scalar or an array as
    long as \code{y}.  See details.
  }
  \item{jacfunc }{if not \code{NULL}, an \R function that computes the
    Jacobian of the system of differential equations dydot(i)/dy(j), or
    a string giving the name of a function or subroutine in
    \file{dllname} that computes the Jacobian (see vignette
    \code{"compiledCode"} for more about this option).

    In some circumstances, supplying
    \code{jacfunc} can speed up the computations, if the system is
    stiff.  The \R calling sequence for \code{jacfunc} is identical to
    that of \code{func}.

    If the Jacobian is a full matrix, \code{jacfunc} should return a
    matrix dydot/dy, where the ith row contains the derivative of
    \eqn{dy_i/dt} with respect to \eqn{y_j}, or a vector containing the
    matrix elements by columns (the way \R and FORTRAN store matrices).

    If the Jacobian is banded, \code{jacfunc} should return a matrix
    containing only the nonzero bands of the Jacobian, rotated
    row-wise. See first example of \code{lsode}.
  }
  \item{jactype }{the structure of the Jacobian, one of
    \code{"fullint"}, \code{"fullusr"}, \code{"bandusr"} or
    \code{"bandint"} - either full or banded and estimated internally or
    by user; overruled if \code{mf} is not \code{NULL}.
  }
  \item{mf }{the "method flag" passed to function vode - overrules
    \code{jactype} - provides more options than \code{jactype} - see
    details.
  }
  \item{verbose }{if TRUE: full output to the screen, e.g. will
    print the \code{diagnostiscs} of the integration - see details.
  }
  \item{tcrit }{if not \code{NULL}, then \code{vode} cannot integrate
    past \code{tcrit}. The FORTRAN routine \code{dvode} overshoots its
    targets (times points in the vector \code{times}), and interpolates
    values for the desired time points.  If there is a time beyond which
    integration should not proceed (perhaps because of a singularity),
    that should be provided in \code{tcrit}.
  }
  \item{hmin }{an optional minimum value of the integration stepsize. In
    special situations this parameter may speed up computations with the
    cost of precision. Don't use hmin if you don't know why!
  }
  \item{hmax }{an optional maximum value of the integration stepsize. If
    not specified, hmax is set to the largest difference in
    \code{times}, to avoid that the simulation possibly ignores
    short-term events. If 0, no maximal size is specified.
  }
  \item{hini }{initial step size to be attempted; if 0, the initial step
    size is determined by the solver.
  }
  \item{ynames }{logical; if \code{FALSE}: names of state variables are not
    passed to function \code{func} ; this may speed up the simulation
    especially for multi-D models.
  }
  \item{maxord }{the maximum order to be allowed. \code{NULL} uses the default,
    i.e. order 12 if implicit Adams method (meth = 1), order 5 if BDF
    method (meth = 2). Reduce maxord to save storage space.
  }
  \item{bandup }{number of non-zero bands above the diagonal, in case
    the Jacobian is banded.
  }
  \item{banddown }{number of non-zero bands below the diagonal, in case
    the Jacobian is banded.
  }
  \item{maxsteps }{maximal number of steps per output interval taken by the
    solver.
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions refered to in \code{func} and \code{jacfunc}.
    See package vignette \code{"compiledCode"}.
  }
  \item{initfunc }{if not \code{NULL}, the name of the initialisation function
    (which initialises values of parameters), as provided in
    \file{dllname}. See package vignette \code{"compiledCode"}.
  }
  \item{initpar }{only when \file{dllname} is specified and an
    initialisation function \code{initfunc} is in the dll: the
    parameters passed to the initialiser, to initialise the common
    blocks (FORTRAN) or global variables (C, C++).
  }
  \item{rpar }{only when \file{dllname} is specified: a vector with
    double precision values passed to the dll-functions whose names are
    specified by \code{func} and \code{jacfunc}.
  }
  \item{ipar }{only when \file{dllname} is specified: a vector with
    integer values passed to the dll-functions whose names are specified
    by \code{func} and \code{jacfunc}.
  }
  \item{nout }{only used if \code{dllname} is specified and the model is
    defined in compiled code: the number of output variables calculated
    in the compiled function \code{func}, present in the shared
    library. Note: it is not automatically checked whether this is
    indeed the number of output variables calculed in the dll - you have
    to perform this check in the code - See package vignette \code{"compiledCode"}.
  }
  \item{outnames }{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{func}, present in the shared library.
    These names will be used to label the output matrix.
  }
  \item{forcings }{only used if \file{dllname} is specified: a list with
    the forcing function data sets, each present as a two-columned matrix,
    with (time,value); interpolation outside the interval
    [min(\code{times}), max(\code{times})] is done by taking the value at
    the closest data extreme.

    See \link{forcings} or package vignette \code{"compiledCode"}.
  }
  \item{initforc }{if not \code{NULL}, the name of the forcing function
    initialisation function, as provided in
    \file{dllname}. It MUST be present if \code{forcings} has been given a
    value.
    See \link{forcings} or package vignette \code{"compiledCode"}.
  }
  \item{fcontrol }{A list of control parameters for the forcing functions.
    \link{forcings} or package vignette \code{"compiledCode"}
  }
  \item{... }{additional arguments passed to \code{func} and
    \code{jacfunc} allowing this to be a generic function.
  }
}
\value{
  A matrix of class \code{deSolve} with up to as many rows as elements
  in \code{times} and as many columns as elements in \code{y} plus the number of "global"
  values returned in the next elements of the return from \code{func},
  plus and additional column for the time value.  There will be a row
  for each element in \code{times} unless the FORTRAN routine `lsoda'
  returns with an unrecoverable error. If \code{y} has a names
  attribute, it will be used to label the columns of the output value.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
## =======================================================================
## ex. 1
## The famous Lorenz equations: chaos in the earth's atmosphere
## Lorenz 1963. J. Atmos. Sci. 20, 130-141.
## =======================================================================

chaos <- function(t, state, parameters)
{
  with(as.list(c(state)),{

    dx     <- -8/3*x+y*z
    dy     <- -10*(y-z)
    dz     <- -x*y+28*y-z

    list(c(dx, dy, dz))
  })
}

state <- c(x = 1, y = 1, z = 1)
times <- seq(0, 100, 0.01)

out   <- vode(state, times, chaos, 0)

plot(out, type="l")   # all versus time
plot(out[,"x"], out[,"y"], type = "l", main = "Lorenz butterfly",
  xlab="x",ylab="y")


## =======================================================================
## ex. 2
## SCOC model, in FORTRAN  - to see the FORTRAN code:
## browseURL(paste(system.file(package="deSolve"),"/dynload/scoc.f",sep=""))
## example from Soetaert and Herman, 2009, chapter 3. (simplified)
## =======================================================================

# Forcing function data
Flux <- matrix(ncol=2,byrow=TRUE,data=c(
  1, 0.654, 11, 0.167,   21, 0.060, 41, 0.070, 73,0.277, 83,0.186,
  93,0.140,103, 0.255,  113, 0.231,123, 0.309,133,1.127,143,1.923,
  153,1.091,163,1.001,  173, 1.691,183, 1.404,194,1.226,204,0.767,
  214, 0.893,224,0.737, 234,0.772,244, 0.726,254,0.624,264,0.439,
  274,0.168,284 ,0.280, 294,0.202,304, 0.193,315,0.286,325,0.599,
  335, 1.889,345, 0.996,355,0.681,365,1.135))

parms <- c(k=0.01)

meanDepo <- mean(approx(Flux[,1],Flux[,2], xout=seq(1,365,by=1))$y)
Yini <- meanDepo/parms

times <- 1:365
out <- as.data.frame( vode(Yini, times, func = "scocder",
    parms = parms, dllname = "deSolve",
    initforc="scocforc",  forcings=Flux,
    initfunc = "scocpar", nout = 2,
    outnames = c("Mineralisation","Depo")))

plot(out$time,out$Depo,type="l",col="red")
lines(out$time,out$Mineralisation,col="blue")

# Constant interpolation of forcing function - left side of interval
fcontrol<-list(method="constant")

out2 <- as.data.frame( vode(Yini, times, func = "scocder",
    parms = parms, dllname = "deSolve",
    initforc="scocforc",  forcings=Flux, fcontrol=fcontrol,
    initfunc = "scocpar", nout = 2,
    outnames = c("Mineralisation","Depo")))

plot(out2$time,out2$Depo,type="l",col="red")
lines(out2$time,out2$Mineralisation,col="blue")

# Constant interpolation of forcing function - middle of interval
fcontrol<-list(method="constant",f=0.5)

out3 <- as.data.frame( vode(Yini, times, func = "scocder",
    parms = parms, dllname = "deSolve",
    initforc="scocforc",  forcings=Flux, fcontrol=fcontrol,
    initfunc = "scocpar", nout = 2,
    outnames = c("Mineralisation","Depo")))

lines(out3$time,out3$Depo,type="l",col="orange",lty=2)

}
\references{
  
  P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, 1989.  VODE: A Variable
  Coefficient ODE Solver, SIAM J. Sci. Stat. Comput., 10, pp. 1038-1051.
  \cr Also, LLNL Report UCRL-98412, June 1988.

  G. D. Byrne and A. C. Hindmarsh, 1975. A Polyalgorithm for the
  Numerical Solution of Ordinary Differential Equations.  ACM
  Trans. Math. Software, 1, pp. 71-96.

  A. C. Hindmarsh and G. D. Byrne, 1977. EPISODE: An Effective Package
  for the Integration of Systems of Ordinary Differential
  Equations. LLNL Report UCID-30112, Rev. 1.

  G. D. Byrne and A. C. Hindmarsh, 1976. EPISODEB: An Experimental
  Package for the Integration of Systems of Ordinary Differential
  Equations with Banded Jacobians. LLNL Report UCID-30132, April 1976.

  A. C. Hindmarsh, 1983. ODEPACK, a Systematized Collection of ODE
  Solvers. in Scientific Computing, R. S. Stepleman et al., eds.,
  North-Holland, Amsterdam, pp. 55-64.
  
  K. R. Jackson and R. Sacks-Davis, 1980. An Alternative Implementation
  of Variable Step-Size Multistep Formulas for Stiff ODEs.  ACM
  Trans. Math. Software, 6, pp. 295-318.
  
  Netlib: \url{http://www.netlib.org}
}
\details{

  Before using the integrator \code{vode}, the user has to decide
  whether or not the problem is stiff.
  
  If the problem is nonstiff, use method flag \code{mf} = 10, which
  selects a nonstiff (Adams) method, no Jacobian used.
  
  If the problem is stiff, there are four standard choices which can be
  specified with \code{jactype} or \code{mf}.
  
  The options for \bold{jactype} are
  \describe{
    \item{jac = "fullint":}{a full Jacobian, calculated internally by
      vode, corresponds to \code{mf} = 22,
    }
    \item{jac = "fullusr":}{a full Jacobian, specified by user function
      \code{jacfunc}, corresponds to \code{mf} = 21,
    }
    \item{jac = "bandusr":}{a banded Jacobian, specified by user
      function \code{jacfunc}; the size of the bands specified by
      \code{bandup} and \code{banddown}, corresponds to \code{mf} = 24,
    }
    \item{jac = "bandint":}{a banded Jacobian, calculated by vode; the
      size of the bands specified by \code{bandup} and \code{banddown},
      corresponds to \code{mf} = 25.
    }
  }
  
  More options are available when specifying \bold{mf} directly.
  
  The legal values of \code{mf} are 10, 11, 12, 13, 14, 15, 20, 21, 22,
  23, 24, 25, -11, -12, -14, -15, -21, -22, -24, -25.
  
  \code{mf} is a signed two-digit integer, \code{mf = JSV*(10*METH +
  MITER)}, where
  \describe{
    \item{JSV = SIGN(mf)}{indicates the Jacobian-saving strategy: JSV =
      1 means a copy of the Jacobian is saved for reuse in the corrector
      iteration algorithm. JSV = -1 means a copy of the Jacobian is not
      saved.
    }
    \item{METH}{indicates the basic linear multistep method: METH = 1
      means the implicit Adams method. METH = 2 means the method based
      on backward differentiation formulas (BDF-s).
    }
    \item{MITER}{indicates the corrector iteration method: MITER = 0
      means functional iteration (no Jacobian matrix is involved).
      
      MITER = 1 means chord iteration with a user-supplied full (NEQ by
      NEQ) Jacobian.
     
      MITER = 2 means chord iteration with an internally generated
      (difference quotient) full Jacobian (using NEQ extra calls to
      \code{func} per df/dy value).
      
      MITER = 3 means chord iteration with an internally generated
      diagonal Jacobian approximation (using 1 extra call to \code{func}
      per df/dy evaluation).
      
      MITER = 4 means chord iteration with a user-supplied banded
      Jacobian.
      
      MITER = 5 means chord iteration with an internally generated
      banded Jacobian (using ML+MU+1 extra calls to \code{func} per
      df/dy evaluation).
    }
  }
  If MITER = 1 or 4, the user must supply a subroutine \code{jacfunc}.
  
  The example for integrator \code{\link{lsode}} demonstrates how to
  specify both a banded and full Jacobian.
  
  The input parameters \code{rtol}, and \code{atol} determine the
  \bold{error control} performed by the solver.  If the request for
  precision exceeds the capabilities of the machine, vode will return an
  error code. See \code{\link{lsoda}} for details.
  

  The diagnostics of the integration can be printed to screen
  by calling \code{\link{diagnostics}}. If \code{verbose} = \code{TRUE},
  the diagnostics will written to the screen at the end of the integration.

  See vignette("deSolve") for an explanation of each element in the vectors
  containing the diagnostic properties and how to directly access them.

  \bold{Models} may be defined in compiled C or FORTRAN code, as well as
  in an R-function. See package vignette \code{"compiledCode"} for details.

  More information about models defined in compiled code is in the package
  vignette ("compiledCode"); information about linking forcing functions
  to compiled code is in \link{forcings}.

  Examples in both C and FORTRAN are in the \file{dynload} subdirectory
  of the \code{deSolve} package directory.


}
\seealso{
  \itemize{
    \item \code{\link{rk}}, \item \code{\link{rk4}} and \code{\link{euler}} for
      Runge-Kutta integrators.
    \item \code{\link{lsoda}}, \code{\link{lsode}},
      \code{\link{lsodes}}, \code{\link{lsodar}},
      \code{\link{daspk}} for other solvers of the Livermore family,
    \item \code{\link{ode}} for a general interface to most of the ODE solvers,
    \item \code{\link{ode.band}} for solving models with a banded
      Jacobian,
    \item \code{\link{ode.1D}} for integrating 1-D models,
    \item \code{\link{ode.2D}} for integrating 2-D models,
    \item \code{\link{ode.3D}} for integrating 3-D models,
  }

  \code{\link{diagnostics}} to print diagnostic messages.
}
\keyword{math}

