\name{daspk}
\alias{daspk}                                                   
\title{Solver for Differential Algebraic Equations (DAE)}
\description{
  Solves either:
  \itemize{
    \item a system of ordinary differential equations (ODE) of the form
      \deqn{y' = f(t, y, ...)} or
    \item a system of differential algebraic equations (DAE) of the form
      \deqn{F(t,y,y') = 0} or
    \item a system of linearly implicit DAES in the 
       form \deqn{M y' = f(t, y)}    
  }
  using a combination of backward differentiation formula (BDF) and a
  direct linear system solution method (dense or banded).
  
  The \R function \code{daspk} provides an interface to the FORTRAN DAE
  solver of the same name, written by Linda R. Petzold, Peter N. Brown,
  Alan C. Hindmarsh and Clement W. Ulrich.
    
  The system of DE's is written as an \R function (which may, of course,
  use \code{\link{.C}}, \code{.Fortran}, \code{\link{.Call}}, etc., to
  call foreign code) or be defined in compiled code that has been
  dynamically loaded.
}
\usage{
daspk(y, times, func = NULL, parms, nind = c(length(y), 0, 0), 
  dy = NULL, res = NULL, nalg = 0, 
  rtol = 1e-6, atol = 1e-6, jacfunc = NULL,
  jacres = NULL, jactype = "fullint", mass = NULL, estini = NULL,
  verbose = FALSE, tcrit = NULL, hmin = 0, hmax = NULL,
  hini = 0, ynames = TRUE, maxord = 5, bandup = NULL,
  banddown = NULL, maxsteps = 5000, dllname = NULL,
  initfunc = dllname, initpar = parms, rpar = NULL,
  ipar = NULL, nout = 0, outnames = NULL,
  forcings=NULL, initforc = NULL, fcontrol=NULL,
  events = NULL, lags = NULL, ...)
}

\arguments{
  \item{y }{the initial (state) values for the DE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time; if only one step is
    to be taken; set \code{times} = \code{NULL}.
  }
  \item{func }{to be used if the model is an ODE, or a DAE written in linearly
    implicit form (M y' = f(t, y)). 
    \code{func} should be an \R-function that computes the
    values of the derivatives in the ODE system (the \emph{model
    definition}) at time t.
  
    \code{func} must be defined as:  \code{func <- function(t, y, parms,...)}.
    \cr \code{t} is the current time point in the
    integration, \code{y} is the current estimate of the variables in
    the ODE system.  If the initial values \code{y} has a \code{names}
    attribute, the names will be available inside \code{func}, unless
    \code{ynames} is FALSE.  \code{parms} is a vector or list of
    parameters.  \code{...} (optional) are any other arguments passed to
    the function.
    
    The return value of \code{func} should be a list,
    whose first element is a vector containing the derivatives of
    \code{y} with respect to \code{time}, and whose next elements are
    global values that are required at each point in \code{times}.
    The derivatives should be specified in the same order as the specification
    of the state variables \code{y}.

    Note that it is not possible to define \code{func} as a compiled
    function in a dynamically loaded shared library. Use \code{res}
    instead.
  }
  \item{parms }{vector or list of parameters used in \code{func},
    \code{jacfunc}, or \code{res}
  }
  \item{nind }{if a DAE system: a three-valued vector with the number of
    variables of index 1, 2, 3 respectively.
    The equations must be defined such that the index 1 variables precede
    the index 2 variables which in turn precede the index 3 variables.
    The sum of the variables of different index should equal N,
    the total number of variables.
	Note that this has been added for consistency with \link{radau}. 
	If used, then the variables are weighed differently than in the 
	original daspk code, i.e. index 2 variables are scaled with 1/h,
	index 3 variables are scaled with 1/h^2. In some cases this allows
	daspk to solve index 2 or index 3 problems.
  }
  \item{dy }{the initial derivatives of the state variables of the DE
    system.  Ignored if an ODE.
  }
  \item{res }{if a DAE system: either an \R-function that computes the
    residual function \eqn{F(t,y,y')} of the DAE system (the model
    defininition) at time \code{t}, or a character string giving the
    name of a compiled function in a dynamically loaded shared library.
      
    If \code{res} is a user-supplied \R-function, it must be defined as:
    \code{res <- function(t, y, dy, parms, ...)}.
    
    Here \code{t} is the current time point in the integration, \code{y}
    is the current estimate of the variables in the ODE system,
    \code{dy} are the corresponding derivatives.  If the initial
    \code{y} or \code{dy} have a \code{names} attribute, the names will be
    available inside \code{res}, unless \code{ynames} is \code{FALSE}.
    \code{parms} is a vector of parameters.
    
    The return value of \code{res} should be a list, whose first element
    is a vector containing the residuals of the DAE system,
    i.e. \eqn{\delta = F(t,y,y')}{delta = F(t,y,y')}, and whose next elements contain output
    variables that are required at each point in \code{times}.
  
    If \code{res} is a string, then \code{dllname} must give the name of
    the shared library (without extension) which must be loaded before
    \code{daspk()} is called (see package vignette \code{"compiledCode"}
    for more information).
  }
  \item{nalg }{if a DAE system: the number of algebraic equations
    (equations not involving derivatives). Algebraic equations should
    always be the last, i.e.  preceeded by the differential equations.
    
    Only used if \code{estini} = 1.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector,
    one value for each y,
  }
  \item{atol }{absolute error tolerance, either a scalar or a vector,
    one value for each y.
  }
  \item{jacfunc }{if not \code{NULL}, an \R function that computes the
    Jacobian of the system of differential equations. Only used in case
    the system is an ODE (\eqn{y' = f(t, y)}), specified by \code{func}. The \R
    calling sequence for \code{jacfunc} is identical to that of
    \code{func}.
    
    If the Jacobian is a full matrix, \code{jacfunc} should return a
    matrix \eqn{\partial\dot{y}/\partial y}{dydot/dy}, where the ith row contains the derivative of
    \eqn{dy_i/dt} with respect to \eqn{y_j}, or a vector containing the
    matrix elements by columns (the way \R and FORTRAN store matrices).
    
    If the Jacobian is banded, \code{jacfunc} should return a matrix
    containing only the nonzero bands of the Jacobian, rotated
    row-wise. See first example of lsode.
  }
  \item{jacres }{ \code{jacres} and not \code{jacfunc} should be used if
    the system is specified by the residual function \eqn{F(t, y, y')},
    i.e. \code{jacres} is used in conjunction with \code{res}.
    
    If \code{jacres} is an \R-function, the calling sequence for
    \code{jacres} is identical to that of \code{res}, but with extra
    parameter \code{cj}.  Thus it should be called as: \code{jacres = 
    func(t, y, dy, parms, cj, ...)}.  Here \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system, \eqn{y'} are the corresponding derivatives
    and \code{cj} is a scalar, which is normally proportional to
    the inverse of the stepsize. If the initial \code{y} or \code{dy}
    have a \code{names} attribute, the names will be available inside
    \code{jacres}, unless
    \code{ynames} is \code{FALSE}. \code{parms} is a vector of
    parameters (which may have a names attribute).
    
    If the Jacobian is a full matrix, \code{jacres} should return the
    matrix \eqn{dG/dy + c_j\cdot dG/dy'}{dG/d y + cj*dG/d y'}, where the \eqn{i}th row is the sum of the
    derivatives of \eqn{G_i} with respect to \eqn{y_j} and the scaled
    derivatives of \eqn{G_i} with respect to \eqn{y'_j}.
    
    If the Jacobian is banded, \code{jacres} should return only the
    nonzero bands of the Jacobian, rotated rowwise. See details for the
    calling sequence when \code{jacres} is a string.
  }
  \item{jactype }{the structure of the Jacobian, one of
    \code{"fullint"}, \code{"fullusr"}, \code{"bandusr"} or
    \code{"bandint"} - either full or banded and estimated internally or
    by the user.
  }
  \item{mass }{the mass matrix. 
      If not \code{NULL}, the problem is a linearly
     implicit DAE and defined as \eqn{M\, dy/dt = f(t,y)}{M dy/dt = f(t,y)}. The mass-matrix \eqn{M}
     should be of dimension \eqn{n*n} where \eqn{n} is the 
     number of \eqn{y}-values.

     If \code{mass=NULL} then the model is either an ODE or a DAE, specified with
     \code{res}
  }
  \item{estini }{only if a DAE system, and if initial values of \code{y}
    and \code{dy} are not consistent (i.e. \eqn{F(t,y,dy) \neq 0}{F(t,
     y, dy) != 0}), setting
    \code{estini} = 1 or 2, will solve for them.  If \code{estini} = 1: dy
    and the algebraic variables are estimated from \code{y}; in this
    case, the number of algebraic equations must be given (\code{nalg}).
    If \code{estini} = 2: \code{y} will be estimated from \code{dy}.
  }
  \item{verbose }{if TRUE: full output to the screen, e.g. will
    print the \code{diagnostiscs} of the integration - see details.
  }
  \item{tcrit }{the FORTRAN routine \code{daspk} overshoots its targets
    (times points in the vector \code{times}), and interpolates values
    for the desired time points.  If there is a time beyond which
    integration should not proceed (perhaps because of a singularity),
    that should be provided in \code{tcrit}.
  }
  \item{hmin }{an optional minimum value of the integration stepsize. In
    special situations this parameter may speed up computations with the
    cost of precision. Don't use \code{hmin} if you don't know why!
  }
  \item{hmax }{an optional maximum value of the integration stepsize. If
    not specified, \code{hmax} is set to the largest difference in
    \code{times}, to avoid that the simulation possibly ignores
    short-term events. If 0, no maximal size is specified.
  }
  \item{hini }{initial step size to be attempted; if 0, the initial step
    size is determined by the solver
  }
  \item{ynames }{logical, if \code{FALSE}, names of state variables are not
    passed to function \code{func}; this may speed up the simulation especially
    for large models.
  }
  \item{maxord }{the maximum order to be allowed. Reduce \code{maxord}
    to save storage space ( <= 5)
  }
  \item{bandup }{number of non-zero bands above the diagonal, in case
    the Jacobian is banded (and \code{jactype} one of
    "bandint", "bandusr")
  }
  \item{banddown }{number of non-zero bands below the diagonal, in case
    the Jacobian is banded (and \code{jactype} one of
    "bandint", "bandusr")
  }
  \item{maxsteps }{maximal number of steps per output interval taken by the
    solver; will be recalculated to be at least 500 and a multiple of
    500; if \code{verbose} is \code{TRUE} 
    the solver will give a warning if more than 500 steps are
    taken, but it will continue till \code{maxsteps} steps. 
    (Note this warning was always given in deSolve versions < 1.10.3).
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions referred to in \code{res} and
    \code{jacres}.  See package vignette \code{"compiledCode"}.
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
    specified by \code{res} and \code{jacres}.
  }
  \item{ipar }{only when \file{dllname} is specified: a vector with
    integer values passed to the dll-functions whose names are specified
    by \code{res} and \code{jacres}.
  }
  \item{nout }{only used if \file{dllname} is specified and the model is
    defined in compiled code: the number of output variables calculated
    in the compiled function \code{res}, present in the shared
    library. Note: it is not automatically checked whether this is
    indeed the number of output variables calculated in the dll - you have
    to perform this check in the code - See package vignette
    \code{"compiledCode"}.
  }
  \item{outnames }{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{res}, present in the shared library.
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
    See \link{forcings} or vignette \code{compiledCode}.
  }
  \item{events }{A list that specifies events, i.e. when the value of a
   state variable is suddenly changed. See \link{events} for more information.
  }
  \item{lags }{A list that specifies timelags, i.e. the number of steps
   that has to be kept. To be used for delay differential equations. 
   See \link{timelags}, \link{dede} for more information.
  }
  \item{... }{additional arguments passed to \code{func},
    \code{jacfunc}, \code{res} and \code{jacres}, allowing this to be a
    generic function.
  }
}
\value{
  A matrix of class \code{deSolve} with up to as many rows as elements in
  \code{times} and as many
  columns as elements in \code{y} plus the number of "global" values
  returned in the next elements of the return from \code{func} or
  \code{res}, plus an additional column (the first) for the time value.
  There will be one row for each element in \code{times} unless the
  FORTRAN routine `daspk' returns with an unrecoverable error.  If
  \code{y} has a names attribute, it will be used to label the columns
  of the output value.
      
}
\author{Karline Soetaert <karline.soetaert@nioz.nl>}
\examples{
## =======================================================================
## Coupled chemical reactions including an equilibrium
## modeled as (1) an ODE and (2) as a DAE
##
## The model describes three chemical species A,B,D:
## subjected to equilibrium reaction D <- > A + B
## D is produced at a constant rate, prod
## B is consumed at 1s-t order rate, r
## Chemical problem formulation 1: ODE
## =======================================================================

## Dissociation constant
K <- 1 

## parameters
pars <- c(
        ka   = 1e6,     # forward rate
        r    = 1,
        prod = 0.1)


Fun_ODE <- function (t, y, pars)
{
  with (as.list(c(y, pars)), {
    ra  <- ka*D        # forward rate
    rb  <- ka/K *A*B   # backward rate

    ## rates of changes
    dD  <- -ra + rb + prod
    dA  <-  ra - rb
    dB  <-  ra - rb - r*B
    return(list(dy = c(dA, dB, dD),
                CONC = A+B+D))
  })
}

## =======================================================================
## Chemical problem formulation 2: DAE
## 1. get rid of the fast reactions ra and rb by taking
## linear combinations   : dD+dA = prod (res1) and
##                         dB-dA = -r*B (res2)
## 2. In addition, the equilibrium condition (eq) reads:
## as ra = rb : ka*D = ka/K*A*B = >      K*D = A*B
## =======================================================================

Res_DAE <- function (t, y, yprime, pars)
{
  with (as.list(c(y, yprime, pars)), {

    ## residuals of lumped rates of changes
    res1 <- -dD - dA + prod
    res2 <- -dB + dA - r*B
    
    ## and the equilibrium equation
    eq   <- K*D - A*B

    return(list(c(res1, res2, eq),
                CONC = A+B+D))
  })
}

## =======================================================================
## Chemical problem formulation 3: Mass * Func
## Based on the DAE formulation
## =======================================================================

Mass_FUN <- function (t, y, pars) {
  with (as.list(c(y, pars)), {

    ## as above, but without the 
    f1 <- prod
    f2 <- - r*B
    
    ## and the equilibrium equation
    f3   <- K*D - A*B

    return(list(c(f1, f2, f3),
                CONC = A+B+D))
  })
}
Mass <- matrix(nrow = 3, ncol = 3, byrow = TRUE, 
  data=c(1,  0, 1,         # dA + 0 + dB
        -1,  1, 0,         # -dA + dB +0
         0,  0, 0))        # algebraic
         
times <- seq(0, 100, by = 2)

## Initial conc; D is in equilibrium with A,B
y     <- c(A = 2, B = 3, D = 2*3/K)

## ODE model solved with daspk
ODE <- daspk(y = y, times = times, func = Fun_ODE,
                     parms = pars, atol = 1e-10, rtol = 1e-10)

## Initial rate of change
dy  <- c(dA = 0, dB = 0, dD = 0) 

## DAE model solved with daspk
DAE <- daspk(y = y, dy = dy, times = times,
         res = Res_DAE, parms = pars, atol = 1e-10, rtol = 1e-10)

MASS<- daspk(y=y, times=times, func = Mass_FUN, parms = pars, mass = Mass)

## ================
## plotting output
## ================

plot(ODE, DAE, xlab = "time", ylab = "conc", type = c("l", "p"),
     pch = c(NA, 1))

legend("bottomright", lty = c(1, NA), pch = c(NA, 1),
  col = c("black", "red"), legend = c("ODE", "DAE"))

# difference between both implementations:
max(abs(ODE-DAE))

## =======================================================================
## same DAE model, now with the Jacobian
## =======================================================================
jacres_DAE <- function (t, y, yprime, pars, cj)
{
    with (as.list(c(y, yprime, pars)), {
##    res1 = -dD - dA + prod
      PD[1,1] <- -1*cj      # d(res1)/d(A)-cj*d(res1)/d(dA)
      PD[1,2] <- 0          # d(res1)/d(B)-cj*d(res1)/d(dB)
      PD[1,3] <- -1*cj      # d(res1)/d(D)-cj*d(res1)/d(dD)
##     res2 = -dB + dA - r*B
      PD[2,1] <- 1*cj
      PD[2,2] <- -r -1*cj
      PD[2,3] <- 0
##    eq = K*D - A*B
      PD[3,1] <- -B
      PD[3,2] <- -A
      PD[3,3] <- K
      return(PD)
   })
}

PD <- matrix(ncol = 3, nrow = 3, 0)

DAE2 <- daspk(y = y, dy = dy, times = times,
          res = Res_DAE, jacres = jacres_DAE, jactype = "fullusr",
          parms = pars, atol = 1e-10, rtol = 1e-10)
         
max(abs(DAE-DAE2))

## See \dynload subdirectory for a FORTRAN implementation of this model

## =======================================================================
## The chemical model as a DLL, with production a forcing function
## =======================================================================
times <- seq(0, 100, by = 2)

pars <- c(K = 1, ka   = 1e6, r    = 1)

## Initial conc; D is in equilibrium with A,B
y     <- c(A = 2, B = 3, D = as.double(2*3/pars["K"]))

## Initial rate of change
dy  <- c(dA = 0, dB = 0, dD = 0)

# production increases with time
prod <- matrix(ncol = 2, 
               data = c(seq(0, 100, by = 10), 0.1*(1+runif(11)*1)))

ODE_dll <- daspk(y = y, dy = dy, times = times, res = "chemres",
          dllname = "deSolve", initfunc = "initparms",
          initforc = "initforcs", parms = pars, forcings = prod,
          atol = 1e-10, rtol = 1e-10, nout = 2, 
          outnames = c("CONC","Prod"))

plot(ODE_dll, which = c("Prod", "D"), xlab = "time",
     ylab = c("/day", "conc"), main = c("production rate","D"))

}
\references{
  L. R. Petzold, A Description of DASSL: A Differential/Algebraic
  System Solver, in Scientific Computing, R. S. Stepleman et al.
  (Eds.), North-Holland, Amsterdam, 1983, pp. 65-68.

  K. E. Brenan, S. L. Campbell, and L. R. Petzold, Numerical
  Solution of Initial-Value Problems in Differential-Algebraic
  Equations, Elsevier, New York, 1989.

  P. N. Brown and A. C. Hindmarsh, Reduced Storage Matrix Methods
  in Stiff ODE Systems, J. Applied Mathematics and Computation,
  31 (1989), pp. 40-91.

  P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Using Krylov
  Methods in the Solution of Large-Scale Differential-Algebraic
  Systems, SIAM J. Sci. Comp., 15 (1994), pp. 1467-1488.

  P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Consistent
  Initial Condition Calculation for Differential-Algebraic
  Systems, LLNL Report UCRL-JC-122175, August 1995; submitted to
  SIAM J. Sci. Comp.
  
  Netlib: \url{http://www.netlib.org}
}
\details{
  The daspk solver uses the backward differentiation formulas of orders
  one through five (specified with \code{maxord}) to solve either:
  \itemize{
    \item an ODE system of the form \deqn{y' = f(t,y,...)}  
    or
    \item a DAE system of the form \deqn{y' = M f(t,y,...)} 
    or
    \item a DAE system of the form \deqn{F(t,y,y') = 0}. 
    The index of the DAE should be preferable <= 1.
  }
  
  ODEs are specified using argument \code{func}, 
  DAEs are specified using argument \code{res}.
  
  If a DAE system, Values for y \emph{and} y' (argument \code{dy}) 
  at the initial time must be given as input. Ideally, these values should be consistent,
  that is, if t, y, y' are the given initial values, they should
  satisfy F(t,y,y') = 0.  \cr However, if consistent values are not
  known, in many cases daspk can solve for them: when \code{estini} = 1,
  y' and algebraic variables (their number specified with \code{nalg})
  will be estimated, when \code{estini} = 2, y will be estimated.
  
  The form of the \bold{Jacobian} can be specified by
  \code{jactype}. This is one of:
  
  \describe{
    \item{jactype = "fullint":}{a full Jacobian, calculated internally
      by \code{daspk}, the default,
    }
    \item{jactype = "fullusr":}{a full Jacobian, specified by user
      function \code{jacfunc} or \code{jacres},
    }
    \item{jactype = "bandusr":}{a banded Jacobian, specified by user
      function \code{jacfunc} or \code{jacres}; the size of the bands
      specified by \code{bandup} and \code{banddown},
    }
    \item{jactype = "bandint":}{a banded Jacobian, calculated by
      \code{daspk}; the size of the bands specified by \code{bandup} and
      \code{banddown}.
    }
  }
  
  If \code{jactype} = "fullusr" or "bandusr" then the user must supply a
  subroutine \code{jacfunc}.

  If jactype = "fullusr" or "bandusr" then the user must supply a
  subroutine \code{jacfunc} or \code{jacres}.
  
  The input parameters \code{rtol}, and \code{atol} determine the
  \bold{error control} performed by the solver.  If the request for
  precision exceeds the capabilities of the machine, \code{daspk} will return
  an error code. See \code{\link{lsoda}} for details.
  
  When the index of the variables is specified (argument \code{nind}), 
  and higher index variables
  are present, then the equations are scaled such that equations corresponding
  to index 2 variables are multiplied with 1/h, for index 3 they are multiplied
  with 1/h^2, where h is the time step. This is not in the standard DASPK code,
  but has been added for consistency with solver \link{radau}. Because of this,
  daspk can solve certain index 2 or index 3 problems.
  
  \bold{res and jacres} may be defined in compiled C or FORTRAN code, as
  well as in an R-function. See package vignette \code{"compiledCode"}
  for details.  Examples
  in FORTRAN are in the \file{dynload} subdirectory of the
  \code{deSolve} package directory.

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
    \item \code{\link{radau}} for integrating DAEs up to index 3,
    \item \code{\link{rk}}, \item \code{\link{rk4}} and \code{\link{euler}} for
      Runge-Kutta integrators.
    \item \code{\link{lsoda}}, \code{\link{lsode}},
      \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}},
      for other solvers of the Livermore family,
    \item \code{\link{ode}} for a general interface to most of the ODE solvers,
    \item \code{\link{ode.band}} for solving models with a banded
      Jacobian,
    \item \code{\link{ode.1D}} for integrating 1-D models,
    \item \code{\link{ode.2D}} for integrating 2-D models,
    \item \code{\link{ode.3D}} for integrating 3-D models,
  }

  \code{\link{diagnostics}} to print diagnostic messages.
}
\note{
  In this version, the Krylov method is not (yet) supported.

  From \code{deSolve} version 1.10.4 and above, the following changes were made
  
  \enumerate{
  \item the argument list to \code{daspk} 
  now also includes \code{nind}, the index of each variable. 
  This is used to scale the variables, such that \code{daspk} in R can also solve
  certain index 2 or index 3 problems, which the original Fortran version 
  may not be able to solve.

  \item the default of \code{atol} was changed from 1e-8 to 1e-6,
    to be consistent with the other solvers.

  \item the multiple warnings from daspk when the number of steps exceed 500
  were toggled off unless \code{verbose} is \code{TRUE}
  }
}
\keyword{math}

