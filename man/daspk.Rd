\name{daspk}
\alias{daspk}
\title{General Solver for Differential Algebraic Equations (DAE)}
\description{
  Solves either:
  \itemize{
    \item a system of ordinary differential equations (ODE) of the form
      \deqn{y' = f(t,y,...)} or
    \item a system of differential algebraic equations (DAE) of the form
      \deqn{F(t,y,y') = 0}
  }
  using a combination of backward differentiation formula (BDF) and a
  direct linear system solution method (dense or banded).
  
  The \R function \code{daspk} provides an interface to the Fortran DAE
  solver of the same name, written by Linda R. Petzold, Peter N. Brown,
  Alan C. Hindmarsh and Clement W. Ulrich.
    
  The system of DE's is written as an \R function (which may, of course,
  use \code{\link{.C}}, \code{.Fortran}, \code{\link{.Call}}, etc., to
  call foreign code) or be defined in compiled code that has been
  dynamically loaded.
}
\usage{
daspk(y, times, func = NULL, parms, dy = NULL, res = NULL, nalg = 0,
  rtol = 1e-6, atol = 1e-8, jacfunc = NULL, jacres = NULL,
  jactype = "fullint", estini = NULL, verbose = FALSE, 
  tcrit = NULL, hmin = 0, hmax = NULL, hini = 0, ynames = TRUE, maxord = 5,
  bandup = NULL, banddown = NULL, maxsteps = 5000,   
  dllname = NULL, initfunc = dllname, initpar = parms, 
  rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, 
  ...)
}

\arguments{
  \item{y }{the initial (state) values for the DE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time; if only one step is
    to be taken; set times = NULL.
  }
  \item{func }{cannot be used if the model is a DAE system. If an ODE
    system, \code{func} should be an \R-function that computes the
    values of the derivatives in the ODE system (the \emph{model
    definition}) at time t.
  
    \code{func} must be defined as: \code{yprime = func(t, y,
    parms,...)}.  \cr \code{t} is the current time point in the
    integration, \code{y} is the current estimate of the variables in
    the ODE system.  If the initial values \code{y} has a names
    attribute, the names will be available inside \code{func}, unless
    \code{ynames} is FALSE.  \code{parms} is a vector or list of
    parameters.  \code{...} (optional) are any other arguments passed to
    the function.  The return value of \code{func} should be a list,
    whose first element is a vector containing the derivatives of
    \code{y} with respect to \code{time}, and whose next elements are
    global values that are required at each point in \code{times}.
  
    Note that it is not possible to define \code{func} as a compiled
    function in a dynamically loaded shared library. Use \code{res}
    instead.
  }
  \item{parms }{vector or list of parameters used in \code{func},
    \code{jacfunc}, or \code{res}
  }
  \item{dy }{the initial derivatives of the state variables of the DE
    system.  Ignored if an ODE.
  }
  \item{res }{if a DAE system: either an \R-function that computes the
    residual function F(t,y,y') of the DAE system (the model
    defininition) at time \code{t}, or a character string giving the
    name of a compiled function in a dynamically loaded shared library.
      
    If \code{res} is a user-supplied R-function, it must be called as:
    \code{F = func(t, y, dy, parms, ...)}.
    
    Here \code{t} is the current time point in the integration, \code{y}
    is the current estimate of the variables in the ODE system,
    \code{dy} are the corresponding rates of change.  If the initial
    \code{y} or \code{dy} have a names attribute, the names will be
    available inside \code{func}, unless \code{ynames} is \code{FALSE}.
    \code{parms} is a vector of parameters.
    
    The return value of \code{res} should be a list, whose first element
    is a vector containing the residuals of the DAE system,
    i.e. delta = F(t,y,y'), and whose next elements contain output
    variables that are required at each point in \code{times}.
  
    If \code{res} is a string, then \code{dllname} must give the name of
    the shared library (without extension) which must be loaded before
    \code{daspk()} is called (see package vignette for more
    information).
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
    jacobian of the system of differential equations. Only used in case
    the system is an ODE (y' = f(t,y)), specified by \code{func}. The \R
    calling sequence for \code{jacfunc} is identical to that of
    \code{func}.
    
    If the jacobian is a full matrix, \code{jacfunc} should return a
    matrix dydot/dy, where the ith row contains the derivative of
    \eqn{dy_i/dt} with respect to \eqn{y_j}, or a vector containing the
    matrix elements by columns (the way \R and Fortran store matrices).
    
    If the jacobian is banded, \code{jacfunc} should return a matrix
    containing only the nonzero bands of the jacobian, rotated
    row-wise. See first example of lsode.
  }
  \item{jacres }{ \code{jacres} and not \code{jacfunc} should be used if
    the system is specified by the residual function F(t,y,y'),
    i.e. \code{jacres} is used in conjunction with \code{res}.
    
    If \code{jacres} is an R-function, the calling sequence for
    \code{jacres} is identical to that of \code{res}, but with extra
    parameter \code{cj}.  Thus it should be called as: \code{jacres = 
    func(t, y, dy, parms, cj, ...)}.  Here \code{t} is the current time
    point in the integration, \code{y} is the current estimate of the
    variables in the ODE system, \eqn{y'} are the corresponding rates of
    change and \code{cj} is a scalar, which is normally proportional to
    the inverse of the stepsize. If the initial y or dy have a names
    attribute, the names will be available inside jacres, unless
    \code{ynames} is \code{FALSE}. \code{parms} is a vector of
    parameters (which may have a names attribute).
    
    If the jacobian is a full matrix, \code{jacres} should return the
    matrix dG/dy + cj*dG/dyprime, where the ith row is the sum of the
    derivatives of \eqn{G_i} with respect to \eqn{y_j} and the scaled
    derivatives of \eqn{G_i} with respect to \eqn{dy_j}.
    
    Ff the jacobian is banded, \code{jacres} should return only the
    nonzero bands of the jacobian, rotated rowwise. See details for the
    calling sequence when \code{jacres} is a string.
  }
  \item{jactype }{the structure of the jacobian, one of
    \code{"fullint"}, \code{"fullusr"}, \code{"bandusr"} or
    \code{"bandint"} - either full or banded and estimated internally or
    by the user.
  }
  \item{estini }{only if a DAE system, and if initial values of \code{y}
    and \code{dy} are not consistent (i.e. F(t,y,dy) is not = 0), setting
    \code{estini} = 1 or 2, will solve for them.  If \code{estini} = 1: dy
    and the algebraic variables are estimated from \code{y}; in this
    case, the number of algebraic equations must be given (\code{nalg}).
    If \code{estini} = 2: \code{y} will be estimated from \code{dy}.
  }
  \item{verbose }{if TRUE: full output to the screen, e.g. will output
    the settings of vector \bold{istate} and \bold{rstate} - see details.
  }
  \item{tcrit }{the Fortran routine \code{daspk} overshoots its targets
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
  \item{ynames }{if FALSE: names of state variables are not passed to
    function \code{func} ; this may speed up the simulation especially
    for large models.
  }
  \item{maxord }{the maximum order to be allowed. Reduce \code{maxord}
    to save storage space ( <= 5)
  }
  \item{bandup }{number of non-zero bands above the diagonal, in case
    the jacobian is banded (and \code{jactype} one of
    "bandint","bandusr")
  }
  \item{banddown }{number of non-zero bands below the diagonal, in case
    the jacobian is banded (and \code{jactype} one of
    "bandint","bandusr")
  }
  \item{maxsteps }{maximal number of steps during one call to the
    solver; will be recalculated to be at least 500 and a multiple of
    500; the solver will give a warning if more than 500 steps are
    taken, but it will continue till \code{maxsteps} steps.
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions referred to in \code{func}/or \code{res} and
    \code{jacfunc}.  See package vignette.
  }
  \item{initfunc }{if not NULL, the name of the initialisation function
    (which initialises values of parameters), as provided in
    \file{dllname}. See package vignette.
  }
  \item{initpar }{only when \file{dllname} is specified and an
    initialisation function \code{initfunc} is in the dll: the
    parameters passed to the initialiser, to initialise the common
    blocks (fortran) or global variables (C, C++).
  }
  \item{rpar }{only when \file{dllname} is specified: a vector with
    double precision values passed to the dll-functions whose names are
    specified by \code{func} and \code{jacfunc}.
  }
  \item{ipar }{only when \file{dllname} is specified: a vector with
    integer values passed to the dll-functions whose names are specified
    by \code{func} and \code{jacfunc}.
  }
  \item{nout }{only used if \file{dllname} is specified and the model is
    defined in compiled code: the number of output variables calculated
    in the compiled function \code{func}, present in the shared
    library. Note: it is not automatically checked whether this is
    indeed the number of output variables calculed in the dll - you have
    to perform this check in the code - See package vignette.
  }
  \item{outnames }{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{func}, present in the shared library.
  }
  \item{... }{additional arguments passed to \code{func},
    \code{jacfunc}, \code{res} and \code{jacres}, allowing this to be a
    generic function.
  }
}
\value{
  A matrix with up to as many rows as elements in times and as many
  columns as elements in \code{y} plus the number of "global" values
  returned in the next elements of the return from \code{func} or
  \code{res}, plus an additional column (the first) for the time value.
  There will be one row for each element in \code{times} unless the
  Fortran routine `daspk' returns with an unrecoverable error.  If
  \code{y} has a names attribute, it will be used to label the columns
  of the output value.
      
  The output will have the attributes \code{istate}, and \code{rstate},
  two vectors with several useful elements.  See details.  The first
  element of istate returns the conditions under which the last call to
  lsoda returned. Normal is \code{istate[1] = 2}.  If \code{verbose} =
  TRUE, the settings of istate and rstate will be written to the screen.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
## ====================================================
## Coupled chemical reactions including an equilibrium
## modeled as (1) an ODE and (2) as a DAE
##
## The model describes three chemical species A,B,D:
## subjected to equilibrium reaction D <- > A + B
## D is produced at a constant rate, prod
## B is consumed at 1s-t order rate, r
## ====================================================

## Dissociation constant
K <- 1 

## parameters
pars <- c(
        ka   = 1e6,     # forward rate
        r    = 1,
        prod = 0.1)

## ====================================
## Chemical problem formulation 1: ODE
## ====================================

Fun_ODE <- function (t, y, pars)
{
  with (as.list(c(y, pars)), {
    ra  = ka*D        # forward rate
    rb  = ka/K *A*B   # backward rate

    ## rates of changes
    dD  = -ra + rb + prod
    dA  =  ra - rb
    dB  =  ra - rb - r*B
    return(list(dy = c(dA, dB, dD),
                CONC = A+B+D))
  })
}

## ======================================================
## Chemical problem formulation 2: DAE
## 1. get rid of the fast reactions ra and rb by taking
## linear combinations   : dD+dA = prod (res1) and
##                         dB-dA = -r*B (res2)
## 2. In addition, the equilibrium condition (eq) reads:
## as ra = rb : ka*D = ka/K*A*B = >      K*D = A*B
## ======================================================

Res_DAE <- function (t, y, yprime, pars)
{
  with (as.list(c(y, yprime, pars)), {

    ## residuals of lumped rates of changes
    res1 = -dD - dA + prod
    res2 = -dB + dA - r*B
    
    ## and the equilibrium equation
    eq = K*D - A*B

    return(list(c(res1, res2, eq),
                CONC = A+B+D))
  })
}

times <- seq(0, 100, by = 2)

## Initial conc; D is in equilibrium with A,B
y     <- c(A = 2, B = 3, D = 2*3/K)

## ODE model solved with daspk
ODE <- as.data.frame(daspk(y = y, times = times, func = Fun_ODE,
                     parms = pars, atol = 1e-10, rtol = 1e-10))

## Initial rate of change
dy  <- c(dA = 0, dB = 0, dD = 0) 

## DAE model solved with daspk
DAE <- as.data.frame(daspk(y = y, dy = dy, times = times,
         res = Res_DAE, parms = pars, atol = 1e-10, rtol = 1e-10))

## ================
## plotting output
## ================
opa <- par(mfrow = c(2,2))
for (i in 2:5) 
{
plot(ODE$time,ODE[,i],xlab = "time",
     ylab = "conc",main = names(ODE)[i],type = "l")
points(DAE$time,DAE[,i],col = "red")
}
legend("bottomright",lty = c(1,NA),pch = c(NA,1),
       col = c("black","red"),legend = c("ODE","DAE"))      
       
# difference between both implementations:
max(abs(ODE-DAE))

## ======================================
## same DAE model, now with the jacobian
## ======================================
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

PD <- matrix(nc = 3, nr = 3, 0)

DAE2 <- as.data.frame(daspk(y = y, dy = dy, times = times,
          res = Res_DAE, jacres = jacres_DAE, jactype = "fullusr",
          parms = pars, atol = 1e-10, rtol = 1e-10))
         
max(abs(DAE-DAE2))

## See \dynload subdirectory for a FORTRAN implementation of this model
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
    \item an ODE system of the form \deqn{y' = f(t,y,...)}  for y = Y,
    or
    \item a DAE system of the form \deqn{F(t,y,y') = 0} for y = Y and y'
    = YPRIME.
  }
  
  ODEs are specified in \code{func}, DAEs are specified in \code{res}.
  
  If a DAE system, Values for Y \emph{and} YPRIME at the initial time
  must be given as input. Ideally,these values should be consistent,
  that is, if T, Y, YPRIME are the given initial values, they should
  satisfy F(T,Y,YPRIME) = 0.  \cr However, if consistent values are not
  known, in many cases daspk can solve for them: when \code{estini} = 1,
  y' and algebraic variables (their number specified with \code{nalg})
  will be estimated, when \code{estini} = 2, y will be estimated.
  
  The form of the \bold{jacobian} can be specified by
  \code{jactype}. This is one of:
  
  \describe{
    \item{jactype = "fullint":}{a full jacobian, calculated internally
      by \code{daspk}, the default,
    }
    \item{jactype = "fullusr":}{a full jacobian, specified by user
      function \code{jacfunc} or \code{jacres},
    }
    \item{jactype = "bandusr":}{a banded jacobian, specified by user
      function \code{jacfunc} or \code{jacres}; the size of the bands
      specified by \code{bandup} and \code{banddown},
    }
    \item{jactype = "bandint":}{a banded jacobian, calculated by
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
  precision exceeds the capabilities of the machine, daspk will return
  an error code. See \code{\link{lsoda}} for details.
  
  \bold{res and jacres} may be defined in compiled C or Fortran code, as
  well as in an R-function. See package vignette for details.  Examples
  in Fortran are in the \file{dynload} subdirectory of the
  \code{deSolve} package directory.

  The output will have the attributes \bold{istate}, and \bold{rstate}, two
  vectors with several useful elements.

  If \code{verbose} = TRUE, the settings of istate and rstate will be
  written to the screen.
  
  The following elements of \bold{istate} are meaningful:
  \itemize{
    \item el 1 : returns the conditions under which the last call to
      daspk returned.
      2  if daspk was successful,
      -1 if excess work done,
      -2 means excess accuracy requested. (Tolerances too small),
      -3 means illegal input detected. (See printed message.),
      -4 means repeated error test failures. (Check all input),
      -5 means repeated convergence failures.
        (Perhaps bad Jacobian supplied or wrong choice of MF or tolerances.),
      -6 means error weight became zero during problem.
	(Solution component i vanished, and atol or atol(i) = 0.)
    \item el 12 : The number of steps taken for the problem so far.
    \item el 13 : The number of function evaluations for the problem so
      far.,
    \item el 14 : The number of Jacobian evaluations so far.,
    \item el 15 : The method order last used (successfully).,
    \item el 16 : The order to be attempted on the next step.,
    \item el 17 : if el 1 = -4,-5: the largest component in the error
      vector,
    \item el 20 : The number of matrix LU decompositions so far.,
    \item el 21 : The number of nonlinear (Newton) iterations so far.,
    \item el 22 : The number of convergence failures of the solver so
      far ,
    \item el 23 : The number of error test failures of the integrator so
      far.
  }
  \bold{rstate} contains the following:
  \itemize{
    \item 1: The step size in t last used (successfully).
    \item 2: The step size to be attempted on the next step.
    \item 3: The current value of the independent variable which the
      solver has actually reached, i.e. the current internal mesh point
      in t.
    \item 4: A tolerance scale factor, greater than 1.0, computed when a
      request for too much accuracy was detected.
  }
  For more information about possible options, see the comments in the
  original code daspk.f.
}
\seealso{
  \code{\link{lsoda}}, \code{\link{lsodar}},
  \code{\link{lsode}}, \code{\link{lsodes}}, \code{\link{vode}}
}
\note{
  In this version, the krylov method is not (yet) supported.
}
\keyword{math}

