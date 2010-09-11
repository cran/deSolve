\name{radau}
\alias{radau}

\title{Implicit Runge-Kutta RADAU IIA}

\description{
  Solves the initial value problem for stiff or nonstiff systems of
  ordinary differential equations (ODE) in the form: \deqn{dy/dt = 
  f(t,y)} or linearly implicit differential algebraic equations in the 
  form: \deqn{M dy/dt = f(t,y)}.
  
  The \R function \code{radau} provides an interface to the Fortran solver 
  RADAU5, written by Ernst Hairer and G. Wanner, which implements the 3-stage
  RADAU IIA method.
  It implements the implicit Runge-Kutta method of order 5 with step size 
  control and continuous output.
  The system of ODEs or DAEs is written as an \R function or can be defined in
  compiled code that has been dynamically loaded.
}
\usage{
radau(y, times, func, parms, nind = c(length(y), 0, 0), 
  rtol = 1e-6, atol = 1e-6, jacfunc = NULL, jactype = "fullint", 
  mass = NULL, massup = NULL, massdown = NULL,
  verbose = FALSE, hmax = NULL, hini = 0, ynames = TRUE, 
  bandup = NULL, banddown = NULL, maxsteps = 5000, 
  dllname = NULL, initfunc = dllname, initpar = parms, 
  rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, 
  forcings = NULL, initforc = NULL, fcontrol = NULL, ...)
}

\arguments{
  \item{y }{the initial (state) values for the ODE system. If \code{y}
    has a name attribute, the names will be used to label the output
    matrix.
  }
  \item{times }{time sequence for which output is wanted; the first
    value of \code{times} must be the initial time; if only one step is
    to be taken; set \code{times} = \code{NULL}.
  }
  \item{func }{either an \R-function that computes the values of the
    derivatives in the ODE system (the \emph{model definition}) at time
    t, or the right-hand side of the equation \deqn{M dy/dt = f(t,y)}
    if a DAE. (if \code{mass} is supplied then the problem is assumed a DAE).
    
    \code{func} can also be a character string giving the name of a compiled 
    function in a dynamically loaded shared library.
  
    If \code{func} is an \R-function, it must be defined as:
    
    \code{func <- function(t, y, parms,...)}.

    \code{t} is the current time
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
    \code{radau()} is called. See deSolve package vignette \code{"compiledCode"}
    for more details.
  }
  \item{parms }{vector or list of parameters used in \code{func} or
    \code{jacfunc}.
  }
  \item{nind }{if a DAE system: a three-valued vector with the number of
    variables of index 1, 2, 3 respectively.
    The equations must be defined such that the index 1 variables precede
    the index 2 variables which in turn precede the index 3 variables.
    The sum of the variables of different index should equal N,
    the total number of variables.
  }
  \item{rtol }{relative error tolerance, either a
    scalar or an array as long as \code{y}. See details.
  }
  \item{atol }{absolute error tolerance, either a scalar or an array as
    long as \code{y}. See details.
  }
  \item{jacfunc }{if not \code{NULL}, an \R function that computes the
    Jacobian of the system of differential equations
    \eqn{\partial\dot{y}_i/\partial y_j}{dydot(i)/dy(j)}, or
    a string giving the name of a function or subroutine in
    \file{dllname} that computes the Jacobian (see vignette
    \code{"compiledCode"} from package deSolve, for more about this option).

    In some circumstances, supplying
    \code{jacfunc} can speed up the computations, if the system is
    stiff.  The \R calling sequence for \code{jacfunc} is identical to
    that of \code{func}.

    If the Jacobian is a full matrix,
    \code{jacfunc} should return a matrix \eqn{\partial\dot{y}/\partial y}{dydot/dy}, where the ith row
    contains the derivative of \eqn{dy_i/dt} with respect to \eqn{y_j},
    or a vector containing the matrix elements by columns (the way \R
    and FORTRAN store matrices).  \cr If the Jacobian is banded,
    \code{jacfunc} should return a matrix containing only the nonzero
    bands of the Jacobian, rotated row-wise. See example.
  }
  \item{jactype }{the structure of the Jacobian, one of
    \code{"fullint"}, \code{"fullusr"}, \code{"bandusr"} or
    \code{"bandint"} - either full or banded and estimated internally or
    by user.
  }
  \item{mass }{the mass matrix. 
      If not \code{NULL}, the problem is a linearly
     implicit DAE and defined as \eqn{M\, dy/dt = f(t,y)}{M dy/dt = f(t,y)}. If the mass-matrix \eqn{M}
     is full, it should be of dimension \eqn{n^2}{n*n} where \eqn{n} is the 
     number of \eqn{y}-values; if banded the number of rows should be less than \eqn{n},
     and the mass-matrix is stored diagonal-wise with element \eqn{(i, j)} stored in 
    \code{mass(i - j + mumas + 1, j)}.

     If \code{mass = NULL} then the model is an ODE (default)  
  }
  \item{massup }{number of non-zero bands above the diagonal of the \code{mass} 
    matrix, in case it is banded.
  }
  \item{massdown }{number of non-zero bands below the diagonal of the \code{mass} 
  matrix, in case it is banded.
  }
  \item{verbose }{if \code{TRUE}: full output to the screen, e.g. will
    print the \code{diagnostiscs} of the integration - see details.
  }
  \item{hmax }{an optional maximum value of the integration stepsize. If
    not specified, \code{hmax} is set to the largest difference in
    \code{times}, to avoid that the simulation possibly ignores
    short-term events. If 0, no maximal size is specified.
  }
  \item{hini }{initial step size to be attempted; if 0, the initial step
    size is set equal to 1e-6. Usually 1e-3 to 1e-5 is good for stiff equations
  }
  \item{ynames }{logical, if \code{FALSE} names of state variables are not
    passed to function \code{func}; this may speed up the simulation especially
    for multi-D models.
  }
  \item{bandup }{number of non-zero bands above the diagonal, in case
    the Jacobian is banded.
  }
  \item{banddown }{number of non-zero bands below the diagonal, in case
    the Jacobian is banded.
  }
  \item{maxsteps }{average maximal number of steps per output interval
    taken by the solver. This argument is defined such as to ensure
    compatibility with the Livermore-solvers. RADAU only accepts the maximal
    number of steps for the entire integration, and this is calculated 
    as \code{length(times) * maxsteps}.
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions refered to in \code{func} and
    \code{jacfunc}. See vignette \code{"compiledCode"}
    from package \code{deSolve}.
  }
  \item{initfunc }{if not \code{NULL}, the name of the initialisation function
    (which initialises values of parameters), as provided in
    \file{dllname}. See vignette \code{"compiledCode"}
    from package \code{deSolve}.
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
    indeed the number of output variables calculed in the DLL - you have
    to perform this check in the code - See vignette \code{"compiledCode"}
    from package \code{deSolve}.
  }
  \item{outnames }{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{func}, present in the shared library.
    These names will be used to label the output matrix.
  }
  \item{forcings }{only used if \file{dllname} is specified: a list with
    the forcing function data sets, each present as a two-columned matrix,
    with (time, value); interpolation outside the interval
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
  \item{... }{additional arguments passed to \code{func} and
    \code{jacfunc} allowing this to be a generic function.
  }
}
\value{
  A matrix of class \code{deSolve} with up to as many rows as elements
  in \code{times} and as many columns as elements in \code{y} plus the number of "global"
  values returned in the next elements of the return from \code{func},
  plus and additional column for the time value.  There will be a row
  for each element in \code{times} unless the FORTRAN routine 
  returns with an unrecoverable error. If \code{y} has a names
  attribute, it will be used to label the columns of the output value.

}
\author{Karline Soetaert 

}
\examples{
## =======================================================================
## Example 1: ODE
##   Various ways to solve the same model.
## =======================================================================

## the model, 5 state variables
f1 <- function  (t, y, parms) {
  ydot <- vector(len = 5)

  ydot[1] <-  0.1*y[1] -0.2*y[2]
  ydot[2] <- -0.3*y[1] +0.1*y[2] -0.2*y[3]
  ydot[3] <-           -0.3*y[2] +0.1*y[3] -0.2*y[4]
  ydot[4] <-                     -0.3*y[3] +0.1*y[4] -0.2*y[5]
  ydot[5] <-                               -0.3*y[4] +0.1*y[5]

  return(list(ydot))
}

## the Jacobian, written as a full matrix
fulljac <- function  (t, y, parms) {
  jac <- matrix(nrow = 5, ncol = 5, byrow = TRUE,
                data = c(0.1, -0.2,  0  ,  0  ,  0  ,
                        -0.3,  0.1, -0.2,  0  ,  0  ,
                         0  , -0.3,  0.1, -0.2,  0  ,
                         0  ,  0  , -0.3,  0.1, -0.2,
                         0  ,  0  ,  0  , -0.3,  0.1))
  return(jac)
}

## the Jacobian, written in banded form
bandjac <- function  (t, y, parms) {
  jac <- matrix(nrow = 3, ncol = 5, byrow = TRUE,
                data = c( 0  , -0.2, -0.2, -0.2, -0.2,
                          0.1,  0.1,  0.1,  0.1,  0.1,
                         -0.3, -0.3, -0.3, -0.3,    0))
  return(jac)
}

## initial conditions and output times
yini  <- 1:5
times <- 1:20

## default: stiff method, internally generated, full Jacobian
out   <- radau(yini, times, f1, parms = 0)
plot(out)

## stiff method, user-generated full Jacobian
out2  <- radau(yini, times, f1, parms = 0, jactype = "fullusr",
               jacfunc = fulljac)

## stiff method, internally-generated banded Jacobian
## one nonzero band above (up) and below(down) the diagonal
out3  <- radau(yini, times, f1, parms = 0, jactype = "bandint",
               bandup = 1, banddown = 1)

## stiff method, user-generated banded Jacobian
out4  <- radau(yini, times, f1, parms = 0, jactype = "bandusr",
               jacfunc = bandjac, bandup = 1, banddown = 1)


## =======================================================================
## Example 2: ODE
##   stiff problem from chemical kinetics
## =======================================================================
Chemistry <- function (t, y, p) {
  dy1 <- -.04*y[1] + 1.e4*y[2]*y[3]
  dy2 <- .04*y[1] - 1.e4*y[2]*y[3] - 3.e7*y[2]^2
  dy3 <- 3.e7*y[2]^2
  list(c(dy1, dy2, dy3))
}

times <- 10^(seq(0, 10, by = 0.1))
yini <- c(y1 = 1.0, y2 = 0, y3 = 0)

out <- radau(func = Chemistry, times = times, y = yini, parms = NULL)
plot(out, log = "x", type = "l", lwd = 2)

## =============================================================================
## Example 3: DAE
## Car axis problem, index 3 DAE, 8 differential, 2 algebraic equations
## from
## F. Mazzia and C. Magherini. Test Set for Initial Value Problem Solvers,
## release 2.4. Department
## of Mathematics, University of Bari and INdAM, Research Unit of Bari,
## February 2008.
## Available at http://www.dm.uniba.it/~testset.
## =============================================================================

## Problem is written as M*y = f(t,y,p).
## caraxisfun implements the right-hand side:

caraxisfun <- function(t, y, parms) {
  with(as.list(y), {
  
    yb <- r * sin(w * t)
    xb <- sqrt(L * L - yb * yb)
    Ll <- sqrt(xl^2 + yl^2)
    Lr <- sqrt((xr - xb)^2 + (yr - yb)^2)
        
    dxl <- ul; dyl <- vl; dxr <- ur; dyr <- vr
        
    dul  <- (L0-Ll) * xl/Ll      + 2 * lam2 * (xl-xr) + lam1*xb
    dvl  <- (L0-Ll) * yl/Ll      + 2 * lam2 * (yl-yr) + lam1*yb - k * g
               
    dur  <- (L0-Lr) * (xr-xb)/Lr - 2 * lam2 * (xl-xr)
    dvr  <- (L0-Lr) * (yr-yb)/Lr - 2 * lam2 * (yl-yr) - k * g
        
    c1   <- xb * xl + yb * yl
    c2   <- (xl - xr)^2 + (yl - yr)^2 - L * L
        
    list(c(dxl, dyl, dxr, dyr, dul, dvl, dur, dvr, c1, c2))
  })
}

eps <- 0.01; M <- 10; k <- M * eps^2/2; 
L <- 1; L0 <- 0.5; r <- 0.1; w <- 10; g <- 1

yini <- c(xl = 0, yl = L0, xr = L, yr = L0,
          ul = -L0/L, vl = 0,
          ur = -L0/L, vr = 0,
          lam1 = 0, lam2 = 0)

# the mass matrix
Mass      <- diag(nrow = 10, 1)
Mass[5,5] <- Mass[6,6] <- Mass[7,7] <- Mass[8,8] <- M * eps * eps/2
Mass[9,9] <- Mass[10,10] <- 0
Mass

# index of the variables: 4 of index 1, 4 of index 2, 2 of index 3
index <- c(4, 4, 2)

times <- seq(0, 3, by = 0.01)
out <- radau(y = yini, mass = Mass, times = times, func = caraxisfun,
        parms = parameter, nind = index)

plot(out, which = 1:4, type = "l", lwd = 2)
}
\references{
  E. Hairer and G. Wanner, 1996. Solving Ordinary Differential Equations II.
  Stiff and Differential-algebraic problems. Springer series in computational
  mathematics 14, Springer-Verlag, second edition.
}
\details{
  The work is done by the FORTRAN subroutine \code{RADAU5}, whose
  documentation should be consulted for details. The implementation
  is based on the Fortran 77 version from January 18, 2002.
    
  There are four standard choices for the Jacobian which can be specified with
  \code{jactype}.

  The options for \bold{jactype} are
  
  \describe{
    \item{jactype = "fullint"}{a full Jacobian, calculated internally by
      the solver.
    }
    \item{jactype = "fullusr"}{a full Jacobian, specified by user
      function \code{jacfunc}.
    }
    \item{jactype = "bandusr"}{a banded Jacobian, specified by user
      function \code{jacfunc}; the size of the bands specified by
      \code{bandup} and \code{banddown}.
    }
    \item{jactype = "bandint"}{a banded Jacobian, calculated by radau;
      the size of the bands specified by \code{bandup} and
      \code{banddown}.
    }
  }

  Inspection of the example below shows how to specify both a banded and
  full Jacobian.
  
  The input parameters \code{rtol}, and \code{atol} determine the
  \bold{error control} performed by the solver, which roughly keeps the 
  local error of \eqn{y(i)} below \eqn{rtol(i)*abs(y(i))+atol(i)}.  
   
  The diagnostics of the integration can be printed to screen
  by calling \code{\link{diagnostics}}. If \code{verbose} = \code{TRUE},
  the diagnostics will be written to the screen at the end of the integration.

  See vignette("deSolve") from the \code{deSolve} package for an 
  explanation of each element in the vectors
  containing the diagnostic properties and how to directly access them.

  \bold{Models} may be defined in compiled C or FORTRAN code, as well as
  in an R-function. See package vignette \code{"compiledCode"} from package 
  \code{deSolve} for details.

  Information about linking forcing functions to compiled code is in 
  \link{forcings} (from package \code{deSolve}).

}
\seealso{
  \itemize{
    \item \code{\link{ode}} for a general interface to most of the ODE solvers ,
    \item \code{\link{ode.1D}} for integrating 1-D models,
    \item \code{\link{ode.2D}} for integrating 2-D models,
    \item \code{\link{ode.3D}} for integrating 3-D models,
    \item \code{\link{daspk}} for integrating DAE models up to index 1
  }

  \code{\link{diagnostics}} to print diagnostic messages.
}
\keyword{math}

