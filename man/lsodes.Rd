\name{lsodes}
\alias{lsodes}                               
\title{Solver for Ordinary Differential Equations (ODE) With
  Sparse Jacobian
}
\description{
  Solves the initial value problem for stiff systems of ordinary
  differential equations (ODE) in the form: \deqn{dy/dt = f(t,y)} and
  where the Jacobian matrix df/dy has an arbitrary sparse structure.
  
  The \R function \code{lsodes} provides an interface to the FORTRAN ODE
  solver of the same name, written by Alan C. Hindmarsh and Andrew
  H. Sherman.

  The system of ODE's is written as an \R function or be defined in
  compiled code that has been dynamically loaded.
}

\usage{
lsodes(y, times, func, parms, rtol = 1e-6, atol = 1e-6, 
  jacvec = NULL, sparsetype = "sparseint", nnz = NULL,
  inz = NULL,  rootfunc = NULL,
  verbose = FALSE, nroot = 0, tcrit = NULL, hmin = 0,
  hmax = NULL, hini = 0, ynames = TRUE, maxord = NULL,
  maxsteps = 5000, lrw = NULL, liw = NULL, dllname = NULL,
  initfunc = dllname, initpar = parms, rpar = NULL,
  ipar = NULL, nout = 0, outnames = NULL, forcings=NULL,
  initforc = NULL, fcontrol=NULL, events=NULL, lags = NULL, 
  ...)
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
    \code{t}, or a character string giving the name of a compiled
    function in a dynamically loaded shared library.
  
    If \code{func} is an \R-function, it must be defined as:
     \code{func <- function(t, y, parms,...)}.  \code{t} is the current
    time point in the integration, \code{y} is the current estimate of
    the variables in the ODE system.  If the initial values \code{y} has
    a \code{names} attribute, the names will be available inside \code{func}.
    \code{parms} is a vector or list of parameters; ... (optional) are
    any other arguments passed to the function.
  
    The return value of \code{func} should be a list, whose first
    element is a vector containing the derivatives of \code{y} with
    respect to \code{time}, and whose next elements are global values
    that are required at each point in \code{times}.  The derivatives
    must be specified in the \bold{same order} as the state variables \code{y}.

    If \code{func} is
    a string, then \code{dllname} must give the name of the shared
    library (without extension) which must be loaded before
    \code{lsodes()} is called. See package vignette \code{"compiledCode"}
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
  \item{jacvec }{if not \code{NULL}, an \R function that computes a
    column of the Jacobian of the system of differential equations
    \eqn{\partial\dot{y}_i/\partial y_j}{dydot(i)/dy(j)}, 
    or a string giving the name of a function or
    subroutine in \file{dllname} that computes the column of the
    Jacobian (see vignette \code{"compiledCode"} for more about this option).

    The \R
    calling sequence for \code{jacvec} is identical to that of
    \code{func}, but with extra parameter \code{j}, denoting the column
    number.  Thus, \code{jacvec} should be called as: \code{jacvec =
    func(t, y, j, parms)} and \code{jacvec} should return a vector
    containing column \code{j} of the Jacobian, i.e.  its i-th value is
    \eqn{\partial\dot{y}_i/\partial y_j}{dydot(i)/dy(j)}.
    If this function is absent, \code{lsodes} will
    generate the Jacobian by differences.
  }
  \item{sparsetype }{the sparsity structure of the Jacobian, one of
    "sparseint" or "sparseusr", "sparsejan", ..., 
    The sparsity can be estimated internally by lsodes (first option)
    or given by the user (last two). See details.
  }
  \item{nnz }{the number of nonzero elements in the sparse Jacobian (if
    this is unknown, use an estimate).
  }
  \item{inz }{if \code{sparsetype} equal to "sparseusr", a two-columned matrix
    with the (row, column) indices to the nonzero elements in the sparse
    Jacobian. If \code{sparsetype} = "sparsejan", a vecotr with the elements 
    ian followed by he elements jan as used in the lsodes code. See details.
    In all other cases, ignored.
  }
  \item{rootfunc }{if not \code{NULL}, an \R function that computes the
    function whose root has to be estimated or a string giving the name
    of a function or subroutine in \file{dllname} that computes the root
    function.  The \R calling sequence for \code{rootfunc} is identical
    to that of \code{func}.  \code{rootfunc} should return a vector with
    the function values whose root is sought.
  }
  \item{verbose }{if \code{TRUE}: full output to the screen, e.g. will
    print the \code{diagnostiscs} of the integration - see details.
  }
  \item{nroot }{only used if \file{dllname}  is specified: the number of
    constraint functions whose roots are desired during the integration;
    if \code{rootfunc} is an R-function, the solver estimates the number
    of roots.
  }
  \item{tcrit }{if not \code{NULL}, then \code{lsodes} cannot integrate
    past \code{tcrit}. The FORTRAN routine \code{lsodes} overshoots its
    targets (times points in the vector \code{times}), and interpolates
    values for the desired time points.  If there is a time beyond which
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
    size is determined by the solver.
  }
  \item{ynames }{logical, if \code{FALSE} names of state variables are
    not passed to function \code{func}; this may speed up the simulation
    especially for multi-D models.
  }
  \item{maxord }{the maximum order to be allowed. \code{NULL} uses the
    default, i.e. order 12 if implicit Adams method (meth = 1), order 5
    if BDF method (meth = 2). Reduce maxord to save storage space.
  }
  \item{maxsteps }{maximal number of steps per output interval taken by the
    solver.
  }
  \item{lrw }{the length of the real work array rwork; due to the
    sparsicity, this cannot be readily predicted. If \code{NULL}, a
    guess will be made, and if not sufficient, \code{lsodes} will return
    with a message indicating the size of rwork actually required.
    Therefore, some experimentation may be necessary to estimate the
    value of \code{lrw}.
    
    For instance, if you get the error:
\preformatted{ 
DLSODES- RWORK length is insufficient to proceed.                               
  Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)                    
  In above message,  I1 =     27627   I2 =     25932 
}    
  set \code{lrw} equal to 27627 or a higher value

  }
  \item{liw }{the length of the integer work array iwork; due to the
    sparsicity, this cannot be readily predicted. If \code{NULL}, a guess will
    be made, and if not sufficient, \code{lsodes} will return with a
    message indicating the size of iwork actually required.  Therefore,
    some experimentation may be necessary to estimate the value of
    \code{liw}.
  }
  \item{dllname }{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions refered to in \code{func} and
    \code{jacfunc}. See package vignette \code{"compiledCode"}.
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
    indeed the number of output variables calculated in the dll - you have
    to perform this check in the code. See package vignette
    \code{"compiledCode"}.
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
    See \link{forcings} or vignette \code{compiledCode}.
  }
  \item{events }{A list that specifies events, i.e. when the value of a 
   state variable is suddenly changed. See \link{events} for more information.
  }
  \item{lags }{A list that specifies timelags, i.e. the number of steps 
   that has to be kept. To be used for delay differential equations. 
   See \link{timelags}, \link{dede} for more information.
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
  for each element in \code{times} unless the FORTRAN routine `lsodes'
  returns with an unrecoverable error. If \code{y} has a names
  attribute, it will be used to label the columns of the output value.
}
\author{Karline Soetaert <karline.soetaert@nioz.nl>}
\examples{
## Various ways to solve the same model.

## =======================================================================
## The example from lsodes source code
## A chemical model
## =======================================================================

n  <- 12
y  <- rep(1, n)
dy <- rep(0, n)

times <- c(0, 0.1*(10^(0:4)))

rtol <- 1.0e-4
atol <- 1.0e-6

parms <- c(rk1  = 0.1,   rk2 = 10.0, rk3 = 50.0,  rk4 = 2.5,  rk5 = 0.1,
           rk6  = 10.0,  rk7 = 50.0, rk8 = 2.5,   rk9 = 50.0, rk10 = 5.0,
           rk11 = 50.0, rk12 = 50.0,rk13 = 50.0, rk14 = 30.0,
           rk15 = 100.0,rk16 = 2.5, rk17 = 100.0,rk18 = 2.5,
           rk19 = 50.0, rk20 = 50.0)

#
chemistry <- function (time, Y, pars) {
  with (as.list(pars), {
    dy[1] <- -rk1 *Y[1]
    dy[2] <-  rk1 *Y[1]        + rk11*rk14*Y[4]  + rk19*rk14*Y[5]  -
              rk3 *Y[2]*Y[3]   - rk15*Y[2]*Y[12] - rk2*Y[2]
    dy[3] <-  rk2 *Y[2]        - rk5 *Y[3]       - rk3*Y[2]*Y[3]   -
              rk7*Y[10]*Y[3]   + rk11*rk14*Y[4]   + rk12*rk14*Y[6]
    dy[4] <-  rk3 *Y[2]*Y[3]   - rk11*rk14*Y[4]  - rk4*Y[4]
    dy[5] <-  rk15*Y[2]*Y[12]  - rk19*rk14*Y[5]  - rk16*Y[5]
    dy[6] <-  rk7 *Y[10]*Y[3]  - rk12*rk14*Y[6]  - rk8*Y[6]
    dy[7] <-  rk17*Y[10]*Y[12] - rk20*rk14*Y[7]  - rk18*Y[7]
    dy[8] <-  rk9 *Y[10]       - rk13*rk14*Y[8]  - rk10*Y[8]
    dy[9] <-  rk4 *Y[4]        + rk16*Y[5]       + rk8*Y[6]         +
              rk18*Y[7]
    dy[10] <- rk5 *Y[3]        + rk12*rk14*Y[6]  + rk20*rk14*Y[7]   +
              rk13*rk14*Y[8]   - rk7 *Y[10]*Y[3] - rk17*Y[10]*Y[12] -
              rk6 *Y[10]       - rk9*Y[10]
    dy[11] <- rk10*Y[8]
    dy[12] <- rk6 *Y[10]       + rk19*rk14*Y[5]  + rk20*rk14*Y[7]   -
              rk15*Y[2]*Y[12]  - rk17*Y[10]*Y[12]
    return(list(dy))
  })
}

## =======================================================================
## application 1. lsodes estimates the structure of the Jacobian
##                and calculates the Jacobian by differences
## =======================================================================
out <- lsodes(func = chemistry, y = y, parms = parms, times = times,
              atol = atol, rtol = rtol, verbose = TRUE)

## =======================================================================
## application 2. the structure of the Jacobian is input
##                lsodes calculates the Jacobian by differences
##                this is not so efficient...
## =======================================================================

## elements of Jacobian that are not zero
nonzero <-  matrix(nc = 2, byrow = TRUE, data = c(
  1, 1,   2, 1,    # influence of sp1 on rate of change of others
  2, 2,   3, 2,   4, 2,   5, 2,  12, 2,
  2, 3,   3, 3,   4, 3,   6, 3,  10, 3,
  2, 4,   3, 4,   4, 4,   9, 4,  # d (dyi)/dy4
  2, 5,   5, 5,   9, 5,  12, 5,
  3, 6,   6, 6,   9, 6,  10, 6,
  7, 7,   9, 7,  10, 7,  12, 7,
  8, 8,  10, 8,  11, 8,
  3,10,   6,10,   7,10,  10,10,  12,10,
  2,12,   5,12,   7,12,  10,12,  12,12)
)

## when run, the default length of rwork is too small
## lsodes will tell the length actually needed
# out2 <- lsodes(func = chemistry, y = y, parms = parms, times = times,
#              inz = nonzero, atol = atol,rtol = rtol)  #gives warning
out2 <- lsodes(func = chemistry, y = y, parms = parms, times = times, 
            sparsetype = "sparseusr", inz = nonzero,   
             atol = atol, rtol = rtol, verbose = TRUE, lrw = 351)
                            
## =======================================================================
## application 3. lsodes estimates the structure of the Jacobian
##                the Jacobian (vector) function is input
## =======================================================================
chemjac <- function (time, Y, j, pars) {
  with (as.list(pars), {
    PDJ <- rep(0,n)

    if (j == 1){
       PDJ[1] <- -rk1
       PDJ[2] <- rk1
    } else if (j == 2) {
       PDJ[2] <- -rk3*Y[3] - rk15*Y[12] - rk2
       PDJ[3] <- rk2 - rk3*Y[3]
       PDJ[4] <- rk3*Y[3]
       PDJ[5] <- rk15*Y[12]
       PDJ[12] <- -rk15*Y[12]
    } else if (j == 3) {
       PDJ[2] <- -rk3*Y[2]
       PDJ[3] <- -rk5 - rk3*Y[2] - rk7*Y[10]
       PDJ[4] <- rk3*Y[2]
       PDJ[6] <- rk7*Y[10]
       PDJ[10] <- rk5 - rk7*Y[10]
    } else if (j == 4) {
       PDJ[2] <- rk11*rk14
       PDJ[3] <- rk11*rk14
       PDJ[4] <- -rk11*rk14 - rk4
       PDJ[9] <- rk4
    } else if (j == 5) {
       PDJ[2] <- rk19*rk14
       PDJ[5] <- -rk19*rk14 - rk16
       PDJ[9] <- rk16
       PDJ[12] <- rk19*rk14
    } else if (j == 6) {
       PDJ[3] <- rk12*rk14
       PDJ[6] <- -rk12*rk14 - rk8
       PDJ[9] <- rk8
       PDJ[10] <- rk12*rk14
    } else if (j == 7) {
       PDJ[7] <- -rk20*rk14 - rk18
       PDJ[9] <- rk18
       PDJ[10] <- rk20*rk14
       PDJ[12] <- rk20*rk14
    } else if (j == 8) {
       PDJ[8] <- -rk13*rk14 - rk10
       PDJ[10] <- rk13*rk14
       PDJ[11] <- rk10
    } else if (j == 10) {
       PDJ[3] <- -rk7*Y[3]
       PDJ[6] <- rk7*Y[3]
       PDJ[7] <- rk17*Y[12]
       PDJ[8] <- rk9
       PDJ[10] <- -rk7*Y[3] - rk17*Y[12] - rk6 - rk9
       PDJ[12] <- rk6 - rk17*Y[12]
    } else if (j == 12) {
       PDJ[2] <- -rk15*Y[2]
       PDJ[5] <- rk15*Y[2]
       PDJ[7] <- rk17*Y[10]
       PDJ[10] <- -rk17*Y[10]
       PDJ[12] <- -rk15*Y[2] - rk17*Y[10]
    }
    return(PDJ)
  })
} 

out3 <- lsodes(func = chemistry, y = y, parms = parms, times = times,
              jacvec = chemjac, atol = atol, rtol = rtol)

## =======================================================================
## application 4. The structure of the Jacobian (nonzero elements) AND
##                the Jacobian (vector) function is input
##                not very efficient...
## =======================================================================
out4 <- lsodes(func = chemistry, y = y, parms = parms, times = times,
               lrw = 351, sparsetype = "sparseusr", inz = nonzero,
               jacvec = chemjac, atol = atol, rtol = rtol,
               verbose = TRUE)
}
\references{
  Alan C. Hindmarsh, ODEPACK, A Systematized Collection of ODE Solvers,
  in Scientific Computing, R. S. Stepleman et al. (Eds.), North-Holland,
  Amsterdam, 1983, pp. 55-64.
  
  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman, Yale
  Sparse Matrix Package: I. The Symmetric Codes,
  Int. J. Num. Meth. Eng., 18 (1982), pp. 1145-1151.
  
  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman, Yale
  Sparse Matrix Package: II. The Nonsymmetric Codes, Research Report
  No. 114, Dept. of Computer Sciences, Yale University, 1977.
}
\details{

  The work is done by the FORTRAN subroutine \code{lsodes}, whose
  documentation should be consulted for details (it is included as
  comments in the source file \file{src/opkdmain.f}). The implementation
  is based on the November, 2003 version of lsodes, from Netlib.
    
  \code{lsodes} is applied for stiff problems, where the Jacobian has a
  sparse structure.
  
  There are several choices depending on whether \code{jacvec}
  is specified and depending on the setting of \code{sparsetype}.
  
  If function \code{jacvec} is present, then it should return the j-th
  column of the Jacobian matrix.
  
  There are also several choices for the sparsity specification, selected by
  argument \code{sparsetype}.
  \itemize{
    \item \code{sparsetype} = \code{"sparseint"}.  The sparsity is estimated
      by the solver, based on numerical differences.
      In this case, it is advisable to provide an estimate of the number
      of non-zero elements in the Jacobian (\code{nnz}).
      This value can be approximate; upon return the number of nonzero
      elements actually required will be known (1st element of attribute
      \code{dims}).
      In this case, \code{inz} need not be specified. 
    \item \code{sparsetype} = \code{"sparseusr"}. The sparsity is determined by
      the user. In this case, \code{inz} should be a \code{matrix}, containing indices 
      (row, column) to the  nonzero elements in the Jacobian matrix.
      The number of nonzeros \code{nnz} will be set equal to the number of rows
      in \code{inz}.
    \item \code{sparsetype} = \code{"sparsejan"}. The sparsity is also determined by
      the user.
      In this case, \code{inz} should be a \code{vector}, containting the \code{ian} and 
      \code{jan} elements of the sparse storage format, as used in the sparse solver.
      Elements of \code{ian} should be the first \code{n+1} elements of this vector, and 
      contain the starting locations in \code{jan} of columns 1.. n. 
      \code{jan} contains the row indices of the nonzero locations of           
      the Jacobian, reading in columnwise order.
      The number of nonzeros \code{nnz} will be set equal to the length of \code{inz} - (n+1).
    \item \code{sparsetype} = \code{"1D"}, \code{"2D"}, \code{"3D"}. 
      The sparsity is estimated by the solver, based on numerical differences.
      Assumes finite differences in a 1D, 2D or 3D regular grid - used by 
      functions \code{ode.1D}, \code{ode.2D}, \code{ode.3D}.
      Similar are \code{"2Dmap"}, and \code{"3Dmap"}, which also include a 
      mapping variable (passed in nnz). 
  }
  
  The input parameters \code{rtol}, and \code{atol} determine the
  \bold{error control} performed by the solver.  See \code{\link{lsoda}}
  for details.
   
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

  Examples in both C and FORTRAN are in the \file{doc/examples/dynload} subdirectory
  of the \code{deSolve} package directory.

  \code{lsodes} can find the root of at least one of a set of constraint functions
  \code{rootfunc} of the independent and dependent variables.  It then returns the
  solution at the root if that occurs sooner than the specified stop
  condition, and otherwise returns the solution according the specified
  stop condition.

  Caution:  Because of numerical errors in the function
  \code{rootfun} due to roundoff and integration error, \code{lsodes} may
  return false roots, or return the same root at two or more
  nearly equal values of \code{time}.

}
\seealso{
  \itemize{
    \item \code{\link{rk}}, \item \code{\link{rk4}} and \code{\link{euler}} for
      Runge-Kutta integrators.
    \item \code{\link{lsoda}}, \code{\link{lsode}},
      \code{\link{lsodar}}, \code{\link{vode}},
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

