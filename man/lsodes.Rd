\name{lsodes}
\alias{lsodes}
\title{General solver for ordinary differential equations (ODE) with sparse Jacobian}
\description{Solves the initial value problem for stiff systems of 
    ordinary differential equations (ODE) in the form:
    \deqn{dy/dt = f(t,y)} and where the Jacobian matrix df/dy has an arbitrary sparse structure.
  
    The \R function \code{lsodes} provides an interface to the
    Fortran ODE solver of the same name, written by Alan C. Hindmarsh and Andrew H. Sherman
   \cr The system of ODE's is written as an \R function or be defined in
    compiled code that has been dynamically loaded.}
\usage{lsodes(y, times, func, parms, rtol=1e-6, atol=1e-6, 
  jacvec=NULL, sparsetype = "sparseint", nnz=NULL, inz=NULL,     
  verbose=FALSE, tcrit=NULL, hmin=0, hmax=NULL, hini=0, ynames=TRUE, 
  maxord=NULL, maxsteps=5000, lrw=NULL, liw=NULL,  
 dllname=NULL, initfunc=dllname, initpar=parms,  
  rpar=NULL, ipar=NULL, nout=0, outnames=NULL, ...)}
\arguments{
  \item{y }{the initial (state) values for the ODE system. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{times }{time sequence for which output is wanted; the first value of \code{times} must be the initial time; if only one step is to be taken; set \code{times} = NULL}
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
      before \code{lsodes()} is called. See package vignette for more details.}
  \item{parms }{vector or list of parameters used in \code{func} or \code{jacfunc}.}
  \item{rtol }{relative error tolerance, either a scalar or an array as
      long as \code{y}.  See details. }
  \item{atol }{absolute error tolerance, either a scalar or an array as
      long as \code{y}.  See details.}
  \item{jacvec }{if not \code{NULL}, an \R function that computes
      a column of the jacobian of the system of differential equations
      dydot(i)/dy(j), or a string giving the name of a function or
      subroutine in \file{dllname} that computes the column of the jacobian (see Details
      below for more about this option).  The \R calling sequence for
      \code{jacvec} is identical to that of \code{func}, but with extra parameter \code{j}, denoting the column number.  
      Thus, \code{jacvec} should be called as: \code{jacvec = func(t, y, j, parms)}
      and \code{jacvec} should return a vector containing column \code{j} of the jacobian, i.e. 
      its i-th value is dydot(i)/dy(j). If this function is absent, \code{lsodes} will generate the 
      jacobian by differences}
  \item{sparsetype }{the sparsity structure of the Jacobian, one of "sparseint" or "sparseusr", sparsity estimated internally by lsodes or given by user}      
  \item{nnz }{the number of nonzero elements in the sparse Jacobian (if this is unknown, use an estimate)}
  \item{inz }{(row,column) indices to the nonzero elements in the sparse Jacobian. Necessary if \code{sparsetype} = "sparseusr"; else ignored }
  \item{verbose }{if TRUE: full output to the screen, e.g. will output the settings of vectors *istate* and *rstate* - see details}
  \item{tcrit }{if not \code{NULL}, then \code{lsodes} cannot integrate past \code{tcrit}. The Fortran routine \code{lsodes} overshoots its targets (times points in the vector \code{times}), and interpolates values
      for the desired time points.  If there is a time beyond which integration should not proceed (perhaps because of a singularity),
      that should be provided in \code{tcrit}.}
  \item{hmin }{an optional minimum value of the integration
      stepsize. In special situations this parameter may speed up computations with
      the cost of precision. Don't use \code{hmin} if you don't know why!}
  \item{hmax }{an optional maximum value of the integration stepsize. If not specified, \code{hmax} is set to the largest difference in \code{times}, to avoid that the simulation possibly ignores short-term events. If 0, no maximal size is specified}
  \item{hini }{initial step size to be attempted; if 0, the initial step size is determined by the solver}
  \item{ynames }{if FALSE: names of state variables are not passed to function \code{func} ; this may speed up the simulation especially for multi-D models}
  \item{maxord }{the maximum order to be allowed. NULL uses the default, i.e. order 12 if implicit Adams method (meth=1), order 5 if BDF method (meth=2). Reduce maxord to save storage space}
  \item{maxsteps }{maximal number of steps during one call to the solver}
  \item{lrw }{the length of the real work array rwork; due to the sparsicity, this cannot be readily predicted. If NULL, a guess will be made, and 
  		if not sufficient, \code{lsodes} will return with a message indicating the size of rwork actually required. 
      Therefore, some experimentation may be necessary to estimate the value of \code{lrw} }
  \item{liw }{the length of the integer work array iwork; due to the sparsicity, this cannot be readily predicted. If NULL, a guess will be made, and 
  		if not sufficient, \code{lsodes} will return with a message indicating the size of iwork actually required. 
      Therefore, some experimentation may be necessary to estimate the value of \code{liw} }
  \item{dllname }{a string giving the name of the shared library (without
      extension) that contains all the compiled function or subroutine
      definitions refered to in \code{func} and \code{jacfunc}. See package vignette.}
  \item{initfunc }{if not NULL, the name of the initialisation function (which initialises values of parameters), as provided in \file{dllname}. See package vignette.}
  \item{initpar }{only when \file{dllname} is specified and an initialisation function \code{initfunc} is in the dll: the parameters passed to the initialiser, to initialise the common blocks (fortran) or global variables (C, C++)}
  \item{rpar }{only when \file{dllname} is specified: a vector with double precision values passed to the dll-functions whose names are specified by \code{func} and \code{jacfunc}}
  \item{ipar }{only when \file{dllname} is specified: a vector with integer values passed to the dll-functions whose names are specified by \code{func} and \code{jacfunc}}
  \item{nout }{only used if \code{dllname} is specified and the model is defined in compiled code: the number of output variables calculated in the compiled function \code{func}, present in the shared library. Note:
      it is not automatically checked whether this is indeed the number of output variables calculed in the dll - you have to perform this check in the code. See package vignette.}
  \item{outnames }{only used if \file{dllname} is specified and \code{nout} > 0: the names of output variables calculated in the compiled function \code{func}, present in the shared library}
  \item{... }{additional arguments passed to \code{func} and \code{jacfunc} allowing this to be a generic function}
}
\value{
  A matrix with up to as many rows as elements in times and as many columns as elements in \code{y} plus the number of "global" values returned 
      in the next elements of the return from \code{func}, plus an additional column (the first) for the time value. 
      There will be one row for each element in \code{times} unless the Fortran routine `lsodes' returns with an unrecoverable error. 
      If \code{y} has a names attribute, it will be used to label the columns of the output value. 
      
      The output will have the attributes \code{istate}, and \code{rstate}, two vectors with several useful elements. 
      See details.
      The first element of istate returns the conditions under which the last call to lsoda returned. Normal is \code{istate[1] = 2}.
      If \code{verbose} = TRUE, the settings of istate and rstate will be written to the screen

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# Various ways to solve the same model.
  
  ###############################
  ## The example from lsodes code
  ## A chemical model
  ###############################
  
  n  <- 12
  y  <- rep(1,n)
  dy <- rep(0,n)
  
  times <- c(0,0.1*(10^(0:4)))
  
  rtol = 1.0e-4
  atol = 1.0e-6
  
  parms <- c(rk1=0.1,   rk2=10.0, rk3=50.0,  rk4=2.5,  rk5=0.1,
             rk6=10.0,  rk7=50.0, rk8=2.5,   rk9=50.0, rk10=5.0,
             rk11=50.0, rk12=50.0,rk13=50.0, rk14=30.0,
             rk15=100.0,rk16=2.5, rk17=100.0,rk18=2.5,
             rk19=50.0, rk20=50.0)
  
  #
  chemistry <- function (time,Y,pars)
  {
  with (as.list(pars),{
  
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
   dy[10]<-  rk5 *Y[3]        + rk12*rk14*Y[6]  + rk20*rk14*Y[7]   +
             rk13*rk14*Y[8]   - rk7 *Y[10]*Y[3] - rk17*Y[10]*Y[12] -
             rk6 *Y[10]       - rk9*Y[10]
   dy[11]<-  rk10*Y[8]
   dy[12]<-  rk6 *Y[10]       + rk19*rk14*Y[5]  + rk20*rk14*Y[7]   -
             rk15*Y[2]*Y[12]  - rk17*Y[10]*Y[12]
   return(list(dy))
  })
  }
  
  #--------------
  # application 1. lsodes estimates the structure of the jacobian 
  #                and calculates the jacobian by differences    
  out <- lsodes(func= chemistry, y = y, parms=parms, times=times,
                atol=atol,rtol=rtol,verbose=TRUE)
  
  #--------------
  # application 2. the structure of the jacobian is input
  #                lsodes calculates the jacobian by differences    
  # this is not so efficient... 
  
  # elements of Jacobian that are not zero
  nonzero <-  matrix(nc=2,byrow=TRUE,data=c(
   1, 1,   2, 1,    #influence of sp1 on rate of change of others
   2, 2,   3, 2,   4, 2,   5, 2,  12, 2,
   2, 3,   3, 3,   4, 3,   6, 3,  10, 3,
   2, 4,   3, 4,   4, 4,   9, 4,  #d (dyi)/dy4
   2, 5,   5, 5,   9, 5,  12, 5,
   3, 6,   6, 6,   9, 6,  10, 6,
   7, 7,   9, 7,  10, 7,  12, 7,
   8, 8,  10, 8,  11, 8,
   3,10,   6,10,   7,10,  10,10,  12,10,
   2,12,   5,12,   7,12,  10,12,  12,12))
  
  # when run, the default length of rwork is too small
  # lsodes will tell the length actually needed
  #out2<- lsodes(func= chemistry, y = y, parms=parms, times=times,
  #             inz=nonzero, atol=atol,rtol=rtol)  #gives warning
  out2<- lsodes(func= chemistry, y = y, parms=parms, times=times, 
              sparsetype="sparseusr", inz=nonzero,   
               atol=atol,rtol=rtol,verbose=TRUE,lrw=351)
                              
  #--------------
  # application 3. lsodes estimates the structure of the jacobian 
  #                the jacobian (vector) function is input
  #
  chemjac <- function (time,Y,j,pars)
  {
   with (as.list(pars),{
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
  
  out3<- lsodes(func= chemistry, y = y, parms=parms, times=times,
                jacvec=chemjac, atol=atol,rtol=rtol)          
  
  #--------------
  # application 4. The structure of the jacobian (nonzero elements) AND
  #                the jacobian (vector) function is input
  # not very efficient...
  
  out4<- lsodes(func= chemistry, y = y, parms=parms, times=times,lrw=351,  
                sparsetype="sparseusr", inz=nonzero, jacvec=chemjac,
                atol=atol, rtol=rtol, verbose=TRUE)}
\references{\item  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
  Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
  North-Holland, Amsterdam, 1983, pp. 55-64.
  
  \item S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
   Yale Sparse Matrix Package: I. The Symmetric Codes,
   Int. J. Num. Meth. Eng., 18 (1982), pp. 1145-1151.
  
  \item S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
   Yale Sparse Matrix Package: II. The Nonsymmetric Codes,
  Research Report No. 114, Dept. of Computer Sciences, Yale
  University, 1977.}
\details{The work is done by the Fortran subroutine \code{lsodes},
    whose documentation should be consulted for details (it is included as
    comments in the source file \file{src/opkdmain.f}).  The implementation is based on the
    November, 2003 version of lsodes, from Netlib. 
    
    \code{lsodes} is applied for stiff problems, where the Jacobian has a sparse structure. 
  \cr There are four choices depending on whether \code{jacvec} and \code{inz} is specified.
  
  If function \code{jacvec} is present, then it should return the j-th column of the Jacobian matrix,
  If matrix \code{inz} is present, then it should contain indices (row, column) to the nonzero elements in the Jacobian matrix.
  
  If  \code{jacvec} and \code{inz} are present, then the jacobian is fully specified by the user
  \cr If \code{jacvec} is present, but not \code{nnz} then the structure of the jacobian will be obtained from NEQ+1 calls to \code{jacvec}
  \cr If \code{nnz} is present, but not \code{jacvec} then the jacobian will be estimated internally, by differences. 
  \cr If neither \code{nnz} nor \code{jacvec} is present, then the jacobian will be generated internally by differences, its structure (indices to nonzero elements) will be obtained from NEQ+1 initial calls to \code{func}
  
  If \code{nnz} is not specified, it is advisable to provide an estimate of the number of non-zero elements in the Jacobian (\code{inz})
  
  
  The input parameters \code{rtol}, and \code{atol} determine the \bold{error
  control} performed by the solver. 
  See \code{\link{lsoda}} for details.
   
  \bold{Models} may be defined in compiled C or Fortran code, as well as in an R-function. See package vignette for details.
  
  
  The output will have the \bold{attributes *istate*, and *rstate*}, two vectors with several useful elements. 
      \cr if \code{verbose} = TRUE, the settings of istate and rstate will be written to the screen. 
      
      the following elements of \bold{istate} are meaningful:    
      \item el 1 : returns the conditions under which the last call to lsodes returned. 
      2  if lsodes was successful, -1 if excess work done, -2 means excess accuracy requested. (Tolerances too small), 
  -3 means illegal input detected. (See printed message.), -4 means repeated error test failures. (Check all input),
  -5 means repeated convergence failures. (Perhaps bad Jacobian supplied or wrong choice of MF or tolerances.),
  -6 means error weight became zero during problem. (Solution component i vanished, and atol or atol(i) = 0.)
      \item el 12 : The number of steps taken for the problem so far.
      \item el 13 : The number of function evaluations for the problem so far, excluding those for structure determination.",
      \item el 14 : The number of Jacobian evaluations  and LU decompositions so far, excluding those for structure determination.",
      \item el 15 : The method order last used (successfully).,
      \item el 16 : The order to be attempted on the next step.,
      \item el 17 : if el 1 =-4,-5: the largest component in the error vector,
      \item el 18 : The length of rwork actually required.,
      \item el 19 : The length of iwork actually required.,
      \item el 20 : The number of nonzero elements in the sparse jacobian,
  
      \bold{rstate} contains the following:
      \item 1: The step size in t last used (successfully).
      \item 2: The step size to be attempted on the next step.
      \item 3: The current value of the independent variable which the solver has actually reached, i.e. the
  current internal mesh point in t.  
      \item 4: A tolerance scale factor, greater than 1.0, computed when a request for too much accuracy was detected.
      
      For more information, see the comments in the original code lsodes.f}
\seealso{\item \code{\link{ode}}, \code{\link{lsoda}}, \code{\link{lsode}}, \code{\link{lsodar}}, \code{\link{vode}}, \code{\link{daspk}}, \code{\link{rk}}.}
\keyword{math}

