\name{checkDLL}
\alias{checkDLL}

\title{
Check shared library (DLL/.so) of a compiled model.
}
\description{
Check shared library (DLL/.so) of a compiled model and create a list of symbols.
}
\usage{
checkDLL(func, jacfunc, dllname, initfunc, verbose, nout, outnames, JT = 1)
}

\arguments{
  \item{func}{
 character: name of the derivative function.
}
  \item{jacfunc}{an \R function, that computes the
    Jacobian of the system of differential equations
    \eqn{\partial\dot{y}_i/\partial y_j}{dydot(i)/dy(j)}, or
    a string giving the name of a function or subroutine in
    \file{dllname} that computes the Jacobian.}

  \item{dllname}{a string giving the name of the shared library
    (without extension) that contains all the compiled function or
    subroutine definitions refered to in \code{func} and
    \code{jacfunc}.}

  \item{initfunc}{ the name of the initialisation function
    (which initialises values of parameters), as provided in
    \file{dllname}. See package vignette \code{"compiledCode"}.
}
  \item{verbose}{reserved for future extensions.
}
  \item{nout}{only used if \code{dllname} is specified and the model is
    defined in compiled code: the number of output variables calculated
    in the compiled function \code{func}, present in the shared
    library. Note: it is not automatically checked whether this is
    indeed the number of output variables calculated in the dll - you have
    to perform this check in the code.
}
  \item{outnames}{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{func}, present in the shared library.
    These names will be used to label the output matrix.
}
  \item{JT}{integer specifying the type of the Jacobian. The default value of 1
    must be set to 2 for solver \code{lsodes} if `func` is specified in a  DLL
    or inline compiled and if  a `jacfunc` is provided.
}
}
\details{
  The function \code{checkDLL} is normally called internally by the solver
  functions. It can be used to avoid overhead, when a small compiled
  model with a low number of integration steps is repeatedly called.
  The feature is currently only available for the \code{lsoda} solver.
}

\value{
List of class \code{deSolve.symbols} with:

\item{ModelInit }{pointer to the init function of the DLL (class "externalptr").}
\item{Func }{pointer to the derivative function in the DLL (class "externalptr").}
\item{JacFunc }{pointer to the Jacobi function in the DLL (class "externalptr").}
\item{Nglobal }{number of output variables calculated in the compiled function.}
\item{Nmtot }{list of names of derivatives and output variables.}
}
\seealso{
\code{\link{lsoda}}
}
\examples{
\dontrun{
symbols <- checkDLL(func = "derivs", jacfunc = NULL, dllname = "lorenzc",
             initfunc = "initmod", verbose = TRUE, nout = 0,
             outnames = NULL, JT = 1)
}

}
\keyword{math}
