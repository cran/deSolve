\name{deSolve-package}
\alias{deSolve-package}
\alias{deSolve}
\docType{package}
\title{
  General Solvers for Initial Value Problems of Ordinary Differential Equations
  (ODE), Partial Differential Equations (PDE), Differential Algebraic
  Equations (DAE) and delay differential equations (DDE).
}
\description{
  Functions that solve initial value problems of a system of
  first-order ordinary differential equations (ODE), of partial
  differential equations (PDE), of differential algebraic equations (DAE)
  and delay differential equations.

  The functions provide an interface to the FORTRAN functions lsoda,
  lsodar, lsode, lsodes of the ODEPACK collection,
  to the FORTRAN functions dvode and daspk and a C-implementation of
  solvers of the Runge-Kutta family with fixed or variable time steps.

  The package contains routines designed for solving ODEs
  resulting from 1-D, 2-D and 3-D
  partial differential equations (PDE) that have been converted to ODEs by
  numerical differencing.
  
  It includes root-finding (or event location) and provides access to lagged
  variables and derivatives.

}
\details{
  \tabular{ll}{
    Package: \tab deSolve\cr
    Type:    \tab Package\cr
    Version: \tab 1.7\cr
    Date:    \tab 2010-02-16\cr
    License: \tab GNU Public License 2 or above\cr
  }

  The system of differential equations is written as an \R function or defined in
  compiled code that has been dynamically loaded, see package vignette
  (\code{vignette(compiledCode)}) for details. The solvers may be used
  as part of a modeling package for differential equations, or for
  parameter estimation using any appropriate modeling tool for
  non-linear models in \R such as \code{\link{optim}},
  \code{\link{nls}}, \code{\link{nlm}} or \code{\link[nlme]{nlme}}.
  
}
\author{
  Karline Soetaert,
  
  Thomas Petzoldt,

  R. Woodrow Setzer (Maintainer)
}
\references{
Karline Soetaert, Thomas Petzoldt, R. Woodrow Setzer (2010).
  Solving Differential Equations in R: Package deSolve
  Journal of Statistical Software, 33(9), 1--25.
  \url{http://www.jstatsoft.org/v33/i09/}

Alan C. Hindmarsh, 1983. ODEPACK, A Systematized Collection of ODE Solvers, in
Scientific Computing, R. S. Stepleman et al. (Eds.), North-Holland, Amsterdam,
pp. 55-64.

L. R. Petzold, 1983.
A Description of DASSL: A Differential/Algebraic System Solver,
in Scientific Computing,
R. S. Stepleman et al. (Eds.), North-Holland, Amsterdam, pp. 65-68.

P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, 1989.
VODE: A Variable Coefficient ODE Solver,
SIAM J. Sci. Stat. Comput., 10, pp. 1038-1051.
}
\seealso{
  \code{\link{ode}} for a general interface to most of the ODE solvers
  
  \code{\link{ode.band}} for solving models with a banded Jacobian

  \code{\link{ode.1D}}, \code{\link{ode.2D}}, \code{\link{ode.3D}},
  for integrating 1-D, 2-D and 3-D models
  
  \code{\link{dede}} for a general interface to the delay differential 
  equation solvers

  \code{\link{lsoda}}, \code{\link{lsode}},
  \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}},
  for ODE solvers of the Livermore family .

  \code{\link{daspk}}, for a DAE solver of the Livermore family.

  \code{\link{rk}}, \code{\link{rk4}}, \code{\link{euler}} for Runge-Kutta
  solvers.

  \code{\link{DLLfunc}}, \code{\link{DLLres}}, for testing model implementations
  in compiled code.

  \code{\link{forcings}},  \code{\link{events}}, for how to implement forcing 
  functions (external variables) and events (sudden changes in state variables)
  
  \code{\link{lagvalue}}, \code{\link{lagderiv}}, for how to get access to 
  lagged values of state variables and derivatives 
  
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(aquaphy)
example(lsoda)
example(ode.band)
example(ode.1D)
example(ode.2D)
example(dede)

## run demos
demo("odedim")     # partial differential equations
demo("CCL4model")  # a model fitting example (this will take some time)

## open the directory with source code of demos
browseURL(paste(system.file(package = "deSolve"), "/demo", sep = ""))

## open the directory with R sourcecode examples
browseURL(paste(system.file(package = "deSolve"), "/doc/examples", sep = ""))
## open the directory with C and FORTRAN sourcecode examples
browseURL(paste(system.file(package = "deSolve"), "/doc/examples/dynload", sep = ""))

## show package vignette with how to use deSolve
## + source code of the vignette
vignette("deSolve")
edit(vignette("deSolve"))

## show package vignette with tutorial about how to use compiled models
## + source code of the vignette
## + directory with C and FORTRAN sources
vignette("compiledCode")
edit(vignette("compiledCode"))
browseURL(paste(system.file(package = "deSolve"), "/doc", sep = ""))
}
}

\keyword{ package }