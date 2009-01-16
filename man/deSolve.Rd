\name{deSolve-package}
\alias{deSolve-package}
\alias{deSolve}
\docType{package}
\title{
  General Solvers for Ordinary Differential Equations
  (ODE) and for Differential Algebraic Equations (DAE)
}
\description{
  Functions that solve a system of first-order ordinary differential
  equations (ODE) and of differential algebraic equations (DAE).

  The functions provide an interface to the FORTRAN functions lsoda,
  lsodar, lsode, lsodes, dvode and daspk. The package also contains
  routines designed for solving uni-and multicomponent 1-D and 2-D
  reactive transport models.
}
\details{
  \tabular{ll}{
    Package: \tab deSolve\cr
    Type:    \tab Package\cr
    Version: \tab 1.1-3\cr
    Date:    \tab 2009-01-12\cr
    License: \tab GNU Public License 2 or above\cr
  }

  The system of ODE's is written as an \R function or be defined in
  compiled code that has been dynamically loaded, see package vignette
  (\code{vignette(compiledCode)}) for details. The solvers may be used
  as part of a modeling package for differential equations, or for
  parameter estimation using any appropriate modeling tool for
  non-linear models in \R such as \code{\link{optim}},
  \code{\link[nls]{nls}}, \code{\link{nlm}} or \code{\link[nlme]{nlme}}.
  
}
\author{
  Karline Soetaert,
  Thomas Petzoldt,
  R. Woodrow Setzer (Maintainer)
}
%\references{
%% to be filled out when our paper is published
%}
\seealso{
  \code{\link{ode}}, \code{\link{lsoda}}, \code{\link{lsode}},
  \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}},
  \code{\link{daspk}}, \code{\link{rk}}.
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(aquaphy)
example(lsoda)
example(ode.band)
example(ode.1D)
example(ode.2D)

## run demos
demo("rk_solvers") # comparison of lsoda with Runge-Kutta-Type Solvers
demo("CCL4model")  # a model fitting example (this will take some time)

## open the directory with source code of demos
browseURL(paste(system.file(package = "deSolve"), "/demo", sep = ""))

## open the directory with R sourcecode examples
browseURL(paste(system.file(package = "deSolve"), "/examples", sep = ""))
## open the directory with C and FORTRAN sourcecode examples
browseURL(paste(system.file(package = "deSolve"), "/dynload", sep = ""))

## show package vignette with tutorial about how to use compiled models
## + source code of the vignette
## + directory with C and FORTRAN sources
vignette("compiledCode")
edit(vignette("compiledCode"))
browseURL(paste(system.file(package = "deSolve"), "/doc", sep = ""))
}
}

\keyword{ package }