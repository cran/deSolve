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
  to the FORTRAN functions dvode, zvode and daspk and a C-implementation of
  solvers of the Runge-Kutta family with fixed or variable time steps.

  The package contains routines designed for solving ODEs
  resulting from 1-D, 2-D and 3-D partial differential equations (PDE)
  that have been converted to ODEs by numerical differencing.
  It includes root-finding (or event location) and provides access to lagged
  variables and derivatives.

  The system of differential equations is written as an \R{} function or
  defined in compiled code that has been dynamically loaded, see
  package vignette \href{../doc/compiledCode.pdf}{compiledCode} for details. The
  solvers may be used as part of a modeling package for differential
  equations, or for parameter estimation using any appropriate
  modeling tool for non-linear models in \R{} such as
  \code{\link{optim}}, \code{\link{nls}}, \code{\link{nlm}} or
  \code{\link[nlme]{nlme}}, or \code{\link[FME]{FME}}.


  \bold{Package Vignettes, Examples, Online Resources}

  \itemize{
    \item Solving Initial Value Differential Equations in R
      (\href{../doc/deSolve.pdf}{pdf}, \href{../doc/deSolve.R}{R code})
    \item Writing Code in Compiled Languages
      (\href{../doc/compiledCode.pdf}{pdf}, \href{../doc/compiledCode.R}{R code})
    \item Examples in R (\url{../doc/examples}), and in Fortran or C (\url{../doc/dynload}, \url{../doc/dynload-dede})
    \item deSolve homepage: \url{http://desolve.r-forge.r-project.org} (Papers, Books, PDFs)
    \item Mailing list: \url{mailto:r-sig-dynamic-models@r-project.org}
   }

}
\author{
  Karline Soetaert, Thomas Petzoldt, R. Woodrow Setzer
}
\references{

  Karline Soetaert, Thomas Petzoldt, R. Woodrow Setzer (2010):  Solving
  Differential Equations in R: Package deSolve Journal of Statistical
  Software, 33(9), 1--25.  \url{https://www.jstatsoft.org/v33/i09/}

  Karline Soetaert, Thomas Petzoldt, R. Woodrow Setzer (2010): Solving
  differential equations in R. The R Journal 2(2), 5-15.
  \href{https://journal.r-project.org/archive/2010-2/RJournal_2010-2_Soetaert~et~al.pdf}{pdf}

  Karline Soetaert, Thomas Petzoldt (2011): Solving ODEs, DAEs, DDEs and
  PDEs in R. Journal of Numerical Analysis, Industrial and Applied
  Mathematics (JNAIAM) 6(1-2),
  51-65. %\href{http://jnaiam.org/uploads/jnaiam_6_4.pdf}{pdf}

  Karline Soetaert, Jeff Cash, Francesca Mazzia, (2012):
  Solving Differential Equations in R. Springer, 248 pp.

  Alan C. Hindmarsh (1983): ODEPACK, A Systematized Collection of ODE
  Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
  North-Holland, Amsterdam, pp. 55-64.

  L. R. Petzold, (1983):  A Description of DASSL: A Differential/Algebraic
  System Solver, in Scientific Computing, R. S. Stepleman et al. (Eds.),
  North-Holland, Amsterdam, pp. 65-68.

  P. N. Brown, G. D. Byrne, A. C. Hindmarsh (1989):  VODE: A Variable
  Coefficient ODE Solver, SIAM J. Sci. Stat. Comput., 10, pp. 1038-1051.

  See also the references given on the specific help pages of the
  different methods.

}
\seealso{
  \code{\link{ode}} for a general interface to most of the ODE solvers,

  \code{\link{ode.band}} for solving models with a banded Jacobian,

  \code{\link{ode.1D}}, \code{\link{ode.2D}}, \code{\link{ode.3D}},
  for integrating 1-D, 2-D and 3-D models,

  \code{\link{dede}} for a general interface to the delay differential
  equation solvers,

  \code{\link{lsoda}}, \code{\link{lsode}},
  \code{\link{lsodes}}, \code{\link{lsodar}}, \code{\link{vode}},
  for ODE solvers of the Livermore family,

  \code{\link{daspk}}, for a DAE solver up to index 1, of the Livermore family,

  \code{\link{radau}} for integrating DAEs up to index 3 using an implicit
    Runge-Kutta,

  \code{\link{rk}}, \code{\link{rkMethod}}, \code{\link{rk4}},
  \code{\link{euler}} for Runge-Kutta solvers,

  \code{\link{DLLfunc}}, \code{\link{DLLres}}, for testing model implementations
  in compiled code,

  \code{\link{forcings}},  \code{\link{events}}, for how to implement forcing
  functions (external variables) and events (sudden changes in state variables),

  \code{\link{lagvalue}}, \code{\link{lagderiv}}, for how to get access to
  lagged values of state variables and derivatives.

}

\examples{
library(deSolve)

## Chaos in the atmosphere
Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <-  a * X + Y * Z
    dY <-  b * (Y - Z)
    dZ <- -X * Y + c * Y - Z
    list(c(dX, dY, dZ))
  })
}

parameters <- c(a = -8/3, b = -10, c = 28)
state      <- c(X = 1, Y = 1, Z = 1)
times      <- seq(0, 100, by = 0.01)

out <- ode(y = state, times = times, func = Lorenz, parms = parameters)

plot(out)

## add a 3D figure if package scatterplot3D is available
if (require(scatterplot3d))
  scatterplot3d(out[,-1], type = "l")

}
\keyword{ package }
