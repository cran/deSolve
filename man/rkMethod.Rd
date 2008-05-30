\name{rkMethod}
\alias{rkMethod}

\title{Collection of Parameter Sets (Butcher Arrays) for 
  the Runge-Kutta Family of ODE Solvers}
\description{
  This function returns a list specifying coefficients and properties of
  ODE solver methods from the Runge-Kutta family.
}
\usage{
rkMethod(method = NULL, ...)
}

\arguments{
  \item{method}{a string constant naming one of the
    pre-defined methods of the Runge-Kutta family of solvers.
    The most common methods are the fixed-step methods "euler", "rk2", "rk4" 
    or the variable step methods "rk23bs", "rk34f", "rk45f" or "rk45dp7".}
  \item{\dots}{specification of user-defined solver, see \emph{Value} and example below.}
}
\details{
  This function supplies \code{method} settings for \code{\link{rk}} or 
  \code{\link{ode}}. If called without arguments, the names of all implemented
  solvers of the Runge-Kutta family is returned.
  
  The following comparison gives an idea how the algorithms of \pkg{deSolve}
  are related to similar algorithms of other simulation languages:
  
  \tabular{lll}{
    \bold{\R}       \tab | \tab \bold{Description} \cr
    "euler"         \tab | \tab Euler's Method\cr
    "rk2"           \tab | \tab 2nd order Runge-Kutta, fixed time step (Heun's method)\cr
    "rk4"           \tab | \tab classical 4th order Runge-Kutta, fixed time step\cr
%    "rk12"         \tab | \tab Euler-Heun, order 1(2)\cr
    "rk23"          \tab | \tab Runge-Kutta, order 2(3), Octave: \bold{ode23}\cr
    "rk23bs", "ode23"   \tab | \tab Bogacki-Shampine, order 2(3), Matlab: \bold{ode23}\cr
    "rk34f"         \tab | \tab Runge-Kutta-Fehlberg, order 3(4)\cr
    "rk45f"         \tab | \tab Runge-Kutta-Fehlberg, order 4(5), Octave: \bold{ode45, pair=1} \cr
    "rk45e"         \tab | \tab Runke-Kutta-England, order 4(5)\cr
    "rk45dp6"       \tab | \tab Dormand-Prince, order 4(5), local order 6\cr
    "rk45dp7", "ode45"  \tab | \tab Dormand-Prince 4(5) coefficients 
        (also known as \bold{dopri5}, Matlab: \bold{ode45}, Octave: \bold{ode45, pair=0})\cr
  }
  
  Note that this table is based on the Runge-Kutta coefficients only, 
  but  the algorithms do also differ in their implementation and in
  their stepsize adaption strategy and interpolation methods.
}

\note{  
  \itemize{
   \item The Runge-Kutta solvers are primarily provided for didactic reasons. 
     For most practical cases, solvers of the Livermore family 
     (\code{\link{lsoda}}, \code{\link{lsode}}, \code{\link{lsodes}}, 
     \code{\link{lsodar}}, \code{\link{vode}}, \code{\link{daspk}})
     are superior because of higher efficiency and faster implementation 
     (FORTRAN and C). 
     In addition to this, some of the Livermore solvers are also suitable for stiff 
     ODEs, differential algebraic equations (DAEs), or partial differential equations 
     (PDEs).
   \item Among the \code{rk} solvers, "rk45dp7" is used by default,
     because of its high order (4), re-use of the last intermediate
     steps (FSAL = first same as last) and built-in polynomial interpolation 
     (dense output). Solver "rk23bs", that supports also FSAL, may be useful
     for slightly stiff systems if demands on precision are low. Classical 
     "rk4" is traditionally used in cases where an adequate stepsize is
     known a-priori or if external forcing data are provided
     for fixed time steps only and interpolation of external data is not desired.
  }
}
\value{
  A list with the following elements:
  \item{ID}{name of the method (character)}
  \item{varstep}{boolean value specifying if the method allows for variable time step (\code{TRUE}) or not (\code{FALSE}).}
    \item{FSAL}{(first step as last) boolean value specifying if the method allows
      re-use of the last function evaluation (\code{TRUE}) or not (\code{FALSE} or \code{NULL}).}
  \item{A}{coefficient matrix of the method. 
    As \code{link{rk}} supports only explicit methods, this matrix must be lower triangular.
    \code{A} can also be a vector if only the subdiagonal values are different from zero}
  \item{b1}{weighting coefficients for averaging the function evaluations of method 1}
  \item{b2}{weighting coefficients for averaging the function evaluations of method 2
    (optional, for embedded methods that allow variable time step)}
  \item{c}{coefficients for calculating the intermediate time steps}  
  \item{d}{coefficients for polynomial interpolation of the outputs from internal
    steps (dense output), currently only available for method \code{rk45dp7} (Dormand-Prince).}
  \item{stage}{number of function evaluations needed (corresponds to number of rows in A)}  
  \item{Qerr}{global error order of the method, important for automatic time-step adjustment.}    
}
\references{
  Runge, C. (1895) Ueber die numerische Aufloesung von Differentialgleichungen,
  Math. Ann. \bold{46}, 167--178.

  Kutta, W. (1901) Beitrag zur naeherungsweisen Integration totaler Differentialgleichungen,
  Z. Math. Phys. \bold{46}, 435--453.

  Dormand, J. R. and Prince, P. J. (1980) A family of embedded Runge-Kutta formulae,
  J. Comput. Appl. Math. \bold{6}(1), 19--26.

  Dormand, J. R. and Prince, P. J. (1981) High order embedded Runge-Kutta formulae,
  J. Comput. Appl. Math. \bold{7}(1), 67--75.
  
  Bogacki, P. and Shampine L.F. (1989) A 3(2) pair of Runge-Kutta formulas, 
    Appl. Math. Lett. \bold{2}, 1--9. 

  Fehlberg, E. (1967) Klassische Runge-Kutta-Formeln fuenfter and siebenter Ordnung mit
  Schrittweiten-Kontrolle, Computing (Arch. Elektron. Rechnen) \bold{4}, 93--106.

  Butcher, J. C. (1987) The numerical analysis of ordinary differential equations,
  Runge-Kutta and general linear methods, Wiley, Chichester and New York.

  Engeln-Muellges, G. and Reutter, F. (1996) Numerik Algorithmen:
  Entscheidungshilfe zur Auswahl und Nutzung. VDI Verlag, Duesseldorf.
  
  Octave-Forge - Extra Packages for GNU Octave, Package OdePkg. 
    \url{http://octave.sourceforge.net/doc/odepkg.html}
    
  MATLAB (R) is a registed property of The Mathworks Inc. \url{http://www.mathworks.com/}  
}
\author{Thomas Petzoldt \email{thomas.petzoldt@tu-dresden.de}}

\seealso{\code{\link{rk}}, \code{\link{ode}}}
\examples{
rkMethod()          # returns the names of all methods
rkMethod("rk45dp7") # parameters of the Dormand-Prince 5(4) method
rkMethod("ode45")   # an alias for the same method

func <- function(t, x, parms) {
  with(as.list(c(parms, x)),{
    dP  <- a * P        - b * K * P
    dK  <- b * P * K  - c * K
    res <- c(dP, dK)
    list(res)
  })
}
times  <- seq(0, 20, length = 21)
parms  <- c(a = 0.1, b = 0.1, c = 0.1)
x <- c(P=2, K=1)

ode(x, times, func, parms, method = rkMethod("rk4"))

ode(x, times, func, parms, method = "ode45")

## disable polynomial interpolation (dense output)
## and fall back to linear approximation
ode(x, times, func, parms, method = rkMethod("rk45dp7", d = NULL))

## define and use a new rk method
ode(x, times, func, parms, 
  method = rkMethod(ID = "midpoint",
    varstep = FALSE,
    #A       = matrix(c(0, 0, 1/2, 0), nrow=2, byrow=TRUE), 
    # or simply, because this A is nonzero only in the subdiagonal
    A      = c(0, 1/2),
    b1      = c(0, 1),
    c       = c(0, 1/2),
    stage   = 2,
    Qerr    = 1
  )
)

}

\keyword{ math }

