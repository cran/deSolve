\name{rkMethod}
\alias{rkMethod}

\title{Collection of Parameter Sets (Butcher Arrays) for 
  the Runge-Kutta Family of ODE Solvers
}
\description{
  This function returns a list specifying coefficients and properties of
  ODE solver methods from the Runge-Kutta family.
}
\usage{
rkMethod(method = NULL, ...)
}
\arguments{
  
  \item{method }{a string constant naming one of the pre-defined methods
    of the Runge-Kutta family of solvers.  The most common methods are
    the fixed-step methods \code{"euler"}, \code{"rk2"}, \code{"rk4"} or
    the variable step methods \code{"rk23bs"} (alias \code{"ode23"}),
    \code{"rk45dp7"} (alias \code{"ode45"}) or \code{"rk78f"}.
  }

  \item{\dots }{specification of a user-defined solver, see \emph{Value}
    and example below.
  }
}
\details{
  This function supplies \code{method} settings for \code{\link{rk}} or
  \code{\link{ode}}. If called without arguments, the names of all
  currently implemented solvers of the Runge-Kutta family are returned.
  
  The following comparison gives an idea how the algorithms of \pkg{deSolve}
  are related to similar algorithms of other simulation languages:
  
  \tabular{lll}{
    \bold{rkMethod}    \tab | \tab \bold{Description} \cr
    "euler"            \tab | \tab Euler's Method\cr
    "rk2"              \tab | \tab 2nd order Runge-Kutta, fixed time step (Heun's method)\cr
    "rk4"              \tab | \tab classical 4th order Runge-Kutta, fixed time step\cr
    "rk23"             \tab | \tab Runge-Kutta, order 2(3); Octave: ode23\cr
    "rk23bs", "ode23"  \tab | \tab Bogacki-Shampine, order 2(3); Matlab: ode23\cr
    "rk34f"            \tab | \tab Runge-Kutta-Fehlberg, order 3(4)\cr
    "rk45ck"           \tab | \tab Runge-Kutta Cash-Karp, order 4(5)\cr
    "rk45f"            \tab | \tab Runge-Kutta-Fehlberg, order 4(5); Octave: ode45, pair=1 \cr
    "rk45e"            \tab | \tab Runge-Kutta-England, order 4(5)\cr
    "rk45dp6"          \tab | \tab Dormand-Prince, order 4(5), local order 6\cr
    "rk45dp7", "ode45" \tab | \tab Dormand-Prince 4(5), local order 7 \cr
                       \tab | \tab (also known as dopri5; MATLAB: ode45; Octave: ode45, pair=0)\cr
    "rk78f"            \tab | \tab Runge-Kutta-Fehlberg, order 7(8)\cr		       
    "rk78dp"           \tab | \tab Dormand-Prince, order 7(8)\cr
  }
  
  Note that this table is based on the Runge-Kutta coefficients only,
  but the algorithms differ also in their implementation, in their
  stepsize adaption strategy and interpolation methods. 

  The table reflects the state at time of writing and it is of course possible
  that implementations change.
  
  Methods \code{"rk45dp7"} (alias \code{"ode45"}) and \code{"rk45ck"} contain
  specific and efficient built-in interpolation schemes (dense output).

  As an alternative, Neville-Aitken polynomials can be used to interpolate between
  time steps. This is available for all RK methods and may be useful to speed
  up computation if no dense-output formula is available. Note however, that 
  this can introduce considerable local error; it is disabled by default 
  (see \code{nknots} below).
}

\note{  
  \itemize{
    \item Adaptive stepsize Runge-Kuttas are preferred if the solution
      contains parts where it changes fast, and parts where nothing much
      happens. They will take small steps over bumpy ground and long steps
      over uninteresting terrain.

    \item As a suggestion, one may use \code{"rk23"} (alias
      \code{"ode23"}) for simple problems and \code{"rk45dp7"} (alias
      \code{"ode45"}) for rough problems. The default solver is
      \code{"rk45dp7"} (alias "ode45"), because of its relatively high
      order (4), re-use of the last intermediate steps (FSAL = first
      same as last) and built-in polynomial interpolation (dense
      output).

    \item Solver \code{"rk23bs"}, that supports also FSAL, may be useful for
      slightly stiff systems if demands on precision are relatively low.

    \item Another good choice, assuring medium accuracy, is the Cash-Karp
      Runge-Kutta method, \code{"rk45ck"}.

    \item Classical \code{"rk4"} is traditionally used in cases where an
      adequate stepsize is known a-priori or if external forcing data
      are provided for fixed time steps only and frequent interpolation
      of external data needs to be avoided.
  
    \item Method \code{"rk45dp7"} (alias \code{"ode45"}) contains an
      efficient built-in interpolation scheme (dense output) based on
      intermediate function evaluations.
    }

    Starting with version 1.8 implicit Runge-Kutta (\code{irk}) methods
    are also supported by the general \code{rk} interface, however their
    implementation is still experimental.  Instead of this you may
    consider \code{\link{radau}} for a specific full implementation of an
    implicit Runge-Kutta method.  
}
   
\value{
  A list with the following elements:
  
  \item{ID}{name of the method (character)}
  \item{varstep}{boolean value specifying if the method allows for
    variable time step (\code{TRUE}) or not (\code{FALSE}).
  }

  \item{FSAL}{(first same as last) optional boolean value specifying if
    the method allows re-use of the last function evaluation
    (\code{TRUE}) or not (\code{FALSE} or \code{NULL}).
  }
  
  \item{A}{coefficient matrix of the method. As \code{link{rk}} supports
    only explicit methods, this matrix must be lower triangular.
    \code{A} must be a vector for fixed step methods where only the
    subdiagonal values are different from zero.
  }
  \item{b1}{coefficients of the lower order Runge-Kutta pair.
  }
  \item{b2}{coefficients of the higher order Runge-Kutta pair
    (optional, for embedded methods that allow variable time step).
  }
  \item{c}{coefficients for calculating the intermediate time steps.}
  
  \item{d}{optional coefficients for built-in polynomial interpolation
    of the outputs from internal steps (dense output), currently only
    available for method \code{rk45dp7} (Dormand-Prince).
  }
  
  \item{densetype}{optional integer value specifying the dense output formula;
    currently only \code{densetype = 1} for \code{rk45dp7} (Dormand-Prince)
    and \code{densetype = 2} for \code{rk45ck} (Cash-Karp) are supported. 
    Undefined values (e.g., \code{densetype = NULL}) disable dense output.
  }
  
  \item{stage}{number of function evaluations needed (corresponds to
    number of rows in A).
  }
  \item{Qerr}{global error order of the method, important for automatic
    time-step adjustment.
  }

  \item{nknots}{integer value specifying the order of interpolation
    polynomials for methods without dense output.  If \code{nknots} < 2
    (the default) then internal interpolation is switched off and
    integration is performed step by step between external time steps.

    If \code{nknots} is between 3 and 8, Neville-Aitken polynomials
    are used, which need at least \code{nknots + 1} internal time steps.
    Interpolation may speed up integration but can lead to local
    errors higher than the tolerance, especially if external and
    internal time steps are very different.
  }

  \item{alpha}{optional tuning parameter for stepsize
    adjustment. If \code{alpha} is omitted, it is set to
    \eqn{1/Qerr - 0.75 beta}. The default value is
    \eqn{1/Qerr} (for \code{beta} = 0).}

  \item{beta}{optional tuning parameter for stepsize adjustment. Typical
    values are \eqn{0} (default) or \eqn{0.4/Qerr}.
  }

}

\references{

  Bogacki, P. and Shampine L.F. (1989) A 3(2) pair of Runge-Kutta
  formulas, Appl. Math. Lett. \bold{2}, 1--9.

  Butcher, J. C. (1987) The numerical analysis of ordinary differential
  equations, Runge-Kutta and general linear methods, Wiley, Chichester
  and New York.

  Cash, J. R. and Karp A. H., 1990. A variable order Runge-Kutta method
  for initial value problems with rapidly varying right-hand sides,
  ACM Transactions on Mathematical Software \bold{16}, 201--222.

  Dormand, J. R. and Prince, P. J. (1980) A family of embedded
  Runge-Kutta formulae, J. Comput. Appl. Math. \bold{6}(1), 19--26.

  Engeln-Muellges, G. and Reutter, F. (1996) Numerik Algorithmen:
  Entscheidungshilfe zur Auswahl und Nutzung. VDI Verlag, Duesseldorf.

  Fehlberg, E. (1967) Klassische Runge-Kutta-Formeln fuenfter and
  siebenter Ordnung mit Schrittweiten-Kontrolle, Computing
  (Arch. Elektron. Rechnen) \bold{4}, 93--106.

  Kutta, W. (1901) Beitrag zur naeherungsweisen Integration totaler
  Differentialgleichungen, Z. Math. Phys. \bold{46}, 435--453.

  Octave-Forge - Extra Packages for GNU Octave, Package OdePkg.
  \url{http://octave.sourceforge.net/odepkg}

  Prince, P. J.  and Dormand, J. R. (1981) High order embedded
  Runge-Kutta formulae, J. Comput. Appl. Math. \bold{7}(1), 67--75.

  Runge, C. (1895) Ueber die numerische Aufloesung von
  Differentialgleichungen, Math. Ann. \bold{46}, 167--178.
    
  MATLAB (R) is a registed property of The Mathworks
  Inc. \url{http://www.mathworks.com/} }

\author{Thomas Petzoldt \email{thomas.petzoldt@tu-dresden.de}}

\seealso{\code{\link{rk}}, \code{\link{ode}}}
\examples{
rkMethod()          # returns the names of all available methods
rkMethod("rk45dp7") # parameters of the Dormand-Prince 5(4) method
rkMethod("ode45")   # an alias for the same method

func <- function(t, x, parms) {
  with(as.list(c(parms, x)),{
    dP  <- a * P      - b * C * P
    dC  <- b * P * C  - c * C
    res <- c(dP, dC)
    list(res)
  })
}
times  <- seq(0, 200, length = 101)
parms  <- c(a = 0.1, b = 0.1, c = 0.1)
x <- c(P = 2, C = 1)

## rk using ode45 as the default method
out <- rk(x, times, func, parms)

## all methods can be called also from 'ode' by using rkMethod
out <- ode(x, times, func, parms, method = rkMethod("rk4"))

## 'ode' has aliases for the most common RK methods
out <- ode(x, times, func, parms, method = "ode45")

##===========================================================================
## Comparison of local error from different interpolation methods
##===========================================================================

## lsoda with lower tolerances (1e-10) used as reference
o0 <- ode(x, times, func, parms, method = "lsoda", atol = 1e-10, rtol = 1e-10)

## rk45dp7 with hmax = 10 > delta_t  = 2
o1 <- ode(x, times, func, parms, method = rkMethod("rk45dp7"), hmax = 10)

## disable dense-output interpolation
## and use only Neville-Aitken polynomials instead
o2 <- ode(x, times, func, parms, 
  method = rkMethod("rk45dp7", densetype = NULL, nknots = 5), hmax = 10)

## stop and go: disable interpolation completely
## and integrate explicitly between external time steps
o3 <- ode(x, times, func, parms, 
  method = rkMethod("rk45dp7", densetype = NULL, nknots = 0, hmax=10))
  
## compare different interpolation methods with lsoda
mf <- par("mfrow" = c(4, 1))
matplot(o1[,1], o1[,-1], type = "l", xlab = "Time", main = "State Variables",
  ylab = "P, C")
matplot(o0[,1], o0[,-1] - o1[,-1], type = "l", xlab = "Time", ylab = "Diff.",
  main="Difference between lsoda and ode45 with dense output")
abline(h = 0, col = "grey")
matplot(o0[,1], o0[,-1] - o2[,-1], type = "l", xlab = "Time", ylab = "Diff.", 
  main="Difference between lsoda and ode45 with Neville-Aitken")
abline(h = 0, col = "grey")
matplot(o0[,1], o0[,-1] - o3[,-1], type = "l", xlab = "Time", ylab = "Diff.",
  main="Difference between lsoda and ode45 in 'stop and go' mode")
abline(h = 0, col = "grey")
par(mf)

##===========================================================================
## rkMethod allows to define user-specified Runge-Kutta methods
##===========================================================================
out <- ode(x, times, func, parms, 
          method = rkMethod(ID = "midpoint",
            varstep = FALSE,
            A      = c(0, 1/2),
            b1      = c(0, 1),
            c       = c(0, 1/2),
            stage   = 2,
            Qerr    = 1
          )
        )
plot(out)

## compare method diagnostics
times  <- seq(0, 200, length = 10)
o1 <- ode(x, times, func, parms, method = rkMethod("rk45ck"))
o2 <- ode(x, times, func, parms, method = rkMethod("rk78dp"))
diagnostics(o1)
diagnostics(o2)

}

\keyword{ math }

