
### rkFixed
### Generalized solver for Runge-Kutta methods with fixed time step
### This function is internal and not intended to be called directly.



rkFixed <- function(y, times, func, parms, tcrit = NULL,
  verbose = FALSE, hini = 0, method = rkMethod("rk4", ... ), ...) {

  stage <- method$stage
  A     <- method$A
  bb1   <- method$b1
  cc    <- method$c
  qerr  <- 1/method$Qerr

  FF    <- matrix(0, nrow=length(y), ncol=stage)

  y0   <- y
  out  <- c(times[1], y0)
  t    <- min(times)
  tmax <- min(max(times), tcrit)     # NULL is handled automatically by min

  # if hini==0: set internal time step size equal to external step size
  if (hini == 0) h <- diff(times) else h <- hini

  if (verbose) {
    cat("method   =", method$ID, "\n")
    cat("stepsize =", h, "\n")
    cat("tmax     =", tmax, "\n")
  }
  i <- 1; iinc <- (length(h) > 1)    # don't increment i if stepsize is constant
  
  if (!is.matrix(A)) {               # "A" coefficients given as subdiagonal
    while (t < tmax) {
      dt  <- min(h[i], tmax - t)
      for (j in 1:stage) {
        if (j == 1) Fj <- 0 else Fj <- A[j] * FF[ ,j - 1]
        FF[, j] <- dt * func(t + dt * cc[j], y0 + Fj, parms)
      }
      dy  <- FF %*% bb1
      y1  <- y0 + dy
      y0  <- y1
      t   <- t + dt
      i   <- i + iinc
      out <- rbind(out, c(t, y1))
    }
  } else {                           # "A" coefficients as matrix
    while (t < tmax) {
      dt  <- min(h[i], tmax - t)
      for (j in 1:stage) {
        k  <- 1
        Fj <- 0
        while (k < j) {
          Fj <- Fj + A[j, k] * FF[ , k]
          k <- k + iinc
        }
        FF[, j] <- dt * func(t + dt * cc[j], y0 + Fj, parms)
      }
      dy  <- FF %*% bb1
      y1  <- y0 + dy
      y0  <- y1
      t   <- t + dt
      i   <- i + iinc
      out <- rbind(out, c(t, y1))      
    }
  } # end if
  out
}
