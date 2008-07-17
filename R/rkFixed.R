## Generalized sover for Runge-Kutta methods with fixed time step
rkFixed <- function( y, times, func, parms, tcrit = NULL,
  verbose=FALSE, hini=0, method = rkMethod("rk4", ... ), ...) {

  stage <- method$stage
  A     <- method$A
  bb1   <- method$b1
  #bb2   <- method$b2
  cc    <- method$c
  qerr  <- 1/method$Qerr

  FF      <- matrix(0, nrow=length(y), ncol=stage)

  y0   <- y
  out  <- c(times[1], y0)
  if (verbose) {
    cat("method=", method$ID, "\n")
    cat("hini=", hini, "\n")
  }
  t    <- min(times)
  tmax <- max(times, tcrit)                  # NULL is handled automatically by max
  dt <- hini
  ## derive internal (!) time step
  times <- unique(c(seq(t, tmax, dt), tmax)) # last step may possibly be shorter
  if (!is.matrix(A)) {                       # "A" coefficients given as subdiagonal
    for (i in 1:(length(times) - 1)) {
      t  <- times[i]
      for (j in 1:stage) {
        if (j == 1) Fj <- 0 else Fj <- A[j] * FF[ ,j - 1]
        FF[, j] <- dt * func(t + dt * cc[j], y0 + Fj, parms)
      }
      dy <- FF %*% bb1
      y1 <- y0 + dy
      out<- rbind(out, c(times[i + 1], y1))
      y0 <- y1
    }
  } else {                                   # "A" coefficients as matrix
    for (i in 1:(length(times) - 1)) {
      t  <- times[i]
      for (j in 1:stage) {
        k  <- 1
        Fj <- 0
        while (k < j) {
          #if (j == 1) Fj <- 0 else Fj <- Fj + A[j, k] * FF[ , k]
          Fj <- Fj + A[j, k] * FF[ , k]
          k <- k + 1
        }
        FF[, j] <- dt * func(t + dt * cc[j], y0 + Fj, parms)
      }
      dy <- FF %*% bb1
      y1 <- y0 + dy
      out<- rbind(out, c(times[i + 1], y1))
      y0 <- y1
    }
  } # end if
  out
}
