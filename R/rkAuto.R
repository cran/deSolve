## Generalized solver for Runge-Kutta methods with variable time step
## This function is internal and not intended for the end user
rkAuto <- function(
  y, times, func, parms, rtol = 1e-6, atol = 1e-6,
	tcrit = NULL, verbose = FALSE,
  hmin = 0, hmax = NULL, hini = 0,
  method = rkMethod("rk45dp7", ... ), maxsteps = 5000, ...) {

  stage <- method$stage
  A     <- method$A
  bb1   <- method$b1
  bb2   <- method$b2
  cc    <- method$c
  dd    <- method$d
  qerr  <- 1/method$Qerr
  FSAL  <- ifelse(is.null(method$FSAL), FALSE, method$FSAL)

  ## track essential internal information
  ## experimental! Do it similar like lsoda
  istate <- numeric(23)

  Nstates <- length(y)
  FF      <- matrix(0, nrow = Nstates, ncol = stage)

  steps <- 0
  y0   <- y
  out  <- c(times[1], y0)
  accept <- FALSE

  if (verbose) {
    cat("method=", method$ID, "\n")
    cat("hini=", hini, "\n")
  }
  t    <- min(times)
  tmax <- max(times, tcrit)   # NULL is handled automatically by max()
  dt   <- min(hmax, hini)
  hmax <- min(hmax, tmax - t) # limit hmax to the time range (esp. if hmax = Inf)

  ## function evaluations
  while ((t + dt) <= tmax) {
    ## save one step if the coefficients allow this (first same as last)
    if (FSAL & accept) {
          j1 <- 2
          FF[, 1] <- FF[, stage]
        } else {
          j1 <- 1
    }
    for (j in j1:stage) {
      Fj <- 0
      k  <- 1
      while (k < j) {
        Fj <- Fj + A[j, k] * FF[ ,k]  * dt
        k <- k + 1
      }
      FF[, j] <- func(t + dt * cc[j], y0 + Fj, parms)
    }
    ## Estimation of new values
    dy1  <- FF %*% bb1
    dy2  <- FF %*% bb2
    y1   <- y0 + dt * dy1
    y2   <- y0 + dt * dy2

    ## stepsize adjustment after Press et. al (2007), Numerical Recipes in C
    yabs  <- pmax.int(abs(y0), abs(y2))
    scal  <- atol + yabs * rtol
    delta <- abs(y2 - y1)

    err <- delta / scal
    err  <- err[is.finite(err)] # remove Inf, in case of atol == rtol == 0
    if (length(err) > 0) {
      ## Press (2007): maximum is fine ...
      err <- max(err)
      ## ... but he takes the euklidean norm
      #err <- sqrt(sum(err)^2 / length(err))
    } else {
      err  <- 1
    }

    if (err == 0) {
       accept <- TRUE
       dtnew <- hmax       # hmax must not be Inf or NULL !!!
       if (verbose) cat("t=", t, " err=", err, " h=", dt, " +++ \n")
    } else if (err < 1) {  # accept
       accept <- TRUE
       ## safety factor = 0.9
       dtnew  <- min(hmax, dt * 0.9 * (err ^ -qerr))
       if (verbose) cat("t=", t, " err=", err, " h=", dt, " + \n")
    } else if (err > 1){  # reject
       accept <- FALSE
       dtnew  <- dt * max(0.9 * (err ^ -qerr), 0.2)
       if (verbose) cat("t=", t, " err=", err, " h=", dt, " - \n")
    } else {             # err == 1
       accept <- TRUE
       dtnew <- dt
       if (verbose) cat("t=", t, " err=", err, " h=", dt, " = \n")
    }

    ## Final check
    if (dt < hmin) {
      accept <- TRUE
      if (verbose) cat("h < hmin \n")
      warning("h < hmin")
      istate[1] <- -2
      dtnew <- hmin
    }
    ## data storage. Store imprecise results too, but warn if h < hmin
    if (accept) {
      ## polynomial interpolation ("dense output")
      if (!is.null(dd)) { # i.e. polynomial constants available
          densr <- denspar(FF, y0, y1, dt, stage, dd)
          tdens <- times[times > t & times <= (t + dt)]
          if (length(tdens) > 0) {
            newout <- densout(densr, t, tdens, dt)
            out <- rbind(out, cbind(tdens, newout))
          }
      }  else {
         out <- rbind(out, c(t + dt, y2))
      }
      t   <- t + dt
      y0  <- y2
    }
    steps <- steps + 1
    if (steps > maxsteps)
      stop("
        An excessive amount of work (> maxsteps ) was done,
        but integration was not successful -
        increase maxsteps, increase atol/rtol, check your equations
        or select an alternative algorithm.
        ")
    dt  <- min(dtnew, tmax - t)
    if (t >= tmax) break
  }

  ## attach essential internal information
  ## experimental! Codes similar like lsoda
  istate[12] <- steps                   # number of steps
  istate[13] <- steps * (stage - FSAL)  # number of function evaluations
  istate[15] <- method$Qerr             # order of the method
  attr(out, "istate") <- istate
  out
}
