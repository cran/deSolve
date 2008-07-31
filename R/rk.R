

### rk   Top level Function for solvers of the Runge-Kutta family
###      See help(rk) for a description of parameters
###
###



rk <- function(y, times, func, parms, rtol=1e-6, atol=1e-6,
	tcrit = NULL, verbose=FALSE, hmin = 0, hmax = NULL, hini = hmax,
  method = rkMethod("rk45dp7", ... ), maxsteps = 5000, ...) {

    ## check input
    if (!is.numeric(y))       stop("`y' must be numeric")
    n <- length(y)
    if (! is.null(times)&&!is.numeric(times))
        stop("`times' must be NULL or numeric")
    if (!is.function(func) && !is.character(func))
      stop("`func' must be a function")
    if (!is.numeric(rtol))   stop("`rtol' must be numeric")
    if (!is.numeric(atol))   stop("`atol' must be numeric")
    if (!is.null(tcrit) & !is.numeric(tcrit)) stop("`tcrit' must be numeric")
    if (length(atol) > 1 && length(atol) != n)
        stop("`atol' must either be a scaler, or as long as `y'")
    if (length(rtol) > 1 && length(rtol) != n)
        stop("`rtol' must either be a scaler, or as long as `y'")
    if (!is.numeric(hmin))   stop("`hmin' must be numeric")
    if (hmin < 0) stop ("`hmin' must be a non-negative value")
    if (is.null(hmax))
       hmax <- ifelse (is.null(times), 0, max(abs(diff(times))))
    if (!is.numeric(hmax))   stop("`hmax' must be numeric")
    if (hmax < 0)            stop ("`hmax' must be a non-negative value")
    if (hini < 0)            stop("`hini' must be a non-negative value")

    if (!is.numeric(y))     stop("`y' must be numeric")
    if (!is.numeric(times)) stop("`times' must be numeric")
    if (!is.function(func)) stop("`func' must be a function")
    
    if (is.character(method)) method <- rkMethod(method)

    ## Pass state names to function
    Ynames <- attr(y, "names")

    Func <- function(time, state, parms) {
      attr(state, "names") <- Ynames
      func(time, state, parms, ...)[[1]]
    }   

    Func2 <- function(time, state, parms) { 
      attr(state, "names") <- Ynames
      func(time, state, parms, ...)
    }    

    ## Call func once to figure out whether and how many "global"
    ## results it wants to return and some other safety checks
    rho <- environment(func)
    tmp <- eval(Func2(times[1], y, parms), rho)
    if (!is.list(tmp)) stop("Model function must return a list\n")
    if (length(tmp[[1]]) != length(y))
      stop(paste("The number of derivatives returned by func() (",
                 length(tmp[[1]]),
                 "must equal the length of the initial conditions vector (",
                 length(y),")", sep="")
      )
    Nglobal <- if (length(tmp) > 1) length(unlist(tmp[-1])) else 0
    Nmtot   <- attr(unlist(tmp[-1]), "names")

    ## -----------------------------------------------------------------------
    varstep <- method$varstep
    if (varstep) { # methods with variable step size
      out <- rkAuto(y, times, Func, parms, rtol = rtol, atol = atol, tcrit = tcrit,
               verbose = verbose, hmin = hmin, hmax = hmax, hini = hini, 
               method = method, maxsteps = maxsteps, ...)
    } else {       # fixed step methods
      out <- rkFixed(y, times, Func, parms, tcrit = tcrit,
         verbose = verbose, hini = hini, method = method, ...)
    }
    istate <- attr(out, "istate") # remember internal information
    ## -----------------------------------------------------------------------
    nm <- c("time",
      if (!is.null(attr(y, "names"))) names(y) else as.character(1:n)
    )
    if (is.null(method$d)) {
      ## simple linear interpolation
      ## for all the methods that have no polynomial coefficients d
      m   <- ncol(out)
      res <- matrix(0, nrow = length(times), ncol = m)
      res[,1] <- times
      for (i in 2:m) {
        res[,i] <- as.vector(approx(out[,1], out[,i], times)$y)
      }
      out <- res
    }
    ## external outputs
    if (Nglobal > 0) {
      out2 <- matrix(nrow = nrow(out), ncol = Nglobal)
      for (i in 1:nrow(out2))
        out2[i,] <- unlist(Func2(out[i, 1], out[i, -1], parms)[-1])
      out <- cbind(out, out2)
      nm  <- c(nm,
        if (!is.null(Nmtot)) Nmtot else as.character((n + 1) : (n + Nglobal))
      )
    }
    dimnames(out) <- list(NULL, nm)
    attr(out, "istate") <- istate
    out
}
