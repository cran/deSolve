### ============================================================================
### Interface to a generalized code for solving explicit variable and fixed
### step ODE solvers of the Runge-Kutta family, see helpfile for details.
### ============================================================================

rk <- function(y, times, func, parms, rtol = 1e-6, atol = 1e-6,
  verbose = FALSE, tcrit = NULL, hmin = 0, hmax = NULL, hini = hmax, ynames=TRUE,
  method = rkMethod("rk45dp7", ... ), maxsteps = 5000,
  dllname = NULL, initfunc=dllname, initpar = parms,
  rpar = NULL,  ipar = NULL, nout = 0, outnames=NULL, ...) {

    ## Check inputs
    if (!is.numeric(y))     stop("`y' must be numeric")
    n <- length(y)
    if (! is.null(times)&&!is.numeric(times))
        stop("`times' must be NULL or numeric")
    if (!is.function(func) && !is.character(func))
      stop("`func' must be a function")
    if (!is.numeric(rtol))  stop("`rtol' must be numeric")
    if (!is.numeric(atol))  stop("`atol' must be numeric")
    if (!is.null(tcrit) & !is.numeric(tcrit)) stop("`tcrit' must be numeric")
    if (length(atol) > 1 && length(atol) != n)
        stop("`atol' must either be a scaler, or as long as `y'")
    if (length(rtol) > 1 && length(rtol) != n)
        stop("`rtol' must either be a scaler, or as long as `y'")
    if (!is.numeric(hmin))  stop("`hmin' must be numeric")
    if (hmin < 0) stop ("`hmin' must be a non-negative value")
    if (is.null(hmax))
       hmax <- ifelse (is.null(times), 0, max(abs(diff(times))))
    if (!is.numeric(hmax))  stop("`hmax' must be numeric")
    if (hmax < 0)           stop ("`hmax' must be a non-negative value")
    if (hmax == 0)          hmax <- .Machine$double.xmax # i.e. practically unlimited
    if (hini < 0)           stop("`hini' must be a non-negative value")
    if (!is.numeric(y))     stop("`y' must be numeric")
    if (!is.numeric(times)) stop("`times' must be numeric")
    if (!is.function(func) && !is.character(func))
      stop("`func' must be a function or character vector")
    if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
      stop("You need to specify the name of the dll or shared library where func can be found (without extension)")
    if (maxsteps < 0)       stop("maxsteps must be positive")
    if (!is.finite(maxsteps)) maxsteps <- .Machine$integer.max
    if (is.character(method)) method <- rkMethod(method)
    if (is.null(tcrit)) tcrit <- max(times)

    ## Check if interpolation is switched off
    if (!is.null(method$interpolation) && !method$interpolation) {
      cat("\nMethod without or with disabled interpolation\n")
    } else {
      trange <- diff(range(times))
      ## ensure that we have at least 4..5 knots
      ## to allow 3rd order polynomial interpolation
      ## for methods without built-in dense output
      if ((is.null(method$d) &                             # has no "dense output"?
        (hmax > 0.2 * trange))) {                          # time steps too large?
        hini <- hmax <- 0.2 * trange
        if (hmin < hini) hmin <- hini
        cat("\nNote: Method ", method$ID,
            " needs intermediate steps for interpolation\n")
        cat("hmax decreased to", hmax, "\n")
      }
    }

    ## Model as shared object (DLL)?
    Ynames <- attr(y,"names")
    Initfunc <- NULL
    if(!is.null(dllname)) {
      if (is.loaded(initfunc, PACKAGE = dllname, type = "") ||
          is.loaded(initfunc, PACKAGE = dllname, type = "Fortran")) {
        Initfunc <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
       } else if (initfunc != dllname && ! is.null(initfunc))
         stop(paste("cannot integrate: initfunc not loaded ", initfunc))
    }

    ## If func is a character vector, then copy its value to funcname
    ## check to make sure it describes a function in a loaded dll
    if (is.character(func)) {
      funcname <- func
      ## get the pointer and put it in func
      if(is.loaded(funcname, PACKAGE = dllname)) {
        Func <- getNativeSymbolInfo(funcname, PACKAGE = dllname)$address
      } else stop(paste("cannot integrate: dyn function not loaded", funcname))

      ## If we go this route, the number of "global" results is in nout
      ## and output variable names are in outnames
      Nglobal <- nout
      if (is.null(outnames))
         { Nmtot   <- NULL} else
      if (length(outnames) == nout)
         { Nmtot   <- outnames} else
      if (length(outnames) > nout)
         Nmtot <- outnames[1:nout] else
         Nmtot <- c(outnames,(length(outnames)+1):nout)

      Nstates <- length(y) # assume length of states is correct
      rho <- NULL
      if (is.null(ipar)) ipar <- 0
      if (is.null(rpar)) rpar <- 0
    } else {
      initpar <- NULL # parameter initialisation not needed if function is not a DLL
      rho <- environment(func)

      ## func is overruled, either including ynames, or not
      ## This allows to pass the "..." arguments and the parameters
        if(ynames) {
         Func   <- function(time,state,parms){
           attr(state,"names") <- Ynames
           func   (time,state,parms,...)}
        } else {                            # no ynames...
         Func   <- function(time,state,parms)
           func   (time,state,parms,...)
        }

      ## Call func once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
      tmp <- eval(Func(times[1], y, parms), rho)
      if (!is.list(tmp)) stop("Model function must return a list\n")

      Nstates <-length(y)
      if (length(tmp[[1]]) != Nstates)
        stop(paste("The number of derivatives returned by func() (",
                   length(tmp[[1]]),
                   "must equal the length of the initial conditions vector (",
                   Nstates,")", sep=""))

      ## use "unlist" here because some output variables are vectors/arrays
      Nglobal <- if (length(tmp) > 1)
          length(unlist(tmp[-1]))  else 0
      Nmtot <- attr(unlist(tmp[-1]),"names")
    }

    ## handle length of atol and rtol
    if (Nstates %% length(atol))
      warning("length of atol does not match number of states")
    if (Nstates %% length(rtol))
      warning("length of rtol does not match number of states")

    atol <- rep(atol, length.out = Nstates)
    rtol <- rep(rtol, length.out = Nstates)

    ## Number of steps until the solver gives up
    nsteps <- min(.Machine$integer.max, maxsteps * length(times))
    varstep <- method$varstep
    vrb <- FALSE # TRUE forces internal debugging output of the C code

    if (varstep) {                        # methods with variable step size
      out <- .Call("call_rkAuto", as.double(y), as.double(times),
        Func, Initfunc, parms,
        as.integer(Nglobal), rho, as.double(atol),
        as.double(rtol), as.double(tcrit), as.integer(vrb),
        as.double(hmin), as.double(hmax), as.double(hini),
        as.double(rpar), as.integer(ipar), method,
        as.integer(nsteps))
    } else {                              # fixed step methods
      out <- .Call("call_rkFixed", as.double(y), as.double(times),
        Func, Initfunc, parms,
        as.integer(Nglobal), rho,
        as.double(tcrit), as.integer(vrb),
        as.double(hini), as.double(rpar), as.integer(ipar), method,
        as.integer(nsteps))
    }

    ## Names for the outputs
    nm <- c("time",
      if (!is.null(attr(y, "names"))) names(y) else as.character(1:n)
    )

    ## Global outputs
    if (Nglobal > 0) {
      nm  <- c(nm,
        if (!is.null(Nmtot)) Nmtot else as.character((n + 1) : (n + Nglobal))
      )
    }

    ## Column names and state information
    dimnames(out) <- list(NULL, nm)
    istate <- attr(out, "istate")
    attr(out, "type") <- "rk"
    if (verbose) diagnostics(out)
    return(out)
}
