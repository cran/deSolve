### ============================================================================
### Interface to a special code for Euler's ODE solver
### with fixed step size and without interpolation, see helpfile for details.
### ============================================================================

euler <- function(y, times, func, parms, verbose = FALSE, ynames=TRUE,
  dllname = NULL, initfunc=dllname, initpar = parms,
  rpar = NULL,  ipar = NULL, nout = 0, outnames=NULL, ...) {

    ## check input
    if (!is.numeric(y))  stop("`y' must be numeric")
    n <- length(y)
    if (! is.null(times)&&!is.numeric(times))
        stop("`times' must be NULL or numeric")
    if (!is.function(func) && !is.character(func))
      stop("`func' must be a function")

    if (!is.numeric(y))     stop("`y' must be numeric")
    if (!is.numeric(times)) stop("`times' must be numeric")

    if (!is.function(func) && !is.character(func))
      stop("`func' must be a function or character vector")
    if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
      stop("You need to specify the name of the dll or shared library where func can be found (without extension)")


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
        } else stop(paste("cannot integrate: dyn function not loaded",funcname))

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
      ## ThPe:
      Nstates <- length(y) # assume length of states is correct
      rho <- NULL
      if (is.null(ipar)) ipar <- 0
      if (is.null(rpar)) rpar <- 0
    } else {
      initpar <- NULL # parameter initialisation not needed if function is not a DLL
      rho <- environment(func)
      # func and jac are overruled, either including ynames, or not
      # This allows to pass the "..." arguments and the parameters
      if(ynames) {
        Func   <- function(time,state,parms) {
          attr(state,"names") <- Ynames
          func   (time,state,parms,...)
        }
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

    ## the CALL to the integrator
    out <- .Call("call_euler", as.double(y), as.double(times),
        Func, Initfunc, parms, as.integer(Nglobal), rho, as.integer(verbose),
        as.double(rpar), as.integer(ipar))

    nm <- c("time",
      if (!is.null(attr(y, "names"))) names(y) else as.character(1:n)
    )

    ## global outputs
    if (Nglobal > 0) {
      nm  <- c(nm,
        if (!is.null(Nmtot)) Nmtot else as.character((n + 1) : (n + Nglobal))
      )
    }
    ## column names and state information
    dimnames(out) <- list(NULL, nm)
    istate <- attr(out, "istate")
    if (!is.null(istate) && istate[1] == -1)

    if (verbose) diagnostics(out)

    attr(out, "type")   <- "rk"
    out
}
