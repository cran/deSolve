### ============================================================================
### Interface to a generalized code for solving explicit variable and fixed
### step ODE solvers of the Runge-Kutta family, see helpfile for details.
### ============================================================================

rk <- function(y, times, func, parms, rtol = 1e-6, atol = 1e-6,
  verbose = FALSE, tcrit = NULL, hmin = 0, hmax = NULL, hini = hmax, ynames=TRUE,
  method = rkMethod("rk45dp7", ... ), maxsteps = 5000,
  dllname = NULL, initfunc=dllname, initpar = parms,
  rpar = NULL,  ipar = NULL, nout = 0, outnames=NULL, forcings=NULL,
  initforc = NULL, fcontrol=NULL, events = NULL, ...) {

    ## Check inputs
    hmax <- checkInput(y, times, func, rtol, atol,
        jacfunc=NULL, tcrit, hmin, hmax, hini, dllname)
    if (hmax == 0) hmax <- .Machine$double.xmax # i.e. practically unlimited

    n <- length(y)

    ## KS -> ThPe: maxsteps/tcrit checks are extra - should they be done in the other?
    if (maxsteps < 0)       stop("maxsteps must be positive")
    if (!is.finite(maxsteps)) maxsteps <- .Machine$integer.max
    if (is.character(method)) method <- rkMethod(method)
    if (is.null(tcrit)) tcrit <- max(times)

    ## Check interpolation order
    if (is.null(method$nknots)) {
      method$nknots <- 5L # fifth order polynomials by default
    } else {
      method$nknots <- as.integer(ceiling(method$nknots))
    }
    nknots <- method$nknots
    if (nknots > 8) {
        warning("Large number of nknots does not make sense.")
    } else if (nknots < 2) {
      # cat("\nMethod without or with disabled interpolation\n")
      method$nknots <- 0L
    } else {
      trange <- diff(range(times))
      ## ensure that we have at least nknots + 2 data points; + 0.5 for safety)
      ## to allow 3rd order polynomial interpolation
      ## for methods without built-in dense output

      if ((is.null(method$d) &                             # has no "dense output"?
        (hmax > 1.0/(nknots + 2.5) * trange))) {           # time steps too large?
        hini <- hmax <- 1.0/(nknots + 2.5) * trange
        if (hmin < hini) hmin <- hini
        cat("\nNote: Method ", method$ID,
            " needs intermediate steps for interpolation\n")
        cat("hmax decreased to", hmax, "\n")
      }
    }

    ## Model as shared object (DLL)?
    Ynames <- attr(y, "names")
    Initfunc <- NULL
    Eventfunc <- NULL
    events <- checkevents(events, times, Ynames, dllname) 
    
    flist    <-list(fmat = 0, tmat = 0, imat = 0, ModelForc = NULL)
    Nstates <- length(y) # assume length of states is correct

    if (is.character(func)) { # function specified in a DLL
      DLL <- checkDLL(func, NULL, dllname,
                      initfunc, verbose, nout, outnames)

      Initfunc <- DLL$ModelInit
      Func     <- DLL$Func
      Nglobal  <- DLL$Nglobal
      Nmtot    <- DLL$Nmtot
      Eventfunc <- events$func
      
      if (! is.null(forcings))
        flist <- checkforcings(forcings, times, dllname, initforc, verbose, fcontrol)

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
          attr(state, "names") <- Ynames
          func(time, state, parms, ...)}
        if (! is.null(events$Type))
          if (events$Type == 2)
            Eventfunc <- function(time, state) {
              attr(state, "names") <- Ynames
              events$func(time, state, parms, ...) 
            }  
      } else {                            # no ynames...
        Func   <- function(time, state, parms)
          func(time, state, parms, ...)
        if (! is.null(events$Type))
          if (events$Type == 2)
            Eventfunc <- function(time, state)  
              events$func(time, state, parms,...) 
      }

      ## Call func once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
      FF <- checkFuncEuler(Func, times, y, parms, rho, Nstates)
      Nglobal <- FF$Nglobal
      Nmtot   <- FF$Nmtot
      
      if (! is.null(events$Type))
        if (events$Type == 2) checkEventFunc(Eventfunc, times, y, rho)
    }

    ## handle length of atol and rtol
    if (Nstates %% length(atol))
      warning("length of atol does not match number of states")
    if (Nstates %% length(rtol))
      warning("length of rtol does not match number of states")

    atol <- rep(atol, length.out = Nstates)
    rtol <- rep(rtol, length.out = Nstates)

    ## Number of steps until the solver gives up
    nsteps  <- min(.Machine$integer.max, maxsteps * length(times))
    varstep <- method$varstep
    vrb <- FALSE # TRUE would force internal debugging output of the C code

    if (varstep) {  # methods with variable step size
      if (is.null(hini)) hini <- hmax
      out <- .Call("call_rkAuto", as.double(y), as.double(times),
        Func, Initfunc, parms, Eventfunc, events,
        as.integer(Nglobal), rho, as.double(atol),
        as.double(rtol), as.double(tcrit), as.integer(vrb),
        as.double(hmin), as.double(hmax), as.double(hini),
        as.double(rpar), as.integer(ipar), method,
        as.integer(nsteps), flist)
    } else {        # fixed step methods
      ## hini=0 for fixed step methods means
      ## that steps in "times" are used as they are
      if (is.null(hini)) hini <- 0
      out <- .Call("call_rkFixed", as.double(y), as.double(times),
        Func, Initfunc, parms, Eventfunc, events,
        as.integer(Nglobal), rho,
        as.double(tcrit), as.integer(vrb),
        as.double(hini), as.double(rpar), as.integer(ipar), method,
        as.integer(nsteps), flist)
    }

    ## saving results
    out <- saveOutrk(out, y, n, Nglobal, Nmtot,
                     iin = c(1,12,13,15), iout = c(1:3, 18))

    attr(out, "type") <- "rk"
    if (verbose) diagnostics(out)
    return(out)
}
