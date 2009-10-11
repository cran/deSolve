### ============================================================================
### Interface to a special code for the classsical Runge-Kutta ODE solver
### with fixed step size and without interpolation, see helpfile for details.
### ============================================================================

rk4 <- function(y, times, func, parms, verbose = FALSE, ynames = TRUE,
  dllname = NULL, initfunc = dllname, initpar = parms,
  rpar = NULL,  ipar = NULL, nout = 0, outnames = NULL, forcings = NULL,
  initforc = NULL, fcontrol = NULL, ...) {

    ## check input
    checkInputEuler(y, times, func, dllname)
    n <- length(y)

    Ynames <- attr(y,"names")
    Initfunc <- NULL
    # KS: moved this upward
    flist    <-list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
    Nstates <- length(y) # assume length of states is correct

    ## Model as shared object (DLL)?
    if (is.character(func)) {
      DLL <- checkDLL(func,NULL,dllname,
                    initfunc,verbose,nout, outnames)

      Initfunc <- DLL$ModelInit
      Func     <- DLL$Func
      Nglobal  <- DLL$Nglobal
      Nmtot    <- DLL$Nmtot

      if (! is.null(forcings))
        flist <- checkforcings(forcings,times,dllname,initforc,verbose,fcontrol)

      rho <- NULL
      if (is.null(ipar)) ipar<-0
      if (is.null(rpar)) rpar<-0

    } else {
      initpar <- NULL # parameter initialisation not needed if function is not a DLL
      rho <- environment(func)
      ## func and jac are overruled, either including ynames, or not
      ## This allows to pass the "..." arguments and the parameters
      if(ynames) {
        Func   <- function(time, state, parms) {
          attr(state, "names") <- Ynames
          func   (time,state,parms,...)
        }
      } else {                            # no ynames...
        Func   <- function(time, state, parms)
          func   (time, state, parms,...)
      }

      ## Call func once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
      FF <- checkFuncEuler(Func,times,y,parms,rho,Nstates)
      Nglobal<-FF$Nglobal
      Nmtot <- FF$Nmtot
    }
    vrb <- FALSE # TRUE forces internal debugging output of the C code

    ## KS -> Thomas: still need to pass flist
    ## the CALL to the integrator
    out <- .Call("call_rk4", as.double(y), as.double(times),
        Func, Initfunc, parms, as.integer(Nglobal), rho, as.integer(vrb),
        as.double(rpar), as.integer(ipar), flist)

    out <- saveOutrk(out, y, n, Nglobal, Nmtot,
                     iin = c(1, 12, 13, 15), iout=c(1:3, 18))

    attr(out, "type") <- "rk"
    if (verbose) diagnostics(out)
    return(out)
}
