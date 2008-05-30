

### lsodar -- solves ordinary differential equation systems 
### This a variant version of the lsode integrator.  It differs from it
### in two ways:
### (a) It switches automatically between stiff and nonstiff methods.
### This means that the user does not have to determine whether the
### problem is stiff or not, and the solver will automatically choose the
### appropriate method.  It always starts with the nonstiff method.
### This is similar to lsoda.
### (b) It finds the root of at least one of a set of constraint
### functions g(i) of the independent and dependent variables.
### It finds only those roots for which some g(i), as a function
### of t, changes sign in the interval of integration.
### It then returns the solution at the root, if that occurs
### sooner than the specified stop condition, and otherwise returns
### the solution according the specified stop condition.


lsodar <- function(y, times, func, parms, rtol=1e-6, atol=1e-6,
	jacfunc=NULL, jactype = "fullint", rootfunc=NULL,
  verbose=FALSE,  nroot = 0, tcrit = NULL, hmin=0, hmax=NULL, hini=0, ynames=TRUE, 
  maxordn = 12, maxords = 5, bandup = NULL, banddown = NULL, 
  maxsteps = 5000, dllname=NULL,initfunc=dllname, initpar=parms,
   rpar=NULL, ipar=NULL, nout=0, outnames=NULL,...)   
{
### check input
    if (!is.numeric(y))        stop("`y' must be numeric")
    n <- length(y)
    if (! is.null(times)&&!is.numeric(times))
        stop("`times' must be NULL or numeric")
    if (!is.function(func) && !is.character(func))
      stop("`func' must be a function or character vector")
    if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
      stop("You need to specify the name of the dll or shared library where func can be found (without extension)")
    if (!is.numeric(rtol))    stop("`rtol' must be numeric")
    if (!is.numeric(atol))    stop("`atol' must be numeric")
    if (!is.null(tcrit) & !is.numeric(tcrit)) stop("`tcrit' must be numeric")
    if (!is.null(jacfunc) && !(is.function(jacfunc) || is.character(jacfunc)))
        stop("`jacfunc' must be a function or character vector")
    if (!is.null(rootfunc) && !(is.function(rootfunc) || is.character(rootfunc)))
        stop("`rootfunc' must be a function or character vector")
    if (length(atol) > 1 && length(atol) != n)
        stop("`atol' must either be a scaler, or as long as `y'")
    if (length(rtol) > 1 && length(rtol) != n)
        stop("`rtol' must either be a scaler, or as long as `y'")
    if (!is.numeric(hmin))    stop("`hmin' must be numeric")
    if (hmin < 0) stop ("`hmin' must be a non-negative value")
    if (is.null(hmax))
       hmax <- ifelse (is.null(times), 0, max(abs(diff(times))))
    if (!is.numeric(hmax))    stop("`hmax' must be numeric")
    if (hmax < 0) stop ("`hmax' must be a non-negative value")
    if (hini < 0)             stop("`hini' must be a non-negative value")
    if (!is.numeric(maxordn)) stop("`maxordn' must be numeric")
    if(maxordn < 1 || maxordn > 12) stop("`maxord' must be >1 and <=12")
    if (!is.numeric(maxords)) stop("`maxords' must be numeric")
    if(maxords < 1 || maxords > 5)stop("`maxords' must be >1 and <=5")

### Jacobian, method flag
       if (jactype == "fullint" ) jt <- 2 # full jacobian, calculated internally
  else if (jactype == "fullusr" ) jt <- 1 # full jacobian, specified by user function
  else if (jactype == "bandusr" ) jt <- 4 # banded jacobian, specified by user function
  else if (jactype == "bandint" ) jt <- 5 # banded jacobian, specified internally
  else stop("jactype must be one of fullint, fullusr, bandusr or bandint")

  # check other specifications depending on jacobian  
    if (jt %in% c(4,5) && is.null(bandup))   
        stop("lsodar: bandup must be specified if banded jacobian")            
    if (jt %in% c(4,5) && is.null(banddown)) 
        stop("lsodar: banddown must be specified if banded jacobian")            
    if (is.null(banddown)) banddown <-1
    if (is.null(bandup  )) bandup   <-1  

      if (jt %in% c(1,4) && is.null(jacfunc))
  stop ("lsoda: cannot perform integration: *jacfunc* NOT specified; either specify *jacfunc* or change *jactype*")

### model and jacobian function

    Ynames <- attr(y,"names")
    JacFunc  <- NULL
    RootFunc <- NULL

    ModelInit <- NULL
    if(!is.null(dllname))
    {
        if (is.loaded(initfunc, PACKAGE = dllname,
           type = "") || is.loaded(initfunc, PACKAGE = dllname,
            type = "Fortran")) 
        { ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
        } else if (initfunc != dllname && ! is.null(initfunc))
            stop(paste("cannot integrate: initfunc not loaded ",initfunc))        
     }

    ## If func is a character vector, then
    ## copy its value to funcname,  
    ## check to make sure it describes
    ## a function in a loaded dll
    if (is.character(func)) {
      funcname <- func
      ## get the pointer and put it in func
      if(is.loaded(funcname, PACKAGE = dllname)) {
        Func <- getNativeSymbolInfo(funcname, PACKAGE = dllname)$address
        } else stop(paste("cannot integrate: dyn function not loaded",funcname))

      ## Is there a jacobian?
      if (!is.null(jacfunc)) {
        if (!is.character(jacfunc))
          stop("If 'func' is dynloaded, so must 'jacfunc' be")
        jacfuncname <- jacfunc
        if(is.loaded(jacfuncname, PACKAGE = dllname))
            {JacFunc <- getNativeSymbolInfo(jacfuncname, PACKAGE = dllname)$address
            } else stop(paste("cannot integrate: jac function not loaded ",jacfunc))
        } 
      ## Is there a root function?
      if (!is.null(rootfunc)) {
        if (!is.character(rootfunc))
          stop("If 'func' is dynloaded, so must 'rootfunc' be")
        rootfuncname <- rootfunc
        if(is.loaded(rootfuncname, PACKAGE = dllname))
            {RootFunc <- getNativeSymbolInfo(rootfuncname, PACKAGE = dllname)$address
            } else stop(paste("cannot integrate: root function not loaded ",rootfunc))
        if (nroot == 0) stop("if rootfunc is specified in a dll, then nroot should be > 0")
        }

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
      rho <- NULL
      if (is.null(ipar)) ipar<-0
      if (is.null(rpar)) rpar<-0

    } else {
      if(is.null(initfunc)) initpar <- NULL # parameter initialisation not needed if function is not a DLL    
      rho <- environment(func)
      # func and jac are overruled, either including ynames, or not
      # This allows to pass the "..." arguments and the parameters

        if(ynames)
        {
         Func    <- function(time,state) 
         { attr(state,"names") <- Ynames 
           func   (time,state,parms,...)[1]}   
         
         Func2   <- function(time,state) 
         { attr(state,"names") <- Ynames
           func   (time,state,parms,...)}    
         
         JacFunc <- function(time,state) 
         { attr(state,"names") <- Ynames
           jacfunc(time,state,parms,...)}    

         RootFunc <- function(time,state) 
         { attr(state,"names") <- Ynames 
           rootfunc(time,state,parms,...)}   

        } else {                            # no ynames...
         Func    <- function(time,state) 
           func   (time,state,parms,...)[1] 
        
         Func2   <- function(time,state) 
           func   (time,state,parms,...)    
         
         JacFunc <- function(time,state) 
           jacfunc(time,state,parms,...)    

         RootFunc <- function(time,state) 
           rootfunc(time,state,parms,...)   

        }
        
      ## Call func once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
      tmp <- eval(Func2(times[1], y), rho) 

      if (!is.list(tmp)) stop("Model function must return a list\n")
      if (length(tmp[[1]]) != length(y))
        stop(paste("The number of derivatives returned by func() (",
                   length(tmp[[1]]),
                   "must equal the length of the initial conditions vector (",
                   length(y),")",sep=""))
                   
      # use "unlist" here because some output variables are vectors/arrays
        Nglobal <- if (length(tmp) > 1)   
            length(unlist(tmp[-1]))  else 0   
        Nmtot <- attr(unlist(tmp[-1]),"names")

      # and for rootfunc
        if (! is.null(rootfunc))
        { tmp2 <- eval(rootfunc(times[1],y,parms,...), rho) 
          if (!is.vector(tmp2)) stop("root function must return a vector\n")        
          nroot <- length(tmp2)
        } else nroot = 0
      if (jt %in% c(1,4))
      {tmp <- eval(JacFunc(times[1], y), rho) 
       if (!is.matrix(tmp)) stop("Jacobian function must return a matrix\n")
       dd <- dim(tmp)
       if((jt ==4 && dd != c(bandup+banddown+1,n)) ||
          (jt ==1 && dd != c(n,n))) stop("Jacobian dimension not ok") 
      } 

    }

### work arrays iwork, rwork
# length of rwork and iwork 

  if(jt %in% c(1,2)) lmat <- n^2+2 else
  if(jt %in% c(4,5)) lmat <- (2*banddown+bandup+1)*n+2
  lrn = 20+n*(maxordn+1)+ 3*n +3*nroot       # length in case non-stiff method
  lrs = 20+n*(maxords+1)+ 3*n +lmat+3*nroot  # length in case stiff method
  lrw = max(lrn,lrs)                         # actual length: max of both
  liw = 20 + n

# only first 20 elements passed to solver; other will be allocated in C-code  
  iwork <- vector("integer",20)
  rwork <- vector("double",20)
  rwork[] <- 0.
  iwork[] <- 0

  iwork[1] <- banddown
  iwork[2] <- bandup
  iwork[6] <- maxsteps
  if (maxordn != 12) iwork[8] <- maxordn
  if (maxords != 5)  iwork[9] <- maxords

  if(! is.null(tcrit)) rwork[1] <- tcrit
  rwork[5] <- hini
  rwork[6] <- hmax
  rwork[7] <- hmin

#  the task to be performed.
  if (! is.null(times))
      itask <- ifelse (is.null (tcrit), 1,4) else      # times specified
      itask <- ifelse (is.null (tcrit), 2,5)           # only one step
  if(is.null(times)) times<-c(0,1e8)

# print to screen...
  if (verbose) 
   {
   print("--------------------")   
   print("time settings")   
   print("--------------------")   
   if (itask==1)print("normal computation of output values of y(t) at t = TOUT") else
   if (itask==2)print("take one step only and return.")                          else
   if (itask==3)print("istop at the first internal mesh point at or beyond t = TOUT and return. ")  else
   if (itask==4)print("normal computation of output values of y(t) at t = TOUT but without overshooting t = TCRIT.") else
   if (itask==5)print("take one step, without passing TCRIT, and return.")
   print("--------------------") 
   print("Integration settings")     
   print("--------------------") 
   if (is.character(func)) print(paste("model function a DLL: ",func)) else
                           print("model function an R-function: ")
   if (is.character(jacfunc)) print(paste ("jacobian specified as a DLL: ",jacfunc)) else
   if (!is.null(jacfunc))     print("jacobian specified as an R-function: ") else
                              print("jacobian not specified")
   print("--------------------")   
   
   }

### calling solver    
    storage.mode(y) <- storage.mode(times) <- "double"
    IN <-4
    # KS: Func and Jacfunc rather than func and jacfunc; initpar rather than par
    # added itask,rwork, iwork,jt,Nglobal    removed hmin,hmax
    out <- .Call("call_lsoda",y,times,Func,as.double(initpar),
                 rtol, atol, rho, tcrit, JacFunc, ModelInit,  
                 as.integer(verbose), as.integer(itask), as.double(rwork),
                 as.integer(iwork), as.integer(jt),as.integer(Nglobal),
                 as.integer(lrw),as.integer(liw),as.integer(IN),RootFunc,
                 as.integer(nroot), as.double (rpar), as.integer(ipar),
                 as.integer(0), PACKAGE="deSolve")   #  

### saving results    

    istate <- attr(out,"istate")
    iroot  <- attr(out, "iroot")   
    rstate <- attr(out, "rstate")    
    nm <- c("time",
            if (!is.null(attr(y,"names"))) names(y) else as.character(1:n))
    if (Nglobal > 0) {
       if (!is.character(func)) {                  # if a DLL: already done...    
        out2 <- matrix( nrow=Nglobal, ncol=ncol(out))
        for (i in 1:ncol(out2)) {
          y <- out[-1,i]
          names(y) <- nm[-1]
          out2[, i] <- unlist(Func2(out[1, i], y)[-1])  # KS: Func2 rather than func
          }
        out <- rbind(out,out2)
        }  # end !is.character func
        nm <- c(nm,
                if (!is.null(Nmtot)) Nmtot else
                                     as.character((n+1) : (n + Nglobal)))
    }
    attr(out,"istate")  <- istate
    attr(out, "rstate") <- rstate        
    attr(out, "iroot")  <- iroot    

    dimnames(out) <- list(nm,NULL)
    
    if (verbose)
    {
      print("--------------------")    
      print("lsodar return code")
      print("--------------------")      
      idid <- istate[1]
      print(paste("istate = ",idid))

      if (idid == 2) print(" LSODAR was successful but no root was found") else
      if (idid == 3) print(" LSODAR was successful and a root was found before reaching the end") else
      if (idid == -1) print(" excess work done on this call. (Perhaps wrong jacobian type)") else
      if (idid == -2) print(" excess accuracy requested. (Tolerances too small.)") else
      if (idid == -3) print(" illegal input detected. (See printed message.)") else
      if (idid == -4) print(" repeated error test failures. (Check all input.)") else
      if (idid == -5) print(" repeated convergence failures. (Perhaps bad Jacobian supplied or wrong choice of jt or tolerances.)") else
      if (idid == -6) print(" error weight became zero during problem. (Solution component i vanished, and ATOL or ATOL(i) = 0.)") else
      if (idid == -7) print(" work space insufficient to finish (see messages)")

      print("--------------------")
      print("ISTATE values")
      print("--------------------")
      df <- c( " istate, the return code",
               " The number of steps taken for the problem so far.",
               " The number of function evaluations for the problem so far.",
               " The number of Jacobian evaluations and LU decompositions so far.",
               " The method order last used (successfully).",
               " The order to be attempted on the next step.",
               " if istate=-4,-5: the index of the component with the largest error vector",
               " The length of rwork actually required.",
               " The length of iwork actually required.",
               " The method indicator for the last succesful step, 1=adams (nonstiff), 2= bdf (stiff)",
               " The current method indicator to be attempted on th next step, 1=adams (nonstiff), 2= bdf (stiff)")

      ii <- c(1,12:21)
       print(data.frame(mess=df, val=istate[ii]))
       print("--------------------")       
       print("RSTATE values")
       print("--------------------")       
      df <- c( " The step size in t last used (successfully).",
    " The step size to be attempted on the next step.",
    " The current value of the independent variable which the solver has actually reached",
    " Tolerance scale factor, greater than 1.0, computed when a request for too much accuracy was detected",
    " the value of t at the time of the last method switch, if any.")
       print(data.frame(mess=df, val=rstate[1:5]))

    }

    t(out)
}
