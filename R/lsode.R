
### lsode -- solves ordinary differential equation systems 
### It differs from lsoda in that the user has to specify whether or not
### the problem is stiff and choose the appropriate method.
### It is very similar to vode, except for some implementation details.
### More specifically, in vode it is possible to choose whether or not a copy
### of the Jacobian is saved for reuse in the corrector iteration algorithm; 
### In lsoda, a copy is not kept; this requires less memory but may be slightly
### slower.


lsode <- function(y, times, func, parms, rtol=1e-6, atol=1e-6, 
  jacfunc=NULL, jactype = "fullint", mf = NULL, verbose=FALSE, 
  tcrit = NULL, hmin=0, hmax=NULL, hini=0, ynames=TRUE, 
  maxord=NULL, bandup=NULL, banddown=NULL, maxsteps=5000, 
  dllname=NULL,initfunc=dllname, initpar=parms, 
  rpar=NULL, ipar=NULL, nout=0, outnames=NULL,...)  
{
### check input
  if (!is.numeric(y))     stop("`y' must be numeric")
  n <- length(y)
  if (! is.null(times)&&!is.numeric(times))
    stop("`times' must be NULL or numeric")
  if (!is.function(func) && !is.character(func))
    stop("`func' must be a function or character vector")
  if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
    stop("You need to specify the name of the dll or shared library where func can be found (without extension)")
  if (!is.numeric(rtol))  stop("`rtol' must be numeric")
  if (!is.numeric(atol))  stop("`atol' must be numeric")
  if (!is.null(tcrit) & !is.numeric(tcrit)) stop("`tcrit' must be numeric")
  if (!is.null(jacfunc) && !(is.function(jacfunc) || is.character(jacfunc)))
    stop("`jacfunc' must be a function or character vector")
  if (length(atol) > 1 && length(atol) != n)
    stop("`atol' must either be a scaler, or as long as `y'")
  if (length(rtol) > 1 && length(rtol) != n)
    stop("`rtol' must either be a scaler, or as long as `y'")
  if (!is.numeric(hmin))   stop("`hmin' must be numeric")
  if (hmin < 0)            stop("`hmin' must be a non-negative value")
  if (is.null(hmax))
    hmax <- if (is.null(times)) 0 else max(abs(diff(times)))
  if (!is.numeric(hmax))   stop("`hmax' must be numeric")
  if (hmax < 0)            stop("`hmax' must be a non-negative value")
  if (hmax == Inf)  hmax <- 0
  if (hini < 0)            stop("`hini' must be a non-negative value")
  if (!is.null(maxord)) 
    if(maxord < 1) stop("`maxord' must be >1")

### Jacobian, method flag
  if (is.null(mf)){ 
    if (jactype == "fullint" ) imp <- 22 # full jacobian, calculated internally
    else if (jactype == "fullusr" ) imp <- 21 # full jacobian, specified by user function
    else if (jactype == "bandusr" ) imp <- 24 # banded jacobian, specified by user function
    else if (jactype == "bandint" ) imp <- 25 # banded jacobian, specified internally
    else stop("jactype must be one of fullint, fullusr, bandusr or bandint if mf not specified")
  } else imp <- mf

  if (! imp %in% c(10:15, 20:25)) 
    stop ("lsode: cannot perform integration: method flag mf not allowed")
  
                                        # check other specifications depending on jacobian
  miter <- imp%%10 
  if (miter %in% c(1,4) & is.null(jacfunc)) 
    stop ("lsode: cannot perform integration: *jacfunc* NOT specified; either specify *jacfunc* or change *jactype* or *mf*")
  meth <- abs(imp)%/%10                # basic linear multistep method

  if (is.null (maxord)) maxord <- if (meth==1) 12 else 5
  if (meth==1 && maxord > 12) stop ("lsode: maxord too large: should be <= 12")
  if (meth==2 && maxord > 5 ) stop ("lsode: maxord too large: should be <= 5")
  if (miter %in% c(4,5) && is.null(bandup))   
    stop("lsode: bandup must be specified if banded jacobian")
  if (miter %in% c(4,5) && is.null(banddown)) 
    stop("lsode: banddown must be specified if banded jacobian")
  if (is.null(banddown)) banddown <-1
  if (is.null(bandup  )) bandup   <-1  

### model and jacobian function  
  JacFunc <- NULL
  Ynames <- attr(y,"names")

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
  ## copy its value to funcname 
  ## check to make sure it describes
  ## a function in a loaded dll
  if (is.character(func)) {
    funcname <- func
    ## get the pointer and put it in func
                                        # KS: changed that...
    if(is.loaded(funcname, PACKAGE = dllname)) {
      Func <- getNativeSymbolInfo(funcname, PACKAGE = dllname)$address
    } else stop(paste("cannot integrate: dyn function not loaded",funcname))

    ## Finally, is there a jacobian?
    if (!is.null(jacfunc)) {
      if (!is.character(jacfunc))
        stop("If 'func' is dynloaded, so must 'jacfunc' be")
      jacfuncname <- jacfunc
      if(is.loaded(jacfuncname, PACKAGE = dllname))
        {JacFunc <- getNativeSymbolInfo(jacfuncname, PACKAGE = dllname)$address
       } else stop(paste("cannot integrate: jac function not loaded ",jacfunc))
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
      } else {                          # no ynames...
        Func    <- function(time,state) 
          func   (time,state,parms,...)[1] 
        
        Func2   <- function(time,state) 
          func   (time,state,parms,...)    
         
        JacFunc <- function(time,state) 
          jacfunc(time,state,parms,...)    
      }
        
    ## Call func once to figure out whether and how many "global"
    ## results it wants to return and some other safety checks
    tmp <- eval(Func2(times[1], y), rho) 
    if (!is.list(tmp)) stop("Model function must return a list\n")
    if (length(tmp[[1]]) != length(y))
      stop(paste("The number of derivatives returned by func() (",
                 length(tmp[[1]]),
                 ") must equal the length of the initial conditions vector (",
                 length(y),")",sep=""))

                                        # use "unlist" here because some output variables are vectors/arrays
    Nglobal <- if (length(tmp) > 1)   
      length(unlist(tmp[-1]))  else 0   
    Nmtot <- attr(unlist(tmp[-1]),"names")

    if (miter %in% c(1,4))
      {tmp <- eval(JacFunc(times[1], y), rho) 
       if (!is.matrix(tmp)) stop("Jacobian function must return a matrix\n")
       dd <- dim(tmp)
       if((miter ==4 && dd != c(bandup+banddown+1,n)) ||
          (miter ==1 && dd != c(n,n))) stop("Jacobian dimension not ok") 
     } 
  }                                                                                


### work arrays iwork, rwork
                                        # length of rwork and iwork 

  lrw = 20+n*(maxord+1)+3*n
  if(miter %in% c(1,2) ) lrw = lrw + 2*n*n+2
  if(miter ==3)          lrw = lrw + n+2
  if(miter %in% c(4,5) ) lrw = lrw + (2*banddown+ bandup+1)*n+2

  liw   <- if (miter %in% c(0,3)) 20 else 20+n

                                        # only first 20 elements passed to solver; other will be allocated in C-code
  iwork <- vector("integer",20)
  rwork <- vector("double",20)
  rwork[] <- 0.
  iwork[] <- 0

  iwork[1] <- banddown
  iwork[2] <- bandup
  iwork[5] <- maxord
  iwork[6] <- maxsteps

  if(!is.null(tcrit)) rwork[1] <- tcrit
  rwork[5] <- hini
  rwork[6] <- hmax
  rwork[7] <- hmin

                                        # the task to be performed.
  itask <- if (! is.null(times)) {
    if (is.null (tcrit)) 1 else 4
  } else  {                             # times specified
    if (is.null (tcrit)) 2 else 5       # only one step
  }
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
      print("integration method")
      print("--------------------")   
      df   <- c("method flag, mf","meth","miter")
      vals <- c(imp,  meth, miter)
      txt  <- "mf= (10*meth + miter)"
   
      if (meth==1)txt<-c(txt,"the basic linear multistep method: the implicit Adams method")                                             else 
      if (meth==2)txt<-c(txt,"the basic linear multistep method: based on backward differentiation formulas")

      if (miter==0)txt<-c(txt,"functional iteration (no Jacobian matrix is involved")                                                    else 
      if (miter==1)txt<-c(txt,"chord iteration with a user-supplied full (NEQ by NEQ) Jacobian")                                         else 
      if (miter==2)txt<-c(txt,"chord iteration with an internally generated full Jacobian, (NEQ extra calls to F per df/dy value)")      else 
      if (miter==3)txt<-c(txt,"chord iteration with an internally generated diagonal Jacobian (1 extra call to F per df/dy evaluation)") else 
      if (miter==4)txt<-c(txt,"chord iteration with a user-supplied banded Jacobian")                                                    else 
      if (miter==5)txt<-c(txt,"chord iteration with an internally generated banded Jacobian (using ML+MU+1 extra calls to F per df/dy evaluation)") 
                                                                                                                                         
      print(data.frame(parameter=df, value=vals,message=txt))
    } 

### calling solver
  storage.mode(y) <- storage.mode(times) <- "double"
  IN <-2
    
  out <- .Call("call_lsoda",y,times,Func,initpar,
               rtol, atol, rho, tcrit, JacFunc, ModelInit,  
               as.integer(verbose), as.integer(itask), as.double(rwork),
               as.integer(iwork), as.integer(imp),as.integer(Nglobal),
               as.integer(lrw),as.integer(liw),as.integer(IN),
               NULL, as.integer(0), as.double (rpar), as.integer(ipar),
               as.integer(0), PACKAGE="deSolve")

### saving results    

  istate <- attr(out,"istate")
  rstate <- attr(out, "rstate")    
  nm <- c("time",
          if (!is.null(attr(y,"names"))) names(y) else as.character(1:n))
  if (Nglobal > 0) {
    if (!is.character(func)) {         # if a DLL: already done...    
      out2 <- matrix( nrow=Nglobal, ncol=ncol(out))
      for (i in 1:ncol(out2)) {
        y <- out[-1,i]
        names(y) <- nm[-1]
        out2[, i] <- unlist(Func2(out[1, i], y)[-1]) # KS: Func2 rather than func
      }
      out <- rbind(out,out2)
    }                                   # end !is.character func
    nm <- c(nm,
            if (!is.null(Nmtot)) Nmtot else
            as.character((n+1) : (n + Nglobal)))
  }
  attr(out,"istate") <- istate
  attr(out, "rstate") <- rstate        

  dimnames(out) <- list(nm,NULL)
    
  if (verbose)
    {
      print("--------------------")    
      print("lsode return code")
      print("--------------------")      
      idid <- istate[1]
      print(paste("istate = ",idid))

      if (idid == 2) print(" lsode was successful") else
      if (idid == -1) print(" excess work done on this call. (Perhaps wrong jacobian type)") else
      if (idid == -2) print(" excess accuracy requested. (Tolerances too small.)") else
      if (idid == -3) print(" illegal input detected. (See printed message.)") else
      if (idid == -4) print(" repeated error test failures. (Check all input.)") else
      if (idid == -5) print(" repeated convergence failures. (Perhaps bad Jacobian supplied or wrong choice of jactype or tolerances.)") else
      if (idid == -6) print(" error weight became zero during problem. (Solution component i vanished, and ATOL or ATOL(i) = 0.)") 

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
              " The length of iwork actually required."              )

      ii <- c(1,12:19)
      print(data.frame(mess=df, val=istate[ii]))

      print("--------------------")
      print("RSTATE values")
      print("--------------------")
      df <- c( " The step size in t last used (successfully).",
              " The step size to be attempted on the next step.",
              " The current value of the independent variable which the solver has actually reached",
              " Tolerance scale factor, greater than 1.0, computed when a request for too much accuracy was detected" 
              )
      print(data.frame(mess=df, val=rstate[1:4]))

    }

  t(out)
}
