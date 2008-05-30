
### daspk -- solves differential algebraic and ordinary differential equation 
###          systems defined in res (DAE) or func (ODE)
###          and outputs values for the times in `times'
###          on input, y and dy contains the initial values of the state 
###          variables and rates of changes for times[1]
###          parms is a vector of parameters for func.  They should not
###          change during the integration. `rtol', and `atol'
###          are, respectively, the relative tolerance parameter, and the
###          absolute tolerance parameter.  `atol' may be scaler or vector.
###          `rtol' is a scalar or a vector.
###
###          The return value is a matrix whose rows correspond to the values
###          in `times', and columns to the elements of `y'.
###
###          'res' may be a string instead of an R function.  If
###          so, then if jacres is not NULL, it must be a character string
###          as well.  In these cases, 'res' is the name
###          of a function to be found in the dll named 'dllname' 
###          (without extension). 'jacres' points to the name of the jacobian.




daspk          <- function(y,               # state variables
                          times,            # time sequence for which output is wanted  
                          func=NULL,        # function that returns rate of change
                          parms,            # other parameters passed to func and jac                        
                          dy=NULL,          # rate of change
                          res=NULL,         # function that returns the residual G(t,y,y') of the differential/algebraic system.
                          nalg=0,           # number of algebraic equations
                          rtol=1e-6,        # relative tolerance  
                          atol=1e-8,        # absolute tolerance  
                          jacfunc=NULL,     # jacobian, if ode specified by func
                          jacres=NULL,      # jacobian, if dae or ode specified by res
                          jactype = "fullint",  # jacobian
                          estini = NULL,    # etsimate iniital conditions
                          verbose=FALSE,    # 
                          tcrit = NULL,
                          hmin=0,           # minimum step size
                          hmax=NULL,        # maximum step size                                           
                          hini=0,           # initial step size
                          ynames=TRUE,      # if false: names of state variables are not passed to function func
                          maxord =5,        # maximal method order; reduce to save storage (<=5)
                          bandup=NULL,         # upper band
                          banddown=NULL,       # lower band
                          maxsteps=5000,    # maximal number of steps in one call to the solver - multipe of 500!
                          dllname=NULL,     # in case initialiser func 
                          initfunc=dllname,                          
                          initpar=parms,    # to initialise common block/global vairaibles
                          rpar=NULL, ipar=NULL,
                          nout=0,           # only if dllname is present: the number of output variables
                          outnames=NULL,    # only if dllnmae is present: the names of output variables...                          
                          ...)              # accessory parameters passed to ??
{
### check input 
    if (!is.numeric(y))
        stop("`y' must be numeric")
    n <- length(y)
    if (! is.null(times)&&!is.numeric(times))
        stop("`times' must be NULL or numeric")
    if (is.null(res) && is.null(func))    
        stop("either `func' or 'res' must be specified")
    if (!is.null(res) && !is.null(func))    
        stop("either `func' OR 'res' must be specified, not both")
    if (!is.null(jacres) && !is.null(jacfunc))    
        stop("either `jacfunc' OR 'jacres' must be specified, not both")
    if (!is.null(func) && !is.function(func))
        stop("`func' must be a function or NULL")
    if (!is.null(res) && !is.function(res) && !is.character(res))
        stop("`res' must be NULL, a function or character vector")
    if (is.character(res) && (is.null(dllname) || !is.character(dllname)))
        stop("You need to specify the name of the dll or shared library where res can be found (without extension)")
    if (!is.numeric(rtol))           
        stop("`rtol' must be numeric")
    if (!is.numeric(atol))           
        stop("`atol' must be numeric")
    if (!is.null(tcrit) & !is.numeric(tcrit))
        stop("`tcrit' must be numeric")
    if (!is.null(jacfunc) && !(is.function(jacfunc) ))
        stop("`jacfunc' must be a function or NULL")
    if (!is.null(jacres) && !(is.function(jacres) || is.character(jacres)))
        stop("`jacres' must be a function or character vector")
    if (length(atol) > 1 && length(atol) != n)
        stop("`atol' must either be a scaler, or as long as `y'")
    if (length(rtol) > 1 && length(rtol) != n)
        stop("`rtol' must either be a scaler, or as long as `y'")
    if (!is.numeric(hmin))        
        stop("`hmin' must be numeric")
    if (hmin < 0)                 
        stop("`hmin' must be a non-negative value")
    if (is.null(hmax))
       hmax <- ifelse (is.null(times), 0, max(abs(diff(times))))
    if (!is.numeric(hmax))      
        stop("`hmax' must be numeric")
    if (hmax < 0)             
        stop("`hmax' must be a non-negative value")
    if (hini < 0)             
        stop("`hini' must be a non-negative value")
    if (!is.numeric(maxord))  
        stop("`maxord' must be numeric") 
    if(maxord < 1 || maxord > 5) 
        stop("`maxord' must be >1 and <=5")
    #max number of iterations ~ maxstep; a multiple of 500
    maxIt <- max(1,(maxsteps+499)%/%500)

### Jacobian, method flag
    if (jactype == "fullint" ) 
        imp <- 22 # full jacobian, calculated internally
    else if (jactype == "fullusr" )
        imp <- 21 # full jacobian, specified by user function
    else if (jactype == "bandusr" ) 
        imp <- 24 # banded jacobian, specified by user function
    else if (jactype == "bandint" ) 
        imp <- 25 # banded jacobian, specified internally
    else stop("jactype must be one of fullint, fullusr, bandusr or bandint")

    if (imp %in% c(24,25) && is.null(bandup))   
        stop("daspk: bandup must be specified if banded jacobian")            
    if (imp %in% c(24,25) && is.null(banddown)) 
        stop("daspk: banddown must be specified if banded jacobian")            

#  if (miter == 4) jacobian should have empty banddown empty rows-vode+daspk only! 
    if (imp == 24) 
        erow<-matrix(nc=n,nr=banddown,0) 
    else erow<-NULL
    
    if (is.null(banddown)) 
       banddown <-1
    if (is.null(bandup  )) 
       bandup   <-1  

    if (is.null(dy))      
       dy <- rep(0,n)
    if (!is.numeric(dy))  
       stop("`dy' must be numeric")

### model and jacobian function
    Ynames  <- attr(y,"names")
    dYnames <- attr(dy,"names")
    Res     <- NULL
    JacRes  <- NULL
    PsolFunc<- NULL
    
    ModelInit <- NULL
    if(!is.null(dllname))
    {
        if (is.loaded(initfunc, PACKAGE = dllname, type = "") || 
           is.loaded(initfunc, PACKAGE = dllname, type = "Fortran"))
        { ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
        } else if (initfunc != dllname && ! is.null(initfunc))
            stop(paste("cannot integrate: initfunc not loaded ",initfunc))        
     }

    ## If res is a character vector, then
    ## check to make sure it describes
    ## a function in a loaded dll
    psolfunc <- NULL  # not yet supported 
    if (is.character(res)) {
        resname <- res
        if(is.loaded(resname, PACKAGE = dllname)) {
        Res <- getNativeSymbolInfo(resname, PACKAGE = dllname)$address  
        } else stop(paste("cannot integrate: res function not loaded",resname))

#        if (is.null(kryltype))
#        {
        if (!is.null(jacres))   {
            if (!is.character(jacres))
                stop("If 'res' is dynloaded, so must 'jacres' be")
            jacname <- jacres
            if(is.loaded(jacname, PACKAGE = dllname))
            {JacRes <- getNativeSymbolInfo(jacname, PACKAGE = dllname)$address
            } else stop(paste("cannot integrate: jacobian function jacres not loaded ",jacres))
                                }
        if (!is.null(psolfunc)) {
            if (!is.character(psolfunc))
                stop("If 'res' is dynloaded, so must 'psolfunc' be")
            if(is.loaded(psolfunc, PACKAGE = dllname))
            {PsolFunc <- getNativeSymbolInfo(psolfunc, PACKAGE = dllname)$address
            } else stop(paste("cannot integrate: psolfunc not loaded ",psolfunc))
                                }
#        } else if (kryltype =="banded")
#        {                        
#        lenpd    <- (2*banddown + bandup +1) * n
#        mband    <-  banddown + bandup +1
#        msave    <- (n/mband) + 1
#        lwp      <- lenpd + 2 * msave
#        lip      <- n
#        if(is.loaded("dbanja",PACKAGE="deSolve"))
#           JacRes   <- getNativeSymbolInfo("dbanja",PACKAGE="deSolve")$address
#        if(is.loaded("dbanps",PACKAGE="deSolve"))
#           PsolFunc <- getNativeSymbolInfo("dbanps",PACKAGE="deSolve")$address
#        ipar     <- c(ipar,banddown,bandup)
#        } else stop(paste("cannot integrate: kryltype not known ",kryltype))
      ## If we go this route, the number of "global" results is in nout
      ## and output variable names are in outnames

      Nglobal <- nout
      rho     <- NULL
      if (is.null(outnames))
         { Nmtot   <- NULL} else
      if (length(outnames) == nout) 
         { Nmtot   <- outnames} else
      if (length(outnames) > nout) 
         Nmtot <- outnames[1:nout] else
         Nmtot <- c(outnames,(length(outnames)+1):nout)
      if (is.null(ipar)) 
         ipar<-0
      if (is.null(rpar)) 
         rpar<-0
  

    }    else {
      if(is.null(initfunc))  
         initpar <- NULL # parameter initialisation not needed if function is not a DLL    
    
        rho <- environment(func)
      # func or res and jac are overruled, either including ynames, or not
      # This allows to pass the "..." arguments and the parameters

        if (is.null(res))              # res is NOT specified, func is
        {                              
         Res    <- function(time,y,dy) 
           {if(ynames)attr(y,"names")  <- Ynames 
           FF <-func   (time,y,parms,...) 
            c(dy-unlist(FF[1]), unlist(FF[-1]))
            }    

         Res2   <- function(time,y,dy) 
            {if(ynames)attr(y,"names") <- Ynames
            func   (time,y,parms,...)}    
          } else {                       # res is specified
         Res    <- function(time,y,dy) 
            {if(ynames){attr(y,"names")  <- Ynames 
                        attr(dy,"names") <- dYnames  }
             unlist(res   (time,y,dy,parms,...))}    

         Res2   <- function(time,y,dy) 
            {if(ynames){attr(y,"names") <- Ynames
                        attr(dy,"names") <- dYnames  }
             res (time,y,dy,parms,...)}    
         }
        # the jacobian
        if (! is.null(jacfunc))        # jacobian associated with func
        {
           tmp <- eval(jacfunc(times[1], y, dy,parms, 1, ...), rho) 
           if (! is.matrix(tmp))stop("jacfunc must return a matrix\n")

          JacRes <- function(Rin,y,dy) 
              {if(ynames) {attr(y,"names")  <- Ynames
                           attr(dy,"names") <- dYnames  }        
              JF <- -1* rbind(jacfunc(Rin[1],y,dy,parms,...),erow)
              if (imp %in% c(24,25)) JF[bandup+1,]<-JF[bandup+1,]+Rin[2] else 
                                     JF           <-JF + diag(nc=n,x=Rin[2])        
              return(JF) }    
         } else if (! is.null(jacres)) { # jacobian associated with res
            tmp <- eval(jacres(times[1], y, dy, parms, 1, ...), rho) 
            if (! is.matrix(tmp))stop("jacres must return a matrix\n")
            dd <- dim(tmp)
            if((imp ==24 && dd != c(bandup+banddown+1,n)) ||
            (imp ==21 && dd != c(n,n))) stop("Jacobian dimension not ok") 

          JacRes <- function(Rin,y,dy) 
             {if(ynames) {attr(y,"names")  <- Ynames
                          attr(dy,"names") <- dYnames  }        
              rbind(jacres(Rin[1],y,dy,parms,Rin[2],...),erow)}     
           } else JacRes <- NULL
         
      ## Call res once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
        tmp <- eval(Res2(times[1], y, dy), rho) 
        if (!is.list(tmp))
            stop("Model function must return a list\n")
        if (length(tmp[[1]]) != length(y))
            stop(paste("The number of derivatives returned by func() (",
                length(tmp[[1]]), "must equal the length of the initial conditions vector (",
                length(y), ")", sep = ""))

        Nglobal <- if (length(tmp) > 1)
            length(unlist(tmp[-1]))  else 0
        # check for NULL? stop("Problem interpreting model output - check for NULL values")
    
        Nmtot <- attr(unlist(tmp[-1]),"names")
        
    }
    
### work arrays INFO, iwork, rwork

# the INFO vector
    info   <- vector("integer",20)
    info[] <- 0
    if (length(atol)==n)
     { if (length(rtol) != n)    rtol <- rep(rtol,len=n) 
     } else if (length(rtol)==n) atol <- rep(atol,len=n)
     
    info[2] <- length(atol)==n
    if (is.null(times)) 
    { info[3]<-1
    times<-c(0,1e8)
    }
#    if (krylov == TRUE) 
#    {if (is.null(kryltype) && is.null(psolfunc))
#        stop ("daspk: cannot perform integration: *psolfunc* NOT specified and krylov method chosen..")        
#     if (is.null(kryltype) && ! is.character (psolfunc)) 
#        stop ("daspk: krylov method in R-functions not yet implemented") 
#     if (is.null(kryltype) && is.null(lwp)) stop("daspk: krylov method chosen, but lwp not defined") 
#     if (is.null(kryltype) && is.null(lip)) stop("daspk: krylov method chosen, but lip not defined") 
#      info[12] <- 1
#      if (is.null(krylpar ))  {
#      krylpar <- c(min(5,n),min(5,n),5,0.05)
#       } else {
#      if (!is.numeric(krylpar)) stop("daspk: krylpar is not numeric")
#      if (length(krylpar)!=4)   stop("daspk: krylpar should contain 4 elements")
#      if (krylpar[1] <1 || krylpar[1]>n) stop("daspk: krylpar[1] MAXL not valid")
#      if (krylpar[2] <1 || krylpar[2]>krylpar[1]) stop("daspk: krylpar[2] KMP not valid")
#      if (krylpar[3] <0 ) stop("daspk: krylpar[3] NRMAX not valid")
#      if (krylpar[4] <0 || krylpar[4]>1) stop("daspk: krylpar[4] EPLI not valid")     
#      info[13] =1
#     }
 #    if (! is.null(JacRes)) info[15] <- 1
 #   } 
    # info[14], [16], [17], [18] not implemented
    if (imp %in% c(22,25)) info[5] <- 0  # internal generation jacobian
    if (imp %in% c(21,24)) info[5] <- 1  # user-defined generation jacobian
    if (imp %in% c(22,21)) info[6] <- 0  # full jacobian
    if (imp %in% c(25,24)) info[6] <- 1  # sparse jacobian
    info[7] <-  hmax != Inf
    info[8] <-  hini != 0
    nrowpd  <- ifelse(info[6]==0, n, 2*banddown+bandup+1) 
    if (info[5]==1 && is.null(jacfunc) && is.null(jacres)) 
  stop ("daspk: cannot perform integration: *jacfunc* or *jacres* NOT specified; either specify *jacfunc* or *jacres* or change *jactype*")  

    info[9] <- maxord!=5
 
    if (! is.null (estini)) info[11] <- estini        # daspk will estimate dy and algebraic equ.
    if (info[11] > 2 || info[11]< 0 ) stop("daspk: illegal value for estini")
    
# length of rwork and iwork 
#    if (info[12]==0) {
    lrw <- 50+max(maxord+4,7)*n
    if (info[6]==0) {lrw <- lrw+ n*n} else {
    if (info[5]==0) lrw <- lrw+ (2*banddown+bandup+1)*n + 2*(n/(bandup+banddown+1)+1) else
                    lrw <- lrw+ (2*banddown+bandup+1)*n  }
    liw <- 40+n 
       
#    } else {
#     maxl <- krylpar[1] 
#     kmp  <- krylpar[2]
#     lrw <- 50+(maxord+5)*n+max(maxl+3+min(1,maxl-kmp))*n + (maxl+3)*maxl+1+lwp
#     liw <- 40+lip 
#    }    

    if (info[10] %in% c(1,3)) liw <- liw+n
    if (info[11] ==1)         liw <- liw+n
    if (info[16] ==1)         liw <- liw+n
    if (info[16] ==1)         lrw <- lrw+n    
    iwork <- vector("integer",liw)
    rwork <- vector("double",lrw)
    
 if(! is.null(tcrit)) {info[4]<-1;rwork[1] <- tcrit}
    
    if(info[6] == 1) {iwork[1]<-banddown;iwork[2]<-bandup}
    if(info[7] == 1) rwork[2] <- hmax 
    if(info[8] == 1) rwork[3] <- hini
    if(info[9] == 1) iwork[3] <- maxord
    # info[10] not implemented
    if (info[11]>0)
    { lid <- ifelse(info[10] %in% c(0,2), 40, 40+n)
      iwork[lid+(1:n)       ]<- - 1
      iwork[lid+(1:(n-nalg))]<-    1
    }
#    if (info[12]==1) 
#     {iwork[27]<-lwp
#     iwork[28]<-lip}
#    if (info[13]==1)
#     {iwork[24:26]<- krylov[1:3]
#     rwork[10]<-krylov[4]}

# print to screen...
#    if (verbose)
#    {
#       if (info[12] == 0) 
#       {print("uses standard direct method") 
#       }else print("uses Krylov iterative method")
#    }

### calling solver
    storage.mode(y) <- storage.mode(dy) <- storage.mode(times) <- "double"
    storage.mode(rtol) <- storage.mode(atol)  <- "double"

    out <- .Call("call_daspk", y, dy, times, Res, as.double(initpar), 
        rtol, atol,rho, tcrit, 
        JacRes, ModelInit, PsolFunc, as.integer(verbose),as.integer(info),
        as.integer(iwork),as.double(rwork), as.integer(Nglobal),as.integer(maxIt),
        as.integer(bandup),as.integer(banddown),as.integer(nrowpd),
        as.double (rpar), as.integer(ipar),PACKAGE = "deSolve")

### saving results    

    out [1,1] <- times[1]                                   
    istate <- attr(out, "istate")
    rstate <- attr(out, "rstate")

    # ordinary output variables already estimated         
    nm <- c("time", if (!is.null(attr(y, "names"))) names(y) else as.character(1:n))
    if (Nglobal > 0) nm <- c(nm, if (!is.null(Nmtot)) Nmtot else as.character((n +
            1):(n + Nglobal)))

    attr(out, "istate") <- istate
    attr(out, "rstate") <- rstate        

    dimnames(out) <- list(nm, NULL)

    if (verbose)
    {
      print("--------------------")    
      print("daspk return code")
      print("--------------------")      

      idid <- istate[1]
      print(paste("idid = ",idid))
      if (idid >0) 
      {
      print (" *** TASK COMPLETED *** ")
      if (idid == 1) print("a step was successfully taken in the intermediate-output mode.  The code has not yet reached TOUT")
      if (idid == 2) print("the integration to TSTOP was successfully completed (T = TSTOP) by stepping exactly to TSTOP")
      if (idid == 3) print("the integration to TOUT was successfullycompleted (T = TOUT) by stepping past TOUT. Y(*) and YPRIME(*) are obtained by interpolation")
      if (idid == 4) print("the initial condition calculation, with INFO(11) > 0, was successful, and INFO(14) = 1. No integration steps were taken, and the solution is not considered to have been started")
      } else if (idid<0& idid>-33) 
      {
      print (" *** TASK INTERRUPTED *** ")
      if (idid == -1) print("a large amount of work has been expended (about 500 steps)") else
      if (idid == -2) print("the error tolerances are too stringent") else
      if (idid == -3) print("the local error test cannot be satisfied because a zero component in ATOL was specified and the corresponding computed solution component is zero.  Thus, a pure relative error test is impossible for this component") else
      if (idid == -5) print("there were repeated failures in the evaluation or processing of the preconditioner (in jacfunc)") else
      if (idid == -6) print("DDASPK had repeated error test failures on the last attempted step") else
      if (idid == -7) print("the nonlinear system solver in the time integration could not converge") else
      if (idid == -8) print("the matrix of partial derivatives appears to be singular (direct method)") else
      if (idid == -9) print("the nonlinear system solver in the time integration failed to achieve convergence, and there were repeated error test failures in this step") else
      if (idid == -10) print("the nonlinear system solver in the time integration failed to achieve convergence because IRES was equal to -1") else
      if (idid == -11) print("IRES = -2 was encountered and control is being returned to the calling program") else
      if (idid == -12) print("DDASPK failed to compute the initial Y, YPRIME") else
      if (idid == -13) print("unrecoverable error encountered inside user's PSOL routine, and control is being returned to the calling program") else
      if (idid == -14) print("the Krylov linear system solver could not achieve convergence") 
      } else if (idid ==-33)  {print (" *** TASK TERMINATED *** ")
                        print("the code has encountered trouble from which it cannot recover.  A message is printed explaining the trouble and control is returned to the calling program.")}

      print("--------------------")
      print("ISTATE values")
      print("--------------------")      

      df <- c( " idid, the return code",
               " The order of the method to be attempted on th next step",
               " The order of the method used on the last step",
               " The number of steps taken for the problem so far.",
               " The number of res evaluations for the problem so far.",
               " The number of jacobian evaluations for the problem so far.",
               " The total number of error test failures so far.",
               " The total number of nonlinear convergence failures so far.",
               " The number of convergence failures of the linear iteration so far.",
               " The length of iwork actually required.",
               " The length of rwork actually required.",
               " The number of nonlinear iterations so far.",
               " The number of linear (Krylov) iterations so far ",
               " The number of psol calls so far.")

      ii <- c(1,8:9,12:22)
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
    return(t(out))
    }
