
### lsodes -- solves ordinary differential equation systems with general
### sparse jacobian matrix. The structure of the jacobian is either specified
### by the user or estimated internally.
### The sparsity is either estimated internally (default), provided by the user
### or of a special type. To date, "1D" and "2D" are supported.
### i.e. sparsity associated with 1- and 2-Dimensional reaction-transport models

lsodes <- function(y, times, func, parms, rtol=1e-6, atol=1e-6, 
  jacvec=NULL, sparsetype="sparseint",nnz = NULL, inz = NULL, verbose=FALSE,   
  tcrit = NULL, hmin=0, hmax=NULL, hini=0, ynames=TRUE, 
  maxord=NULL, maxsteps=5000, lrw=NULL, liw=NULL, 
  dllname=NULL, initfunc=dllname,
  initpar=parms, rpar=NULL, ipar=NULL, 
  nout=0, outnames=NULL,...)  
{

### check input
###########################
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
    if (!is.null(jacvec) && !(is.function(jacvec) || is.character(jacvec)))
        stop("`jacvec' must be a function or character vector")
    if (length(atol) > 1 && length(atol) != n)
        stop("`atol' must either be a scaler, or as long as `y'")
    if (length(rtol) > 1 && length(rtol) != n)
        stop("`rtol' must either be a scaler, or as long as `y'")
    if (!is.numeric(hmin))   stop("`hmin' must be numeric")
    if (hmin < 0)            stop("`hmin' must be a non-negative value")
    if (is.null(hmax))
       hmax <- ifelse (is.null(times), 0, max(abs(diff(times))))
    if (!is.numeric(hmax))   stop("`hmax' must be numeric")
    if (hmax < 0)            stop("`hmax' must be a non-negative value")
    if (hmax == Inf)  hmax <- 0
    if (hini < 0)            stop("`hini' must be a non-negative value")
    if (is.null (maxord))       maxord <- 5
    if (maxord > 5 )         stop ("maxord too large: should be <= 5")
    if (maxord < 1 )         stop ("`maxord' must be >1")

### Sparsity type and Jacobian method flag imp
##############################################
    if(sparsetype=="sparseusr" && is.null(inz)) 
                             stop("inz must be specified if sparsetype==sparseusr")  
    if (sparsetype=="1D" && ! is.null(jacvec))
                             stop("cannot combine sparsetype=1D and jacvec")
    if (sparsetype=="2D" && ! is.null(jacvec))
                             stop("cannot combine sparsetype=2D and jacvec")

  # imp = method flag as used in lsodes
  # inz supplied,jac supplied
   if (! is.null(jacvec) &&  sparsetype=="sparseusr") imp <- 21 else 
  # inz internally generated,jac supplied
   if (! is.null(jacvec) && !sparsetype=="sparseusr") imp <- 121 else       
  # inz supplied,jac not supplied
   if (is.null(jacvec) &&  sparsetype%in%c("sparseusr","1D","2D")) imp <- 22  else
  # sparse jacobian, calculated internally
                                                      imp <- 222

### Special-purpose sparsity structures: 1-D and 2-D reaction-transport problems
  # Typically these applications are called via ode.1D and ode.2D
  # Here the sparsity is determined in the c-code; this needs extra input:
  # the number of components *Nspec* and the number of boxes. 
  # This information is passed by ode.1D and ode.2D in parameter nnz (which is a vector)
  # nnz is altered to include the number of nonzero elements (element 1).
  # Type contains the type of sparsity + nspec + num boxes

     if (sparsetype=="1D")                  {
       nspec  <- nnz[1] 
       Type   <- c(2,nnz)
       nnz    <- n*(2+nspec)-2*nspec   
      } else if (sparsetype=="2D")          {
      nspec  <- nnz[1] 
      dimens <- nnz[2:3]
      Type   <- c(3,nnz)
      nnz    <- n*(4+nspec)-2*nspec*(sum(dimens))
      } else if (sparsetype == "sparseusr") { 
      Type <- 0 
      nnz  <- nrow(inz) 
      } else                                { 
      Type <- 1
      if (is.null(nnz))   nnz <- n*n         }

     if (nnz<1) stop ("Jacobian should at least contain one non-zero value")



### model and jacobian function 
############################### 
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
      if (!is.null(jacvec)) {
        if (!is.character(jacvec))
          stop("If 'func' is dynloaded, so must 'jacvec' be")
        jacvecname <- jacvec
        if(is.loaded(jacvecname, PACKAGE = dllname))
            {JacFunc <- getNativeSymbolInfo(jacvecname, PACKAGE = dllname)$address
            } else stop(paste("cannot integrate: jac function not loaded ",jacvec))
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
         
         JacFunc <- function(time,state,J) 
         { attr(state,"names") <- Ynames
           jacvec(time,state,J,parms,...)}    
        } else {                            # no ynames...
         Func    <- function(time,state) 
           func   (time,state,parms,...)[1] 
        
         Func2   <- function(time,state) 
           func   (time,state,parms,...)    
         
         JacFunc <- function(time,state,J) 
           jacvec(time,state,J,parms,...)    
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
      
    }


### work arrays iwork, rwork
############################
# 1. Estimate length of rwork and iwork if not provided via arguments lrw, liw
  moss  <- imp%/%100         # method to be used to obtain sparsity
  meth  <- imp%%100%/%10     # basic linear multistep method
  miter <- imp%%10           # corrector iteration method
  lenr = 2     # real to integer wordlength ratio (2 due to double precision)

  if (is.null(lrw))          # make a guess of real work space needed
  {

  lrw = 20+n*(maxord+1)+3*n +20  #extra 20 to make sure

  if(miter == 1) lrw = lrw + 2*nnz + 2*n + (nnz+9*n)/lenr
  if(miter == 2) lrw = lrw + 2*nnz + 2*n + (nnz+10*n)/lenr  
  if(miter == 3) lrw = lrw + n + 2
  
  if (sparsetype == "1D") lrw <- lrw*1.1 # increase to be sure it is enough...
  }
  
  if (is.null(liw))          # make a guess of itneger work space needed
  { 
  if (moss == 0 && miter %in% c(1,2)) liw <- 31+n+nnz +30 else  # extra 30
                                      liw <- 30
  }

# 2. Allocate and set values
# only first 20 elements of rwork passed to solver; 
# other elements will be allocated in C-code
# for iwork: only first 30 elements, except when sparsity imposed

  rwork <- vector("double",20)
  rwork[] <- 0.

# iwork will contain sparsity structure (ian,jan)
# See documentation of DLSODES how this is done
  if(sparsetype=="sparseusr")   
  {
   iwork   <- vector("integer",liw)   
   iwork[] <- 0

   iw       <- 32+n
   iwork[31]<- iw

   # input = 2-columned matrix inz; converted to ian,jan and put in iwork
   # column indices should be sorted...
   rr  <- inz[,2]
   if (min(rr[2:nnz]-rr[1:(nnz-1)])<0) stop ("cannot proceed: row indices in inz should be sorted")

   for(i in 1:n)
   {
    ii <- which (rr==i)
    il <- length(ii)
    i1 <- iwork[i+30]
    i2 <- iwork[i+30]+il-1
    iwork[i+31] <- i2+1
    if (il>0) iwork[i1:i2] <- inz[ii,1]
   }
   iwork[31:(31+n)] <- iwork[31:(31+n)]-31-n
  } else   {   # sparsity not imposed; only 30 element of iwork allocated.
   iwork <- vector("integer",30)
   iwork[] <- 0
  }

# other elements of iwork, rwork  
  iwork[5] <- maxord
  iwork[6] <- maxsteps
  
  if(! is.null(tcrit)) rwork[1] <- tcrit
  rwork[5] <- hini
  rwork[6] <- hmax
  rwork[7] <- hmin

# the task to be performed.
  if (! is.null(times))
      itask <- ifelse (is.null (tcrit), 1,4) else      # times specified
      itask <- ifelse (is.null (tcrit), 2,5)           # only one step
  if(is.null(times)) times <- c(0,1e8)

### print to screen...
#######################
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
                           print(paste("model function an R-function: "))
   if (is.character(jacvec)) print(paste ("jacobian specified as a DLL: ",jacvec)) else
   if (!is.null(jacvec)) print(paste ("jacobian specified as an R-function: ")) else
                        print("jacobian not specified")
   print("--------------------")   
   print("integration method")
   print("--------------------")   
   if (imp == 21)  txt <-" the user has supplied indices to nonzero elements of jacobian, and a jacobian function"
   if (imp == 22)  {
     if (sparsetype=="sparseusr") txt <-" the user has supplied indices to nonzero elements of jacobian, the jacobian will be estimated internally, by differences"
     if (sparsetype=="1D")        txt <-" the nonzero elements are according to a 1-D model, the jacobian will be estimated internally, by differences"
     if (sparsetype=="2D")        txt <-" the nonzero elements are according to a 2-D model, the jacobian will be estimated internally, by differences"
                   }
   if (imp == 122) txt <-" the user has supplied the jacobian, its structure (indices to nonzero elements) will be obtained from NEQ+1 initial calls to jacvec"
   if (imp == 222) txt <-" the jacobian will be generated internally, its structure (indices to nonzero elements) will be obtained from NEQ+1 initial calls to func"
   
   print(txt)
  } 


### calling solver
###################
    storage.mode(y) <- storage.mode(times) <- "double"
    IN <-3

    out <- .Call("call_lsoda",y,times,Func,as.double(initpar),
                 rtol, atol, rho, tcrit, JacFunc, ModelInit,  
                 as.integer(verbose), as.integer(itask), as.double(rwork),
                 as.integer(iwork), as.integer(imp),as.integer(Nglobal),
                 as.integer(lrw),as.integer(liw),as.integer(IN),
                 NULL, as.integer(0), as.double (rpar), as.integer(ipar),
                 as.integer(Type),PACKAGE="deSolve")

### saving results    
###################
    istate <- attr(out,"istate")
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
    attr(out,"istate") <- istate
    attr(out, "rstate")<- rstate    

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
      if (idid == -6) print(" error weight became zero during problem. (Solution component i vanished, and ATOL or ATOL(i) = 0.)") else
      if (idid == -7) print(" a fatal error came from sparse solver CDRV by way of DPRJS or DSOLSS") 
            
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
               " The number of nonzero elements in the sparse jacobian")

      ii <- c(1,12:20)
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
