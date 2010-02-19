# ks 21-12-09: Func <- ... unlist() ... output variables now set in C-code
#              lags added

### ============================================================================
### lsodes -- solves ordinary differential equation systems with general
### sparse Jacobian matrix.
### The sparsity structure of the Jacobian is either specified
### by the user, estimated internally (default), or of a special type.
### To date, "1D", "2D", "3D" are supported as special types.
### These are the sparsity associated with 1- 2- and 3-Dimensional PDE models
### ============================================================================

lsodes <- function(y, times, func, parms, rtol=1e-6, atol=1e-6, 
  jacvec=NULL, sparsetype="sparseint", nnz = NULL, inz = NULL,
  verbose=FALSE, tcrit = NULL, hmin=0, hmax=NULL, hini=0, ynames=TRUE,
  maxord=NULL, maxsteps=5000, lrw=NULL, liw=NULL, 
  dllname=NULL, initfunc=dllname,
  initpar=parms, rpar=NULL, ipar=NULL, nout=0, outnames=NULL, forcings=NULL,
  initforc = NULL, fcontrol=NULL, events=NULL, lags = NULL, ...)  {

### check input
  hmax <- checkInput (y, times, func, rtol, atol,
    jacvec, tcrit, hmin, hmax, hini, dllname,"jacvec")

  n <- length(y)

  if (is.null (maxord))
    maxord <- 5
  if (maxord > 5 )
    stop ("'maxord' too large: should be <= 5")
  if (maxord < 1 )
    stop ("`maxord' must be >1")

### Sparsity type and Jacobian method flag imp

  if (sparsetype=="sparseusr" && is.null(inz))
    stop("'inz' must be specified if 'sparsetype' = 'sparseusr'")
  if (sparsetype=="1D" && ! is.null(jacvec))
    stop("cannot combine 'sparsetype=1D' and 'jacvec'")
  if (sparsetype=="2D" && ! is.null(jacvec))
    stop("cannot combine 'sparsetype=2D' and 'jacvec'")
  if (sparsetype=="3D" && ! is.null(jacvec))
    stop("cannot combine 'sparsetype=3D' and 'jacvec'")

  # imp = method flag as used in lsodes
  if (! is.null(jacvec) &&  sparsetype=="sparseusr")
    imp <- 21   # inz supplied,jac supplied
  else if (! is.null(jacvec) && !sparsetype=="sparseusr")
    imp <- 121  # inz internally generated,jac supplied
  else if (is.null(jacvec) &&  sparsetype%in%c("sparseusr","1D","2D","3D"))
    imp <- 22   # inz supplied,jac not supplied
  else
    imp <- 222  # sparse Jacobian, calculated internally

## Special-purpose sparsity structures: 1-D and 2-D reaction-transport problems
## Typically these applications are called via ode.1D, ode.2D and ode.3D
## Here the sparsity is specified in the C-code; this needs extra input:
## the number of components *nspec* and the dimensionality of the problem
## (number of boxes in each direction).
## This information is passed by ode.1D, ode.2D and ode.3D in parameter
## nnz (a vector).
## nnz is altered to include the number of nonzero elements (element 1).
## 'Type' contains the type of sparsity + nspec + num boxes + cyclicBnd

  if (sparsetype=="1D") {
    nspec  <- nnz[1]
    Type   <- c(2,nnz)    #type=2
    nnz    <- n*(2+nspec)-2*nspec
  } else if (sparsetype=="2D")  {
    nspec  <- nnz[1]
    dimens <- nnz[2:3]
    Type   <- c(3,nnz)   #type=3
    nnz    <- n*(4+nspec)-2*nspec*(sum(dimens))

    if (Type[5]==1) { # cyclic boundary in x-direction
      nnz <- nnz + 2*dimens[1]*nspec
    }
    if (Type[6] ==1) {# cyclic boundary in y-direction
      nnz <- nnz + 2*dimens[2]*nspec
    }
  } else if (sparsetype=="3D")  {
    nspec  <- nnz[1]
    dimens <- nnz[2:4]    #type=4
    Type   <- c(4,nnz)
    nnz    <- n*(6+nspec)-2*nspec*(sum(dimens))

    if (Type[6]==1) { # cyclic boundary in x-direction
      nnz <- nnz + 2*dimens[1]*nspec
    }
    if (Type[7] ==1) {# cyclic boundary in y-direction
      nnz <- nnz + 2*dimens[2]*nspec
    }
    if (Type[8] ==1) {# cyclic boundary in y-direction
      nnz <- nnz + 2*dimens[3]*nspec
    }
  } else if (sparsetype == "sparseusr") {
    Type <- 0
    nnz  <- nrow(inz)
  } else  {
    Type <- 1
    if (is.null(nnz))   nnz <- n*n
  }
  if (nnz<1)
    stop ("Jacobian should at least contain one non-zero value")

### model and Jacobian function 
  JacFunc   <- NULL
  Ynames    <- attr(y,"names")
  flist     <- list(fmat=0,tmat=0,imat=0,ModelForc=NULL)
  ModelInit <- NULL
  Eventfunc <- NULL
  
  events <- checkevents(events, times, Ynames, dllname) 

  if (is.character(func)) {   # function specified in a DLL
    DLL <- checkDLL(func,jacvec,dllname,
                    initfunc,verbose,nout, outnames, JT=2)

    ModelInit <- DLL$ModelInit
    Func    <- DLL$Func
    JacFunc <- DLL$JacFunc
    Nglobal <- DLL$Nglobal
    Nmtot   <- DLL$Nmtot

    if (! is.null(forcings))
      flist <- checkforcings(forcings,times,dllname,initforc,verbose,fcontrol)

    rho <- NULL
    if (is.null(ipar)) ipar<-0
    if (is.null(rpar)) rpar<-0
    Eventfunc <- events$func
  } else {

    if(is.null(initfunc))
      initpar <- NULL # parameter initialisation not needed if function is not a DLL
    rho <- environment(func)
    # func and jac are overruled, either including ynames, or not
    # This allows to pass the "..." arguments and the parameters

    if (ynames) {
      Func    <- function(time,state) {
        attr(state,"names") <- Ynames
        unlist(func   (time,state,parms,...))
      }

      Func2   <- function(time,state) {
        attr(state,"names") <- Ynames
        func   (time,state,parms,...)
      }
         
      JacFunc <- function(time,state,J){
        attr(state,"names") <- Ynames
        jacvec(time,state,J,parms,...)
      }
       if (! is.null(events$Type))
       if (events$Type == 2)
         Eventfunc <- function(time,state) {
           attr(state,"names") <- Ynames
           events$func(time,state,parms,...) 
         }
    } else {                            # no ynames...
      Func    <- function(time,state)
        unlist(func   (time,state,parms,...))
        
      Func2   <- function(time,state)
        func   (time,state,parms,...)
         
      JacFunc <- function(time,state,J)
        jacvec(time,state,J,parms,...)
     
      if (! is.null(events$Type))
      if (events$Type == 2)
         Eventfunc <- function(time,state)  
           events$func(time,state,parms,...) 
    }
        
    ## Check function and return the number of output variables +name
    FF <- checkFunc(Func2,times,y,rho)
    Nglobal<-FF$Nglobal
    Nmtot <- FF$Nmtot

    if (! is.null(events$Type))
      if (events$Type == 2) 
        checkEventFunc(Eventfunc,times,y,rho)
    
  }

### work arrays iwork, rwork
  # 1. Estimate length of rwork and iwork if not provided via arguments lrw, liw
  moss  <- imp%/%100         # method to be used to obtain sparsity
  meth  <- imp%%100%/%10     # basic linear multistep method
  miter <- imp%%10           # corrector iteration method
  lenr = 2     # real to integer wordlength ratio (2 due to double precision)

  if (is.null(lrw)) {         # make a guess of real work space needed
    lrw = 20+n*(maxord+1)+3*n +20  #extra 20 to make sure

    if(miter == 1) lrw = lrw + 2*nnz + 2*n + (nnz+9*n)/lenr
    if(miter == 2) lrw = lrw + 2*nnz + 2*n + (nnz+10*n)/lenr
    if(miter == 3) lrw = lrw + n + 2
  
    if (sparsetype == "1D") lrw <- lrw*1.1 # increase to be sure it is enough...
  }
  
  if (is.null(liw)) {         # make a guess of itneger work space needed
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
  if(sparsetype=="sparseusr")  {
    iwork   <- vector("integer",liw)
    iwork[] <- 0

    iw       <- 32+n
    iwork[31]<- iw

    # input = 2-columned matrix inz; converted to ian,jan and put in iwork
    # column indices should be sorted...
    rr  <- inz[,2]
    if (min(rr[2:nnz]-rr[1:(nnz-1)])<0)
      stop ("cannot proceed: row indices in inz should be sorted")

    for(i in 1:n)  {
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
  if (verbose)   {
    printtask(itask,func,jacvec)      
    printM("\n--------------------")
    printM("Integration method")
    printM("--------------------\n")
    if (imp == 21)
      txt <- "  The user has supplied indices to nonzero elements of Jacobian,
      and a Jacobian function"  else
    if (imp == 22)  {
      if (sparsetype=="sparseusr")
        txt <-"  The user has supplied indices to nonzero elements of Jacobian,
        the Jacobian will be estimated internally, by differences"
      if (sparsetype=="1D")
        txt <-"  The nonzero elements are according to a 1-D model,
        the Jacobian will be estimated internally, by differences"
      if (sparsetype=="2D")
        txt <-"  The nonzero elements are according to a 2-D model,
        the Jacobian will be estimated internally, by differences"
      if (sparsetype=="3D")
        txt <-"  The nonzero elements are according to a 3-D model,
        the Jacobian will be estimated internally, by differences"
                   } else
    if (imp == 122)
      txt <-"  The user has supplied the Jacobian,
      its structure (indices to nonzero elements) will be obtained from NEQ+1 calls to jacvec"  else
    if (imp == 222)
      txt <-"  The Jacobian will be generated internally,
      its structure (indices to nonzero elements) will be obtained from NEQ+1 calls to func"
   
    printM(txt)
  } 

### calling solver
  storage.mode(y) <- storage.mode(times) <- "double"
  IN <-3

  lags <- checklags(lags) 

  out <- .Call("call_lsoda",y,times,Func,initpar,
               rtol, atol, rho, tcrit, JacFunc, ModelInit, Eventfunc,
               as.integer(verbose), as.integer(itask), as.double(rwork),
               as.integer(iwork), as.integer(imp),as.integer(Nglobal),
               as.integer(lrw),as.integer(liw),as.integer(IN),
               NULL, as.integer(0), as.double (rpar), as.integer(ipar),
               as.integer(Type),flist, events, lags, PACKAGE="deSolve")

### saving results    
  out <- saveOut(out, y, n, Nglobal, Nmtot, func, Func2,
                 iin=c(1,12:20), iout=c(1:3,14,5:9,17))

  attr(out, "type") <- "lsodes"
  if (verbose) diagnostics(out)
  out
}
