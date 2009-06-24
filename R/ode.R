### ============================================================================
###
### ode.1D, ode.2D ode.band: special-purpose integration routines
### ode.1D is designed for solving multi-component 1-D reaction-transport models
### ode.2D is designed for solving multi-component 2-D reaction-transport models
### ode.band is designed for solving single-component 1-D reaction-transport models
### ode.1D,ode.band offer the choice between the integrators vode,
###                  lsode, lsoda, lsodar and lsodes.
### ode.2D uses lsodes.
###
### ============================================================================

ode    <- function (y, times, func, parms,
                    method= c("lsoda","lsode","lsodes","lsodar","vode","daspk",
                              "euler", "rk4", "ode23", "ode45"),
                    ...)  {
  if (is.null(method)) method <- "lsoda"
  if (is.list(method)) {
#  is() should work from R 2.7 on ...
#   if (!is(method, "rkMethod"))
    if (!"rkMethod" %in% class(method))
      stop("'method' should be given as string or as a list of class 'rkMethod'")
    out <- rk(y, times, func, parms, method = method, ...)
  } else if (is.function(method))
    out <- method(y, times, func, parms,...)
  else
    out <- switch(match.arg(method),
      lsoda = lsoda(y, times, func, parms, ...),
      vode  = vode(y, times, func, parms, ...),
      lsode = lsode(y, times, func, parms, ...),
      lsodes= lsodes(y, times, func, parms, ...),
      lsodar= lsodar(y, times, func, parms, ...),
      daspk = daspk(y, times, func, parms, ...),
      euler = rk(y, times, func, parms, method = "euler",  ...),
      rk4   = rk(y, times, func, parms, method = "rk4", ...),
      ode23 = rk(y, times, func, parms, method = "ode23", ...),
      ode45 = rk(y, times, func, parms, method = "ode45", ...)
    )

  return(out)
}

### ============================================================================

ode.1D    <- function (y, times, func, parms, nspec = NULL, dimens = NULL,
                       method= "lsode", ...)   {
# check input
  if (any(!is.na(pmatch(names(list(...)), "jacfunc"))))
    stop ("cannot run ode.1D with jacfunc specified - remove jacfunc from call list")

  if (is.null(nspec) && is.null(dimens))
    stop ("cannot run ode.1D: nspec OR dimens should be specified")

  if (nspec == 1 & method != "lsodes") {
    out <- ode.band(y, times, func, parms, nspec = nspec, method = method, ...)
    return(out)
  }

  N     <- length(y)
  if (is.null(nspec)  )
    nspec = N/dimens
  if (N%%nspec !=0    )
    stop ("cannot run ode.1D: nspec is not an integer fraction of number of state variables")

# if lsodes is used
  if (is.character(func) || method=="lsodes") {
    if ( method != "lsodes") warning("ode.1D: R-function specified in a DLL-> integrating with lsodes")
    if (is.null(dimens) ) dimens    <- N/nspec
    out <- lsodes(y=y,times=times,func=func,parms,sparsetype="1D",nnz=c(nspec,dimens),...)

  } else {
  # internal function #
    bmodel <- function (time,state,pars,model,...) {
      Modconc <-  model(time,state[ij],pars,...)   # ij: reorder state variables
      list(c(Modconc[[1]][ii]),Modconc[-1])        # ii: reorder rate of change
    }

    if (is.character(func))
      stop ("cannot run ode.1D with R-function specified in a DLL")

    ii    <- as.vector(t(matrix(ncol=nspec,1:N)))   # from ordering per slice -> per spec
    ij    <- as.vector(t(matrix(nrow=nspec,1:N)))   # from ordering per spec -> per slice

    bmod  <- function(time,state,pars,...)
      bmodel(time,state,pars,func,...)

    if (is.null(method))
      method <- "lsode"
    if (method == "vode")
      out <- vode(y[ii], times, func=bmod, parms=parms, bandup=nspec,
                  banddown=nspec, jactype="bandint", ...)
    else if (method == "lsode")
      out <- lsode(y[ii], times, func=bmod, parms=parms, bandup=nspec,
                   banddown=nspec, jactype="bandint", ...)
    else if (method == "lsoda")
      out <- lsoda(y[ii], times, func=bmod, parms=parms, bandup=nspec,
                   banddown=nspec, jactype="bandint", ...)
    else if (method == "lsodar")
      out <- lsodar(y[ii], times, func=bmod, parms=parms, bandup=nspec,
                   banddown=nspec, jactype="bandint", ...)
    else
      stop ("cannot run ode.1D: method should be one of vode, lsoda, lsodar, lsode")
      out[,(ii+1)] <- out[,2:(N+1)]
  }
  return(out)
}

### ============================================================================

ode.2D    <- function (y, times, func, parms, nspec=NULL, dimens,
   cyclicBnd = NULL, ...)  {

 # check input
  if (any(!is.na(pmatch(names(list(...)), "jacfunc"))))
    stop ("cannot run ode.2D with jacfunc specified - remove jacfunc from call list")
  if (is.null(dimens))
     stop ("cannot run ode.2D: dimens should be specified")
  if (length(dimens)!=2)
     stop ("cannot run ode.2D: dimens should contain 2 values")

  N     <- length(y)
  if (N%%prod(dimens) !=0    )
    stop ("cannot run ode.2D: dimensions are not an integer fraction of number of state variables")

  if (is.null (nspec))
    nspec = N/prod(dimens) else
  if (nspec*prod(dimens) != N)
    stop ("cannot run ode.2D: dimens[1]*dimens[2]*nspec is not equal to number of state variables")

  Bnd <- c(0,0)
  if (! is.null(cyclicBnd)) {
    if (max(cyclicBnd) > 2 )
      stop ("cannot run steady.2D: cyclicBnd should be a vector or number not exceeding 2")
    Bnd[cyclicBnd[cyclicBnd>0]]<-1
  }

  # use lsodes - note:expects rev(dimens)...
   out <- lsodes(y=y, times=times, func=func, parms, sparsetype="2D",
          nnz=c(nspec,rev(dimens), rev(Bnd)), ...)

  return(out)
}

### ============================================================================

ode.3D    <- function (y, times, func, parms, nspec=NULL, dimens, ...){
 # check input
  if (any(!is.na(pmatch(names(list(...)), "jacfunc"))))
    stop ("cannot run ode.3D with jacfunc specified - remove jacfunc from call list")
  if (is.null(dimens))
     stop ("cannot run ode.3D: dimens should be specified")
  if (length(dimens)!=3)
     stop ("cannot run ode.3D: dimens should contain 3 values")

  N     <- length(y)
  if (N%%prod(dimens) !=0    )
    stop ("cannot run ode.3D: dimensions are not an integer fraction of number of state variables")

  if (is.null (nspec))
    nspec = N/prod(dimens) else
  if (nspec*prod(dimens) != N)
    stop ("cannot run ode.3D: dimens[1]*dimens[2]*dimens[3]*nspec is not equal to number of state variables")

  Bnd <- c(0,0,0)    #  cyclicBnd not included

  # use lsodes - note:expects rev(dimens)...
   out <- lsodes(y=y, times=times, func=func, parms, sparsetype="3D",
          nnz=c(nspec,rev(dimens), rev(Bnd)), ...)

  return(out)
}

### ============================================================================

ode.band  <- function (y, times, func, parms, nspec=NULL, bandup=nspec,
                       banddown=nspec, method = "lsode", ...)  {

  if (is.null(bandup)  )
    stop ("cannot run ode.band: bandup is not specified")
  if (is.null(banddown))
    stop ("cannot run ode.band: banddown is not specified")
  if (is.null(method))
    method <- "lsode"
  if (method == "vode")
   vode(y, times, func, parms=parms, bandup=bandup, banddown=banddown,
        jactype="bandint", ...)
  else if (method == "lsode")
   lsode(y, times, func, parms=parms, bandup=bandup, banddown=banddown,
         jactype="bandint", ...)
  else if (method == "lsoda")
   lsoda(y, times, func, parms=parms, bandup=bandup, banddown=banddown,
         jactype="bandint", ...)
  else if (method == "lsodar")
   lsodar(y, times, func, parms=parms, bandup=bandup, banddown=banddown,
          jactype="bandint", ...)
  else
   stop ("cannot run ode.band: method should be one of vode, lsoda, lsodar or lsode")

}

