### ============================================================================
### ============================================================================
### S3 methods
### karline+Thomas: from version 1.9, also possible to plot multiple
###                 outputs and to add observations.
### ============================================================================
### ============================================================================


### ============================================================================
### first some common functions
### ============================================================================
# Update range, taking into account neg values for log transformed values
Range <- function(Range, x, log) {
   if (log)
      x[x <= 0] <- min(x[x>0])  # remove zeros
   return( range(Range, x, na.rm = TRUE) )
}

## =============================================================================
## function for checking and expanding arguments in dots (...) with default
## =============================================================================

expanddots <- function (dots, default, n) {
  dots <- if (is.null(dots)) default else dots
  rep(dots, length.out = n)
}

# ks->Th: for xlim and ylim....
expanddotslist <- function (dots, n) {
  if (is.null(dots)) return(dots)
  dd <- if (!is.list(dots )) list(dots) else dots
  rep(dd, length.out = n)
}

## =============================================================================
## functions for expanding arguments in dots  (...)
## =============================================================================

repdots <- function(dots, n)
  if (is.function(dots)) dots else rep(dots, length.out = n)

setdots <- function(dots, n) lapply(dots, repdots, n)

## =============================================================================
## function for extracting element 'index' from dots  (...)
## =============================================================================

extractdots <- function(dots, index) {
  ret <- lapply(dots, "[", index)
  ret <- lapply(ret, unlist) # flatten list
  return(ret)
}

### ============================================================================
### Merge two observed data files; assumed that first column = 'x' and ignored
### ============================================================================

# from 3-columned format (what, where, value) to wide format...
convert2wide <- function(Data) {
    cnames   <- as.character(unique(Data[,1]))

    MAT      <- Data[Data[,1] == cnames[1], 2:3]
    colnames.MAT <- c("x", cnames[1])

    for ( ivar in cnames[-1]) {
      sel <- Data[Data[,1] == ivar, 2:3]
      nt  <- cbind(sel[,1],matrix(nrow = nrow(sel), ncol = ncol(MAT)-1, data = NA),sel[,2])
      MAT <- cbind(MAT, NA)
      colnames(nt) <- colnames(MAT)
      MAT <- rbind(MAT, nt)
      colnames.MAT <- c(colnames.MAT, ivar)
    }
  colnames(MAT) <- colnames.MAT
  return(MAT)
}


mergeObs <- function(obs, Newobs) {

  if (! class(Newobs) %in% c("data.frame","matrix"))
    stop ("the elements in 'obs' should be either a 'data.frame' or a 'matrix'")

  if (is.character(Newobs[,1]) | is.factor(Newobs[,1]))
    Newobs <- convert2wide(Newobs)

  obsname <- colnames(obs)

## check if some observed variables in NewObs are already in obs
  newname <- colnames(Newobs)[-1]    # 1st column = x-var and ignored
  ii <- which (newname %in% obsname)
  if (length(ii) > 0)
    obsname <- c(obsname, newname[-ii] )
  else
    obsname <- c(obsname, newname)

## padding with NA of the two datasets
  O1 <- matrix(nrow = nrow(Newobs), ncol = ncol(obs), data = NA)
  O1[ ,1] <- Newobs[,1]
  for (j in ii) {   # obseerved data in common are put in correct position
    jj <- which (obsname == newname[j])
    O1[,jj] <- Newobs[,j+1]
  }
  O1 <- cbind(O1, Newobs[,-c(1,ii+1)] )
  colnames(O1) <- obsname

  nnewcol <- ncol(Newobs)-1 - length (ii)  # number of new columns
  if (nnewcol > 0) {
     O2 <- matrix(nrow = nrow(obs), ncol = nnewcol, data = NA)
     O2 <- cbind(obs, O2)
     colnames(O2) <- obsname
  } else O2 <- obs

  obs <- rbind(O2, O1)
  return(obs)
}

## =============================================================================
## Set the mfrow parameters and whether to "ask" for opening a new device
## =============================================================================

setplotpar <- function(nmdots, dots, nv, ask) {
    if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
      nc <- min(ceiling(sqrt(nv)), 3)
      nr <- min(ceiling(nv/nc), 3)
      mfrow <- c(nr, nc)
    } else if ("mfcol" %in% nmdots)
      mfrow <- rev(dots$mfcol)
    else mfrow <- dots$mfrow

    if (! is.null(mfrow))  mf <- par(mfrow = mfrow)

   ## interactively wait if there are remaining figures
    if (is.null(ask))
      ask <- prod(par("mfrow")) < nv && dev.interactive()

    return(ask)
}

## =============================================================================
## find a variable
## =============================================================================

selectvar <- function (which, var, NAallowed = FALSE) {
  if (!is.numeric(which)) {
    ln <- length(which)
    ## keep ordering...
    Select <- NULL
    for ( i in 1:ln) {
      ss <- which(which[i]==var)
      if (length(ss) ==0 & ! NAallowed)
        stop("variable ", which[i], " not in variable names")
      else if (length(ss) == 0)
        Select <- c(Select,NA)
      else
        Select <- c(Select,ss)
    }
  } else {
    Select <- which + 1  # "Select" now refers to the column number
    if (max(Select) > length(var))
        stop("index in 'which' too large: ", max(Select)-1)
    if (min(Select) < 1)
        stop("index in 'which' should be > 0")
  }
  return(Select)
}

### ============================================================================
### print a deSolve object
### ============================================================================

print.deSolve <- function(x, ...)
  print(as.data.frame(x), ...)

### ============================================================================
### Create a histogram for a list of variables
### ============================================================================

hist.deSolve <- function (x, which = 1:(ncol(x)-1), ask = NULL, ...) {

    t        <- 1     # column with "times"
    varnames <- colnames(x)
    which    <- selectvar(which, varnames)

    np     <- length(which)
    dots   <- list(...)
    nmdots <- names(dots)

## Set par mfrow and ask.
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

## expand all dots to np values (no defaults
    Dots  <- setdots(dots, np)

    # different from default settings
    Dots$main <- expanddots (dots$main, varnames[which], np)
    Dots$xlab <- expanddots (dots$xlab, varnames[t],     np)
    Dots$xlab <- expanddots (dots$xlabb, ""        ,     np)

## xlim and ylim are special: they are vectors or lists
    xxlim <- expanddotslist(dots$xlim, np)
    yylim <- expanddotslist(dots$ylim, np)

## plotting
    for (ii in 1:np) {
        i <- which[ii]
        dots <- extractdots(Dots, ii)
        if (! is.null(xxlim[[ii]])) dots$xlim <- xxlim[[ii]]
        if (! is.null(yylim[[ii]])) dots$ylim <- yylim[[ii]]
        do.call("hist", c(alist(x[, i]), dots))
    }
}
### ============================================================================
### Image, filled.contour and persp plots
### ============================================================================
image.deSolve <- function (x, which = NULL, ask = NULL,
  add.contour = FALSE, grid = NULL, method = "image", legend = FALSE, ...) {

    dimens <- attributes(x)$dimens
    if (is.null(dimens))
      stop("cannot make an image from deSolve output which is 0-dimensional")
    else if (length(dimens) ==1)  # 1-D
      plot.ode1D(x, which, ask, add.contour, grid, method=method,
        legend = legend, ...)
    else if (length(dimens) ==2)  # 2-D
      plot.ode2D(x, which, ask, add.contour, grid, method=method,
        legend = legend, ...)
    else
      stop("cannot make an image from deSolve output with more than 2 dimensions")
}

### ============================================================================
### Plot utilities for the S3 plot method, 0-D, 1-D, 2-D
### Plotting 0-D variables
### ============================================================================

plot.deSolve <- function (x, ..., which = NULL, ask = NULL, obs = NULL,
    obspar = list()) {

## check observed data - can be a list
    nobs <- 0

    if (! is.null(obs)) {

      if (!is.data.frame(obs) & is.list(obs)) { # a list with different data sets
       Obs <- obs
       obs <- Obs[[1]]
       obs.pos <- matrix(nrow = 1, c(1, nrow(obs)))
       if (! class(obs) %in% c("data.frame", "matrix"))
         stop ("'obs' should be either a 'data.frame' or a 'matrix'")
       if (length(Obs) > 1)
         for ( i in 2 : length(Obs)) {
           obs <- mergeObs(obs, Obs[[i]])
           obs.pos <- rbind(obs.pos, c(obs.pos[nrow(obs.pos),2] +1, nrow(obs)))
         }
       obsname <- colnames(obs)
      } else {
       if (is.character(obs[,1]) | is.factor(obs[,1]))   # long format - convert
          obs <- convert2wide(obs)
       obsname <- colnames(obs)
       if (! class(obs) %in% c("data.frame", "matrix"))
         stop ("'obs' should be either a 'data.frame' or a 'matrix'")
       obs.pos <- matrix(nrow = 1, c(1, nrow(obs)))
      }
    DD <- duplicated(obsname)
    if (sum(DD) > 0)
      obs <- mergeObs(obs[,!DD], cbind(obs[,1],obs[,DD]))
    nobs <- nrow(obs.pos)
    }

## variables to be plotted
    varnames <- colnames(x)
    Which   <- which

    if (is.null(Which) & is.null(obs))  # All variables plotted
      Which <- 1 : (ncol(x)-1)

    else if (is.null(Which)) {          # All common variables in x and obs plotted
     Which <- which(varnames %in% obsname)
     Which <- Which[Which != 1]         # remove first element (x-value)
     Which <- varnames[Which]           # names rather than numbers
    }

## Position of variables to be plotted in "x"
    t       <- 1     # column with "times"
    xWhich  <- selectvar(Which, varnames)
    np      <- length(xWhich)

    ldots   <- list(...)
    ndots   <- names(ldots)

## number of figures in a row and interactively wait if remaining figures
    ask <- setplotpar(ndots, ldots, np, ask)
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }

## create two lists: x2:   other deSolve objects,
##                   dots: remaining (plotting) parameters
    x2     <- list()
    dots   <- list()
    nd     <- 0
    nother <- 0

    if (length(ldots) > 0)
     for ( i in 1:length(ldots))
      if ("deSolve" %in% class(ldots[[i]])) { # a deSolve object
        x2[[nother <- nother + 1]] <- ldots[[i]]
        names(x2)[nother] <- ndots[i]
        # a list of deSolve objects
       } else if (is.list(ldots[[i]]) & "deSolve" %in% class(ldots[[i]][[1]])) {
        for (j in 1:length(ldots[[i]])) {
          x2[[nother <- nother+1]] <- ldots[[i]][[j]]
          names(x2)[nother] <- names(ldots[[i]])[[j]]
        }
       } else if (! is.null(ldots[[i]])) {  # a graphical parameter
        dots[[nd <- nd+1]] <- ldots[[i]]
        names(dots)[nd] <- ndots[i]
      }
    nmdots <- names(dots)

## check compatibility of all deSolve objects
    if (nother > 0) {
      for ( i in 1:nother) {
        if (min(colnames(x2[[i]]) == varnames) == 0)
          stop("'x' is not compatible with other deSolve objects - colnames not the same")
      }
    }

    nx <- nother + 1 # total number of deSolve objects to be plotted

## Position of variables in "obs" (NA = not observed)
    if (nobs > 0) {
      ObsWhich <- selectvar(varnames[xWhich], obsname, NAallowed = TRUE)
      ObsWhich [ ObsWhich > ncol(obs)] <- NA
      Obspar <- setdots(obspar, nobs)
    } else
      ObsWhich <- rep(NA, np)

## plotting parameters : split in plot parameters and point parameters
    plotnames <- c("xlab","ylab","xlim","ylim","main","sub","log","asp",
                   "ann","axes","frame.plot","panel.first","panel.last",
                   "cex.lab","cex.axis","cex.main")

    # plot.default parameters
    ii <- names(dots) %in% plotnames
    dotmain <- dots[ii]
    dotmain <- setdots(dotmain, np)  # expand to np for each plot

    # these are different from the default
    dotmain$xlab <- expanddots(dots$xlab, varnames[t]     , np)
    dotmain$ylab <- expanddots(dots$ylab, ""              , np)
    dotmain$main <- expanddots(dots$main, varnames[xWhich], np)

    # ylim and xlim can be lists and are at least two values
    yylim  <- expanddotslist(dots$ylim, np)
    xxlim  <- expanddotslist(dots$xlim, np)

    # point parameters
    ip <- !names(dots) %in% plotnames
    dotpoints <- dots[ip]
    dotpoints <- setdots(dotpoints, nx)   # expand all dots to nx values

    # these are different from default
    dotpoints$type <- expanddots(dots$type, "l", nx)
    dotpoints$lty  <- expanddots(dots$lty, 1:nx, nx)
    dotpoints$pch  <- expanddots(dots$pch, 1:nx, nx)
    dotpoints$col  <- expanddots(dots$col, 1:nx, nx)
    dotpoints$bg   <- expanddots(dots$bg,  1:nx, nx)


## for each output variable (plot)
    for (i in 1 : np) {
      ii <- xWhich[i]     # position of variable in 'x'
      io <- ObsWhich[i]   # position of variable in 'obs'

      # plotting parameters for deSolve output 1 (opens a plot)
      Dotmain   <- extractdots(dotmain, i)
      Dotpoints <- extractdots(dotpoints, 1)

      Xlog <- Ylog <- FALSE
      if (! is.null(Dotmain$log)) {
        Ylog  <- length(grep("y",Dotmain$log))
        Xlog  <- length(grep("x",Dotmain$log))
      }

      ## ranges
      if ( is.null (yylim[[i]])) {
        yrange <- Range(NULL, x[, ii], Ylog)
        if (nother>0)
         for (j in 1:nother)
           yrange <- Range(yrange, x2[[j]][,ii], Ylog)
        if (! is.na(io)) yrange <- Range(yrange, obs[,io], Ylog)
          Dotmain$ylim <- yrange
      } else
        Dotmain$ylim  <- yylim[[i]]


      if ( is.null (xxlim[[i]])) {
        xrange <- Range(NULL, x[, t], Xlog)
        if (nother>0)
         for (j in 1:nother)
           xrange <- Range(xrange, x2[[j]][,t], Xlog)
        if (! is.na(io)) xrange <- Range(xrange, obs[,1], Xlog)
          Dotmain$xlim <- xrange
      } else
        Dotmain$xlim  <- xxlim[[i]]


      ## first deSolve object plotted (new plot created)
      do.call("plot", c(alist(x[, t], x[, ii]), Dotmain, Dotpoints))

      if (nother > 0)        # if other deSolve outputs
        for (j in 2:nx)
          do.call("lines", c(alist(x2[[j-1]][, t], x2[[j-1]][, ii]),
                  extractdots(dotpoints, j)) )

      ## if observed variables: select correct pars
      if (! is.na(io))
           for (j in 1: nobs)
              if (length (i.obs <- obs.pos[j, 1]:obs.pos[j, 2]) > 0)
                do.call("points", c(alist(obs[i.obs, 1], obs[i.obs, io]),
                         extractdots(Obspar, j) ))

    }
}

### ============================================================================
## to draw a legend
### ============================================================================

drawlegend <- function (parleg, dots) {
        Plt <- par(plt = parleg)
        par(new = TRUE)
        ix <- 1
        minz <- dots$zlim[1]
        maxz <- dots$zlim[2]
        binwidth <- (maxz - minz)/64
        iy <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
        iz <- matrix(iy, nrow = 1, ncol = length(iy))

        image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "",
              ylab = "", col = dots$col)

        do.call("axis", list(side = 4, mgp = c(3,1,0), las=2))

        par(plt = Plt)
}

### ============================================================================
## to drape a color over a persp plot.
### ============================================================================

drapecol <- function (A,
          col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100),
              NAcol = "white", Range = NULL)
{
    nr <- nrow(A)
    nc <- ncol(A)
    ncol <- length(col)

    AA <- 0.25 * (A[1:(nr - 1), 1:(nc - 1)] + A[1:(nr - 1), 2:nc] +
        A[2:nr, 1:(nc - 1)] + A[2:nr, 2:nc])
    if (is.null(Range))
      Range <- range(A, na.rm = TRUE)
    else {
      AA[AA > Range[2]] <- Range[2]
      AA[AA < Range[1]] <- Range[1]
    }
    Ar <- Range
    rn <- Ar[2] - Ar[1]
    ifelse(rn != 0, drape <- col[1 + trunc((AA - Ar[1])/rn *
        (ncol - 1))], drape <- rep(col[1], ncol))
    drape[is.na(drape)] <- NAcol
    return(drape)
}

### ============================================================================
### Finding 1-D variables
### ============================================================================

select1dvar <- function (which, var, att) {

    proddim <- prod(att$dimens)
    ln   <- length(which)
    csum <- cumsum(att$lengthvar) + 2

    if (!is.numeric(which)) {
        # keep ordering...
        Select <- NULL
        for ( i in 1 : ln) {
          ss <- which(which[i] == var)
          if (length(ss) == 0)
            stop("variable ", which[i], " not in variable names")
          Select <- c(Select, ss)
        }
    }
    else {
        Select <- which  # "Select now refers to the column number
        if (max(Select) > length(var))
            stop("index in 'which' too large")
        if (min(Select) < 1)
            stop("index in 'which' should be > 0")
    }

    istart <- numeric(ln)
    istop  <- numeric(ln)
    for ( i in 1 : ln) {
      if (Select[i] <= att$nspec) {
        ii <- Select[i]
        istart[i] <- (ii-1)*proddim + 2
        istop[i]  <- istart[i] + proddim-1
      } else {
        ii <- Select[i] - att$nspec
        istart[i] <- csum[ii]
        istop[i]  <- csum[ii+1]-1
      }
      if (istart[i] == istop[i])
        stop ("variable ",which[i], " is not a 1-D variable")

    }
  return(list(which = Select, istart = istart, istop = istop))
 }
### ============================================================================
### Finding 2-D variables
### ============================================================================

select2dvar <- function (which, var, att) {

    proddim <- prod(att$dimens)
    ln   <- length(which)
    csum <- cumsum(att$lengthvar) + 2

    if (!is.numeric(which)) {
        # keep ordering...
        Select <- NULL
        for ( i in 1 : ln) {
          ss <- which(which[i] == var)
          if (length(ss) == 0)
            stop("variable ", which[i], " not in variable names")
          Select <- c(Select, ss)
        }
    }
    else {
        Select <- which  # "Select now refers to the column number
        if (max(Select) > length(var))
            stop("index in 'which' too large")
        if (min(Select) < 1)
            stop("index in 'which' should be > 0")
    }

    istart <- numeric(ln)
    istop  <- numeric(ln)
    dimens <- list()
    for ( i in 1 : ln) {
      if (Select[i] <= att$nspec) {  # a state variable
        ii <- Select[i]
        istart[i] <- (ii-1)*proddim + 2
        istop[i]  <- istart[i] + proddim-1
        dimens[[i]] <- att$dimens
      } else {
        ii <- Select[i] - att$nspec
        istart[i] <- csum[ii]
        istop[i]  <- csum[ii+1]-1
        ij <- which(names(att$dimvar) == var[Select[i]])
        if (length(ij) == 0)
          stop("variable ",var[Select]," is not two-dimensional")
        dimens[[i]] <- att$dimvar[[ij]]

      }
    }

  return(list(which = Select, istart = istart, istop = istop, dim = dimens))
 }

### ============================================================================
### plotting 1-D variables as line plot
### ============================================================================

plot.1D <- function (x, which = NULL, ask = NULL, grid = NULL,
  xyswap = FALSE, delay=0, ...) {

## Check settings of x
    att <- attributes(x)
    nspec <- att$nspec
    dimens <- att$dimens
    proddim <- prod(dimens)
    if (length(dimens) != 1)
      stop ("plot.1D only works for models solved with 'ode.1D'")

    if ((ncol(x)- nspec*proddim) < 1)
      stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

    if (is.null(which))
      which <- 1:nspec
    np <- length(which)

    varnames <-  if (! is.null(att$ynames))
      att$ynames else 1:nspec

    if (! is.null(att$lengthvar))
      varnames <- c(varnames, names(att$lengthvar)[-1])

    Select <- select1dvar(which, varnames, att)
    which <- Select$which

    dots <- list(...)
    nmdots <- names(dots)

## number of figures in a row and interactively wait if remaining figures
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    Dots <- setdots(dots, np) # expand all dots to np values (no defaults)

    # These are different from defaulst
    Dots$xlab <- expanddots(dots$xlab,  "x", np)
    Dots$ylab <- expanddots(dots$ylab,  varnames[which], np)


    # allow individual xlab and ylab (vectorized)
    times <- x[,1]
    Dots$main <- expanddots(Dots$main, paste("time",times), np)
    Dotsmain <- expanddots(dots$main, paste("time",times), length(times))

    ## xlim and ylim are special:
    xxlim <- expanddotslist(dots$xlim, np)
    yylim <- expanddotslist(dots$ylim, np)

    xyswap <- rep(xyswap, length = np)

    for (j in 1:length(times)) {
      for (i in 1:np) {
        istart <- Select$istart[i]
        istop  <- Select$istop[i]

        out <- x[j,istart:istop]

        if (! is.null(grid))
          Grid <- grid
        else
          Grid <- 1:length(out)

        dots      <- extractdots(Dots, i)
        dots$main <- Dotsmain[j]
        if (! is.null(xxlim[[i]])) dots$xlim <- xxlim[[i]]
        if (! is.null(yylim[[i]])) dots$ylim <- yylim[[i]]
        if (is.null(yylim[[i]]) & xyswap[i])
           dots$ylim <- rev(range(Grid))    # y-axis


        if (! xyswap[i]) {
          do.call("plot", c(alist(Grid, out), dots))
        } else {
          if (is.null(Dots$xlab[i]) | is.null(Dots$ylab[i])) {
            dots$ylab <- Dots$xlab[i]
            dots$xlab <- Dots$ylab[i]
          }

          do.call("plot", c(alist(out, Grid), dots))
        }
       }
       if (delay > 0) Sys.sleep(0.001 * delay)
     }
}

### ============================================================================

plot.ode1D <- function (x, which, ask, add.contour, grid,
   method = "image", legend, ...) {

## Default color scheme
   BlueRed <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

## if x is vector, check if there are enough columns ...
   att <- attributes(x)
   nspec <- att$nspec
   dimens <- att$dimens
   proddim <- prod(dimens)

    if ((ncol(x)- nspec * proddim) < 1)
      stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

## variables to be plotted
   if (is.null(which))
      which <- 1 : nspec

   np <- length(which)

   varnames <-  if (! is.null(att$ynames))
      att$ynames else 1:nspec

   if (! is.null(att$lengthvar))
      varnames <- c(varnames, names(att$lengthvar)[-1])

   Select <- select1dvar(which, varnames, att)
   which <- Select$which

   dots <- list(...)
   nmdots <- names(dots)

## number of figures in a row and interactively wait if remaining figures
   ask <- setplotpar(nmdots, dots, np, ask)
   if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

   Dots  <- setdots(dots, np)   # expand dots to np values (no defaults)

    # different from the default
   Dots$main  <- expanddots(dots$main, varnames[which], np)
   Dots$xlab  <- expanddots(dots$xlab, "times",  np)
   Dots$ylab  <- expanddots(dots$ylab, "",       np)

   # colors - different if persp, image or filled.contour

   if (method == "persp")
      dotscol <- dots$col

   else if (method == "filled.contour")  {
      dotscolorpalette <- if (is.null(dots$color.palette))
        BlueRed else dots$color.palette
      dotscol <- dotscolorpalette(100)
      add.contour <- FALSE
      legend <- FALSE
   } else
     if (is.null(dots$col))
       dotscol <- BlueRed(100) else dotscol <- dots$col

   Addcontour <- rep(add.contour, length = np)

## xlim, ylim and zlim are special:
   xxlim <- expanddotslist(dots$xlim, np)
   yylim <- expanddotslist(dots$ylim, np)
   zzlim <- expanddotslist(dots$zlim, np)

   times <- x[,1]

   if (legend) {
      parplt <- par("plt") - c(0,0.07,0,0)
      parleg <- c(parplt[2]+0.02, parplt[2]+0.05, parplt[3], parplt[4])
      plt.or <- par(plt = parplt)
      on.exit(par(plt = plt.or))
   }
# Check if grid is increasing...
   if (! is.null(grid))
      gridOK <- min(diff (grid)) >0
   else
      gridOK <- TRUE

   if (! gridOK) grid <- rev(grid)

## for each output variable (plot)
   for (i in 1:np) {

      ii     <- which[i]
      istart <- Select$istart[i]
      istop  <- Select$istop[i]
      if (gridOK)
        out    <- x[ ,istart:istop]
      else
        out    <- x[ ,istop:istart]

      dots      <- extractdots(Dots, i)
      if (! is.null(xxlim)) dots$xlim <- xxlim[[i]]
      if (! is.null(yylim)) dots$ylim <- yylim[[i]]
      if (! is.null(zzlim))
        dots$zlim <- zzlim[[i]]
      else
        dots$zlim <- range(out, na.rm=TRUE)
      List <- alist(z = out, x = times)
      if (! is.null(grid)) List$y = grid

      if (method == "persp") {
          if (is.null(dots$zlim))  # this to prevent error when range = 0
            if (diff(range(out, na.rm=TRUE)) == 0)
              dots$zlim <- c(0, 1)

          if(is.null(dotscol))
             dots$col <- drapecol(out, col = BlueRed (100), Range = dots$zlim)
          else
             dots$col <- drapecol(out, col = dotscol, Range = dots$zlim)

        } else if (method == "filled.contour")
          dots$color.palette <- dotscolorpalette
        else
          dots$col <- dotscol

        do.call(method, c(List, dots))
        if (Addcontour[i]) do.call("contour", c(List, add = TRUE))

      if (legend) {
        if (method == "persp")
           if (is.null(dotscol))
             dots$col <- BlueRed(100)
           else
              dots$col <- dotscol
        if (is.null(dots$zlim)) dots$zlim <- range(out, na.rm=TRUE)

        drawlegend(parleg, dots)
      }
  }
}
### ============================================================================
### plotting 2-D variables
### ============================================================================

plot.ode2D <- function (x, which, ask, add.contour, grid, method = "image",
   legend = TRUE, ...) {

## Default color scheme
   BlueRed <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

## if x is vector, check if there are enough columns ...
    att <- attributes(x)
    nspec <- att$nspec
    dimens <- att$dimens
    proddim <- prod(dimens)

    if ((ncol(x)- nspec*proddim) < 1)
      stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

## variables to be plotted
    if (is.null(which))
      which <- 1:nspec
    np <- length(which)

    varnames <-  if (! is.null(att$ynames))
      att$ynames else 1:nspec
    if (! is.null(att$lengthvar))
      varnames <- c(varnames, names(att$lengthvar)[-1])


    Select <- select2dvar(which,varnames,att)
    which <- Select$which

    dots <- list(...)
    nmdots <- names(dots)

## number of figures in a row and interactively wait if remaining figures
    Ask <- setplotpar(nmdots, dots, np, ask)

# here ask is always true by default...
    if (is.null(ask)) ask <- TRUE
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }

#    Main <-  if (is.null(dots$main)) varnames else rep(dots$main, length.out =np)
    N <- np * nrow(x)

    if (method == "filled.contour") {
      add.contour <- FALSE
      legend <- FALSE
    }
    Dots  <- setdots(dots, N)  # expand dots to np values (no defaults)

    # different from the default
    Dots$main <- expanddots(dots$main, varnames[which], N)
    Dots$xlab <- expanddots(dots$xlab,  "x"  , N)
    Dots$ylab <- expanddots(dots$ylab,  "y"  , N)

    if (method == "persp")
      dotscol <- dots$col

    else if (method == "filled.contour") {
      dotscolorpalette <- if (is.null(dots$color.palette))
        BlueRed else dots$color.palette
      dotscol <- dotscolorpalette(100)
      add.contour <- FALSE
      legend <- FALSE
    }
    else
    if (is.null(dots$col))
      dotscol <- BlueRed(100) else  dotscol <-  dots$col

    dotslim <- dots$zlim

    xxlim <- expanddotslist(dots$xlim, np)
    yylim <- expanddotslist(dots$ylim, np)
    zzlim <- expanddotslist(dots$zlim, np)

    Addcontour <- rep(add.contour, length = np)

    i <- 0
    if (legend) {
      parplt <- par("plt") - c(0,0.05,0,0)
      parleg <- c(parplt[2]+0.02, parplt[2]+0.05, parplt[3], parplt[4])
      plt.or <- par(plt = parplt)
      on.exit(par(plt = plt.or))
    }
    for (nt in 1:nrow(x)) {
      for (ii in 1:np) {
        i       <- i+1
        istart <- Select$istart[ii]
        istop  <- Select$istop[ii]
        out <- x[nt,istart:istop]
        dim(out) <- Select$dim[[ii]]

        dots    <- extractdots(Dots, i)
        if (! is.null(xxlim)) dots$xlim <- xxlim[[ii]]
        if (! is.null(yylim)) dots$ylim <- yylim[[ii]]
        if (! is.null(zzlim)) dots$zlim <- zzlim[[ii]]
      else  {
        dots$zlim <- range(out, na.rm=TRUE)
        if (diff(dots$zlim ) == 0 )
          dots$zlim[2] <- dots$zlim[2] +1
      }
        List <- alist(z = out)
        if (! is.null(grid)) {
          List$x <- grid$x
          List$y <- grid$y
        }

        if (method == "persp") {
          if (is.null(dots$zlim))
            if (diff(range(out, na.rm = TRUE)) == 0)
              dots$zlim <- c(0, 1)

          if(is.null(dotscol))
            dots$col <- drapecol(out, col = BlueRed(100), Range = dots$zlim)
          else
            dots$col <- drapecol(out, col = dotscol, Range = dots$zlim)

        } else if (method == "image")
          dots$col <- dotscol
        else if (method == "filled.contour")
          dots$color.palette <- dotscolorpalette

        do.call(method, c(List, dots))
        if (add.contour) do.call("contour", c(List, add = TRUE))

        if (legend) {
          if (method == "persp")
            if (is.null(dotscol))
             dots$col <- BlueRed(100)
            else
              dots$col <- dotscol
          if (is.null(dots$zlim)) dots$zlim <- range(out, na.rm=TRUE)

          drawlegend(parleg, dots)
        }
       }
     }
     if (sum(par("mfrow") - c(1, 1)) == 0)
       mtext(outer = TRUE, side = 3, paste("time ", x[nt, 1]),
         cex = 1.5, line = -1.5)
}

### ============================================================================
### Summaries of ode variables
### ============================================================================

summary.deSolve <- function(object, ...){

  att <- attributes(object)
  svar   <- att$lengthvar[1]   # number of state variables
  lvar   <- att$lengthvar[-1]  # length of other variables
  nspec  <- att$nspec          # for models solved with ode.1D, ode.2D
  dimens <- att$dimens

  if (is.null(svar)) svar <- att$dim[2]-1  # models solved as DLL

  # variable names: information for state and ordinary variables is different
  if (is.null(att$ynames))
    if (is.null(dimens))
      varnames <- colnames(object)[2:(svar+1)]
    else
      varnames <- 1:nspec
  else
    varnames <- att$ynames   # this gives one name for multi-dimensional var.

  if (length(lvar) > 0) {
    lvarnames <- names(lvar)
    if (is.null(lvarnames)) lvarnames <- (length(varnames)+1):(length(varnames)+length(lvar))
    varnames <- c(varnames, lvarnames)
  }
  # length of state AND other variables
  if (is.null(dimens))                 # all 0-D state variables
    lvar <- c(rep(1, len = svar), lvar)
  else
    lvar <- c(rep(prod(dimens), nspec), lvar) # multi-D state variables

  # summaries for all variables
  Summ <- NULL
  for (i in 1:length(lvar)) {
    if (lvar[i] > 1)
      { Select <- select1dvar(i, varnames, att)
        out <- as.vector(object[,Select$istart:Select$istop])
      }
    else {
      Select <- selectvar(varnames[i], colnames(object), NAallowed = TRUE)
      if (is.na(Select))   # trick for composite names, e.g. "A.x" rather than "A"
         Select <- cumsum(lvar)[i]
      out <- object[ ,Select]
    }
  Summ <- rbind(Summ, c(summary(out, ...), N = length(out), sd = sd(out)))
  }
  rownames(Summ) <- varnames  # rownames or an extra column?
  data.frame(t(Summ))         # like this or not transposed?
}

### ============================================================================
### Subsets of ode variables
### ============================================================================

subset.deSolve  <- function(x, subset = NULL, select = NULL, which = select, ...) {

  Which <- which # for compatibility between plot.deSolve and subset

  if (missing(subset))
        r <- TRUE
  else {
        e <- substitute(subset)
        r <- eval(e, as.data.frame(x), parent.frame())
        if (!is.logical(r))
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
  }

  if (is.numeric(Which))
    return(x[r ,Which+1])

  if (is.null(Which))
    return(x[r , -1])         # Default: all variables, except time

  att   <- attributes(x)
  svar  <- att$lengthvar[1]   # number of state variables
  lvar  <- att$lengthvar[-1]  # length of other variables
  nspec <- att$nspec          # for models solved with ode.1D, ode.2D

  if (is.null(svar)) svar <- att$dim[2]-1  # models solved as DLL
  if(is.null(nspec)) nspec <- svar
  dimens <- att$dimens

  # variable names: information for state and ordinary variables is different
  if (is.null(att$ynames))
    if (is.null(dimens))
      varnames <- colnames(x)[2:(svar+1)]
    else
      varnames <- 1:nspec
  else
    varnames <- att$ynames   # this gives one name for multi-dimensional var.
  varnames <- c("time",varnames)
  if (length(lvar) > 0) {
    lvarnames <- names(lvar)
    if (is.null(lvarnames))
      lvarnames <- (length(varnames)+1):(length(varnames)+length(lvar))
    varnames <- c(varnames, lvarnames)
  }

  # length of state AND other variables
  if (is.null(dimens))                 # all 0-D state variables
    lvar <- c(rep(1, len = svar), lvar)
  else
    lvar <- c(rep(prod(dimens), nspec), lvar) # multi-D state variables

  cvar <- cumsum(c(1,lvar))

  # Add selected variables to Out
  Out <- NULL
  for (iw in 1:length(Which)) {
    i <- which (varnames == Which[iw])
    if (length(i) == 0) {
      i <- which (colnames(x) == Which[iw])
      if (length(i) == 0)
         stop ("cannot find variable ", Which[iw], " in output")
      Out <- cbind(Out, x[,i])
    } else {
      if (is.null(i))
         stop ("cannot find variable ", Which[iw], " in output")
      istart <- 1
      if (i > 1) istart <- cvar[i-1]+1
      istop <- cvar[i]
       Out <- cbind(Out, x[ ,istart:istop])
    }
  }
  if (length(Which) == ncol(Out)) colnames(Out) <- Which
  return(Out[r,])
}
