### ============================================================================
### first some common functions
### ============================================================================

## =============================================================================
## Set the mfrow parameters and whether to "ask" for opening a new device
## =============================================================================

setplotpar <- function(nmdots, dots, nv, ask) {
    if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
      nc <- min(ceiling(sqrt(nv)),3)
      nr <- min(ceiling(nv/nc),3)
      mfrow <- c(nr, nc)
    }
    else if ("mfcol" %in% nmdots)
        mfrow <- rev(dots$mfcol)
    else mfrow <- dots$mfrow

    if (! is.null(mfrow)) {
      mf <- par(mfrow=mfrow)
    }

   ## interactively wait if there are remaining figures
    if (is.null(ask))
      ask <- prod(par("mfrow")) < nv && dev.interactive()

    return(ask)
}

## =============================================================================
## find a variable
## =============================================================================
## KS->ThPe: CHANGED to use a loop, such as to keep ordering of inputted vars.

selectvar <- function (which,var) {

    if (!is.numeric(which)) {
        ln <- length(which)    
        # keep ordering...
        Select <- NULL
        for ( i in 1:ln) {
          ss <- which(which[i]==var)
          if (length(ss) ==0)
            stop("variable ", which[i], " not in var")
          Select <- c(Select,ss)
        }        
    }
    else {
        Select <- which + 1  # "Select now refers to the column number
        if (max(Select) > length(var))
            stop("index in 'which' too large")
        if (min(Select) < 1)
            stop("index in 'which' should be > 0")
    }
  return(Select)
}

### ============================================================================
### S3 methods
### ============================================================================

print.deSolve <- function(x, ...)
  print(as.data.frame(x), ... )

### ============================================================================

hist.deSolve <- function (x, which = 1:(ncol(x)-1), ask = NULL, ...) {
    t     <- 1     # column with "times"
    var   <- colnames(x)
    which <- selectvar(which,var)

    np <- length(which)

    dots   <- list(...)
    nmdots <- names(dots)

    ## Set par mfrow and ask.
    ask <- setplotpar(nmdots, dots, np, ask)

    ## interactively wait if there are remaining figures
    if (is.null(ask))
        ask <- prod(par("mfrow")) < length(which) && dev.interactive()
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    Main <- is.null(dots$main)

    xxlab <- if (is.null(dots$xlab))  colnames(x)[t]  else dots$xlab
    yylab <- if (is.null(dots$ylab))  ""              else dots$ylab

    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)

    for (i in which) {
        if (Main)
            dots$main <- colnames(x)[i]
            dots$xlab <- xxlab[i-1]
            dots$ylab <- yylab[i-1]
        do.call("hist", c(alist(x[, i]), dots))
    }
}
### ============================================================================

## KS->ThPe: CHANGED  main: can also be inputted as vector; changed selectvar to keep the 
## ordering; xlab did not work...
## ideal would be to also be able to have erything in the ... a vector.

image.deSolve <- function (x, which = NULL, ask = NULL, 
  add.contour = FALSE, grid=NULL, ...) {

    dimens <- attributes(x)$dimens
    if (is.null(dimens))
      stop("cannot make an image from deSolve output which is 0-dimensional")
    else if (length(dimens) ==1)  # 1-D
      plot.ode1D(x, which, ask, add.contour, grid, ...)
    else if (length(dimens) ==2)  # 2-D
      plot.ode2D(x, which, ask, add.contour, grid, ...)
    else 
      stop("cannot make an image from deSolve output with more than 2 dimensions")
}

### ============================================================================
### Plot utilities for the S3 plot method, 0-D, 1-D, 2-D 
### ============================================================================


plot.deSolve <- function (x, which = 1:(ncol(x)-1), ask = NULL, ...) {
    if (is.null(which))
      which <- 1:(ncol(x)-1)
    t <- 1     # column with "times"
    var <- colnames(x)
    which <- selectvar(which,var)
    
    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and 
    # interactively wait if there are remaining figures
   
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    

    xxlab <- if (is.null(dots$xlab))  colnames(x)[t]  else dots$xlab
    yylab <- if (is.null(dots$ylab))  ""              else dots$ylab
    Main <-  if (is.null(dots$main))  colnames(x)[which] else dots$main

    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)
    Main <- rep(Main,length.out=np)
     
    for (i in 1:np) {
        ii <- which[i]        
        dots$main <- Main[i]
        dots$xlab <- xxlab[i]
        dots$ylab <- yylab[i]
        do.call("plot", c(alist(x[, t], x[, ii]), dots))
    }
}


### ============================================================================
 select1dvar <- function (which,var) {

    if (!is.numeric(which)) {
        ln <- length(which)    
        # keep ordering...
        Select <- NULL
        for ( i in 1:ln) {
          ss <- which(which[i]==var)
          if (length(ss) ==0)
            stop("variable ", which[i], " not in var")
          Select <- c(Select,ss)
        }        
    }
    else {
        Select <- which  # "Select now refers to the column number
        if (max(Select) > length(var))
            stop("index in 'which' too large")
        if (min(Select) < 1)
            stop("index in 'which' should be > 0")
    }
  return(Select)
 }


### ============================================================================

plot.ode1D <- function (x, which, ask, add.contour, grid, ...) {

# if x is vector, check if there are enough columns ...
    nspec <- attributes(x)$nspec
    dimens <- attributes(x)$dimens
    proddim <- prod(dimens)

    if ((ncol(x)- nspec*proddim) < 1)
      stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

    if (is.null(which))
      which <- 1:nspec

    var <-  if (! is.null(attributes(x)$ynames))
      attributes(x)$ynames else 1:nspec

    which <- select1dvar(which,var)

    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and
    # interactively wait if there are remaining figures

    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

    Main <-  if (is.null(dots$main)) var else rep(dots$main, length.out =np)

    labs <- (is.null(dots$xlab) && is.null(dots$ylab))
    xxlab <- if (is.null(dots$xlab))  "times"  else dots$xlab
    yylab <- if (is.null(dots$ylab))  ""   else dots$ylab
    dots$col <- if (is.null(dots$col))
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100) else dots$col

    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)
    times <- x[,1]
    for (i in which) {
        dots$main <- Main[i]
        istart <- (i-1)*proddim + 1
        out <- x[,(istart+1):(istart+prod(dimens))]

        dots$xlab <- xxlab[i]
        dots$ylab <- yylab[i]
        
        List <- alist(z=out,x=times)
        if (! is.null(grid)) List$y = grid
        
        do.call("image", c(List, dots)) 
        if (add.contour) do.call("contour", c(List, add=TRUE))
    }
}

### ============================================================================

plot.ode2D <- function (x, which, ask, add.contour, grid, ...) {

# if x is vector, check if there are enough columns ...
    nspec <- attributes(x)$nspec
    dimens <- attributes(x)$dimens
    proddim <- prod(dimens)

    if ((ncol(x)- nspec*proddim) < 1)
      stop("ncol of 'x' should be > 'nspec' * dimens if x is a vector")

    if (is.null(which))
      which <- 1:nspec

    var <-  if (! is.null(attributes(x)$ynames))
      attributes(x)$ynames else 1:nspec

    which <- select1dvar(which,var)

    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and
    # interactively wait if there are remaining figures

    Ask <- setplotpar(nmdots, dots, np, ask)

# here ask is always true...
    if (is.null(ask)) ask <- TRUE
    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
    }
    
    Main <-  if (is.null(dots$main)) var else rep(dots$main, length.out =np)

    labs <- (is.null(dots$xlab) && is.null(dots$ylab))
    xxlab <- if (is.null(dots$xlab))  "times"  else dots$xlab
    yylab <- if (is.null(dots$ylab))  ""   else dots$ylab
    dots$col <- if (is.null(dots$col))
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100) else dots$col

    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)

    for (nt in 1:nrow(x)) {
      for (i in which) {
        dots$main <- Main[i]
        istart <- (i-1)*proddim + 1
        out <- x[nt,(istart+1):(istart+prod(dimens))]
        dim(out) <- dimens
        dots$xlab <- xxlab[i]
        dots$ylab <- yylab[i]

        List <- alist(z=out)
        if (! is.null(grid)) {
          List$x <- grid$x
          List$y <- grid$y
        }  

        do.call("image", c(List, dots)) 
        if (add.contour) do.call("contour", c(List, add=TRUE))
     }
     mtext(outer=TRUE, side=3,paste("time ",x[nt,1]), cex=1.5, line=-1.5)

   } 
}

### ============================================================================

