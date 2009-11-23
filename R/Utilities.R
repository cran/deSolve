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

selectvar <- function (which,var) {

    if (!is.numeric(which)) {
        ln <- length(which)
        Select <- which(var %in% which)
        if (length(Select) != ln)
            stop("not all variables in 'which' are in 'x'")
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

plot.deSolve <- function (x, which = 1:(ncol(x)-1), ask = NULL, ...) {
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
        do.call("plot", c(alist(x[, t], x[, i]), dots))
    }
}

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
