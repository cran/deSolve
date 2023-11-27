Changes version 1.40
================================
* compacted vignettes

Changes version 1.39
================================
* small documentation fix
* printf fixes in forcings.c and rprintf.c

Changes version 1.37
================================
* changed encoding of DESCIPTION to UTF-8
* fixed also a NOTE about documentation, 
* format of citations

Changes version 1.37
================================
* remove GNU extensions and improve compatibility by 
  replacing KIND=8 with KIND=0.0d0 in Fortran code

Changes version 1.36
================================
* remove GNU extensions to improve Fortran compatibility
  for complex numbers in file zvode.f 

Changes version 1.35
================================
* speedup of multiple `lsoda` solver calls with a pre-identified
  symbol table for compiled models.
  Feature request and idea from Johannes Ranke.
* use new CITATION format
* convert NEWS to NEWS.md in markdown format
* documentation fixes

Changes version 1.34
================================
* fix deprecated function declaration without prototype
* fix LaTeX style issue
* fix netlib.org URL

Changes version 1.33
================================
* back-compatibility with R <= 3.6.2
  (patch provided by Davis Vaughan) 

Changes version 1.32
================================
* `USE_FC_LEN_T` for Fortran BLAS/LAPACK routines
  to ensure compatibility with R > 4.2

Changes version 1.31
================================
* replace deprecated S constant `DOUBLE_XMAX` with `DBL_MAX`

Changes version 1.30
================================
* add doi to references

Changes version 1.29
================================
* replace legacy S macros `PROBLEM` and `ERROR` in C code with `Rf_error()`
* replace http:// with https:// if possible

Changes version 1.28
================================
* minor: avoid an implicit type conversion at the C level

Changes version 1.27.1
================================
* C level change to fix symbol clash on MacOS
  solution kindly provided by Brian Ripley

Changes version 1.27
================================
* C level changes to ensure compatibility with gcc 10
  solution kindly provided by Brian Ripley

Changes version 1.26
================================
* fix outdated class checks to ensure R 4.0.0 compatibility
  (Karline, Thomas)
* Fortran modernization: initialization of variables (Karline)

Changes version 1.25
================================
* add original authors of LINPACK, ODEPACK and SPARSKIT to Authors@R
  (Thomas, Karline, Woody)

Changes version 1.24
================================
* fix compiler warnings to improve Fortran compatibility 
  (Thomas, Karline, Woody)
  thanks to Brian Ripley and Kurt Hornik
* `iteration`: set attribute before calling diagnostics (Karline)
* add open researcher id, ORCID (Thomas)

Changes version 1.23
================================
* add FME to 'Suggests'

Changes version 1.22
================================
* small updates of examples (Thomas)
* improve Fortran compiler compatibility (Thomas)

Changes version 1.21
================================
* change way how PROTECT/UNPROTECT is handled (Thomas)
  (suggested by Tomas Kalibera / R-Core)
* fixed inconsistency in the aphid model example
  (suggested by Sarah Kintner)

Changes version 1.20
================================
* register native routines (Thomas)
* check if event data frame has ordered time (and if not, order)
* change 'event list' to event matrix or data frame in docs
* intentional version jump to indicate chances at the C level

Changes version 1.14
================================
* `matplot.deSolve` is not anymore exported as `matplot` to avoid the 
  respective startup message
* please use `matplot.deSolve` or the alias `matplot.0D` instead  (Thomas)
* small fix that allows parameters in list format for `DLLfunc` and `DLLres`
* a little bit Fortran modification (e.g. avoid `real*8` and `complex*16` types)

Changes version 1.13
================================
* observed data and `plot.deSolve` / `matplot` for multiple outputs  (Karline)
* combining compiled code function with R code event function (Karline)
* check sorting of event times (Karline)
* fix bug related to negative event time (patch supplied by J. Stott)
* relax setting of `tcrit` to make integration with events slightly faster
  (patch from J. Stott)
* adapt maxstep calculation for rk methods,
print a warning if maxsteps is exceeded, fix diagnostics (Thomas)
* more argument checking for rk solvers (Thomas)
* add reference to book of Soetaert, Cash and Mazzia (2012)

Changes version 1.12
================================
* new functions `matplot.deSolve` and `matplot.1D`
* fix valgrind issue (detected by new compilers)
* small improvments of plotting functions
* import standard packages as required by upcoming R versions

Changes version 1.11
================================
* compiledCode vignette now with dede example
* warning and error bug resolved
* Time SEXP incompatibility with R 3.1.1 resolved
* CFunc compatibility (compiled code)

Changes version 1.10.9
================================
* documentation updates, hyperlinks to examples and vignettes
* moved example directories

Changes version 1.10.8
================================
* remove redundant .R files from inst/doc
* fixed bug in event code (patch contributed by Jonathan Stott)

Changes version 1.10.7
================================
* Fortran examples of compiled dede models (Woody)
* vignettes moved to /vignettes
* roles of authors (Authors@R)
* function timestep is now internal
* small documentation updates

Changes version 1.10.6 (Thomas) 
================================
* change declaration of variable dimensions from `(1)` to `(*)` in legacy
  Fortran code to pass automatic bounds check
* remove the Jacobian examples from `?ode` because `banddown=0` can
  lead to problems on some systems; examples will come back in a
  next release
* fixed bug in the `iteration` solver
* small documentation updates

Changes version 1.10.5 (Karline, Thomas)
================================
* extended subset.deSolve with argument `arr`, when TRUE returns an
  array for >2-D output
* fixed the R compiler notes
* `plot.ode.2D` now has an `mtext` argument, via the `...`, to label multiple 
  figures in margin... CHECK - see ode.2D
* `subset` can also be a vector with indices in addition to logical
* image with `legend = TRUE` changed size of plot in different layouts - now solved
  (by adding `par(mar = par("mar"))`)
* new method to output warnings and error messages 
* add data type check for external outputs in rk_util.c
* add interface for compiled dede models
* emphasize consistent order of states in y and return value of func
* changes of Fortran error messages (to be continued)

Changes version 1.10-4 (Thomas, Karline)
================================
* allow reverted time vector for fixed step solvers
    - todo: find solution for dense output methods, and Livermore solvers
* all solvers now have default `atol = 1e-6`; before this `daspk` and `vode` had 1e-8.
* multiple warnings from `daspk` if num steps = 500 toggled off.
* added input argument `nind` to `daspk`, to make it compatible with `radau`.
  this also changes the way the variables are weighed, 
  hence this differs from the original daspk 2.0 code.
* improved warning printing in `daspk` and `vode` 
* extended sparse Jacobian input in `lsodes`. (2-D and 3-D sparsity
  with mapping var and arbitrary sparsity in ian/jan format).

Changes version 1.10-3 (Karline)
================================
* rwork and iwork in lsodes from Fortran -> C (to remove compiler warnings)
* roots + events:  now certain roots can stop simulation + fixed bug in radau root
* improved events\roots help file
* diagnostics(out) gave error in case method=iteration (no rstate) now fixed
* the package authors agreed to assign the maintainer role to T.P.,
but the order of authorship and credits remain unchanged.

Changes version 1.10-2 (Karline)
================================
* remove NAs from forcing functions - when used in DLL (file forcings.R)
* new argument `restructure` in ode.1D, for use with implicit solvers not in deSolve
* removed requirement to have eventfunc in compiled code when func is in compiled code
* subsetting on summary.deSolve

Changes version 1.10-1 (Thomas)
===============================
* remove several redundant variables from C code
* add NEWS file

Changes version 1.10 (Karline, Thomas)
======================================
* compiled code using mass in daspk
* cleanEventTimes

Changes version 1.9+ (Karline)
==============================
* roots, events, lags in radau
* roots in lsodes
* lags in daspk
* ode (method = "iteration")

Changes version 1.9 (Karline, Thomas)
=====================================
* summary.deSolve
* subset.deSolve
* plotting deSolve objects improved:
    - plot more than one output in same figures (scenarios), 
    - add observations
* vignette improved
* fixed bug in 'timesteps'

Changes version 1.8.1 (Thomas, Woody, Karline)
==============================================
* fixed compiler warnings using valgrind
* fixed compiler warning C-code

Changes version 1.8 (Thomas)
============================
* Dormand-Prince 8(7) coefficients use now common 
  instead of decimal fractions


Changes version 1.8 (Karline)
=============================
* Runge-Kuttas:
    - extra output: number of failed steps (see also 2)
    - number of function evaluations + 1 for initial condition
    - dense output for cash-karp
    - dopri8(7) added
    - radau added!! implicit runge kutta, solves also DAE up to index 3!
* other:
    - image function for ode.2-D added.
    - changed warning printing in FORTRAN code
    - common interface for radau and daspk:
      both can solve systems written as M*dy = f(x,y).
      daspk can also solve systems written as 0 = g(x,y,dy) (=default for daspk)



