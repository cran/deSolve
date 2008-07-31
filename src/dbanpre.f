C-----------------------------------------------------------------------
C 
C             Preconditioner Routines for Banded Problems
C                          14 September 1995
C
C The following pair of subroutines -- DBANJA and DBANPS -- provides a
C general-purpose banded preconditioner matrix for use with the DDASPK
C solver, with the Krylov linear system method.  When using DDASPK to
C solve a problem G(t,y,y') = 0, whose iteration matrix (Jacobian)
C    J = dG/dy + c * dG/dy'  (c = scalar)
C is either banded or approximately equal to a banded matrix, these
C routines can be used to generate a banded approximation to J as the
C preconditioner and to solve the resulting banded linear system, in 
C conjunction with the Krylov method option (INFO(12) = 1) in DDASPK.
C
C Other than the user-supplied residual routine RES defining G(t,y,y'),
C the only other inputs required by these routines are the
C half-bandwidth parameters ML and MU of the approximate banded
C Jacobian.  If the system size is NEQ, the half-bandwidths are
C defined as integers between 0 and NEQ - 1 such that only elements
C with indices (i,j) satisfying
C    -ML .le. j - i .le. MU
C are to be retained in the preconditioner.  E.g., if ML = MU = 0, a
C diagonal matrix will be generated as the preconditioner.  The banded
C preconditioner is obtained by difference quotient approximations.  If
C the true problem Jacobian is not banded but is approximately equal to
C a matrix that is banded, the procedure used here will have the effect
C of lumping the elements outside of the band onto the elements within
C the band.
C
C To use these routines in conjunction with DDASPK, the user's calling
C program should include the following, in addition to setting the other
C DDASPK input parameters.
C
C (a) Dimension the array IPAR to have length at least 2, and load the
C     half-bandwidths into IPAR as
C       IPAR(1) = ML   and   IPAR(2) = MU
C     IPAR is used to communicate these parameters to DBANJA and DBANPS.
C     If the user program also uses IPAR for communication with RES,
C     that data should be located beyond the first 2 words of IPAR.
C
C (b) Include the names DBANJA and DBANPS in an EXTERNAL statement.
C     Set INFO(15) = 1 to indicate that a JAC routine exists.
C     Then in the call to DDASPK, pass the names DBANJA and DBANPS as
C     the arguments JAC and PSOL, respectively.
C
C (c) The DDASPK work arrays RWORK and IWORK must include segments WP
C     and IWP for use by DBANJA/DBANPS.  The lengths of these depend on
C     the problem size and half-bandwidths, as follows:
C       LWP =  length of RWORK segment WP = 
C                     (2*ML + MU + 1)*NEQ + 2*( (NEQ/(ML+MU+1)) + 1)
C       LIWP = length of IWORK segment IWP = NEQ
C     (Note the integer divide in LWP.)  Load these lengths in IWORK as
C       IWORK(27) = LWP 
C       IWORK(28) = LIWP
C     and include these values in the declared size of RWORK and IWORK.
C
C
C The DBANJA and DBANPS routines generate and solve the banded
C preconditioner matrix P within the preconditioned Krylov algorithm
C used by DDASPK when INFO(12) = 1.  P is generated and LU-factored 
C periodically during the integration, and the factors are used to 
C solve systems Px = b as needed.
c Karline: changed the rpar and ipar input
c in deSolve, ipar contains : ipar(1)= nout, ipar(2) = length (rpar)
c ipar(3) = length(ipar); ipar(ipar(3)-2), and ipar(ipar(3)-1) contain
c ml and mu
C-----------------------------------------------------------------------
