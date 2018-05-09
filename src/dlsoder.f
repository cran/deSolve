*DECK DLSODER
C DLSODER was created by merging DLSODE with DLSODAR - Karline Soetaert
      SUBROUTINE DLSODER (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1                  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF,
     2                  G, NG, JROOT, rpar, ipar)
      IMPLICIT NONE
      EXTERNAL F, JAC, G
CKS: added rpar, ipar, and G
      INTEGER ipar(*)
      DOUBLE PRECISION rpar(*)

      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      INTEGER NG, JROOT
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),
     1   JROOT(NG)

C-----------------------------------------------------------------------
C***BEGIN PROLOGUE  DLSODER
C***PURPOSE  Livermore Solver for Ordinary Differential Equations.
C            DLSODER solves the initial-value problem for stiff or
C            nonstiff systems of first-order ODE's,
C               dy/dt = f(t,y),   or, in component form,
C               dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(N)),  i=1,...,N.
C            and with Root-finding.

C***CATEGORY  I1A
C***TYPE      DOUBLE PRECISION (SLSODE-S, DLSODE-D)
C***KEYWORDS  ORDINARY DIFFERENTIAL EQUATIONS, INITIAL VALUE PROBLEM,
C             STIFF, NONSTIFF
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C             Center for Applied Scientific Computing, L-561
C             Lawrence Livermore National Laboratory
C             Livermore, CA 94551.
C           Root function added by Karline Soetaert
C***DESCRIPTION - see DLSODE and DLSODAR
C  Note: length of RWORK array = 20 + NYH*(MAXORD + 1) + 3*NEQ + LWM + 3*NG.
C
C-----------------------------------------------------------------------


C  Declare externals.
      EXTERNAL DPREPJ, DSOLSY
      DOUBLE PRECISION DUMACH, DVNORM
C
C  Declare all other variables.
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C KS: added next line
      INTEGER LG0, LG1, LGX, IOWNR3, IRFND, ITASKC, NGC, NGE
      INTEGER I, I1, I2, IFLAG, IMXER, KGO, LF0,
     1   LENIW, LENRW, LENWM, ML, MORD, MU, MXHNL0, MXSTP0
      INTEGER IRFP, IRT, LENYH, LYHNEW
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION ROWNR3, T0, TLAST, TOUTC
      DOUBLE PRECISION ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER(LEN=80) MSG
      SAVE MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
C The following internal Common block contains
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine ("own" variables), and
C (b) variables which are communicated between subroutines.
C The block DLS001 is declared in subroutines DLSODER, DINTDY, DSTODE,
C DPREPJ, and DSOLSY.
C The block DLSR01 is declared in subroutines DLSODAR, DRCHEK, DROOTS.
C Groups of variables are replaced by dummy arrays in the Common
C declarations in routines where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C karline: added next common block
      COMMON /DLSR01/ ROWNR3(2), T0, TLAST, TOUTC,
     1   LG0, LG1, LGX, IOWNR3(2), IRFND, ITASKC, NGC, NGE
C
      DATA  MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .GT. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT  DLSODER
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
C Karline: added nest sentence
      ITASKC = ITASK
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT,
C MF, ML, MU, and NG..
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0 .OR. MITER .GT. 5) GO TO 608
      IF (MITER .LE. 3) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
C karline: added next four lines
      IF (NG .LT. 0) GO TO 630
      IF (ISTATE .EQ. 1) GO TO 35
      IF (IRFND .EQ. 0 .AND. NG .NE. NGC) GO TO 631
 35   NGC = NG
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = 0.0D0
      HMXI = 0.0D0
      HMIN = 0.0D0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0D0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0D0) GO TO 615
      HMXI = 0.0D0
      IF (HMAX .GT. 0.0D0) HMXI = 1.0D0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0D0) GO TO 616
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted  G0, G1, GX, YH, WM,
C EWT, SAVF, ACOR.
C-----------------------------------------------------------------------
CKS: init changes 60   LYH = 21
60    IF (ISTATE .EQ. 1) NYH = N
      LG0 = 21
      LG1 = LG0 + NG
      LGX = LG1 + NG
      LYHNEW = LGX + NG
      IF (ISTATE .EQ. 1) LYH = LYHNEW
      IF (LYHNEW .EQ. LYH) GO TO 62
C If ISTATE = 3 and NG was changed, shift YH to its new location. ------
      LENYH = L*NYH
      IF (LRW .LT. LYHNEW-1+LENYH) GO TO 62
      I1 = 1
      IF (LYHNEW .GT. LYH) I1 = -1
      CALL DCOPY (LENYH, RWORK(LYH), I1, RWORK(LYHNEW), I1)
      LYH = LYHNEW
62    CONTINUE
CKS end of changes
      LWM = LYH + (MAXORD + 1)*NYH
      IF (MITER .EQ. 0) LENWM = 0
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) LENWM = N*N + 2
      IF (MITER .EQ. 3) LENWM = N + 2
      IF (MITER .GE. 4) LENWM = (2*ML + MU + 1)*N + 2
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LACOR = LSAVF + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 20 + N
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) LENIW = 20
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0D0) GO TO 619
        IF (ATOLI .LT. 0.0D0) GO TO 620
 70     CONTINUE
      IF (ISTATE .EQ. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to DSTODE. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
C MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      DO 80 I = 1,N
        RWORK(I+LSAVF-1) = RWORK(I+LWM-1)
 80   CONTINUE
C Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
 90   IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
        RWORK(I) = 0.0D0
 95   CONTINUE
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = DUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0D0) GO TO 625
      IF (H0 .NE. 0.0D0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0D0)
     1   H0 = TCRIT - T
 110  JSTART = 0
      IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      NHNIL = 0
      NST = 0
      NJE = 0
      NSLAST = 0
      HU = 0.0D0
      NQU = 0
      CCMAX = 0.3D0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (NEQ, T, Y, RWORK(LF0), rpar, ipar)
      NFE = 1
C Load the initial value vector in YH. ---------------------------------
      DO 115 I = 1,N
        RWORK(I+LYH-1) = Y(I)
 115  CONTINUE
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = 1.0D0
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 621
        RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
 120  CONTINUE
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(I))
C if this is positive, or MAX(ATOL(I)/ABS(Y(I))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by..
C                                      NEQ
C   H0**2 = TOL / ( w0**-2 + (1/NEQ) * SUM ( f(i)/ywt(i) )**2  )
C                                       1
C where   w0     = MAX ( ABS(T), ABS(TOUT) ),
C         f(i)   = i-th component of initial value of f,
C         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
C The sign of H0 is inferred from the initial values of TOUT and T.
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0D0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0D0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 140
      DO 130 I = 1,N
        TOL = MAX(TOL,RTOL(I))
 130  CONTINUE
 140  IF (TOL .GT. 0.0D0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0D0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0D0*UROUND)
      TOL = MIN(TOL,0.001D0)
      SUM = DVNORM (N, RWORK(LF0), RWORK(LEWT))
      SUM = 1.0D0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0D0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0D0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
        RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
 190  CONTINUE
CKS: start changes      GO TO 270
C
C Check for a zero of g at T. ------------------------------------------
      IRFND = 0
      TOUTC = TOUT
      IF (NGC .EQ. 0) GO TO 270
      CALL DRCHEK (1, G, NEQ, Y, RWORK(LYH), NYH,
     1   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT, rpar, ipar)
      IF (IRT .EQ. 0) GO TO 270
      GO TO 632
CKS: end changes
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C First, DRCHEK is called to check for a root within the last step
C taken, other than the last root found there, if any.
C If ITASK = 2 or 5, and y(TN) has not yet been returned to the user
C because of an intervening root, return through Block G.
C-----------------------------------------------------------------------
 200  NSLAST = NST
C karline: added from here
      IRFP = IRFND
      IF (NGC .EQ. 0) GO TO 205
      IF (ITASK .EQ. 1 .OR. ITASK .EQ. 4) TOUTC = TOUT
      CALL DRCHEK (2, G, NEQ, Y, RWORK(LYH), NYH,
     1   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT, rpar, ipar)
      IF (IRT .NE. 1) GO TO 205
      IRFND = 1
      ISTATE = 3
      T = T0
      GO TO 425
 205  CONTINUE
      IRFND = 0
      IF (IRFP .EQ. 1 .AND. TLAST .NE. TN .AND. ITASK .EQ. 2) GO TO 400
C karline: till here
      IF (ITASK .EQ. 1) THEN
        GOTO 210
      ELSE IF (ITASK .EQ. 2) THEN
        GOTO 250
      ELSE IF (ITASK .EQ. 3) THEN
        GOTO 220
      ELSE IF (ITASK .EQ. 4) THEN
        GOTO 230
      ELSE IF (ITASK .EQ. 5) THEN
        GOTO 240
      ENDIF
C      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0D0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0D0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 245
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
C this in lsoda, not in lsode...       IF (IHIT) T = TCRIT + karline added next line
      IF (IRFP .EQ. 1 .AND. TLAST .NE. TN .AND. ITASK .EQ. 5) GO TO 400
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator DSTODE.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 510
        RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
 260    CONTINUE
 270  TOLSF = UROUND*DVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0D0) GO TO 280
      TOLSF = TOLSF*2.0D0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'DLSODER-  Warning..internal T (=R1) and H (=R2) are'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWD (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      (H = step size). Solver will continue anyway'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'DLSODER-  Above warning has been issued I1 times.  '
      CALL XERRWD (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      It will not be issued again for this problem'
      CALL XERRWD (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
 290  CONTINUE
C-----------------------------------------------------------------------
C  CALL DSTODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPREPJ,DSOLSY)
C-----------------------------------------------------------------------
      CALL DSTODE (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), IWORK(LIWM),
     2   F, JAC, DPREPJ, DSOLSY, rpar,ipar)
      KGO = 1 - KFLAG
      IF (KGO .EQ. 1) THEN
        GOTO 300
      ELSE IF (KGO .EQ. 2) THEN
        GOTO 530
      ELSE IF (KGO .EQ. 3) THEN
        GOTO 540
      ENDIF
      
C      GO TO (300, 530, 540), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C Then call DRCHEK to check for a root within the last step.
C Then, if no root was found, check for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
Ckarline: changed this
      IF (NGC .EQ. 0) GO TO 315
      CALL DRCHEK (3, G, NEQ, Y, RWORK(LYH), NYH,
     1   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT, rpar, ipar)
      IF (IRT .NE. 1) GO TO 315
      IRFND = 1
      ISTATE = 3
      T = T0
      GO TO 425
 315  CONTINUE
C karline: end of changes
      IF (ITASK .EQ. 1) THEN
        GOTO 320
      ELSE IF (ITASK .EQ. 2) THEN
        GOTO 400
      ELSE IF (ITASK .EQ. 3) THEN
        GOTO 330
      ELSE IF (ITASK .EQ. 4) THEN
        GOTO 340
      ELSE IF (ITASK .EQ. 5) THEN
        GOTO 350
      ENDIF
C      GO TO (320, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 320  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0D0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 345
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from DLSODER.
C If ITASK .NE. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
        Y(I) = RWORK(I+LYH-1)
 410  CONTINUE
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
 425  CONTINUE
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
C karline: added next two lines
      IWORK(10) = NGE
      TLAST = T
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.  The optional outputs
C are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'DLSODER-  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWD (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWD (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0D0)
      ISTATE = -1
      GO TO 580
C EWT(I) .LE. 0.0 for some I (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSODER-  At T (=R1), EWT(I1) has become R2 .LE. 0.'
      CALL XERRWD (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'DLSODER-  At T (=R1), too much accuracy requested  '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      for precision of machine..  see TOLSF (=R2) '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'DLSODER-  At T(=R1) and step size H(=R2), the error'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'DLSODER-  At T (=R1) and step size H (=R2), the    '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWD (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = 0.0D0
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional outputs. -------------------------------
 580  DO 590 I = 1,N
        Y(I) = RWORK(I+LYH-1)
 590  CONTINUE
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
C Karline: added next two lines
      IWORK(10) = NGE
      TLAST = T
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'DLSODER-  ISTATE (=I1) illegal '
      CALL XERRWD (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0D0, 0.0D0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'DLSODER-  ITASK (=I1) illegal  '
      CALL XERRWD (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 603  MSG = 'DLSODER-  ISTATE .GT. 1 but DLSODER not initialized '
      CALL XERRWD (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 604  MSG = 'DLSODER-  NEQ (=I1) .LT. 1     '
      CALL XERRWD (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 605  MSG = 'DLSODER-  ISTATE = 3 and NEQ increased (I1 to I2)  '
      CALL XERRWD (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 606  MSG = 'DLSODER-  ITOL (=I1) illegal   '
      CALL XERRWD (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 607  MSG = 'DLSODER-  IOPT (=I1) illegal   '
      CALL XERRWD (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 608  MSG = 'DLSODER-  MF (=I1) illegal     '
      CALL XERRWD (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 609  MSG = 'DLSODER-  ML (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)'
      CALL XERRWD (MSG, 50, 9, 0, 2, ML, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 610  MSG = 'DLSODER-  MU (=I1) illegal.. .LT.0 or .GE.NEQ (=I2)'
      CALL XERRWD (MSG, 50, 10, 0, 2, MU, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 611  MSG = 'DLSODER-  MAXORD (=I1) .LT. 0  '
      CALL XERRWD (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 612  MSG = 'DLSODER-  MXSTEP (=I1) .LT. 0  '
      CALL XERRWD (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 613  MSG = 'DLSODER-  MXHNIL (=I1) .LT. 0  '
      CALL XERRWD (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 614  MSG = 'DLSODER-  TOUT (=R1) behind T (=R2)      '
      CALL XERRWD (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWD (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0D0)
      GO TO 700
 615  MSG = 'DLSODER-  HMAX (=R1) .LT. 0.0  '
      CALL XERRWD (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0D0)
      GO TO 700
 616  MSG = 'DLSODER-  HMIN (=R1) .LT. 0.0  '
      CALL XERRWD (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0D0)
      GO TO 700
 617  CONTINUE
      MSG='DLSODER- RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
      GO TO 700
 618  CONTINUE
      MSG='DLSODER- IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWD (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)
      GO TO 700
 619  MSG = 'DLSODER-  RTOL(I1) is R1 .LT. 0.0        '
      CALL XERRWD (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0D0)
      GO TO 700
 620  MSG = 'DLSODER-  ATOL(I1) is R1 .LT. 0.0        '
      CALL XERRWD (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0D0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSODER-  EWT(I1) is R1 .LE. 0.0         '
      CALL XERRWD (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0D0)
      GO TO 700
 622  CONTINUE
      MSG='DLSODER- TOUT (=R1) too close to T(=R2) to start integration'
      CALL XERRWD (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  CONTINUE
      MSG='DLSODER- ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWD (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  CONTINUE
      MSG='DLSODER- ITASK = 4 OR 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWD (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  CONTINUE
      MSG='DLSODER- ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWD (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'DLSODER-  At start of problem, too much accuracy   '
      CALL XERRWD (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWD (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0D0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'DLSODER-  Trouble in DINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWD (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0D0)
C Karline: added next error messages
 630  MSG = 'DLSODER-  NG (=I1) .lt. 0     '
      CALL XERRWD (MSG, 30, 30, 0, 1, NG, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 631  MSG = 'DLSODER-  NG changed (from I1 to I2) illegally,   '
      CALL XERRWD (MSG, 50, 31, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      i.e. not immediately after a root was found.'
      CALL XERRWD (MSG, 50, 31, 0, 2, NGC, NG, 0, 0.0D0, 0.0D0)
      GO TO 700
 632  MSG = 'DLSODER-  One or more components of g has a root  '
      CALL XERRWD (MSG, 50, 32, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      too near to the initial point.    '
      CALL XERRWD (MSG, 40, 32, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
C
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'DLSODER-  Run aborted.. apparent infinite loop     '
      CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)
      RETURN
C----------------------- END OF SUBROUTINE DLSODER ----------------------
      END


*DECK DLSODESR
C DLSODESR was created by merging DLSODES with DLSODAR - Karline Soetaert
      SUBROUTINE DLSODESR (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF, 
     2            G, NG, JROOT, rpar, ipar)
      IMPLICIT NONE
      EXTERNAL F, JAC, G
CKS: added rpar, ipar, and G
      INTEGER ipar(*)
      DOUBLE PRECISION rpar(*)      
      
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF
      INTEGER NG, JROOT
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),
     &   JROOT(NG)

C-----------------------------------------------------------------------
C DLSODES solves the initial value problem for stiff or nonstiff
C systems of first order ODEs,
C     dy/dt = f(t,y) ,  or, in component form,
C     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
C DLSODES is a variant of the DLSODE package, and is intended for
C problems in which the Jacobian matrix df/dy has an arbitrary
C sparse structure (when the problem is stiff).
C
C Authors:       Alan C. Hindmarsh
C                Center for Applied Scientific Computing, L-561
C                Lawrence Livermore National Laboratory
C                Livermore, CA 94551
C and
C                Andrew H. Sherman
C                J. S. Nolen and Associates
C                Houston, TX 77084
C
C           Root function added by Karline Soetaert
C
C***DESCRIPTION - see DLSODES
C  Note: length of RWORK array = 20 + NYH*(MAXORD + 1) + 3*NEQ + LWM + 3*NG.
C-----------------------------------------------------------------------
      EXTERNAL DPRJS, DSOLSS
      DOUBLE PRECISION DUMACH, DVNORM
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     1   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     2   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     3   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
C KS: added next lines
      INTEGER LG0, LG1, LGX, IOWNR3, IRFND, ITASKC, NGC, NGE
      INTEGER IRFP, IRT, LYHNEW
      DOUBLE PRECISION ROWNR3, T0, TLAST, TOUTC
      INTEGER I, I1, I2, IFLAG, IMAX, IMUL, IMXER, IPFLAG, IPGO, IREM,
     1   J, KGO, LENRAT, LENYHT, LENIW, LENRW, LF0, LIA, LJA,
     2   LRTEM, LWTEM, LYHD, LYHN, MF1, MORD, MXHNL0, MXSTP0, NCOLM
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH
      DOUBLE PRECISION ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER(LEN=60) MSG
      SAVE LENRAT, MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
C The following two internal Common blocks contain
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine ("own" variables), and
C (b) variables which are communicated between subroutines.
C The block DLS001 is declared in subroutines DLSODES, DIPREP, DPREP,
C DINTDY, DSTODE, DPRJS, and DSOLSS.
C The block DLSS01 is declared in subroutines DLSODES, DIPREP, DPREP,
C DPRJS, and DSOLSS.
C Groups of variables are replaced by dummy arrays in the Common
C declarations in routines where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      COMMON /DLSS01/ CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH,
     1   IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP,
     2   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,
     3   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,
     4   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU
C karline:
      COMMON /DLSR01/ ROWNR3(2), T0, TLAST, TOUTC,
     1   LG0, LG1, LGX, IOWNR3(2), IRFND, ITASKC, NGC, NGE
C
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
C-----------------------------------------------------------------------
C In the Data statement below, set LENRAT equal to the ratio of
C the wordlength for a real number to that for an integer.  Usually,
C LENRAT = 1 for single precision and 2 for double precision.  If the
C true ratio is not an integer, use the next smaller integer (.ge. 1).
C-----------------------------------------------------------------------
      DATA LENRAT/2/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
C Karline: added nest sentence
      ITASKC = ITASK
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C If ISTATE = 1, the final setting of work space pointers, the matrix
C preprocessing, and other initializations are done in Block C.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT,
C MF, ML, MU, and NG.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      MOSS = MF/100
      MF1 = MF - 100*MOSS
      METH = MF1/10
      MITER = MF1 - 10*METH
      IF (MOSS .LT. 0 .OR. MOSS .GT. 2) GO TO 608
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0 .OR. MITER .GT. 3) GO TO 608
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) MOSS = 0
C Karline: start add
      IF (NG .LT. 0) GO TO 680
      IF (ISTATE .EQ. 1) GO TO 35
      IF (IRFND .EQ. 0 .AND. NG .NE. NGC) GO TO 681
 35   NGC = NG
C Karline: end added
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = 0.0D0
      HMXI = 0.0D0
      HMIN = 0.0D0
      SETH = 0.0D0
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. 0.0D0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0D0) GO TO 615
      HMXI = 0.0D0
      IF (HMAX .GT. 0.0D0) HMXI = 1.0D0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0D0) GO TO 616
      SETH = RWORK(8)
      IF (SETH .LT. 0.0D0) GO TO 609
C Check RTOL and ATOL for legality. ------------------------------------
 60   RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 65 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0D0) GO TO 619
        IF (ATOLI .LT. 0.0D0) GO TO 620
 65     CONTINUE
C-----------------------------------------------------------------------
C Compute required work array lengths, as far as possible, and test
C these against LRW and LIW.  Then set tentative pointers for work
C arrays.  Pointers to RWORK/IWORK segments are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted  G0, G1, GX, WM, YH, SAVF, EWT, ACOR.
C If MITER = 1 or 2, the required length of the matrix work space WM
C is not yet known, and so a crude minimum value is used for the
C initial tests of LRW and LIW, and YH is temporarily stored as far
C to the right in RWORK as possible, to leave the maximum amount
C of space for WM for matrix preprocessing.  Thus if MITER = 1 or 2
C and MOSS .ne. 2, some of the segments of RWORK are temporarily
C omitted, as they are not needed in the preprocessing.  These
C omitted segments are: ACOR if ISTATE = 1, EWT and ACOR if ISTATE = 3
C and MOSS = 1, and SAVF, EWT, and ACOR if ISTATE = 3 and MOSS = 0.
C-----------------------------------------------------------------------
      LRAT = LENRAT
      IF (ISTATE .EQ. 1) NYH = N
      LWMIN = 0
      IF (MITER .EQ. 1) LWMIN = 4*N + 10*N/LRAT
      IF (MITER .EQ. 2) LWMIN = 4*N + 11*N/LRAT
      IF (MITER .EQ. 3) LWMIN = N + 2
      LENYH = (MAXORD+1)*NYH
      LREST = LENYH + 3*N
      LENRW = 20 + LWMIN + LREST
      IWORK(17) = LENRW
      LENIW = 30
      IF (MOSS .EQ. 0 .AND. MITER .NE. 0 .AND. MITER .NE. 3)
     1   LENIW = LENIW + N + 1
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
      LIA = 31
      IF (MOSS .EQ. 0 .AND. MITER .NE. 0 .AND. MITER .NE. 3)
     1   LENIW = LENIW + IWORK(LIA+N) - 1
      IWORK(18) = LENIW
      IF (LENIW .GT. LIW) GO TO 618
      LJA = LIA + N + 1
      LIA = MIN(LIA,LIW)
      LJA = MIN(LJA,LIW)

C      LWM = 21
C Karline: start changes
      IF (ISTATE .EQ. 1) NYH = N
      LG0 = 21
      LG1 = LG0 + NG
      LGX = LG1 + NG
      LYHNEW = LGX + NG
      IF (ISTATE .EQ. 1) LYH = LYHNEW
      IF (LYHNEW .EQ. LYH) GO TO 67
C If ISTATE = 3 and NG was changed, shift YH to its new location. ------
      LENYH = L*NYH
      IF (LRW .LT. LYHNEW-1+LENYH) GO TO 67
      I1 = 1
      IF (LYHNEW .GT. LYH) I1 = -1
      CALL DCOPY (LENYH, RWORK(LYH), I1, RWORK(LYHNEW), I1)
      LYH = LYHNEW
67    CONTINUE
CKS end of changes
      LWM = LYHNEW
      IF (ISTATE .EQ. 1) NQ = 1
      NCOLM = MIN(NQ+1,MAXORD+2)
      LENYHM = NCOLM*NYH
      LENYHT = LENYH
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) LENYHT = LENYHM
      IMUL = 2
      IF (ISTATE .EQ. 3) IMUL = MOSS
      IF (MOSS .EQ. 2) IMUL = 3
      LRTEM = LENYHT + IMUL*N
      LWTEM = LWMIN
C      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) LWTEM = LRW - 20 - LRTEM
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) LWTEM = LRW -(LWM-1)-LRTEM
      LENWK = LWTEM
      LYHN = LWM + LWTEM
      LSAVF = LYHN + LENYHT
      LEWT = LSAVF + N
      LACOR = LEWT + N
      ISTATC = ISTATE
      IF (ISTATE .EQ. 1) GO TO 100
C-----------------------------------------------------------------------
C ISTATE = 3.  Move YH to its new location.
C Note that only the part of YH needed for the next step, namely
C MIN(NQ+1,MAXORD+2) columns, is actually moved.
C A temporary error weight array EWT is loaded if MOSS = 2.
C Sparse matrix processing is done in DIPREP/DPREP if MITER = 1 or 2.
C If MAXORD was reduced below NQ, then the pointers are finally set
C so that SAVF is identical to YH(*,MAXORD+2).
C-----------------------------------------------------------------------

      LYHD = LYH - LYHN
      IMAX = LYHN - 1 + LENYHM
C Move YH.  Move right if LYHD < 0; move left if LYHD > 0. -------------
      IF (LYHD .LT. 0) THEN
        DO 72 I = LYHN,IMAX
          J = IMAX + LYHN - I
          RWORK(J) = RWORK(J+LYHD)
 72     CONTINUE
      ENDIF
      IF (LYHD .GT. 0) THEN
        DO 76 I = LYHN,IMAX
          RWORK(I) = RWORK(I+LYHD)
 76     CONTINUE
      ENDIF
      LYH = LYHN
      IWORK(22) = LYH
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) GO TO 92
      IF (MOSS .NE. 2) GO TO 85
C Temporarily load EWT if MITER = 1 or 2 and MOSS = 2. -----------------
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 82 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 621
        RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
 82   CONTINUE
 85   CONTINUE
C DIPREP and DPREP do sparse matrix preprocessing if MITER = 1 or 2. ---
      LSAVF = MIN(LSAVF,LRW)
      LEWT = MIN(LEWT,LRW)
      LACOR = MIN(LACOR,LRW)
CKS
      CALL DIPREP (NEQ, Y, RWORK, IWORK(LIA),IWORK(LJA), IPFLAG, F, JAC,
     & rpar, ipar )
      LENRW = LWM - 1 + LENWK + LREST
      IWORK(17) = LENRW
      IF (IPFLAG .NE. -1) IWORK(23) = IPIAN
      IF (IPFLAG .NE. -1) IWORK(24) = IPJAN
      IPGO = -IPFLAG + 1

      IF (IPGO .EQ. 1) THEN
        GOTO 90 
      ELSE IF (IPGO .EQ. 2) THEN
        GOTO 628
      ELSE IF (IPGO .EQ. 3) THEN
        GOTO 629
      ELSE IF (IPGO .EQ. 4) THEN
        GOTO 630
      ELSE IF (IPGO .EQ. 5) THEN
        GOTO 631
      ELSE IF (IPGO .EQ. 6) THEN
        GOTO 632
      ELSE IF (IPGO .EQ. 7) THEN
        GOTO 633
      ENDIF           
C      GO TO (90, 628, 629, 630, 631, 632, 633), IPGO
      
 90   IWORK(22) = LYH
      IF (LENRW .GT. LRW) GO TO 617
C Set flag to signal parameter changes to DSTODE. ----------------------
 92   JSTART = -1
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
        RWORK(I) = 0.0D0
 95   CONTINUE
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C the sparse matrix preprocessing (MITER = 1 or 2), and the
C calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  CONTINUE
      LYH = LYHN
      IWORK(22) = LYH
      TN = T
      NST = 0
      H = 1.0D0
      NNZ = 0
      NGP = 0
      NZL = 0
      NZU = 0
C Load the initial value vector in YH. ---------------------------------
      DO 105 I = 1,N
        RWORK(I+LYH-1) = Y(I)
 105  CONTINUE
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (NEQ, T, Y, RWORK(LF0), rpar, ipar)
      NFE = 1
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 110 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 621
        RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
 110  CONTINUE
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) GO TO 120
C DIPREP and DPREP do sparse matrix preprocessing if MITER = 1 or 2. ---
      LACOR = MIN(LACOR,LRW)
      CALL DIPREP (NEQ, Y, RWORK, IWORK(LIA),IWORK(LJA), IPFLAG, F, JAC,
     & rpar, ipar)
      LENRW = LWM - 1 + LENWK + LREST
      IWORK(17) = LENRW
      IF (IPFLAG .NE. -1) IWORK(23) = IPIAN
      IF (IPFLAG .NE. -1) IWORK(24) = IPJAN
      IPGO = -IPFLAG + 1
      IF (IPGO .EQ. 1) THEN
        GOTO 115
      ELSE IF (IPGO .EQ. 2) THEN
        GOTO 628
      ELSE IF (IPGO .EQ. 3) THEN
        GOTO 629
      ELSE IF (IPGO .EQ. 4) THEN
        GOTO 630
      ELSE IF (IPGO .EQ. 5) THEN
        GOTO 631
      ELSE IF (IPGO .EQ. 6) THEN
        GOTO 632
      ELSE IF (IPGO .EQ. 7) THEN
        GOTO 633
      ENDIF
C      GO TO (115, 628, 629, 630, 631, 632, 633), IPGO
 115  IWORK(22) = LYH
      IF (LENRW .GT. LRW) GO TO 617
C Check TCRIT for legality (ITASK = 4 or 5). ---------------------------
 120  CONTINUE
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 125
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0D0) GO TO 625
      IF (H0 .NE. 0.0D0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0D0)
     1   H0 = TCRIT - T
C Initialize all remaining parameters. ---------------------------------
 125  UROUND = DUMACH()
      JSTART = 0
      IF (MITER .NE. 0) RWORK(LWM) = SQRT(UROUND)
      MSBJ = 50
      NSLJ = 0
      CCMXJ = 0.2D0
      PSMALL = 1000.0D0*UROUND
      RBIG = 0.01D0/PSMALL
      NHNIL = 0
      NJE = 0
      NLU = 0
      NSLAST = 0
      HU = 0.0D0
      NQU = 0
      CCMAX = 0.3D0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
C if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by..
C                                      NEQ
C   H0**2 = TOL / ( w0**-2 + (1/NEQ) * Sum ( f(i)/ywt(i) )**2  )
C                                       1
C where   w0     = MAX ( ABS(T), ABS(TOUT) ),
C         f(i)   = i-th component of initial value of f,
C         ywt(i) = EWT(i)/TOL  (a weight for y(i)).
C The sign of H0 is inferred from the initial values of TOUT and T.
C ABS(H0) is made .le. ABS(TOUT-T) in any case.
C-----------------------------------------------------------------------
      LF0 = LYH + NYH
      IF (H0 .NE. 0.0D0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0D0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 140
      DO 130 I = 1,N
        TOL = MAX(TOL,RTOL(I))
 130  CONTINUE
 140  IF (TOL .GT. 0.0D0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0D0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0D0*UROUND)
      TOL = MIN(TOL,0.001D0)
      SUM = DVNORM (N, RWORK(LF0), RWORK(LEWT))
      SUM = 1.0D0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0D0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0D0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
        RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
 190  CONTINUE
CKS: start changes      GO TO 270
C
C Check for a zero of g at T. ------------------------------------------
      IRFND = 0
      TOUTC = TOUT
      IF (NGC .EQ. 0) GO TO 270
      CALL DRCHEK (1, G, NEQ, Y, RWORK(LYH), NYH,
     1   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT, rpar, ipar)
      IF (IRT .EQ. 0) GO TO 270
      GO TO 682
CKS: end changes
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C First, DRCHEK is called to check for a root within the last step
C taken, other than the last root found there, if any.
C If ITASK = 2 or 5, and y(TN) has not yet been returned to the user
C because of an intervening root, return through Block G.
C-----------------------------------------------------------------------
 200  NSLAST = NST
C karline: added from here
      IRFP = IRFND
      IF (NGC .EQ. 0) GO TO 205
      IF (ITASK .EQ. 1 .OR. ITASK .EQ. 4) TOUTC = TOUT
      CALL DRCHEK (2, G, NEQ, Y, RWORK(LYH), NYH,
     1   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT, rpar, ipar)
      IF (IRT .NE. 1) GO TO 205
      IRFND = 1
      ISTATE = 3
      T = T0
      GO TO 425
 205  CONTINUE
      IRFND = 0
      IF (IRFP .EQ. 1 .AND. TLAST .NE. TN .AND. ITASK .EQ. 2) GO TO 400
C karline: till here
      IF (ITASK .EQ. 1) THEN
        GOTO 210
      ELSE IF (ITASK .EQ. 2) THEN
        GOTO 250
      ELSE IF (ITASK .EQ. 3) THEN
        GOTO 220
      ELSE IF (ITASK .EQ. 4) THEN
        GOTO 230
      ELSE IF (ITASK .EQ. 5) THEN
        GOTO 240
      ENDIF

C      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0D0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0D0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 245
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
C karline:added next line
      IF (IRFP .EQ. 1 .AND. TLAST .NE. TN .AND. ITASK .EQ. 5) GO TO 400
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      IF (ISTATE .EQ. 2) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator DSTODE.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 510
        RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
 260  CONTINUE
 270  TOLSF = UROUND*DVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0D0) GO TO 280
      TOLSF = TOLSF*2.0D0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'DLSODES- Warning..Internal T (=R1) and H (=R2) are'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWD (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'DLSODES- Above warning has been issued I1 times.  '
      CALL XERRWD (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWD (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
 290  CONTINUE
C-----------------------------------------------------------------------
C    CALL DSTODE(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,WM,F,JAC,DPRJS,DSOLSS)
C-----------------------------------------------------------------------
      CALL DSTODE (NEQ, Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), RWORK(LACOR), RWORK(LWM), RWORK(LWM),
     2   F, JAC, DPRJS, DSOLSS, rpar,ipar)
      KGO = 1 - KFLAG
      IF (KGO .EQ. 1) THEN
        GOTO 300
      ELSE IF (KGO .EQ. 2) THEN
        GOTO 530
      ELSE IF (KGO .EQ. 3) THEN
        GOTO 540
      ELSE IF (KGO .EQ. 4) THEN
        GOTO 550
      ENDIF
C      GO TO (300, 530, 540, 550), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).  Test for stop conditions.
C Then call DRCHEK to check for a root within the last step.
C Then, if no root was found, check for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
Ckarline: changed this
      IF (NGC .EQ. 0) GO TO 305
      CALL DRCHEK (3, G, NEQ, Y, RWORK(LYH), NYH,
     1   RWORK(LG0), RWORK(LG1), RWORK(LGX), JROOT, IRT, rpar, ipar)
      IF (IRT .NE. 1) GO TO 305
      IRFND = 1
      ISTATE = 3
      T = T0
      GO TO 425
 305  CONTINUE
      IF (ITASK .EQ. 1) THEN
        GOTO 310
      ELSE IF (ITASK .EQ. 2) THEN
        GOTO 400
      ELSE IF (ITASK .EQ. 3) THEN
        GOTO 330
      ELSE IF (ITASK .EQ. 4) THEN
        GOTO 340
      ELSE IF (ITASK .EQ. 5) THEN
        GOTO 350
      ENDIF
C GO TO (310, 400, 330, 340, 350), ITASK
      
C ITASK = 1.  if TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0D0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 345
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from DLSODES.
C If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
        Y(I) = RWORK(I+LYH-1)
 410  CONTINUE
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
 425  CONTINUE
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NNZ
      IWORK(20) = NGP
      IWORK(21) = NLU
      IWORK(25) = NZL
      IWORK(26) = NZU
C karline: added next two lines
      IWORK(10) = NGE
      TLAST = T
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.
C The optional outputs are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'DLSODES- At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWD (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWD (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0D0)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSODES- At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWD (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'DLSODES- At T (=R1), too much accuracy requested  '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'DLSODES- At T(=R1) and step size H(=R2), the error'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'DLSODES- At T (=R1) and step size H (=R2), the    '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWD (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 560
C KFLAG = -3.  Fatal error flag returned by DPRJS or DSOLSS (CDRV). ----
 550  MSG = 'DLSODES- At T (=R1) and step size H (=R2), a fatal'
      CALL XERRWD (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      error flag was returned by CDRV (by way of  '
      CALL XERRWD (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      Subroutine DPRJS or DSOLSS)       '
      CALL XERRWD (MSG, 40, 207, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -7
      GO TO 580
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = 0.0D0
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional outputs. -------------------------------
 580  DO 590 I = 1,N
        Y(I) = RWORK(I+LYH-1)
 590  CONTINUE
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NNZ
      IWORK(20) = NGP
      IWORK(21) = NLU
      IWORK(25) = NZL
      IWORK(26) = NZU
C karline: added next two lines
      IWORK(10) = NGE
      TLAST = T
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'DLSODES- ISTATE (=I1) illegal.'
      CALL XERRWD (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0D0, 0.0D0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'DLSODES- ITASK (=I1) illegal. '
      CALL XERRWD (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 603  MSG = 'DLSODES- ISTATE.gt.1 but DLSODES not initialized. '
      CALL XERRWD (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 604  MSG = 'DLSODES- NEQ (=I1) .lt. 1     '
      CALL XERRWD (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 605  MSG = 'DLSODES- ISTATE = 3 and NEQ increased (I1 to I2). '
      CALL XERRWD (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 606  MSG = 'DLSODES- ITOL (=I1) illegal.  '
      CALL XERRWD (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 607  MSG = 'DLSODES- IOPT (=I1) illegal.  '
      CALL XERRWD (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 608  MSG = 'DLSODES- MF (=I1) illegal.    '
      CALL XERRWD (MSG, 30, 8, 0, 1, MF, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 609  MSG = 'DLSODES- SETH (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 9, 0, 0, 0, 0, 1, SETH, 0.0D0)
      GO TO 700
 611  MSG = 'DLSODES- MAXORD (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 11, 0, 1, MAXORD, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 612  MSG = 'DLSODES- MXSTEP (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 613  MSG = 'DLSODES- MXHNIL (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 614  MSG = 'DLSODES- TOUT (=R1) behind T (=R2)      '
      CALL XERRWD (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWD (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0D0)
      GO TO 700
 615  MSG = 'DLSODES- HMAX (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0D0)
      GO TO 700
 616  MSG = 'DLSODES- HMIN (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0D0)
      GO TO 700
 617  MSG = 'DLSODES- RWORK length is insufficient to proceed. '
      CALL XERRWD (MSG, 50, 17, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
      GO TO 700
 618  MSG = 'DLSODES- IWORK length is insufficient to proceed. '
      CALL XERRWD (MSG, 50, 18, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='        Length needed is .ge. LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWD (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)
      GO TO 700
 619  MSG = 'DLSODES- RTOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0D0)
      GO TO 700
 620  MSG = 'DLSODES- ATOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0D0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSODES- EWT(I1) is R1 .le. 0.0         '
      CALL XERRWD (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0D0)
      GO TO 700
 622  MSG='DLSODES- TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWD (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='DLSODES- ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWD (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='DLSODES- ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWD (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='DLSODES- ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWD (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'DLSODES- At start of problem, too much accuracy   '
      CALL XERRWD (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWD (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0D0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'DLSODES- Trouble in DINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWD (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0D0)
      GO TO 700
 628  MSG='DLSODES- RWORK length insufficient (for Subroutine DPREP).  '
      CALL XERRWD (MSG, 60, 28, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 28, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
      GO TO 700
 629  MSG='DLSODES- RWORK length insufficient (for Subroutine JGROUP). '
      CALL XERRWD (MSG, 60, 29, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 29, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
      GO TO 700
 630  MSG='DLSODES- RWORK length insufficient (for Subroutine ODRV).   '
      CALL XERRWD (MSG, 60, 30, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 30, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
      GO TO 700
 631  MSG='DLSODES- Error from ODRV in Yale Sparse Matrix Package.     '
      CALL XERRWD (MSG, 60, 31, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      IMUL = (IYS - 1)/N
      IREM = IYS - IMUL*N
      MSG='      At T (=R1), ODRV returned error flag = I1*NEQ + I2.   '
      CALL XERRWD (MSG, 60, 31, 0, 2, IMUL, IREM, 1, TN, 0.0D0)
      GO TO 700
 632  MSG='DLSODES- RWORK length insufficient (for Subroutine CDRV).   '
      CALL XERRWD (MSG, 60, 32, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='        Length needed is .ge. LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 32, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
      GO TO 700
 633  MSG='DLSODES- Error from CDRV in Yale Sparse Matrix Package.     '
      CALL XERRWD (MSG, 60, 33, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      IMUL = (IYS - 1)/N
      IREM = IYS - IMUL*N
      MSG='      At T (=R1), CDRV returned error flag = I1*NEQ + I2.   '
      CALL XERRWD (MSG, 60, 33, 0, 2, IMUL, IREM, 1, TN, 0.0D0)
      IF (IMUL .EQ. 2) THEN
      MSG='        Duplicate entry in sparsity structure descriptors.  '
      CALL XERRWD (MSG, 60, 33, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      ENDIF
      IF (IMUL .EQ. 3 .OR. IMUL .EQ. 6) THEN
      MSG='        Insufficient storage for NSFC (called by CDRV).     '
      CALL XERRWD (MSG, 60, 33, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      ENDIF

C Karline: added next error messages
 680  MSG = 'DLSODES-  NG (=I1) .lt. 0     '
      CALL XERRWD (MSG, 30, 30, 0, 1, NG, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 681  MSG = 'DLSODES-  NG changed (from I1 to I2) illegally,   '
      CALL XERRWD (MSG, 50, 31, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      i.e. not immediately after a root was found.'
      CALL XERRWD (MSG, 50, 31, 0, 2, NGC, NG, 0, 0.0D0, 0.0D0)
      GO TO 700
 682  MSG = 'DLSODES-  One or more components of g has a root  '
      CALL XERRWD (MSG, 50, 32, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      too near to the initial point.    '
      CALL XERRWD (MSG, 40, 32, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
C
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'DLSODES- Run aborted.. apparent infinite loop.    '
      CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)
      RETURN
C----------------------- End of Subroutine DLSODES ---------------------
      END
