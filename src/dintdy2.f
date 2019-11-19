C  The code in this file is based of function DINTDY from file
C    opdka1.f from https://www.netlib.org/odepack/
C 
C  Original Author: Hindmarsh, Alan C., (LLNL)  
C  Adapted for use in R package deSolve by the deSolve authors.

      SUBROUTINE INTERPOLY(T, K, I, YH, NYH, DKY, nq, tn, h)

C***PURPOSE  Interpolate solution derivatives to be used in C-code.
C  computes interpolated values of the K-th derivative of the i-th
C  dependent variable vector y, and stores it in DKY.  This routine
C  is called within the package with K = 0 and T = TOUT, but may
C  also be called by the user for any K up to the current order.
C  (See detailed instructions in the usage documentation.)
C
C  The computed values in DKY are gotten by interpolation using the
C  Nordsieck history array YH.  This array corresponds uniquely to a
C  vector-valued polynomial of degree NQCUR or less, and DKY is set
C  to the K-th derivative of this polynomial at T.
C  The formula for DKY is:
C               q
C   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
C              j=K
C  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
C The above sum is done in reverse order.
C  IFLAG is returned negative if either K or T is out of bounds.
C
C***BASED ON  DINTDY
      IMPLICIT NONE
      INTEGER K, NYH, NQ, I, IC, J, JB, JB2, JJ, JJ1, JP1
      DOUBLE PRECISION T, DKY, H, C, R, S, Tn
      DOUBLE PRECISION YH(NYH,*) 
C
C***FIRST EXECUTABLE STATEMENT  
      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = nq+1 - K
      DO 10 JJ = JJ1,NQ
        IC = IC*JJ
 10   CONTINUE
 15   C = IC
      DKY = C*YH(I,nq+1)
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1,JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1,J
          IC = IC*JJ
 30     CONTINUE
 35     C = IC
        DKY = C*YH(I,JP1) + S*DKY
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      DKY = R*DKY
      RETURN
C----------------------- END OF SUBROUTINE InterpolY ----------------------
      END
