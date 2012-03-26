      FUNCTION BESSY (N,X)
! ------------------------------------------------------------------
!     This subroutine calculates the second kind Bessel Function of
!     integer order N, for any real X. We use here the classical
!     recursive formula.
! ------------------------------------------------------------------
      IMPLICIT NONE
      REAL *8 X,BESSY,BESSY0,BESSY1,TOX,BY,BYM,BYP
      INTEGER N,J
      IF (N.EQ.0) THEN
      BESSY = BESSY0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSY = BESSY1(X)
      RETURN
      ENDIF
      IF (X.EQ.0.) THEN
      BESSY = -1.E30
      RETURN
      ENDIF
      TOX = 2./X
      BY  = BESSY1(X)
      BYM = BESSY0(X)
      DO 11 J = 1,N-1
      BYP = J*TOX*BY-BYM
      BYM = BY
      BY  = BYP
   11 CONTINUE
      BESSY = BY
      RETURN
      END
