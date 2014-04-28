      FUNCTION BESSY0 (X)
      IMPLICIT NONE
      REAL *8 X,BESSY0,BESSJ0,FS,FR,Z,FP,FQ,XX
! ---------------------------------------------------------------------
!     This subroutine calculates the Second Kind Bessel Function of
!     order 0, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
! ---------------------------------------------------------------------
      REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
      -.2073370639D-5,.2093887211D-6 /
      DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3,  &
      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
      DATA R1,R2,R3,R4,R5,R6 /-2957821389.D0,7062834065.D0, &
      -512359803.6D0,10879881.29D0,-86327.92757D0,228.4622733D0 /
      DATA S1,S2,S3,S4,S5,S6 /40076544269.D0,745249964.8D0, &
      7189466.438D0,47447.26470D0,226.1030244D0,1.D0 /
      IF (X.EQ.0.D0) THEN
      BESSY0 = -1.E30
      RETURN
      ENDIF
      IF (X.LT.8.D0) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSY0 = FR/FS+.636619772D0*BESSJ0(X)*LOG(X)
      ELSE
      Z = 8.D0/X
      Y = Z*Z
      XX = X-.785398164D0
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSY0 = SQRT(.636619772D0/X)*(FP*SIN(XX)+Z*FQ*COS(XX))
      ENDIF
      RETURN
      END
