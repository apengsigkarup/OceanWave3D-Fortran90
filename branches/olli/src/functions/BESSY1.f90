      FUNCTION BESSY1 (X)
      IMPLICIT NONE
      REAL *8 X,BESSY1,BESSJ1,FR,FS,Z,FP,FQ,XX
! ----------------------------------------------------------------------
!     This subroutine calculates the Second Kind Bessel Function of
!     order 1, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
! ----------------------------------------------------------------------
      REAL *8 Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6,S7
      DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4, &
      .2457520174D-5,-.240337019D-6 /
      DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,  &
      .8449199096D-5,-.88228987D-6,.105787412D-6 /
      DATA R1,R2,R3,R4,R5,R6 /-.4900604943D13,.1275274390D13,   &
      -.5153438139D11,.7349264551D9,-.4237922726D7,.8511937935D4 /
      DATA S1,S2,S3,S4,S5,S6,S7 /.2499580570D14,.4244419664D12, &
      .3733650367D10,.2245904002D8,.1020426050D6,.3549632885D3,1.D0 /
      IF (X.EQ.0.) THEN
      BESSY1 = -1.E30
      RETURN
      ENDIF
      IF (X.LT.8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*(S6+Y*S7)))))
      BESSY1 = X*(FR/FS)+.636619772*(BESSJ1(X)*LOG(X)-1./X)
      ELSE
      Z = 8./X
      Y = Z*Z
      XX = X-2.356194491
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSY1 = SQRT(.636619772/X)*(SIN(XX)*FP+Z*COS(XX)*FQ)
      ENDIF
      RETURN
      END
