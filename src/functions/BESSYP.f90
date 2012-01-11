      FUNCTION BESSYP (N,X)
      IMPLICIT NONE
      REAL *8 X,BESSYP,BESSY
      INTEGER N
      IF (N.EQ.0) THEN
        BESSYP=-BESSY(1,X)
      ELSE IF(X.EQ.0.D0) THEN
        X=1.D-30
      ELSE
        BESSYP=BESSY(N-1,X)-(DFLOAT(N)/X)*BESSY(N,X)
      ENDIF
      RETURN
      END
