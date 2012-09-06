      FUNCTION BESSJP (N,X)
      IMPLICIT NONE
      REAL *8 X,BESSJP,BESSJ
      INTEGER N
      IF (N.EQ.0) THEN
        BESSJP=-BESSJ(1,X)
      ELSE IF(X.EQ.0.D0) THEN
        X=1.D-30
      ELSE
        BESSJP=BESSJ(N-1,X)-(DFLOAT(N)/X)*BESSJ(N,X)
      ENDIF
      RETURN
      END
