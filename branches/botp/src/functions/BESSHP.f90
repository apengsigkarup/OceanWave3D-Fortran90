      FUNCTION BESSHP (N,X)
      IMPLICIT NONE
      REAL *8 X
      COMPLEX *16 BESSHP,BESSH
      INTEGER N
      IF (N.EQ.0) THEN
        BESSHP=-BESSH(1,X)
      ELSE IF(X.EQ.0.D0) THEN
        X=1.D-30
      ELSE
        BESSHP=BESSH(N-1,X)-(DFLOAT(N)/X)*BESSH(N,X)
      ENDIF
      RETURN
      END
