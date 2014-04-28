SUBROUTINE LUBKSB(A,N,NP,INDX,B)
  USE Precision
  IMPLICIT NONE
  INTEGER :: i,j, ii, ll, n, np, indx(n)
  REAL(kind=long) :: A(NP,NP), B(N), SUM
  II=0
  DO  I=1,N
     LL=INDX(I)
     SUM=B(LL)
     B(LL)=B(I)
     IF (II.NE.0)THEN
        DO  J=II,I-1
           SUM=SUM-A(I,J)*B(J)
        END DO
     ELSE IF (SUM.NE.0.) THEN
        II=I
     ENDIF
     B(I)=SUM
  END DO
  DO  I=N,1,-1
     SUM=B(I)
     IF(I.LT.N)THEN
        DO  J=I+1,N
           SUM=SUM-A(I,J)*B(J)
        END DO
     ENDIF
     B(I)=SUM/A(I,I)
  END DO
  RETURN
END SUBROUTINE LUBKSB
