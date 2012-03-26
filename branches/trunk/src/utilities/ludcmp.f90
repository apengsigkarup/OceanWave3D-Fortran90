SUBROUTINE LUDCMP(A,N,NP,INDX,D)
  USE Precision
  IMPLICIT NONE
  REAL(kind=long), PARAMETER :: TINY=1.0E-20
  INTEGER i, j, k, imax, n, np, indx(n)
  REAL(kind=long) ::  A(NP,NP), VV(N), AAMAX, SUM, D, DUM, one=1._long, &
       zero=0._long
  D=one
  DO  I=1,N
     AAMAX=zero
     DO  J=1,N
        IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
     END DO
     IF (AAMAX.EQ.0.) then
	 write(*,'(A)') 'Singular matrix.'
	read*
! pause
     ENDIF
     VV(I)=1./AAMAX
  END DO
  DO  J=1,N
     DO  I=1,J-1
        SUM=A(I,J)
        DO  K=1,I-1
           SUM=SUM-A(I,K)*A(K,J)
        END DO
        A(I,J)=SUM
     END DO
     AAMAX=zero
     DO  I=J,N
        SUM=A(I,J)
        DO  K=1,J-1
           SUM=SUM-A(I,K)*A(K,J)
        END DO
        A(I,J)=SUM
        DUM=VV(I)*ABS(SUM)
        IF (DUM.GE.AAMAX) THEN
           IMAX=I
           AAMAX=DUM
        ENDIF
     END DO
     IF (J.NE.IMAX)THEN
        DO  K=1,N
           DUM=A(IMAX,K)
           A(IMAX,K)=A(J,K)
           A(J,K)=DUM
        END DO
        D=-D
        VV(IMAX)=VV(J)
     ENDIF
     INDX(J)=IMAX
     IF(A(J,J).EQ.0.)A(J,J)=TINY
     IF(J.NE.N)THEN
        DUM=1./A(J,J)
        DO  I=J+1,N
           A(I,J)=A(I,J)*DUM
        END DO
     ENDIF
  END DO
  RETURN
END SUBROUTINE LUDCMP
