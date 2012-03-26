SUBROUTINE FOUR1(DATA,NN,ISIGN)
!
! Replaces data(1:2nn) by its discrete Fourier transform, if isign is 1; or
! with nn times the inverse discrete tranform if isign is -1.  data is a complex
! array of length nn or a real array of length 2*nn.  nn must be an integer
! power of 2.
!
  REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
  DIMENSION DATA(*)
  N=2*NN
  J=1
  DO 11 I=1,N,2
     IF(J.GT.I)THEN
        TEMPR=DATA(J)
        TEMPI=DATA(J+1)
        DATA(J)=DATA(I)
        DATA(J+1)=DATA(I+1)
        DATA(I)=TEMPR
        DATA(I+1)=TEMPI
     ENDIF
     M=N/2
1    IF ((M.GE.2).AND.(J.GT.M)) THEN
        J=J-M
        M=M/2
        GO TO 1
     ENDIF
     J=J+M
11 CONTINUE
  MMAX=2
2 IF (N.GT.MMAX) THEN
     ISTEP=2*MMAX
     THETA=6.28318530717959D0/(ISIGN*MMAX)
     WPR=-2.D0*DSIN(0.5D0*THETA)**2
     WPI=DSIN(THETA)
     WR=1.D0
     WI=0.D0
     DO 13 M=1,MMAX,2
        DO 12 I=M,N,ISTEP
           J=I+MMAX
           TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
           TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
           DATA(J)=DATA(I)-TEMPR
           DATA(J+1)=DATA(I+1)-TEMPI
           DATA(I)=DATA(I)+TEMPR
           DATA(I+1)=DATA(I+1)+TEMPI
12      CONTINUE
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
13     CONTINUE
       MMAX=ISTEP
       GO TO 2
    ENDIF
    RETURN
  END SUBROUTINE FOUR1
