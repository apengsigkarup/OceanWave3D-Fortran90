      SUBROUTINE DSINFT(Y,N)
! Double-precision version of Numerical Recipes fast sine transform.
      implicit none
      integer n, j
      real*8 y(n)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,y1,y2,sum

      THETA=3.141592653589793D0/DBLE(N)
      WR=1.0D0
      WI=0.0D0
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      Y(1)=0.0d0
      DO 11 J=1,N/2
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
        Y1=WI*(Y(J+1)+Y(N-J+1))
        Y2=0.5d0*(Y(J+1)-Y(N-J+1))
        Y(J+1)=Y1+Y2
        Y(N-J+1)=Y1-Y2
11    CONTINUE
      CALL DREALFT(Y,N,+1)
      SUM=0.0d0
      Y(1)=0.5d0*Y(1)
      Y(2)=0.0d0
      DO 12 J=1,N-1,2
        SUM=SUM+Y(J)
        Y(J)=Y(J+1)
        Y(J+1)=SUM
12    CONTINUE
      RETURN
      END
