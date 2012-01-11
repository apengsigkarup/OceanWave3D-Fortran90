      SUBROUTINE dcosft(Y,N)
!
! This is a double-precision version of a fast cosine transform routine
! from Numerical Recipes.  Computes the transform of Y(1:N+1) where N is a
! power of two.  This is also the inverse transform but in that case
! multiply Y by 2/N.
!
      IMPLICIT NONE
      INTEGER n, j
      REAL*8 Y(N+1)
      REAL*8 WR, WI, WPR, WPI, WTEMP, THETA, SUM, y1, y2
      THETA=3.141592653589793D0/DBLE(n)
      WR=1.0D0
      WI=0.0D0
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      SUM=.5d0*(Y(1)-y(n+1))
      y(1)=.5d0*(y(1)+y(n+1))
      DO 11 J=1,n/2 -1
         WTEMP=WR
         WR=WR*WPR-WI*WPI+WR
         WI=WI*WPR+WTEMP*WPI+WI
         Y1=0.5d0*(Y(J+1)+Y(n-J+1))
         Y2=(Y(J+1)-Y(n-J+1))
         Y(J+1)=Y1-WI*Y2
         Y(n-J+1)=Y1+WI*Y2
         SUM=SUM+WR*Y2
 11   CONTINUE
      CALL drealft(Y,n,+1)
      y(n+1)=y(2)
      Y(2)=SUM
      DO 12 J=4,N,2
         SUM=SUM+Y(J)
         Y(J)=SUM
 12   CONTINUE
      RETURN
      END SUBROUTINE DCOSFT
