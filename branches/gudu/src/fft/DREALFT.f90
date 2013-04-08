      SUBROUTINE DREALFT(DATA,N,ISIGN)
!
! Double-precision version of Numerical recipes REALFT routine.
!
      implicit none
      integer n, isign
      real*8 DATA(n)
!local variables
      integer i,i1,i2,i3,i4,n2p3
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA,c1,c2,h1i,h1r,h2i,h2r,wis,wrs

      THETA=3.141592653589793D0/DBLE(N/2)
      C1=0.5d00
      IF (ISIGN.EQ.1) THEN
         C2=-0.5d00
         CALL dFOUR1(DATA,N/2,+1)
      ELSE
         C2=0.5d00
         THETA=-THETA
      ENDIF
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.0D0+WPR
      WI=WPI
      N2P3=N+3
      DO 11 I=2,N/4
         I1=2*I-1
         I2=I1+1
         I3=N2P3-I2
         I4=I3+1
         WRS= (WR)
         WIS= (WI)
         H1R=C1*(DATA(I1)+DATA(I3))
         H1I=C1*(DATA(I2)-DATA(I4))
         H2R=-C2*(DATA(I2)+DATA(I4))
         H2I=C2*(DATA(I1)-DATA(I3))
         DATA(I1)=H1R+WRS*H2R-WIS*H2I
         DATA(I2)=H1I+WRS*H2I+WIS*H2R
         DATA(I3)=H1R-WRS*H2R+WIS*H2I
         DATA(I4)=-H1I+WRS*H2I+WIS*H2R
         WTEMP=WR
         WR=WR*WPR-WI*WPI+WR
         WI=WI*WPR+WTEMP*WPI+WI
 11   CONTINUE
      IF (ISIGN.EQ.1) THEN
         H1R=DATA(1)
         DATA(1)=H1R+DATA(2)
         DATA(2)=H1R-DATA(2)
      ELSE
         H1R=DATA(1)
         DATA(1)=C1*(H1R+DATA(2))
         DATA(2)=C1*(H1R-DATA(2))
         CALL dFOUR1(DATA,N/2,-1)
      ENDIF
      RETURN
      END SUBROUTINE DREALFT
