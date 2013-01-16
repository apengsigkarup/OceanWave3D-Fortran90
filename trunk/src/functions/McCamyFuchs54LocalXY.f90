SUBROUTINE McCamyFuchs54LocalXY(x,y,z,a,Ainc,k,d,g,time,Etol,MaxIter,eta,phi,w)
! By Allan P. Engsig-Karup.
!
USE Precision
USE DataTypes
USE Constants
!
IMPLICIT NONE
INTEGER :: m,MaxIter
REAL(KIND=long) :: x,y,z,a,Ainc,k,d,g,time,Etol
REAL(KIND=long) :: r,theta,eta,phi,w,J0P,JmP,bessj
COMPLEX(KIND=long) :: nextterm,oldterm,eta_c,phi_c,w_c
COMPLEX(KIND=long) :: HmP,H0P,bessh
!
! r   : radius
! phi : angle with abcissa
r = SQRT(x**2+y**2)
!
IF (y.GT.0) THEN
    theta = ACOS(x/r)/pi*180.0_long
ELSE
    theta = (-ACOS(x/r)+two*pi)/pi*180.0_long
ENDIF
!
J0P = -bessj(1,k*a)
H0P = -bessh(1,k*a)
eta_c = bessj(0,k*r) - bessh(0,k*r)* J0P/H0P
nextterm = Etol*10
oldterm = Etol*10
m=1
IF (-(r-a).GT.1e-3) THEN
    eta_c = 0
    !print*, r, a
    !print*, 'Point taken inside cylinder (r<a).'
ELSEIF ((r-a).GT.1e-3) THEN
	!print*,r,a
    !Do nothing
    DO WHILE ((m.LE.MaxIter).AND.(((abs(nextterm).GT.Etol).OR.(abs(oldterm).GT.Etol))))
      oldterm = nextterm
      !
      JmP = half * ( bessj(m-1,k*a) - bessj(m+1,k*a) )
      HmP = half * ( bessh(m-1,k*a) - bessh(m+1,k*a) )
      nextterm = two*(0.0_long,1.0_long)**m*( bessj(m,k*r) - bessh(m,k*r)* JmP/HmP)*cos(m*theta/360.0_long*two*pi)
      !IF (isnan(nextterm)) THEN
      !    print*,'NaN errors. (Increase N or lower Etol.) '
      !    stop
      !ENDIF
      eta_c = eta_c + nextterm
      !IF ((abs(nextterm).LT.Etol).AND.(abs(oldterm).LT.Etol)) THEN
          IF (m.EQ.MaxIter) THEN
              print*,'Max. number of terms reached. (Increase N or lower Etol.), ', MaxIter, Etol
              print*,'test',abs(nextterm),abs(oldterm)
              print*,'x=',x
              print*,'y=',y
read*
!              pause
          ENDIF
          !stop
      !ENDIF
      m=m+1
      !print*,m-1,abs(nextterm),abs(oldterm)
      !print*,'x=',x,'y=',y
  ENDDO
ELSE
  ! assume point is on cylinder.
  r = a
  !print*,r,a
  !DO m = 1,MaxIter
  DO WHILE ((m.LE.MaxIter).AND.(((abs(nextterm).GT.Etol).OR.(abs(oldterm).GT.Etol))))
      oldterm = nextterm
      !
      JmP = half * ( bessj(m-1,k*a) - bessj(m+1,k*a) )
      HmP = half * ( bessh(m-1,k*a) - bessh(m+1,k*a) )
      nextterm = two*(0.0_long,1.0_long)**m*( bessj(m,k*r) - bessh(m,k*r)* JmP/HmP)*cos(m*theta/360.0_long*two*pi)
      !IF (isnan(nextterm)) THEN
      !    print*,'NaN errors. (Increase N or lower Etol.) '
      !    stop
      !ENDIF
      eta_c = eta_c + nextterm
      !IF ((abs(nextterm).LT.Etol).AND.(abs(oldterm).LT.Etol)) THEN
          IF (m.EQ.MaxIter) THEN
              print*,'Max. number of terms reached. (Increase N or lower Etol.), ', MaxIter, Etol, abs(nextterm)
              print*,'x=',x
              print*,'y=',y
read*
!              pause
          ENDIF
          !stop
      !ENDIF
      m=m+1
      !print*,m-1,abs(nextterm),abs(oldterm)
      !pause
  ENDDO
ENDIF
!
eta_c = Ainc*EXP(-(0.0_long,1.0_long)*sqrt(g*k*TANH(k*d))*time)*eta_c
phi_c = -((0.0_long,1.0_long))*g/(sqrt(g*k*TANH(k*d)))*eta_c*COSH(k*(z+d))/COSH(k*d)
w_c   = -((0.0_long,1.0_long))*(sqrt(g*k*TANH(k*d)))*eta_c
! Check this...
eta = REAL(eta_c)
phi = REAL(phi_c)
w   = REAL(w_c)

END SUBROUTINE McCamyFuchs54LocalXY
