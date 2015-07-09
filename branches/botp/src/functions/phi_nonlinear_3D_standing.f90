SUBROUTINE phi_nonlinear_3D_standing(H,angle,Lx,Ly,xc,yc,z,d,eta,omega,time,grav,phix,phiy,phiz,phi,phit)
!
USE Precision
USE Constants
IMPLICIT NONE
!
REAL(KIND=long), INTENT(IN)  :: H,Lx,Ly,xc,yc,z,d,omega,time,grav,eta
REAL(KIND=long), INTENT(OUT) :: phix,phiy,phiz
REAL(KIND=long), INTENT(OUT), OPTIONAL :: phi,phit
! REAL(KIND=long), OPTIONAL :: y,phi,phiyy !Boundary fitted ?
! Local variables
REAL(KIND=long) :: kx, ky, k_omegat, sink_omegat, cosk_omegat, temp1, temp2
!
REAL(KIND=long) :: HH,k,w,g,c,angle,x,y,tmp,c2
!
x = xc*COS(-angle)-yc*SIN(-angle)
y = xc*SIN(-angle)+yc*COS(-angle)

g = grav
kx    = two*pi/Lx/two; ky = zero*two*pi/Ly
k = SQRT(kx**2+ky**2)
w = omega
!w = SQRT(g*k*TANH(k*d))
c = SQRT(g/k*TANH(k*d))
c2 = SQRT(g/(k*TANH(k*d)))
!c2 = c
HH = H

! we use just at time=0 ...
IF (time/=0.d0) THEN
  print*, 'check for time /= 0'
  stop
ENDIF

IF (time==0.d0) THEN
!IF (0.d0==0.d0) THEN
phix = half*c2*COSH(k*(z+d))/COSH(k*(eta+d))*(0.8867_long*10.0**(-1)*pi*SIN(pi*X)+&
            SQRT(SQRT(two))*0.5243*10.0**(-2)*(two*pi)*SIN(two*pi*X)+&
            SQRT(SQRT(three))*0.4978*10.0**(-3)*(three*pi)*SIN(three*pi*X)+&
            SQRT(SQRT(four))*0.6542*10.0**(-4)*(four*pi)*SIN(four*pi*X)+&
            SQRT(SQRT(five))*0.1007*10.0**(-4)*(five*pi)*SIN(five*pi*X)+&
            SQRT(SQRT(six))*0.1653*10.0**(-5)*(six*pi)*SIN(six*pi*X)+&
            SQRT(SQRT(seven))*0.2753*10.0**(-6)*(seven*pi)*SIN(seven*pi*X)+&
            SQRT(SQRT(eight))*0.4522*10.0**(-7)*(eight*pi)*SIN(eight*pi*X))*SIN(w*time)
phiy  = zero
phiz = -half*c2*k*SINH(k*(z+d))/SINH(k*(eta+d))*(0.8867_long*10.0**(-1)*COS(pi*X)+&
            SQRT(SQRT(two))*0.5243*10.0**(-2)*COS(two*pi*X)+&
            SQRT(SQRT(three))*0.4978*10.0**(-3)*COS(three*pi*X)+&
            SQRT(SQRT(four))*0.6542*10.0**(-4)*COS(four*pi*X)+&
            SQRT(SQRT(five))*0.1007*10.0**(-4)*COS(five*pi*X)+&
            SQRT(SQRT(six))*0.1653*10.0**(-5)*COS(six*pi*X)+&
            SQRT(SQRT(seven))*0.2753*10.0**(-6)*COS(seven*pi*X)+&
            SQRT(SQRT(eight))*0.4522*10.0**(-7)*COS(eight*pi*X))*SIN(w*time)

!
IF(present(phit)) THEN
    phit = -half*c2*w*COSH(k*(z+d))/COSH(k*(eta+d))*(0.8867_long*10.0**(-1)*COS(pi*X)+&
            SQRT(SQRT(two))*0.5243*10.0**(-2)*COS(two*pi*X)+&
            SQRT(SQRT(three))*0.4978*10.0**(-3)*COS(three*pi*X)+&
            SQRT(SQRT(four))*0.6542*10.0**(-4)*COS(four*pi*X)+&
            SQRT(SQRT(five))*0.1007*10.0**(-4)*COS(five*pi*X)+&
            SQRT(SQRT(six))*0.1653*10.0**(-5)*COS(six*pi*X)+&
            SQRT(SQRT(seven))*0.2753*10.0**(-6)*COS(seven*pi*X)+&
            SQRT(SQRT(eight))*0.4522*10.0**(-7)*COS(eight*pi*X))*COS(w*time)
ENDIF
IF(present(phi)) THEN
    phi = -half*c2*COSH(k*(z+d))/COSH(k*(eta+d))*(0.8867_long*10.0**(-1)*COS(pi*X)+&
            SQRT(SQRT(two))*0.5243*10.0**(-2)*COS(two*pi*X)+&
            SQRT(SQRT(three))*0.4978*10.0**(-3)*COS(three*pi*X)+&
            SQRT(SQRT(four))*0.6542*10.0**(-4)*COS(four*pi*X)+&
            SQRT(SQRT(five))*0.1007*10.0**(-4)*COS(five*pi*X)+&
            SQRT(SQRT(six))*0.1653*10.0**(-5)*COS(six*pi*X)+&
            SQRT(SQRT(seven))*0.2753*10.0**(-6)*COS(seven*pi*X)+&
            SQRT(SQRT(eight))*0.4522*10.0**(-7)*COS(eight*pi*X))*SIN(w*time)
ENDIF
ELSE
  phix    =  pi*half*0.08934064760240*c*COSH(pi*(z+d))/COSH(pi*d)*SIN(w*time)*SIN(pi*X)
  phiy    =  zero
  phiz    = -half*pi*0.08934064760240*c*SINH(pi*(z+d))/COSH(pi*d)*SIN(w*time)*COS(pi*X)
  !
  IF(present(phit)) THEN
      phit   = -half*w*0.08934064760240*c*COSH(pi*(z+d))/COSH(pi*d)*COS(w*time)*COS(pi*X)
  ENDIF
  IF(present(phi)) THEN
      phi    = -half*0.08934064760240*c*COSH(pi*(z+d))/COSH(pi*d)*SIN(w*time)*COS(pi*X)
  ENDIF
ENDIF
!
tmp = phix
phix = tmp*COS(angle)-phiy*SIN(angle)
phiy = tmp*SIN(angle)+phiy*COS(angle)
!
!xc = x*COS(angle)-y*SIN(angle)
!yc = x*SIN(angle)+y*COS(angle)

END SUBROUTINE phi_nonlinear_3D_standing
