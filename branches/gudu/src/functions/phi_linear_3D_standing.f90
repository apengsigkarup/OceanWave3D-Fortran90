SUBROUTINE phi_linear_3D_standing(H,angle,Lx,Ly,xc,yc,z,d,omega,time,grav,phix,phiy,phiz,phi,phit)
!
USE Precision
USE Constants
IMPLICIT NONE
!
REAL(KIND=long), INTENT(IN)  :: H,Lx,Ly,xc,yc,z,d,omega,time,grav
REAL(KIND=long), INTENT(OUT) :: phix,phiy,phiz
REAL(KIND=long), INTENT(OUT), OPTIONAL :: phi,phit
! REAL(KIND=long), OPTIONAL :: y,phi,phiyy !Boundary fitted ?
! Local variables
REAL(KIND=long) :: kx, ky, k_omegat, sink_omegat, cosk_omegat, temp1, temp2
!
REAL(KIND=long) :: HH,k,w,g,c,angle,x,y,tmp
!
x = xc*COS(-angle)-yc*SIN(-angle)
y = xc*SIN(-angle)+yc*COS(-angle)

g = grav
kx    = two*pi/Lx; ky = two*pi/Ly
k = SQRT(kx**2+ky**2)
w = omega
c = SQRT(g/k*TANH(k*d))
HH = H

phix    =  kx*half*HH*c*COSH(k*(z+d))/SINH(k*d)*SIN(w*time)*SIN(kx*X)*COS(ky*Y)
phiy    =  ky*half*HH*c*COSH(k*(z+d))/SINH(k*d)*SIN(w*time)*COS(kx*X)*SIN(ky*Y)
phiz    = -half*k*HH*c*SINH(k*(z+d))/SINH(k*d)*SIN(w*time)*COS(kx*X)*COS(ky*Y)
!
IF(present(phit)) THEN
    phit   = -half*w*HH*c*COSH(k*(z+d))/SINH(k*d)*COS(w*time)*COS(kx*X)*COS(ky*Y)
ENDIF
IF(present(phi)) THEN
    phi    = -half*HH*c*COSH(k*(z+d))/SINH(k*d)*SIN(w*time)*COS(kx*X)*COS(ky*Y)
ENDIF
!
tmp = phix
phix = tmp*COS(angle)-phiy*SIN(angle)
phiy = tmp*SIN(angle)+phiy*COS(angle)
!
!xc = x*COS(angle)-y*SIN(angle)
!yc = x*SIN(angle)+y*COS(angle)

END SUBROUTINE phi_linear_3D_standing
