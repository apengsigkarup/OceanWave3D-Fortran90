SUBROUTINE phi_linear_standing(H,Lx,x,z,d,omega,time,grav,phix,phiz,phi,phit)
!
USE PRECISION
USE Constants
IMPLICIT NONE
REAL(KIND=long), INTENT(IN)  :: H,Lx,x,z,d,omega,time,grav
REAL(KIND=long), INTENT(OUT) :: phix,phiz
REAL(KIND=long), INTENT(OUT), OPTIONAL :: phi,phit
   ! Local variables
REAL(KIND=long) :: kx, sinkx, temp1, temp2
!
REAL(KIND=long) :: HH,k,w,ky,c,g,y,Ly
!
y = zero ; Ly = zero ; ky = zero
kx    = two*pi/Lx
g = grav
k = SQRT(kx**2+ky**2)
w = omega
c = SQRT(g/k*TANH(k*d))
HH = H

phix    =  kx*half*HH*c*COSH(k*(z+d))/SINH(k*d)*SIN(w*time)*SIN(kx*X)*COS(ky*Y)
phiz    = -half*k*HH*c*SINH(k*(z+d))/SINH(k*d)*SIN(w*time)*COS(kx*X)*COS(ky*Y)
!
IF(present(phit)) THEN
    phit   = -half*w*HH*c*COSH(k*(z+d))/SINH(k*d)*COS(w*time)*COS(kx*X)*COS(ky*Y)
ENDIF
IF(present(phi)) THEN
    phi    = -half*HH*c*COSH(k*(z+d))/SINH(k*d)*SIN(w*time)*COS(kx*X)*COS(ky*Y)
ENDIF

END SUBROUTINE phi_linear_standing
