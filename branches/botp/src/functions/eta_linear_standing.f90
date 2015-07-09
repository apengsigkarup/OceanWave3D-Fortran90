SUBROUTINE eta_linear_standing(H,Lx,x,omega,time,eta,etax,etaxx,etat)
!
USE Precision
USE Constants
IMPLICIT NONE
REAL(KIND=long), INTENT(IN)  :: H,Lx,x,omega,time
REAL(KIND=long), INTENT(OUT) :: eta,etax
REAL(KIND=long), INTENT(OUT) :: etaxx,etat
! REAL(KIND=long), OPTIONAL :: y,etay,etayy !Boundary fitted ?
! Local variables
REAL(KIND=long) :: kx, half_H, sinkx
!
REAL(KIND=long) :: HH,k,w,ky,y,Ly
!
y = zero ; Ly = zero ; ky = zero
kx    = two*pi/Lx
k = SQRT(kx**2+ky**2)
w = omega
HH = H
!
eta   =  half*HH*COS(w*time)*COS(kx*X)*COS(ky*Y)
etax  = -kx*half*HH*COS(w*time)*SIN(kx*X)*COS(ky*Y)
etaxx = -kx**2*half*HH*COS(w*time)*COS(kx*X)*COS(ky*Y)
etat  = -half*w*HH*SIN(w*time)*COS(kx*X)*COS(ky*Y)


END SUBROUTINE eta_linear_standing
