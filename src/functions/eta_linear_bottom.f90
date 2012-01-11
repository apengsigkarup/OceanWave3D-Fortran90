SUBROUTINE eta_linear_bottom(H,k,x,int_kx,omega,time,eta,etax,etaxx,etat)
!
USE Precision
USE Constants
IMPLICIT NONE
REAL(KIND=long), INTENT(IN)  :: H,k,x,omega,time,int_kx
REAL(KIND=long), INTENT(OUT) :: eta,etax
REAL(KIND=long), INTENT(OUT) :: etaxx,etat
! REAL(KIND=long), OPTIONAL :: y,etay,etayy !Boundary fitted ?
! Local variables
REAL(KIND=long) :: kx, half_H, sinkx
!
kx=int_kx-omega*time
half_H = H*half
sinkx = sin(kx)
!
eta   =  half_H*cos(kx)
etax  = -half_H*k*sinkx
etat  =  omega*half_H*sinkx
etaxx = -eta*(k)**2

END SUBROUTINE eta_linear_bottom
