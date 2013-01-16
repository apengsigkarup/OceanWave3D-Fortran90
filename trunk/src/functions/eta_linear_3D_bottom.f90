SUBROUTINE eta_linear_3D_bottom(H,k,angle,x,y,int_kx,omega,time,eta,etax,etay,etaxx,etayy,etat)
!
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
IMPLICIT NONE
REAL(KIND=long), INTENT(IN)  :: H,k,angle,x,y,omega,time,int_kx
REAL(KIND=long), INTENT(OUT) :: eta,etax,etay
REAL(KIND=long), INTENT(OUT) :: etaxx,etayy,etat
! REAL(KIND=long), OPTIONAL :: y,etay,etayy !Boundary fitted ?
! Local variables
REAL(KIND=long) :: kx, ky, half_H, k_omegat, sink_omegat
!
kx=k*cos(angle)
ky=k*sin(angle)
!
k_omegat=int_kx+ky*y-omega*time
!k_omegat=kx*x+ky*y-omega*time
half_H = H*half
sink_omegat = sin(k_omegat)
!
eta   =  half_H*cos(k_omegat)
etax  = -half_H*kx*sink_omegat
etay  = -half_H*ky*sink_omegat

etat  =  omega*half_H*sink_omegat
etayy  = -eta*(ky)**2
etaxx = -eta*(kx)**2


END SUBROUTINE eta_linear_3D_bottom
