SUBROUTINE eta_nonlinear_3D(k,angle,x,y,omega,time,n_four_modes,yy,eta,etax,etay,etaxx,etayy,etat)
!
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: n_four_modes
REAL(KIND=long), INTENT(IN)  :: k,angle,x,y,omega,time,yy(n_four_modes)
REAL(KIND=long), INTENT(OUT) :: eta,etax,etay
REAL(KIND=long), INTENT(OUT) :: etaxx,etayy,etat
! REAL(KIND=long), OPTIONAL :: y,etay,etayy !Boundary fitted ?
! Local variables
REAL(KIND=long) :: kx, ky, k_omegat, cosk_omegat, sink_omegat, km
INTEGER    :: m
!
!
kx=cos(angle)
ky=sin(angle)
!
!kx = k*x-omega*time
k_omegat=kx*k*x+ky*k*y-omega*time
km=DBLE(n_four_modes)
cosk_omegat=half*yy(n_four_modes)*cos(km*k_omegat)
sink_omegat=half*yy(n_four_modes)*sin(km*k_omegat)
!
eta   = cosk_omegat
!
etax  = -km*kx*sink_omegat
etaxx = -(km*kx)**2*cosk_omegat
etay  = -km*ky*sink_omegat
etayy = -(km*ky)**2*cosk_omegat
!
etat  = omega*km*sink_omegat
!FIXME optimizaton of loops (GD)
DO  m=1,n_four_modes-1
   km=DBLE(m)
   cosk_omegat=yy(m)*cos(km*k_omegat)
   sink_omegat=yy(m)*sin(km*k_omegat)
   ! Non-dimensional eta and derivatives
   eta  = eta   + cosk_omegat
   !
   etax = etax  - km*kx*sink_omegat
   etaxx= etaxx - (km*kx)**2*cosk_omegat
   etay = etay  - km*ky*sink_omegat
   etayy= etayy - (km*ky)**2*cosk_omegat
   !
   etat = etat  + omega*km*sink_omegat
ENDDO

END SUBROUTINE eta_nonlinear_3D
