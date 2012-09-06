SUBROUTINE eta_nonlinear(k,x,omega,time,n_four_modes,yy,eta,etax,etaxx,etat)
!
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: n_four_modes
REAL(KIND=long), INTENT(IN)  :: k,x,omega,time, yy(n_four_modes) ! GD: change yy(*)
REAL(KIND=long), INTENT(OUT) :: eta,etax
REAL(KIND=long), INTENT(OUT) :: etaxx,etat
! REAL(long), OPTIONAL :: y,etay,etayy !Boundary fitted ?
! Local variables
REAL(KIND=long) :: kx, sinkx, coskx, km
INTEGER    :: m
!
kx = k*x-omega*time
km=DBLE(n_four_modes)
coskx=half*yy(n_four_modes)*cos(km*kx)
sinkx=half*yy(n_four_modes)*sin(km*kx)
!
eta   = coskx
etax  = -km*sinkx
etaxx = -(km)**2*coskx
etat  = omega*km*sinkx
!FIXME optimizaton of loops (GD)
DO  m=1,n_four_modes-1
   km=DBLE(m)
   coskx=yy(m)*cos(km*kx)
   sinkx=yy(m)*sin(km*kx)
   ! Non-dimensional eta and derivatives
   eta  = eta   + coskx
   etax = etax  - km*sinkx
   etaxx= etaxx - km**2*coskx
   etat = etat  + omega*km*sinkx
ENDDO

END SUBROUTINE eta_nonlinear
