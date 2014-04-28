SUBROUTINE phi_nonlinear(k,x,z,omega,time,n_four_modes,zz,phix,phiz,phi,phit)
!
USE PRECISION
USE Constants
IMPLICIT NONE
INTEGER :: n_four_modes
REAL(KIND=long), INTENT(IN)  :: k,x,z,omega,time,zz(2*n_four_modes+10) ! GD: change zz(*)
REAL(KIND=long), INTENT(OUT) :: phix,phiz
REAL(KIND=long), INTENT(OUT), OPTIONAL :: phi,phit
   ! Local variables
REAL(KIND=long) :: kx, km, sinkmkx, sinhzh, temp1, e, b, e_inv, dir
INTEGER    :: m
!
dir = omega/ABS(omega) !direction of the wavefield
kx=k*x-omega*time
!
phix = zero
phiz = zero
IF(present(phit)) phit = zero
IF(present(phi))  phi  = zero
!
! FIXME: optimization of loops (remove tests inside?)
DO m=1,n_four_modes
   km=DBLE(m)
   sinkmkx = SIN(km*kx)
   e=EXP(km*(zz(1)+z))
   temp1=zz(n_four_modes+m+10)*half/COSH(km*zz(1))
   e_inv=one/e
   b=(e+e_inv)*temp1
   sinhzh=km*(e-e_inv)*temp1
   phiz=phiz+sinhzh*sinkmkx
   !
   temp1=b*km*COS(km*kx)
   phix=phix+temp1
   IF(present(phit)) phit=phit-omega*temp1
   IF(present(phi))  phi=phi+b*sinkmkx
ENDDO
!
phiz= dir*phiz
phix= dir*(phix + zz(5))

IF(present(phi))  phi = dir*(phi + kx*zz(5)) !it's k*x-omega*time      !check the time dependence...
IF(present(phit)) phit= dir*(phit)+(-zz(9)+zz(4)**2/2.d0)-dir*omega*zz(5) !check the time dependence...

END SUBROUTINE phi_nonlinear
