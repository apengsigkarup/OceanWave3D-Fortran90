SUBROUTINE phi_nonlinear_3D(k,angle,x,y,z,omega,time,n_four_modes,zz,phix,phiy,phiz,phi,phit)
!
USE PRECISION
USE Constants
IMPLICIT NONE
INTEGER :: n_four_modes
REAL(KIND=long), INTENT(IN)  :: k,angle,x,y,z,omega,time,zz(2*n_four_modes+10)
REAL(KIND=long), INTENT(OUT) :: phix,phiy,phiz
REAL(KIND=long), INTENT(OUT), OPTIONAL :: phi,phit
   ! Local variables
REAL(KIND=long) :: k_omegat, km, sinkmk_omegat, sinhzh, temp1, e, b, e_inv, dir
REAL(KIND=long) :: kx, ky
INTEGER    :: m
!
dir = omega/ABS(omega) !direction of the wavefield
!
kx=cos(angle)
ky=sin(angle)
!
!kx = k*x-omega*time
k_omegat=kx*k*x+ky*k*y-omega*time
km=DBLE(n_four_modes)
!
phix = zero
phiy = zero
phiz = zero
IF(present(phit)) phit = zero
IF(present(phi))  phi  = zero
!
! FIXME: optimization of loops (remove tests inside?)
DO m=1,n_four_modes
   km=DBLE(m)
   sinkmk_omegat = SIN(km*k_omegat)
   !sinkmkx = SIN(km*kx)
   e=EXP(km*(zz(1)+z))
   temp1=zz(n_four_modes+m+10)*half/COSH(km*zz(1))
   e_inv=one/e
   b=(e+e_inv)*temp1
   sinhzh=km*(e-e_inv)*temp1
   phiz=phiz+sinhzh*sinkmk_omegat
   !
   temp1=b*km*COS(km*k_omegat)
   phix=phix+temp1*kx
   phiy=phiy+temp1*ky
   IF(present(phit)) phit=phit-omega*temp1
   IF(present(phi))  phi=phi+b*sinkmk_omegat
ENDDO
!
phiz= dir*phiz
phix= dir*(phix + kx*zz(5))
phiy= dir*(phiy + ky*zz(5))

IF(present(phi))  phi = dir*(phi + k_omegat*zz(5)) !it's k*x-omega*time      !check the time dependence...
IF(present(phit)) phit= dir*(phit)+(-zz(9)+zz(4)**2/2.d0)-dir*omega*zz(5) !check the time dependence...

END SUBROUTINE phi_nonlinear_3D
