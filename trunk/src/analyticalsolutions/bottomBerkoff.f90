SUBROUTINE BottomBerkoff(FineGrid,GX,GY)
!
! Define bottom profile for Berkhoff experiment (3D)
!
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: FineGrid
REAL(KIND=long) :: y, x, xr, yr, k, angle, dxrdx, dxrdy, dyrdx, dyrdy
INTEGER :: Nx, Ny, i, j, GX, GY

!hbb Added the ghost points here.
Nx = FineGrid%Nx+2*GX
Ny = FineGrid%Ny+2*GY
angle = 20
k = angle/180.0_long*pi
DO j = 1, Ny
!hbb x and y are 2D arrays in the FineGrid. 
   DO i = 1 , Nx
      x = FineGrid%x(i,1) - half*FineGrid%x(Nx,1) ! position shoal in middle relative to y-axis
      y = FineGrid%y(1,j) - 11.0_long ! FIXME: translation relative to standard y-coordinates
      xr = cos(k)*x - sin(k)*y
      yr = sin(k)*x + cos(k)*y
      dxrdx = cos(k)
      dyrdx = sin(k)
      dxrdy = -sin(k)
      dyrdy = cos(k)
      IF (yr<=-5.82_long) THEN
         FineGrid%h(i,j)   = 0.45_long
         FineGrid%hx(i,j)  = zero
         FineGrid%hy(i,j)  = zero
         FineGrid%hxx(i,j) = zero
         FineGrid%hyy(i,j) = zero
      ELSE ! IF (yr>-5.82_long) THEN
         FineGrid%h(i,j)   = 0.45_long-0.02_long*(5.82_long+yr)
         FineGrid%hx(i,j)  = -0.02_long*dyrdx
         FineGrid%hy(i,j)  = -0.02_long*dyrdy
         FineGrid%hxx(i,j) = zero
         FineGrid%hyy(i,j) = zero
      END IF
      IF ((xr/four)**2 + (yr/three)**2<1) THEN
         FineGrid%h(i,j)   = FineGrid%h(i,j) + 0.3_long - half*sqrt(one-(xr/five)**2-(yr/3.75_long)**2)
!WRITE(*,*) '  WARNING: The bottom slopes for the Berkhoff experiment is not defined analytically (i.e. not implemented).'
!WRITE(*,*) '  WARNING: Therefore the derivatives of the bottom profile function is determined numerically.'
         FineGrid%hx(i,j)  = FineGrid%hx(i,j)  + (0.02_long*xr*dxrdx + &
             0.0355556_long*yr*dyrdx)/sqrt(one-0.04_long*xr**2-0.0711111_long*yr**2) 
         FineGrid%hy(i,j)  = FineGrid%hy(i,j)  + (0.02_long*xr*dxrdy + & 
             0.0355556_long*yr*dyrdy)/sqrt(one-0.04_long*xr**2-0.0711111_long*yr**2) 
!         FineGrid%hxx(i,j) = FineGrid%hxx(i,j) + 
!         FineGrid%hyy(i,j) = FineGrid%hyy(i,j) + 
      END IF
      IF (FineGrid%h(i,j)<0.1_long) THEN
         ! flat bottom profile after shoal 
         ! - this is very absorption zone should effectively damp outgoing waves
         FineGrid%h(i,j)   = 0.1_long
         FineGrid%hx(i,j)  = zero
         FineGrid%hy(i,j)  = zero
         FineGrid%hxx(i,j) = zero
         FineGrid%hyy(i,j) = zero
      END IF
   END DO
END DO
WRITE(*,*) '  WARNING: The bottom slopes for the Berkhoff experiment is not defined analytically (i.e. not implemented).'
WRITE(*,*) '  WARNING: Therefore the derivatives of the bottom profile function is determined numerically.'
END SUBROUTINE BottomBerkoff
