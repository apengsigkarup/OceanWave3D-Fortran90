SUBROUTINE BottomVincent(FineGrid,GX,GY)
!
! Define bottom profile for Vincent and Briggs experiment (3D)
!
! By Maite Gouin.
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
angle = 0
k = angle/180.0_long*pi
DO j = 1, Ny
!hbb x and y are 2D arrays in the FineGrid. 
   DO i = 1 , Nx
      x = FineGrid%x(i,1) - 6.1_long -4.51_long ! FIXME: translation relative to standard y-coordinates
      y = FineGrid%y(1,j) - 14.54_long !half*FineGrid%x(Nx,1) ! position shoal in middle relative to y-axis
      yr = cos(k)*y - sin(k)*x
      xr = sin(k)*y + cos(k)*x
      dyrdy = cos(k)
      dxrdy = sin(k)
      dyrdx = -sin(k)
      dxrdx = cos(k)
      !
      FineGrid%h(i,j)   = 0.4572_long
      !
      IF ((yr/3.96_long)**2 + (xr/3.05_long)**2<1) THEN
         FineGrid%h(i,j)   = FineGrid%h(i,j) + 0.4572_long - 0.7620_long*sqrt(one-(yr/4.95_long)**2-(xr/3.81_long)**2)
!WRITE(*,*) '  WARNING: The bottom slopes for the Vincent and Briggs experiment is not defined analytically (i.e. not implemented).'
!WRITE(*,*) '  WARNING: Therefore the derivatives of the bottom profile function is determined numerically.'
         !FineGrid%hx(i,j)  = FineGrid%hx(i,j)  + (0.02_long*xr*dxrdx + &
            ! 0.0355556_long*yr*dyrdx)/sqrt(one-0.04_long*xr**2-0.0711111_long*yr**2) 
         !FineGrid%hy(i,j)  = FineGrid%hy(i,j)  + (0.02_long*xr*dxrdy + & 
            ! 0.0355556_long*yr*dyrdy)/sqrt(one-0.04_long*xr**2-0.0711111_long*yr**2) 
!         FineGrid%hxx(i,j) = FineGrid%hxx(i,j) + 
!         FineGrid%hyy(i,j) = FineGrid%hyy(i,j) + 
      END IF
   END DO
END DO
WRITE(*,*) '  WARNING: The bottom slopes for the Vincent and Briggs experiment is not defined analytically (i.e. not implemented).'
WRITE(*,*) '  WARNING: Therefore the derivatives of the bottom profile function is determined numerically.'
END SUBROUTINE BottomVincent
