SUBROUTINE BottomVincent(FineGrid,GX,GY)
!
! Define bottom profile for Vincent and Briggs experiment (3D), waves are assumed to 
! be in the x-direction and the shoal center is placed at x=10.61m, y=13.72. 
!
! By Maite Gouin. Modified by H.B. Bingham.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: FineGrid
REAL(KIND=long) :: y, x, xr, yr, k, angle, dy
INTEGER :: Nx, Ny, i, j, GX, GY


Nx = FineGrid%Nx+2*GX
Ny = FineGrid%Ny+2*GY
dy = FineGrid%y(1,Ny)-FineGrid%y(1,Ny-1)
angle = 0              ! We assume waves in the x-direction.
k = angle/180.0_long*pi
DO j = 1, Ny
   DO i = 1 , Nx
      x = FineGrid%x(i,1) - 6.1_long -4.51_long ! shoal is centered at 6.1 + 2 lambda in x.
      y = FineGrid%y(1,j) - 13.72_long+GY*dy ! shoal at 13.72m in y.
      yr = cos(k)*y - sin(k)*x
      xr = sin(k)*y + cos(k)*x
      !
      FineGrid%h(i,j)   = 0.4572_long
      !
      IF ((yr/3.96_long)**2 + (xr/3.05_long)**2<1) THEN
         FineGrid%h(i,j)   = FineGrid%h(i,j) + 0.4572_long - 0.7620_long*sqrt(one-(yr/4.95_long)**2-(xr/3.81_long)**2)
      END IF
   END DO
END DO
WRITE(*,*) '  WARNING: The bottom slopes for the Vincent and Briggs experiment are not defined analytically.,'
WRITE(*,*) '  WARNING: Therefore the derivatives of the bottom profile function are determined numerically.'
WRITE(*,*) '  WARNING: Shoal center is at x=10.61m, y=13.72m. Waves should be in the x-direction'
END SUBROUTINE BottomVincent
