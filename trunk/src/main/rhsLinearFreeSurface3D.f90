SUBROUTINE rhsLinearFreeSurface3D(t,Wavefield,g,rhsE,rhsP,Nx,Ny)
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx,Ny
REAL(KIND=long) :: g,t
REAL(KIND=long), DIMENSION(Nx,Ny) :: rhsE, rhsP
TYPE(Wavefield_FS) :: Wavefield
rhsE = Wavefield%W     ! kinematic FS BC
rhsP = -g*Wavefield%E  + Wavefield%P0 ! dynamic FS BC!
END SUBROUTINE rhsLinearFreeSurface3D
