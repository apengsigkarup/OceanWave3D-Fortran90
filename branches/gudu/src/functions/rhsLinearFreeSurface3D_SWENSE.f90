SUBROUTINE rhsLinearFreeSurface3D_SWENSE(t,Wavefield,g,rhsE,rhsP,Nx,Ny)
!
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx,Ny
REAL(KIND=long) :: g,t
REAL(KIND=long), DIMENSION(Nx,Ny) :: rhsE, rhsP ! GD remove: E, Ex, Ey, P, Px, Py, W,
!REAL(KIND=long), DIMENSION(Nx,Ny) :: E_I,Ex_I,Ey_I,Et_I,Px_I,Py_I,W_I,Pt_I
TYPE(Wavefield_FS) :: Wavefield

rhsE = Wavefield%W + Wavefield%Pz_I_s - Wavefield%Et_I     ! kinematic FS BC
rhsP = -g*(Wavefield%E+Wavefield%E_I) - Wavefield%Pt_I_s   ! dynamic FS BC!

END SUBROUTINE rhsLinearFreeSurface3D_SWENSE
