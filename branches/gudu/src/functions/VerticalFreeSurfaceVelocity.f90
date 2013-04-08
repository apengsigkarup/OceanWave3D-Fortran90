SUBROUTINE VerticalFreeSurfaceVelocity(W,Nx,Ny,Nz,sol,DiffStencils,dsdz,gamma)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Diff_def), INTENT(IN) :: DiffStencils
INTEGER :: Nx, Ny, Nz, gamma, i, j
REAL(KIND=long), DIMENSION(Nz,Nx,Ny), INTENT(IN) :: dsdz, sol
REAL(KIND=long), DIMENSION(Nx,Ny), INTENT(OUT) :: W
INTEGER, DIMENSION(2*gamma+1) :: idx
DO i = 1, 2*gamma+1
  idx(i) = Nz-2*gamma + i - 1
END DO
DO j = 1, Ny
  DO i = 1, Nx
    ! Determine for free surface points only
	W(i,j) = dsdz(Nz,i,j)*DOT_PRODUCT(DiffStencils%StencilZ(Nz,:,1),sol(idx,i,j))
  END DO
END DO
END SUBROUTINE VerticalFreeSurfaceVelocity
