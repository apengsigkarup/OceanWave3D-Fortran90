SUBROUTINE DiffZArbitrary(var,varZ,order,Nx,Ny,Nz,DiffStencils,gamma)
!
! Explicit order differentiation in Z-direction.
!
! Full rank stencils always applied.
!
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Diff_def) :: DiffStencils
INTEGER :: Nx, Ny, Nz, rank, order, gamma,i,j,k
REAL(KIND=long), DIMENSION(Nz,Nx,Ny) :: var, varZ
rank = 2*gamma+1
DO j = 1, Ny
	DO i = 1, Nx
		DO k = 1, gamma
			varZ(k,i,j) = DOT_PRODUCT(DiffStencils%StencilZ(k,1:rank,order),var(1:rank,i,j))
		END DO
	END DO
END DO
DO j = 1, Ny
	DO i = 1, Nx
		DO k = gamma+1, Nz-gamma
			varZ(k,i,j) = DOT_PRODUCT(DiffStencils%StencilZ(k,1:rank,order),var(k-gamma:k+gamma,i,j))
		END DO
	END DO
END DO
DO j = 1, Ny
	DO i = 1, Nx
		DO k = Nz-gamma+1, Nz
			varZ(k,i,j) = DOT_PRODUCT(DiffStencils%StencilZ(k,1:rank,order),var(Nz+1-rank:Nz,i,j))
		END DO
	END DO
END DO
END SUBROUTINE DiffZArbitrary
