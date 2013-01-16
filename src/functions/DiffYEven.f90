SUBROUTINE DiffYEven(var,varY,order,Nx,Ny,Nz,DiffStencils,beta)
!
! Explicit order differentiation in Y-direction.
!
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Diff_def) :: DiffStencils
INTEGER :: Nx, Ny, Nz, rank, order, beta,i,j,k
REAL(KIND=long), DIMENSION(Nz,Nx,Ny) :: var, varY
REAL(KIND=long), DIMENSION(2*beta+1) :: Stencil
rank = 2*beta+1
DO j = 1, beta
	Stencil = DiffStencils%StencilY(j,:,order)
	DO i = 1, Nx
		DO k = 1, Nz
			varY(k,i,j) = DOT_PRODUCT(Stencil,var(k,i,1:rank))
		END DO
	END DO
END DO
DO j = beta+1, Ny-beta
	Stencil = DiffStencils%StencilY(j,:,order)
	DO i = 1, Nx
		DO k = 1, Nz
			varY(k,i,j) = DOT_PRODUCT(Stencil,var(k,i,j-beta:j+beta))
		END DO
	END DO
END DO
DO j = Ny-beta+1, Ny
	Stencil = DiffStencils%StencilY(j,:,order)
	DO i = 1, Nx
		DO k = 1, Nz
			varY(k,i,j) = DOT_PRODUCT(Stencil,var(k,i,Ny+1-rank:Ny))
		END DO
	END DO
END DO
END SUBROUTINE DiffYEven
