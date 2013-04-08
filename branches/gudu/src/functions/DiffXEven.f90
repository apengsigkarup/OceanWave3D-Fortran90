SUBROUTINE DiffXEven(var,varX,order,Nx,Ny,Nz,DiffStencils,alpha)
!
! Explicit order differentiation in X-direction.
!
! No conditional loops.
!
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Diff_def) :: DiffStencils
INTEGER :: Nx, Ny, Nz, rank, order, alpha, i,j,k
REAL(KIND=long), DIMENSION(Nz,Nx,Ny) :: var, varX
REAL(KIND=long), DIMENSION(2*alpha+1) :: Stencil
rank = 2*alpha+1
DO i = 1, alpha
	Stencil = DiffStencils%StencilX(i,1:rank,order)
	DO j = 1, Ny
		DO k = 1, Nz
			varX(k,i,j) = DOT_PRODUCT(Stencil,var(k,1:rank,j))
		END DO
	END DO
END DO
DO i = alpha+1, Nx-alpha
	Stencil = DiffStencils%StencilX(i,1:rank,order)
	DO j = 1, Ny
		DO k = 1, Nz
			varX(k,i,j) = DOT_PRODUCT(Stencil,var(k,i-alpha:i+alpha,j))
		END DO
	END DO
END DO
DO i = Nx-alpha+1, Nx
	Stencil = DiffStencils%StencilX(i,1:rank,order)
	DO j = 1, Ny
		DO k = 1, Nz
			varX(k,i,j) = DOT_PRODUCT(Stencil,var(k,Nx+1-rank:Nx,j))
		END DO
	END DO
END DO
END SUBROUTINE DiffXEven
