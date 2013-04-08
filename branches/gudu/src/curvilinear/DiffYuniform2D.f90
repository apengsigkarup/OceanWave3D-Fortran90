SUBROUTINE DiffYuniform2D(var,varY,order,DiffStencilArray,Nx,Ny,beta)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx, Ny, rank, order, beta, i,j
REAL(KIND=long), DIMENSION(2*beta+1,2*beta+1,2) :: DiffStencilArray
REAL(KIND=long), DIMENSION(Nx,Ny) :: var, varY
REAL(KIND=long), DIMENSION(2*beta+1) :: Stencil
rank = 2*beta+1
! one-sided
DO j = 1, beta
	Stencil = DiffStencilArray(1:rank,j,order)
	DO i = 1, Nx
		varY(i,j) = DOT_PRODUCT(Stencil,var(i,1:rank))
	END DO
END DO
! centered
Stencil = DiffStencilArray(1:rank,beta+1,order)
DO i = 1, Nx
	DO j = beta+1, Ny-beta
		varY(i,j) = DOT_PRODUCT(Stencil,var(i,j-beta:j+beta))
	END DO
END DO
! one-sided
DO j = Ny-beta+1, Ny
	Stencil = DiffStencilArray(1:rank,2*beta+1-(Ny-j),order)
	DO i = 1,Nx
		varY(i,j) = DOT_PRODUCT(Stencil,var(i,Ny+1-rank:Ny))
	END DO
END DO

END SUBROUTINE DiffYuniform2D
