SUBROUTINE DiffYuniform3D(var,varY,order,DiffStencilArray,Nx,Ny,Nz,beta)
! By Allan P. Engsig-Karup.
!FIXME: this script has not been tested for bugs
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx, Ny, Nz, rank, order, beta, i,j, k
REAL(KIND=long), DIMENSION(2*beta+1,2*beta+1,2) :: DiffStencilArray
REAL(KIND=long), DIMENSION(Nz,Nx,Ny) :: var, varY
REAL(KIND=long), DIMENSION(2*beta+1) :: Stencil
rank = 2*beta+1
! one-sided
DO j = 1, beta
	Stencil = DiffStencilArray(1:rank,j,order)
	DO i = 1, Nx
		DO k = 1, Nz
			varY(k,i,j) = DOT_PRODUCT(Stencil,var(k,i,1:rank))
		END DO
	END DO
END DO
! centered
Stencil = DiffStencilArray(1:rank,beta+1,order)
DO i = 1, Nx
	DO j = beta+1, Ny-beta
		DO k = 1, Nz
			varY(k,i,j) = DOT_PRODUCT(Stencil,var(k,i,j-beta:j+beta))
		END DO
	END DO
END DO
! one-sided
DO j = Ny-beta+1, Ny
	Stencil = DiffStencilArray(1:rank,2*beta+1-(Ny-j),order)
	DO i = 1,Nx
		DO k = 1, Nz
			varY(k,i,j) = DOT_PRODUCT(Stencil,var(k,i,Ny+1-rank:Ny))
		END  DO
	END DO
END DO

END SUBROUTINE DiffYuniform3D
