SUBROUTINE DiffZuniform3D(var,varZ,order,DiffStencilArray,Nx,Ny,Nz,alpha)
! By Allan P. Engsig-Karup.
!FIXME: this script has not been tested for bugs
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx, Ny, Nz, rank, order, alpha, i,j,k
REAL(KIND=long), DIMENSION(2*alpha+1,2*alpha+1,2) :: DiffStencilArray
REAL(KIND=long), DIMENSION(Nz,Nx,Ny) :: var, varZ
REAL(KIND=long), DIMENSION(2*alpha+1) :: Stencil
rank = 2*alpha+1
! one-sided
DO k = 1, alpha
	Stencil = DiffStencilArray(1:rank,k,order)
	DO j = 1, Ny
		DO i = 1, Nx
			varZ(k,i,j) = DOT_PRODUCT(Stencil,var(1:rank,i,j))
		END DO
	END DO
END DO
! centered
Stencil = DiffStencilArray(1:rank,alpha+1,order)
DO j = 1, Ny
	DO i = 1, Nx
		DO k = alpha+1, Nz-alpha
			varZ(k,i,j) = DOT_PRODUCT(Stencil,var(k-alpha:k+alpha,i,j))
		END DO
	END DO
END DO
! one-sided
DO k = Nz-alpha+1, Nz
	Stencil = DiffStencilArray(1:rank,2*alpha+1-(Nz-k),order)
	DO j = 1, Ny
		DO i = 1, Nx
			varZ(k,i,j) = DOT_PRODUCT(Stencil,var(Nz+1-rank:Nz,i,j))
		END DO
	END DO
END DO

END SUBROUTINE DiffZuniform3D
