SUBROUTINE DiffXuniform3D(var,varX,order,DiffStencilArray,Nx,Ny,Nz,alpha)
! By Allan P. Engsig-Karup.
!FIXME: this script has not been tested for bugs
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx, Ny, Nz, rank, order, alpha, i,j,k
REAL(KIND=long), DIMENSION(2*alpha+1,2*alpha+1,2) :: DiffStencilArray
REAL(KIND=long), DIMENSION(Nz,Nx,Ny) :: var, varX
REAL(KIND=long), DIMENSION(2*alpha+1) :: Stencil
rank = 2*alpha+1
! one-sided
DO i = 1, alpha
	Stencil = DiffStencilArray(1:rank,i,order)
	DO j = 1, Ny
		DO k = 1, Nz
			varX(k,i,j) = DOT_PRODUCT(Stencil,var(k,1:rank,j))
		END DO
	END DO
END DO
! centered
Stencil = DiffStencilArray(1:rank,alpha+1,order)
DO i = alpha+1, Nx-alpha
	DO j = 1, Ny
		DO k = 1, Nz
			varX(k,i,j) = DOT_PRODUCT(Stencil,var(k,i-alpha:i+alpha,j))
		END DO
	END DO
END DO
! one-sided
DO i = Nx-alpha+1, Nx
	Stencil = DiffStencilArray(1:rank,2*alpha+1-(Nx-i),order)
	DO j = 1, Ny
		DO k = 1, Nz
			varX(k,i,j) = DOT_PRODUCT(Stencil,var(k,Nx+1-rank:Nx,j))
		END DO
	END DO
END DO

END SUBROUTINE DiffXuniform3D
