SUBROUTINE DiffXuniform2D(var,varX,order,DiffStencilArray,Nx,Ny,alpha)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx, Ny, rank, order, alpha, i,j
REAL(KIND=long), DIMENSION(2*alpha+1,2*alpha+1,2) :: DiffStencilArray
REAL(KIND=long), DIMENSION(Nx,Ny) :: var, varX
REAL(KIND=long), DIMENSION(2*alpha+1) :: Stencil
rank = 2*alpha+1
! one-sided
DO i = 1, alpha
	Stencil = DiffStencilArray(1:rank,i,order)
	DO j = 1, Ny
		varX(i,j) = DOT_PRODUCT(Stencil,var(1:rank,j))
	END DO
END DO
! centered
Stencil = DiffStencilArray(1:rank,alpha+1,order)
DO i = alpha+1, Nx-alpha
	DO j = 1, Ny
		varX(i,j) = DOT_PRODUCT(Stencil,var(i-alpha:i+alpha,j))
	END DO
END DO
! one-sided
DO i = Nx-alpha+1, Nx
	Stencil = DiffStencilArray(1:rank,2*alpha+1-(Nx-i),order)
	DO j = 1, Ny
		varX(i,j) = DOT_PRODUCT(Stencil,var(Nx+1-rank:Nx,j))
	END DO
END DO

END SUBROUTINE DiffXuniform2D
