SUBROUTINE DetermineGenericStencilsUniform(StencilsG,alpha)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: alpha, i, rank, order
REAL(KIND=long), DIMENSION(2*alpha+1,2*alpha+1,2) :: StencilsG
REAL(KIND=long) :: x0, grid(2*alpha+1)
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: c
rank = 2*alpha+1
ALLOCATE( c(rank,3) )
rank = 2*alpha+1
DO i = 1, rank
	grid(i) = REAL(i,long)
END DO
DO order = 1, 2
	DO i = 1, rank
		! Expansion point where the approximation is to be defined
		x0 = REAL(i,long)
		! Determine local finite difference stencil coefficients
		CALL weights(x0,grid,rank-1,rank-1,order,c)
		StencilsG(1:rank,i,order) = c(1:rank,order+1)
	END DO
END DO
DEALLOCATE(c)
END SUBROUTINE DetermineGenericStencilsUniform
