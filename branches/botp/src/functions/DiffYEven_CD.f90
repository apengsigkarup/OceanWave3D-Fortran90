SUBROUTINE DiffYEven_CD(var,varY,IndexesY,StencilsY,Nx,Ny,Nz,beta)
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
! GD: addition for correct cross derivatives near corners...
INTEGER, DIMENSION(Nz,Nx,Ny,2*beta+1) :: IndexesY
REAL(KIND=long), DIMENSION(Nz,Nx,Ny,2*beta+1) :: StencilsY
INTEGER :: ii
rank = 2*beta+1

DO j = 1, Ny
	DO i = 1, Nx
		DO k = 1, Nz
			varY(k,i,j) = DOT_PRODUCT(StencilsY(k,i,j,1:rank),var(k,i,IndexesY(k,i,j,1:rank)))
		END DO
	END DO
END DO

END SUBROUTINE DiffYEven_CD
