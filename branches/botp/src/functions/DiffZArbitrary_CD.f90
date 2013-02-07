SUBROUTINE DiffZArbitrary_CD(var,varZ,IndexesZ,StencilsZ,Nx,Ny,Nz,gamma)
!
! Explicit order differentiation in Z-direction.
!
! Full rank stencils always applied.
!
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx, Ny, Nz, rank, order, gamma,i, j, k
REAL(KIND=long), DIMENSION(Nz,Nx,Ny) :: var, varZ
!GD: addition for correct cross derivatives near corners...
INTEGER, DIMENSION(Nz,Nx,Ny,2*gamma+1) :: IndexesZ
REAL(KIND=long), DIMENSION(Nz,Nx,Ny,2*gamma+1) :: StencilsZ

rank = 2*gamma+1

DO j = 1, Ny
	DO i = 1, Nx
		DO k = 1, Nz
			varZ(k,i,j) = DOT_PRODUCT(StencilsZ(k,i,j,1:rank),var(IndexesZ(k,i,j,1:rank),i,j))
		END DO
	END DO
END DO

END SUBROUTINE DiffZArbitrary_CD
