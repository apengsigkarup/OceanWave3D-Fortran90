SUBROUTINE DiffYuniform3D_CD(var,varY,IndexesY,StencilsY,Nx,Ny,Nz,beta)
! By Allan P. Engsig-Karup.
!FIXME: this script has not been tested for bugs
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx, Ny, Nz, rank, beta, i, j, k
REAL(KIND=long), DIMENSION(Nz,Nx,Ny) :: var, varY
! GD: addition for correct cross derivatives near corners...
INTEGER, DIMENSION(Nz,Nx,Ny,2*beta+1) :: IndexesY
REAL(KIND=long), DIMENSION(Nz,Nx,Ny,2*beta+1) :: StencilsY
!
rank = 2*beta+1
!
DO i=1,Nx
  DO j = 1, Ny
		DO k = 1, Nz
            varY(k,i,j) = DOT_PRODUCT(StencilsY(k,i,j,1:rank),var(k,i,IndexesY(k,i,j,1:rank)))
		END DO
	END DO
ENDDO

END SUBROUTINE DiffYuniform3D_CD
