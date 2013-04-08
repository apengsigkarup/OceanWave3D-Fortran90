SUBROUTINE DiffZuniform3D_CD(var,varZ,order,IndexesZ,StencilsZ,Nx,Ny,Nz,alpha)
! By Allan P. Engsig-Karup.
!FIXME: this script has not been tested for bugs
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx, Ny, Nz, rank, order, alpha, i, j, k
REAL(KIND=long), DIMENSION(Nz,Nx,Ny) :: var, varZ
! GD: addition for correct cross derivatives near corners...
INTEGER, DIMENSION(Nz,Nx,Ny,2*alpha+1) :: IndexesZ
REAL(KIND=long), DIMENSION(Nz,Nx,Ny,2*alpha+1) :: StencilsZ
!
rank = 2*alpha+1
!
DO i=1,Nx
  DO j = 1, Ny
		DO k = 1, Nz
            varZ(k,i,j) = DOT_PRODUCT(StencilsZ(k,i,j,1:rank),var(IndexesZ(k,i,j,1:rank),i,j))
		END DO
	END DO
ENDDO

END SUBROUTINE DiffZuniform3D_CD
