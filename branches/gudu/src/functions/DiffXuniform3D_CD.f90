SUBROUTINE DiffXuniform3D_CD(var,varX,IndexesX,StencilsX,Nx,Ny,Nz,alpha)
! By Allan P. Engsig-Karup.
!FIXME: this script has not been tested for bugs
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx, Ny, Nz, rank, alpha, i, j, k
REAL(KIND=long), DIMENSION(Nz,Nx,Ny) :: var, varX
! GD: addition for correct cross derivatives near corners...
INTEGER, DIMENSION(Nz,Nx,Ny,2*alpha+1) :: IndexesX
REAL(KIND=long), DIMENSION(Nz,Nx,Ny,2*alpha+1) :: StencilsX
!
rank = 2*alpha+1
!
DO i=1,Nx
  DO j = 1, Ny
		DO k = 1, Nz
            varX(k,i,j) = DOT_PRODUCT(StencilsX(k,i,j,1:rank),var(k,IndexesX(k,i,j,1:rank),j))
		END DO
	END DO
ENDDO

END SUBROUTINE DiffXuniform3D_CD
