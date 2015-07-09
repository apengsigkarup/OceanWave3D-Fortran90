SUBROUTINE DiffXEvenCurvilinear(var,varX,Nxg,Nyg,Nzg,CompGrid,kappa)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: CompGrid
INTEGER :: Nxg, Nyg, Nzg, order, kappa, i,j,k
REAL(KIND=long), DIMENSION(Nzg,Nxg,Nyg) :: var, varX, vare, varn
REAL(KIND=long), DIMENSION(2*kappa+1,2*kappa+1,2) :: DiffStencilArray
REAL(KIND=long) :: ex, nx

DiffStencilArray = CompGrid%CurvilinearStuff%DiffStencils%StencilG

CALL DiffXuniform3D(var,vare,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
CALL DiffYuniform3D(var,varn,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?

DO i = 1, Nxg
	DO j = 1, Nyg
		ex = CompGrid%CurvilinearStuff%ex(i,j)
		nx = CompGrid%CurvilinearStuff%nx(i,j)
		DO k = 1, Nzg
			varX(k,i,j) = ex*vare(k,i,j)+nx*varn(k,i,j)
		END DO
	END DO
END DO
END SUBROUTINE DiffXEvenCurvilinear
