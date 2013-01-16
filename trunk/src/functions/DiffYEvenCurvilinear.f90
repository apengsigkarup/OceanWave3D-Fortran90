SUBROUTINE DiffYEvenCurvilinear(var,varY,Nxg,Nyg,Nzg,CompGrid,kappa)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: CompGrid
INTEGER :: Nxg, Nyg, Nzg, order, kappa, i,j,k
REAL(KIND=long), DIMENSION(Nzg,Nxg,Nyg) :: var, varY, vare, varn
REAL(KIND=long), DIMENSION(2*kappa+1,2*kappa+1,2) :: DiffStencilArray
REAL(KIND=long) :: ey, ny

DiffStencilArray = CompGrid%CurvilinearStuff%DiffStencils%StencilG

CALL DiffXuniform3D(var,vare,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
CALL DiffYuniform3D(var,varn,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?

DO i = 1, Nxg
	DO j = 1, Nyg
		ey = CompGrid%CurvilinearStuff%ey(i,j)
		ny = CompGrid%CurvilinearStuff%ny(i,j)
		DO k = 1, Nzg
			varY(k,i,j) = ey*vare(k,i,j)+ny*varn(k,i,j)
		END DO
	END DO
END DO
END SUBROUTINE DiffYEvenCurvilinear
