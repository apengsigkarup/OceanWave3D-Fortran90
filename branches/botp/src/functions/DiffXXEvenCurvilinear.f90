SUBROUTINE DiffXXEvenCurvilinear(var,varXX,Nxg,Nyg,Nzg,CompGrid,kappa)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: CompGrid
INTEGER :: Nxg, Nyg, Nzg, order, kappa, i,j,k
REAL(KIND=long), DIMENSION(Nzg,Nxg,Nyg) :: var, varXX, vare, varn, varee, varnn, varen
REAL(KIND=long), DIMENSION(2*kappa+1,2*kappa+1,2) :: DiffStencilArray
REAL(KIND=long) :: ex, exx, nx, nxx ,DMe, DMn, DM2en, DM2e, DM2n

DiffStencilArray = CompGrid%CurvilinearStuff%DiffStencils%StencilG

CALL DiffXuniform3D(var,vare,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
CALL DiffYuniform3D(var,varn,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
CALL DiffXuniform3D(var,varee,2,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
CALL DiffYuniform3D(var,varnn,2,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
! mixed
CALL DiffYuniform3D(vare,varen,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?

DO i = 1, Nxg
	DO j = 1, Nyg
		ex   = CompGrid%CurvilinearStuff%ex(i,j)
		nx   = CompGrid%CurvilinearStuff%nx(i,j)
		exx  = CompGrid%CurvilinearStuff%exx(i,j)
		nxx  = CompGrid%CurvilinearStuff%nxx(i,j)
		DO k = 1, Nzg
			DMe  = vare(k,i,j)
			DM2e = varee(k,i,j)
			DMn  = varn(k,i,j)
			DM2n = varnn(k,i,j)
			DM2en = varen(k,i,j)
			varXX(k,i,j) = exx*DMe+nxx*DMn+two*ex*nx*DM2en+ex**2*DM2e+nx**2*DM2n
		END DO
	END DO
END DO
END SUBROUTINE DiffXXEvenCurvilinear
