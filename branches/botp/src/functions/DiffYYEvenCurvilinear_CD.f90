SUBROUTINE DiffYYEvenCurvilinear_CD(var,varYY,Nxg,Nyg,Nzg,CompGrid,kappa,IndexX,StencilX,IndexY,StencilY)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: CompGrid
INTEGER :: Nxg, Nyg, Nzg, order, kappa, i, j, k
REAL(KIND=long), DIMENSION(Nzg,Nxg,Nyg) :: var, varYY, vare, varn, varee, varnn, varen, vare2
REAL(KIND=long), DIMENSION(2*kappa+1,2*kappa+1,2) :: DiffStencilArray
REAL(KIND=long) :: ey, eyy, ny, nyy ,DMe, DMn, DM2en, DM2e, DM2n
! GD: addition for correct cross derivatives near corners...
INTEGER, DIMENSION(Nzg,Nxg,Nyg,2*kappa+1) :: IndexX, IndexY
REAL(KIND=long), DIMENSION(Nzg,Nxg,Nyg,2*kappa+1) :: StencilX, StencilY

DiffStencilArray = CompGrid%CurvilinearStuff%DiffStencils%StencilG

CALL DiffXuniform3D(var,vare,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
CALL DiffYuniform3D(var,varn,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
CALL DiffXuniform3D(var,varee,2,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
CALL DiffYuniform3D(var,varnn,2,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
!! mixed
!CALL DiffYuniform3D(vare,varen,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
! GD: modif for correct cross derivatives near corners...
! mixed
CALL DiffXuniform3D_CD(var,vare2,IndexX(:,:,:,:),StencilX(:,:,:,:),Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
CALL DiffYuniform3D_CD(vare2,varen,IndexY(:,:,:,:),StencilY(:,:,:,:),Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?

DO i = 1, Nxg
	DO j = 1, Nyg
		ey   = CompGrid%CurvilinearStuff%ey(i,j)
		ny   = CompGrid%CurvilinearStuff%ny(i,j)
		eyy  = CompGrid%CurvilinearStuff%eyy(i,j)
		nyy  = CompGrid%CurvilinearStuff%nyy(i,j)

		DO k = 1, Nzg
			DMe  = vare(k,i,j)
			DM2e = varee(k,i,j)
			DMn  = varn(k,i,j)
			DM2n = varnn(k,i,j)
			DM2en = varen(k,i,j)
			varYY(k,i,j) = eyy*DMe+nyy*DMn+two*ey*ny*DM2en+ey**2*DM2e+ny**2*DM2n
		END DO
	END DO
END DO
END SUBROUTINE DiffYYEvenCurvilinear_CD
