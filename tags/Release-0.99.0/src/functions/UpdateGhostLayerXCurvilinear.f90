SUBROUTINE UpdateGhostLayerXCurvilinear(var,S,NNx,NNy,CompGrid,alpha,beta,GhostGridX,GhostGridY)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: CompGrid
INTEGER :: NNx, NNy, rank, alpha, beta, i,j, GhostGridX, GhostGridY, diffb
REAL(KIND=long), DIMENSION(NNx,NNy) :: var, S
REAL(KIND=long), DIMENSION(2*alpha+1) :: Stencile, Stenciln
REAL(KIND=long) :: ex, ey, nx, ny, tmp, diage
INTEGER, DIMENSION(2*beta+1) :: idx
IF (alpha/=beta) THEN
	print*,'ERROR: alpha/=beta. (UpdateGhostLayerXCurvilinear)'
	STOP
END IF
rank = 2*alpha+1
DO i = -beta, beta
	idx(i+beta+1) = i
END DO
! two boundaries
i=1
DO j = 1+GhostGridY, NNy-GhostGridY
	if (j-1<beta+1) then
		diffb = (beta+1)-(j-1)  ! substract one to ensure that it is only dependent on one ghost point.
	else if (j+1>NNy-beta) then
		diffb = (NNy-beta)-(j+1)  ! substract one to ensure that it is only dependent on one ghost point.
	else
		diffb = 0
	end if

	! use geometric factors for wall boundary
	ex = CompGrid%CurvilinearStuff%ex(i+GhostGridX,j)
!	ey = CompGrid%CurvilinearStuff%ey(i+GhostGridX,j)
	nx = CompGrid%CurvilinearStuff%nx(i+GhostGridX,j)
!	ny = CompGrid%CurvilinearStuff%ny(i+GhostGridX,j)
	
	Stencile = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rank,1+GhostGridX,1)
	Stenciln = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rank,1+beta-diffb,1)

	diage = ex*Stencile(1)
	tmp = ex * DOT_PRODUCT(Stencile(2:rank),var(2:rank,j) ) + nx * DOT_PRODUCT( Stenciln,var(i+GhostGridX,j+idx+diffb) )

    var(i,j) = (S(i+GhostGridX,j)- tmp) / diage

END DO
i = NNx
DO j = 1+GhostGridY, NNy-GhostGridY
	if (j-1<beta+1) then
		diffb = (beta+1)-(j-1)
	else if (j+1>NNy-beta) then
		diffb = (NNy-beta)-(j+1)
	else
		diffb = 0
	end if

	! use geometric factors for wall boundary
	ex = CompGrid%CurvilinearStuff%ex(i-GhostGridX,j)
!	ey = CompGrid%CurvilinearStuff%ey(i-GhostGridX,j)
	nx = CompGrid%CurvilinearStuff%nx(i-GhostGridX,j)
!	ny = CompGrid%CurvilinearStuff%ny(i-GhostGridX,j)
	
	Stencile = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rank,rank-GhostGridX,1)
	Stenciln = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rank,1+beta-diffb,1)

	diage = ex*Stencile(rank)

	tmp = ex * DOT_PRODUCT(Stencile(1:rank-1),var(i-rank+1:i-1,j) ) + nx * DOT_PRODUCT( Stenciln,var(i-GhostGridX,j+idx+diffb) )

    var(i,j) = (S(i-GhostGridX,j)- tmp) / diage
END DO
END SUBROUTINE UpdateGhostLayerXCurvilinear
