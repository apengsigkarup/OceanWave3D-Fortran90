SUBROUTINE UpdateGhostLayerYCurvilinear(var,S,NNx,NNy,CompGrid,alpha,beta,GhostGridX,GhostGridY)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: CompGrid
INTEGER :: NNx, NNy, rank, alpha, beta, i,j, GhostGridX, GhostGridY, diffa
REAL(KIND=long), DIMENSION(NNx,NNy) :: var, S
REAL(KIND=long), DIMENSION(2*alpha+1) :: Stencile, Stenciln
REAL(KIND=long) :: ex, ey, nx, ny, tmp, diagn
INTEGER, DIMENSION(2*alpha+1) :: idx
IF (alpha/=beta) THEN
	print*,'ERROR: alpha/=beta. (UpdateGhostLayerYCurvilinear)'
END IF
rank = 2*beta+1 !GD change: 2*alpha+1
DO i = -alpha,alpha
	idx(i+alpha+1) = i
END DO
! two boundaries
j=1
DO i = 1+GhostGridX, NNx-GhostGridX
	if (i-1<alpha+1) then
		diffa = (alpha+1)-(i-1)
	else if (i+1>NNx-alpha) then
		diffa = (NNx-alpha)-(i+1) !GD: change (NNx-beta)-(i+1)
	else
		diffa = 0
	end if

	! use geometric factors for wall boundary
!	ex = CompGrid%CurvilinearStuff%ex(i,j+GhostGridY)
	ey = CompGrid%CurvilinearStuff%ey(i,j+GhostGridY)
!	nx = CompGrid%CurvilinearStuff%nx(i,j+GhostGridY)
	ny = CompGrid%CurvilinearStuff%ny(i,j+GhostGridY)
	
	Stencile = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rank,1+alpha-diffa,1)
	Stenciln = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rank,1+GhostGridY,1)
	diagn = ny*Stenciln(1)

	tmp = ey * DOT_PRODUCT(Stencile,var(i+idx+diffa,j+GhostGridY) ) + ny * DOT_PRODUCT( Stenciln(2:rank),var(i,2:rank) )

    var(i,j) = (S(i,j+GhostGridY) - tmp) / diagn

END DO
j = NNy
DO i = 1+GhostGridX, NNx-GhostGridX
	if (i-1<alpha+1) then
		diffa = (alpha+1)-(i-1)
	else if (i+1>NNx-alpha) then
		diffa = (NNx-alpha)-(i+1)
	else
		diffa = 0
	end if

	! use geometric factors for wall boundary
!	ex = CompGrid%CurvilinearStuff%ex(i,j-GhostGridY)
	ey = CompGrid%CurvilinearStuff%ey(i,j-GhostGridY)
!	nx = CompGrid%CurvilinearStuff%nx(i,j-GhostGridY)
	ny = CompGrid%CurvilinearStuff%ny(i,j-GhostGridY)
	
	Stencile = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rank,1+alpha-diffa,1)
	Stenciln = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rank,rank-GhostGridY,1)

	diagn = ny*Stenciln(rank)

	tmp = ey * DOT_PRODUCT(Stencile,var(i+idx+diffa,j-GhostGridY) ) + ny * DOT_PRODUCT( Stenciln(1:rank-1),var(i,j-rank+1:j-1) )

    var(i,j) = (S(i,j-GhostGridY)- tmp) / diagn
END DO
END SUBROUTINE UpdateGhostLayerYCurvilinear
