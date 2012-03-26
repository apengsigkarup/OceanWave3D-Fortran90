SUBROUTINE UpdateGhostLayerECurvilinear(var,Sx,Sy,NNx,NNy,NNz,CompGrid,alpha,beta,GhostGridX,GhostGridY)
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: CompGrid
INTEGER :: NNx, NNy, NNz, rank, alpha, beta, i,j,k, GhostGridX, GhostGridY, diffb
REAL(KIND=long), DIMENSION(NNx,NNy) :: var, Sx, Sy
REAL(KIND=long), DIMENSION(2*alpha+1) :: Stencile
REAL(KIND=long) :: xe, ye, Source, NormalX, NormalY, ex, ey
IF (alpha/=beta) THEN
	print*,'ERROR: alpha/=beta. (UpdateGhostLayerECurvilinear)'
	STOP
END IF
rank = 2*alpha+1
! On the free surface
k = NNz
! two boundaries
i=1
DO j = 1+GhostGridY, NNy-GhostGridY
  	Stencile = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rank,1+GhostGridX,1)
    !! Source terms on the boundary! FIXME: check the formula...
    !xe = CompGrid%CurvilinearStuff%xe(i+GhostGridX,j)
    !ye = CompGrid%CurvilinearStuff%ye(i+GhostGridX,j)
    ! use geometric factors for wall boundary
	ex = CompGrid%CurvilinearStuff%ex(i+GhostGridX,j)
	ey = CompGrid%CurvilinearStuff%ey(i+GhostGridX,j)
    !Source = xe*Sx(i+GhostGridX,j) + ye*Sy(i+GhostGridX,j)
    NormalX = CompGrid%CurvilinearStuff%NormalX(k,i,j)
	NormalY = CompGrid%CurvilinearStuff%NormalY(k,i,j)
    Source = (NormalX*Sx(i+GhostGridX,j) + NormalY*Sy(i+GhostGridX,j))/(NormalX*ex+NormalY*ey)
    !
    var(i,j) = (Source-DOT_PRODUCT(Stencile(2:rank),var(2:rank,j))) / Stencile(1)
END DO
i = NNx
DO j = 1+GhostGridY, NNy-GhostGridY
   	Stencile = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rank,rank-GhostGridX,1)
    !! Source terms on the boundary! FIXME: check the formula...
    !xe = CompGrid%CurvilinearStuff%xe(i-GhostGridX,j)
    !ye = CompGrid%CurvilinearStuff%ye(i-GhostGridX,j)
    ! use geometric factors for wall boundary
	ex = CompGrid%CurvilinearStuff%ex(i-GhostGridX,j)
	ey = CompGrid%CurvilinearStuff%ey(i-GhostGridX,j)
    !Source = xe*Sx(i-GhostGridX,j) + ye*Sy(i-GhostGridX,j)
    NormalX = CompGrid%CurvilinearStuff%NormalX(k,i,j)
	NormalY = CompGrid%CurvilinearStuff%NormalY(k,i,j)
    Source = (NormalX*Sx(i-GhostGridX,j) + NormalY*Sy(i-GhostGridX,j))/(NormalX*ex+NormalY*ey)
    !
    var(i,j) = (Source-DOT_PRODUCT(Stencile(1:rank-1),var(i-rank+1:i-1,j))) / Stencile(rank)
END DO

END SUBROUTINE UpdateGhostLayerECurvilinear
