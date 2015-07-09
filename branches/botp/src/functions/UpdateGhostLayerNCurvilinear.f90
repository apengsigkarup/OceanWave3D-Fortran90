SUBROUTINE UpdateGhostLayerNCurvilinear(var,Sx,Sy,NNx,NNy,NNz,CompGrid,alpha,beta,GhostGridX,GhostGridY)
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: CompGrid
INTEGER :: NNx, NNy, NNz, rank, alpha, beta, i,j,k, GhostGridX, GhostGridY, diffb
REAL(KIND=long), DIMENSION(NNx,NNy) :: var, Sx, Sy
REAL(KIND=long), DIMENSION(2*alpha+1) :: Stenciln
REAL(KIND=long) :: xn, yn, Source, NormalX, NormalY, nx, ny
INTEGER, DIMENSION(2*beta+1) :: idx
IF (alpha/=beta) THEN
	print*,'ERROR: alpha/=beta. (UpdateGhostLayerECurvilinear)'
	STOP
END IF
rank = 2*beta+1
! On the free surface
k = NNz
! two boundaries
j=1
DO i = 1+GhostGridX, NNx-GhostGridX
  	Stenciln = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rank,1+GhostGridY,1)
    !! Source terms on the boundary ! FIXME: check the formula...
    !xn = CompGrid%CurvilinearStuff%xn(i,j+GhostGridY)
    !yn = CompGrid%CurvilinearStuff%yn(i,j+GhostGridY)
    !Source = xn*Sx(i,j+GhostGridY) + yn*Sy(i,j+GhostGridY)
    ! use geometric factors for wall boundary
	nx = CompGrid%CurvilinearStuff%nx(i,j+GhostGridY)
	ny = CompGrid%CurvilinearStuff%ny(i,j+GhostGridY)
   	! Normal vectors are defined at the ghost points used to impose the kinematic boundary conditions
	NormalX = CompGrid%CurvilinearStuff%NormalX(k,i,j)
	NormalY = CompGrid%CurvilinearStuff%NormalY(k,i,j)
    Source = (NormalX*Sx(i,j+GhostGridY) + NormalY*Sy(i,j+GhostGridY))/(NormalX*nx+NormalY*ny)
    !
    var(i,j) = (Source-DOT_PRODUCT(Stenciln(2:rank),var(i,2:rank))) / Stenciln(1)
END DO
j = NNy
DO i = 1+GhostGridX, NNx-GhostGridX
   	Stenciln = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rank,rank-GhostGridY,1)
    !! Source terms on the boundary! FIXME: check the formula...
    !xn = CompGrid%CurvilinearStuff%xn(i,j-GhostGridY)
    !yn = CompGrid%CurvilinearStuff%yn(i,j-GhostGridY)
    !Source = xn*Sx(i,j-GhostGridY) + yn*Sy(i,j-GhostGridY)
    ! use geometric factors for wall boundary
	nx = CompGrid%CurvilinearStuff%nx(i,j-GhostGridY)
	ny = CompGrid%CurvilinearStuff%ny(i,j-GhostGridY)
    NormalX = CompGrid%CurvilinearStuff%NormalX(k,i,j)
	NormalY = CompGrid%CurvilinearStuff%NormalY(k,i,j)
    Source = (NormalX*Sx(i,j-GhostGridY) + NormalY*Sy(i,j-GhostGridY))/(NormalX*nx+NormalY*ny)
    !
    var(i,j) = (Source-DOT_PRODUCT(Stenciln(1:rank-1),var(i,j-rank+1:j-1))) / Stenciln(rank)
END DO

END SUBROUTINE UpdateGhostLayerNCurvilinear
