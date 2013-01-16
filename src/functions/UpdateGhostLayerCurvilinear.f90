SUBROUTINE UpdateGhostLayerCurvilinear(var,Sx,Sy,NNx,NNy,NNz,CompGrid,alpha,beta,GhostGridX,GhostGridY)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: CompGrid
INTEGER :: NNx, NNy, NNz, ranka, rankb,alpha,beta,i,j,k,GhostGridX, GhostGridY, diffa, diffb
REAL(KIND=long), DIMENSION(NNx,NNy) :: var, Sx, Sy
REAL(KIND=long), DIMENSION(2*alpha+1) :: Stencile
REAL(KIND=long), DIMENSION(2*beta+1)  :: Stenciln
REAL(KIND=long) :: ex, ey, nx, ny, tmp, diage, diagn
REAL(KIND=long) :: NormalX, NormalY, NormalZ
INTEGER, DIMENSION(2*alpha+1) :: idxa
INTEGER, DIMENSION(2*beta+1)  :: idxb
IF (alpha/=beta) THEN
	print*,'ERROR: alpha/=beta. (UpdateGhostLayerCurvilinear)'
	STOP
END IF
ranka = 2*alpha+1
DO i = -alpha, alpha
	idxa(i+alpha+1) = i !GD: modif, idxa(i+beta+1) = i
END DO
rankb = 2*beta+1
DO i = -beta, beta
	idxb(i+beta+1) = i
END DO
! four boundaries
k = NNz
! WEST
i=1
DO j = 1+GhostGridY, NNy-GhostGridY
	if (j-1<beta+1) then
		diffb = (beta+1)-(j-1)  ! substract one to ensure that it is only dependent on one ghost point.
	else if (j+1>NNy-beta) then
		diffb = (NNy-beta)-(j+1)  ! substract one to ensure that it is only dependent on one ghost point.
	else
		diffb = 0
	end if

	! Normal vectors are defined at the ghost points used to impose the kinematic boundary conditions
	NormalX = CompGrid%CurvilinearStuff%NormalX(k,i,j)
	NormalY = CompGrid%CurvilinearStuff%NormalY(k,i,j)

	! use geometric factors for wall boundary
	ex = CompGrid%CurvilinearStuff%ex(i+GhostGridX,j)
	nx = CompGrid%CurvilinearStuff%nx(i+GhostGridX,j)
	ey = CompGrid%CurvilinearStuff%ey(i+GhostGridX,j)
	ny = CompGrid%CurvilinearStuff%ny(i+GhostGridX,j)
	
	Stencile = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:ranka,1+GhostGridX,1)
	Stenciln = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rankb,1+beta-diffb,1)

	diage = (NormalX*ex+NormalY*ey)*Stencile(1)
	tmp = (NormalX*ex+NormalY*ey) * DOT_PRODUCT(Stencile(2:ranka),var(2:ranka,j) ) + &
	      (NormalX*nx+NormalY*ny) * DOT_PRODUCT(Stenciln,var(i+GhostGridX,j+idxb+diffb) )

    var(i,j) = (NormalX*Sx(i+GhostGridX,j)+NormalY*Sy(i+GhostGridX,j)- tmp) / diage
END DO
! EAST
i = NNx
DO j = 1+GhostGridY, NNy-GhostGridY
	if (j-1<beta+1) then
		diffb = (beta+1)-(j-1)
	else if (j+1>NNy-beta) then
		diffb = (NNy-beta)-(j+1)
	else
		diffb = 0
	end if

	! Normal vectors are defined at the ghost points used to impose the kinematic boundary conditions
	NormalX = CompGrid%CurvilinearStuff%NormalX(k,i,j)
	NormalY = CompGrid%CurvilinearStuff%NormalY(k,i,j)

	! use geometric factors for wall boundary
	ex = CompGrid%CurvilinearStuff%ex(i-GhostGridX,j)
	ey = CompGrid%CurvilinearStuff%ey(i-GhostGridX,j)
	nx = CompGrid%CurvilinearStuff%nx(i-GhostGridX,j)
	ny = CompGrid%CurvilinearStuff%ny(i-GhostGridX,j)
	
	Stencile = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:ranka,ranka-GhostGridX,1)
	Stenciln = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rankb,1+beta-diffb,1)

	diage = (NormalX*ex+NormalY*ey)*Stencile(ranka)

	tmp = (NormalX*ex+NormalY*ey) * DOT_PRODUCT(Stencile(1:ranka-1),var(i-ranka+1:i-1,j) ) +&
		  (NormalX*nx+NormalY*ny) * DOT_PRODUCT( Stenciln,var(i-GhostGridX,j+idxb+diffb) ) ! GD: error corrected ey->ny

    var(i,j) = (NormalX*Sx(i-GhostGridX,j)+NormalY*Sy(i-GhostGridX,j)- tmp) / diage
END DO
! SOUTH
j=1
DO i = 1+GhostGridX, NNx-GhostGridX
	if (i-1<alpha+1) then
		diffa = (alpha+1)-(i-1)
	else if (i+1>NNx-alpha) then
		diffa = (NNx-alpha)-(i+1) !GD: change (NNx-beta)-(i+1)
	else
		diffa = 0
	end if

	! Normal vectors are defined at the ghost points used to impose the kinematic boundary conditions
	NormalX = CompGrid%CurvilinearStuff%NormalX(k,i,j)
	NormalY = CompGrid%CurvilinearStuff%NormalY(k,i,j)

	! use geometric factors for wall boundary
	ex = CompGrid%CurvilinearStuff%ex(i,j+GhostGridY)
	ey = CompGrid%CurvilinearStuff%ey(i,j+GhostGridY)
	nx = CompGrid%CurvilinearStuff%nx(i,j+GhostGridY)
	ny = CompGrid%CurvilinearStuff%ny(i,j+GhostGridY)
	
	Stencile = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:ranka,1+alpha-diffa,1)
	Stenciln = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rankb,1+GhostGridY,1)
	diagn = (NormalX*nx + NormalY*ny)*Stenciln(1)

	tmp = (NormalX*ex+NormalY*ey) * DOT_PRODUCT(Stencile,var(i+idxa+diffa,j+GhostGridY) ) +&
		  (NormalX*nx+NormalY*ny) * DOT_PRODUCT( Stenciln(2:rankb),var(i,2:rankb) )

    var(i,j) = (NormalX*Sx(i,j+GhostGridY)+NormalY*Sy(i,j+GhostGridY) - tmp) / diagn
END DO
! NORTH
j = NNy
DO i = 1+GhostGridX, NNx-GhostGridX
	if (i-1<alpha+1) then
		diffa = (alpha+1)-(i-1)
	else if (i+1>NNx-alpha) then
		diffa = (NNx-alpha)-(i+1)
	else
		diffa = 0
	end if

	! Normal vectors are defined at the ghost points used to impose the kinematic boundary conditions
	NormalX = CompGrid%CurvilinearStuff%NormalX(k,i,j)
	NormalY = CompGrid%CurvilinearStuff%NormalY(k,i,j)

	! use geometric factors for wall boundary
	ex = CompGrid%CurvilinearStuff%ex(i,j-GhostGridY)
	ey = CompGrid%CurvilinearStuff%ey(i,j-GhostGridY)
	nx = CompGrid%CurvilinearStuff%nx(i,j-GhostGridY)
	ny = CompGrid%CurvilinearStuff%ny(i,j-GhostGridY)
	
	Stencile = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:ranka,1+alpha-diffa,1)
	Stenciln = CompGrid%CurvilinearStuff%DiffStencils%StencilG(1:rankb,rankb-GhostGridY,1)
	diagn = (NormalX*nx + NormalY*ny)*Stenciln(rankb)

	tmp = (NormalX*ex+NormalY*ey) * DOT_PRODUCT(Stencile,var(i+idxa+diffa,j-GhostGridY) ) +&
	      (NormalX*nx+NormalY*ny) * DOT_PRODUCT(Stenciln(1:rankb-1),var(i,j-rankb+1:j-1) )

    var(i,j) = (NormalX*Sx(i,j-GhostGridY)+NormalY*Sy(i,j-GhostGridY)- tmp) / diagn
END DO

END SUBROUTINE UpdateGhostLayerCurvilinear
