SUBROUTINE ComputeNormalVectors(CompGrid,GhostGridX,GhostGridY,GhostGridZ)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def), INTENT(INOUT) :: CompGrid
INTEGER :: i,j,k,GhostGridX,GhostGridY,GhostGridZ, Nxg, Nyg, Nzg
REAL(KIND=long) :: s, ex, ey, nx, ny, hx, hy

Nxg = CompGrid%Nx + 2*GhostGridX
Nyg = CompGrid%Ny + 2*GhostGridY
Nzg = CompGrid%Nz + GhostGridZ

ALLOCATE( CompGrid%CurvilinearStuff%NormalX(Nzg,Nxg,Nyg) )
ALLOCATE( CompGrid%CurvilinearStuff%NormalY(Nzg,Nxg,Nyg) )
ALLOCATE( CompGrid%CurvilinearStuff%NormalZ(Nzg,Nxg,Nyg) )
CompGrid%CurvilinearStuff%NormalX = zero
CompGrid%CurvilinearStuff%NormalY = zero
CompGrid%CurvilinearStuff%NormalZ = zero

IF (CompGrid%Nx>1) THEN
! West boundary
i = 1
DO j = 1+GhostGridY, CompGrid%Ny+GhostGridY ! omit ghost boundary corners
	DO k = 1+GhostGridZ, CompGrid%Nz+GhostGridZ ! omit bottom layer
		ex = CompGrid%CurvilinearStuff%ex(i+1,j)
		ey = CompGrid%CurvilinearStuff%ey(i+1,j)
		s = SQRT(ex**2+ey**2)
		CompGrid%CurvilinearStuff%NormalX(k,i,j) = -ex/s !west
		CompGrid%CurvilinearStuff%NormalY(k,i,j) = -ey/s !west
!		CompGrid%CurvilinearStuff%NormalZ(k,i,j) = zero 		
	END DO
END DO

! East boundary
i = CompGrid%Nx+2*GhostGridX
DO j = 1+GhostGridY, CompGrid%Ny+GhostGridY ! omit ghost boundary corners
	DO k = 1+GhostGridZ, CompGrid%Nz+GhostGridZ ! omit bottom layer
		ex = CompGrid%CurvilinearStuff%ex(i-1,j)
		ey = CompGrid%CurvilinearStuff%ey(i-1,j)
		s = SQRT(ex**2+ey**2)
		CompGrid%CurvilinearStuff%NormalX(k,i,j) = ex/s !east
		CompGrid%CurvilinearStuff%NormalY(k,i,j) = ey/s !east
!		CompGrid%CurvilinearStuff%NormalZ(k,i,j) = zero 		
	END DO
END DO
END IF

IF (CompGrid%Ny>1) THEN
! South boundary
j = 1
DO i = 1+GhostGridX, CompGrid%Nx+GhostGridX ! omit ghost boundary corners
	DO k = 1+GhostGridZ, CompGrid%Nz+GhostGridZ ! omit bottom layer
		nx = CompGrid%CurvilinearStuff%nx(i,j+1)
		ny = CompGrid%CurvilinearStuff%ny(i,j+1)
		s = SQRT(nx**2+ny**2)
		CompGrid%CurvilinearStuff%NormalX(k,i,j) = -nx/s !south
		CompGrid%CurvilinearStuff%NormalY(k,i,j) = -ny/s !south
!		CompGrid%CurvilinearStuff%NormalZ(k,i,j) = zero 		
	END DO
END DO

! North boundary
j = CompGrid%Ny+2*GhostGridY
DO i = 1+GhostGridX, CompGrid%Nx+GhostGridX ! omit ghost boundary corners
	DO k = 1+GhostGridZ, CompGrid%Nz+GhostGridZ ! omit bottom layer
		nx = CompGrid%CurvilinearStuff%nx(i,j-1)
		ny = CompGrid%CurvilinearStuff%ny(i,j-1)
		s = SQRT(nx**2+ny**2)
		CompGrid%CurvilinearStuff%NormalX(k,i,j) = nx/s !north
		CompGrid%CurvilinearStuff%NormalY(k,i,j) = ny/s !north
!		CompGrid%CurvilinearStuff%NormalZ(k,i,j) = zero 		
	END DO
END DO
END IF

! Down boundary
k = 1
DO i = 1+GhostGridX, CompGrid%Nx+GhostGridX ! omit ghost boundary corners
	DO j = 1+GhostGridY, CompGrid%Ny+GhostGridY ! omit ghost boundary corners
		IF (CompGrid%Nx>1 .AND. CompGrid%Ny>1) THEN
			hx = CompGrid%hx(i,j)
			hy = CompGrid%hy(i,j)
			s = SQRT(hx**2+hy**2+one)
			CompGrid%CurvilinearStuff%NormalX(k,i,j) = -hx/s  !down
			CompGrid%CurvilinearStuff%NormalY(k,i,j) = -hy/s  !down
			CompGrid%CurvilinearStuff%NormalZ(k,i,j) = -one/s !down	
		ELSE IF (CompGrid%Nx>1) THEN
			hx = CompGrid%hx(i,j)
			s = SQRT(hx**2+one)
			CompGrid%CurvilinearStuff%NormalX(k,i,j) = -hx/s  !down
			CompGrid%CurvilinearStuff%NormalZ(k,i,j) = -one/s !down	
		ELSE IF (CompGrid%Ny>1) THEN
			hy = CompGrid%hy(i,j)
			s = SQRT(hy**2+one)
			CompGrid%CurvilinearStuff%NormalY(k,i,j) = -hy/s  !down
			CompGrid%CurvilinearStuff%NormalZ(k,i,j) = -one/s !down	
		ENDIF
	END DO
END DO

! APEK: NOT NEEDED
! Up boundary
!k = CompGrid%Nz + GhostGridZ
!DO i = 1+GhostGridX, CompGrid%Nx+GhostGridX ! omit ghost boundary corners
!	DO j = 1+GhostGridY, CompGrid%Ny+GhostGridY ! omit ghost boundary corners
!		etax = CompGrid%Wavefield_FS%Ex(i,j)
!		etay = CompGrid%Wavefield_FS%Ey(i,j)
!		s = SQRT(etax**2+etay**2)
!		CompGrid%CurvilinearStuff%NormalX(k,i,j) = etax/s  !up
!		CompGrid%CurvilinearStuff%NormalY(k,i,j) = etay/s  !up
!		CompGrid%CurvilinearStuff%NormalZ(k,i,j) = one/s   !up
!	END DO
!END DO

END SUBROUTINE ComputeNormalVectors
