SUBROUTINE SetupCompDomain
USE GlobalVariables
IMPLICIT NONE
! By Allan P. Engsig-Karup.
REAL(KIND=long), DIMENSION(:,:), POINTER :: x, y
REAL(KIND=long), DIMENSION(:), POINTER   :: z
INTEGER :: i, j, k, Nx, Ny, Nz

Nx = FineGrid%Nx
Ny = FineGrid%Ny
Nz = FineGrid%Nz

! Allocate space for physical grids and ghost point layers
ALLOCATE(FineGrid%x(Nx+2*GhostGridX,Ny+2*GhostGridY))
ALLOCATE(FineGrid%y(Nx+2*GhostGridX,Ny+2*GhostGridY))
ALLOCATE(FineGrid%z(Nz+GhostGridZ))

! TEMPORARY POINTERS/VARIABLES
x => FineGrid%x
y => FineGrid%y
z => FineGrid%z

! ALLOCATE BATHYMETRY DATA
ALLOCATE(FineGrid%h(Nx+2*GhostGridX,Ny+2*GhostGridY))
!IF (Nx>1) THEN
  ALLOCATE(FineGrid%hx(Nx+2*GhostGridX,Ny+2*GhostGridY))
  ALLOCATE(FineGrid%hxx(Nx+2*GhostGridX,Ny+2*GhostGridY))
!ENDIF
!IF (Ny>1) THEN
  ALLOCATE(FineGrid%hy(Nx+2*GhostGridX,Ny+2*GhostGridY))
  ALLOCATE(FineGrid%hyy(Nx+2*GhostGridX,Ny+2*GhostGridY))
!ENDIF

! Define fine grid
IF (Nx==1) THEN
	x(1,1:Ny+2*GhostGridY) = Lx;
ELSE
	SELECT CASE (GridX)
		CASE (0) ! Even Grid
			DO i = -GhostGridX,Nx-1+GhostGridX
			   x(i+1+GhostGridX,1:Ny+2*GhostGridY) = REAL(i,long)/REAL((Nx-1),long)*Lx;
			END DO
		CASE DEFAULT
			PRINT *,'Error: Specified grid for x-direction invalid.'
			STOP
	END SELECT
ENDIF
IF (Ny==1) THEN
	y(1:Nx+2*GhostGridX,1) = Ly;
ELSE
	SELECT CASE (GridY)
		CASE (0) ! Even Grid
			DO i = -GhostGridY,Ny-1+GhostGridY
 			   y(1:Nx+2*GhostGridX,i+1+GhostGridY) = REAL(i,long)/REAL((Ny-1),long)*Ly;
			END DO
		CASE DEFAULT
			PRINT *,'Error: Specified grid for y-direction invalid.'
			STOP
	END SELECT
END IF
! Define fine grid
SELECT CASE (GridZ)
	CASE (0) ! Even Grid
		DO i = 0,Nz-1
		   z(i+1+GhostGridZ) = REAL(i,long)/REAL((Nz-1),long)*Lz;
		END DO
	CASE (1) ! Uneven Grid
		DO i = 0,Nz-1
			z(i+1+GhostGridZ) = SIN(pi*REAL(i,long)/REAL(2*(Nz-1),long))
		END DO
	CASE DEFAULT
		PRINT *,'Error: Specified grid for z-direction invalid.'
		STOP
END SELECT
IF (GhostGridZ==1) THEN
	! include ghost point below bottom. Reflect point just above bottom symmetrically
	! to ghost point by reflection about z=0
	z(1)=-z(3)
ENDIF
! Allocate array for coefficients involving sigma and its derivatives
! to be used in defining the transformed laplace problem.
ALLOCATE(FineGrid%dsigmanew(Nz+GhostGridZ,Nx+2*GhostGridX,Ny+2*GhostGridY,5))
FineGrid%dsigmanew=zero
END SUBROUTINE SetupCompDomain
