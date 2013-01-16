SUBROUTINE SetupCompDomain_test
! By Allan P. Engsig-Karup.
USE GlobalVariables
IMPLICIT NONE

REAL(KIND=long), DIMENSION(:,:), POINTER :: x, y
REAL(KIND=long), DIMENSION(:), POINTER   :: z
INTEGER :: Nx, Ny, Nz
! Local variable...
REAL(KIND=long), DIMENSION(:), POINTER :: R, angle ! Radius and angle...
REAL(KIND=long), DIMENSION(:,:), POINTER :: e, n
INTEGER :: i, j, k, Nxstart, Nxcircle

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
IF (Nx>1) THEN
  ALLOCATE(FineGrid%hx(Nx+2*GhostGridX,Ny+2*GhostGridY))
  ALLOCATE(FineGrid%hxx(Nx+2*GhostGridX,Ny+2*GhostGridY))
ENDIF
IF (Ny>1) THEN
  ALLOCATE(FineGrid%hy(Nx+2*GhostGridX,Ny+2*GhostGridY))
  ALLOCATE(FineGrid%hyy(Nx+2*GhostGridX,Ny+2*GhostGridY))
ENDIF

!$$$$$$ ! Semi-circular channel... (with a zone to generate and a zone to absorb...)
!$$$$$$ ! FIXME it is probably not the correct way to do this... (define grid as before, do this in the initialization...)
!$$$$$$ IF ((Nx>1).AND.(Ny>1)) THEN
!$$$$$$
!$$$$$$     ! Index of the beginning of the circle
!$$$$$$     Nxstart = Nx/4
!$$$$$$     ! Index of the end of the circle
!$$$$$$     !Nxstop = NINT(2*Nx/3)
!$$$$$$     ! Number of points in the circle part
!$$$$$$     Nxcircle = Nx/2
!$$$$$$
!$$$$$$     ! Allocate the temporary radius and angle...
!$$$$$$     ALLOCATE(angle(1:Nxcircle))
!$$$$$$     ALLOCATE(R(-GhostGridY:Ny-1+GhostGridY))
!$$$$$$     !
!$$$$$$     DO i = 1,Nxcircle
!$$$$$$         angle(i) = PI-(REAL(i-1,long)/REAL((Nxcircle-1),long)*PI)
!$$$$$$     ENDDO
!$$$$$$     DO j = -GhostGridY,Ny-1+GhostGridY
!$$$$$$         R(j) = Lx + REAL(j,long)/REAL((Ny-1),long)*(Ly)
!$$$$$$     END DO
!$$$$$$     ! Define fine grid
!$$$$$$     DO i=1,Nxcircle
!$$$$$$         DO j=-GhostGridY,Ny-1+GhostGridY
!$$$$$$             x(i+Nxstart-1,j+1+GhostGridY) = R(j)*COS(angle(i))
!$$$$$$             y(i+Nxstart-1,j+1+GhostGridY) = R(j)*SIN(angle(i))
!$$$$$$         ENDDO
!$$$$$$     ENDDO
!$$$$$$     !
!$$$$$$     DEALLOCATE(angle, R)
!$$$$$$     !
!$$$$$$     DO i=1,Nxstart-1
!$$$$$$       DO j=1,Ny+2*GhostGridY
!$$$$$$         x(i,j)=x(Nxstart,j)
!$$$$$$         y(i,j)=y(Nxstart,j)-REAL(Nxstart-i,long)/REAL(Nxstart-1-GhostGridX,long)*2*SFsol%L!4*Lx !FIXME approx...
!$$$$$$       ENDDO
!$$$$$$     ENDDO
!$$$$$$     !
!$$$$$$     DO i=Nxstart+Nxcircle,Nx+2*GhostGridX
!$$$$$$       DO j=1,Ny+2*GhostGridY
!$$$$$$         x(i,j)= x(Nxstart+Nxcircle-1,j)
!$$$$$$         y(i,j)= y(Nxstart+Nxcircle-1,j)-REAL(i-(Nxstart+Nxcircle)+1,long)/&
!$$$$$$             REAL(Nx+2*GhostGridX-(Nxstart+Nxcircle)+1-GhostGridX,long)*2*SFsol%L!4*Lx !FIXME approx...*4*Lx
!$$$$$$       ENDDO
!$$$$$$     ENDDO
!$$$$$$ ELSE
!$$$$$$   PRINT *,'Error this is not done...'
!$$$$$$   STOP
!$$$$$$ ENDIF
! Cylinder case...
! Define the radius and angle
IF ((Nx>1).AND.(Ny>1)) THEN
    ! Allocate the temporary radius and angle...
    ALLOCATE(angle(-GhostGridX:Nx-1+GhostGridX))
    ALLOCATE(R(-GhostGridY:Ny-1+GhostGridY))
    !
    DO i = -GhostGridX,Nx-1+GhostGridX
        angle(i) = PI-(REAL(i,long)/REAL((Nx-1),long)*PI) !PI-0.5d0/REAL((Nx-1),long)*PI-(REAL(i,long)/REAL((Nx-1),long)*PI) !
!$$$$$$         print*,'Initial',angle(i)
!$$$$$$         angle(i) = (PI+(COS(PI*REAL(i,long)/REAL((Nx-1),long)))*PI)/two
!$$$$$$         angle(i) = (PI+TANH(1.5_long*REAL(Nx-1-2*i,long)/REAL(Nx-1))/TANH(1.5_long)*PI)/two
    ENDDO
    i=0
    DO j = -GhostGridY,Ny-1+GhostGridY
      i=i+1
        R(j) = Lx + REAL(j,long)/REAL((Ny-1),long)*(Ly)
        !angle(j) = REAL(j,long)/REAL((Ny-1),long)*PI
        ! cosine grid (seems too strng)
        !R(j) = Lx + (COS(pi*REAL(j,long)/REAL(2*(Ny-1),long)))*(Ly)
        ! stretch grid
!$$$$$$         R(Ny-1+GhostGridY+1-i) = Lx + (1.0_long-TANH(1.3_long*REAL(j,long)/REAL(Ny-1))/TANH(1.3_long))*(Ly)
    END DO
    ! Define symetrical points...
    angle(-1)=angle(0)-(angle(1)-angle(0))
    angle(Nx-1+GhostGridX)=angle(Nx-1)-(angle(Nx-2)-angle(Nx-1))
    !
    R(-1)=R(0)-(R(1)-R(0))
    R(Ny-1+GhostGridY)=R(Ny-1)-(R(Ny-2)-R(Ny-1))
    ! Define fine grid
    DO i=-GhostGridX,Nx-1+GhostGridX
        DO j=-GhostGridY,Ny-1+GhostGridY
            x(i+1+GhostGridX,j+1+GhostGridY) = R(j)*COS(angle(i))
            y(i+1+GhostGridX,j+1+GhostGridY) = R(j)*SIN(angle(i))
        ENDDO
    ENDDO
    ! Prevent from symetrical ghost points...
    !x(1,:)=x(2,:)
    !y(1,:)=y(2,:)-(y(3,2)-y(2,2))
    !x(Nx+2*GhostGridX,:)=x(Nx+GhostGridX,:)
    !y(Nx+2*GhostGridX,:)=y(Nx+GhostGridX,:)-(y(Nx,2)-y(Nx+GhostGridX,2))
    !
    DEALLOCATE(angle, R)
ELSE
  PRINT *,'Error this is not done...'
  STOP
ENDIF
! Grid with straight boundaries
!$$$$$$ ! Define fine grid
!$$$$$$ IF (Nx==1) THEN
!$$$$$$     x(1,1:Ny+2*GhostGridY) = Lx;
!$$$$$$ ELSE
!$$$$$$     SELECT CASE (GridX)
!$$$$$$         CASE (0) ! Even Grid
!$$$$$$             DO i = -GhostGridX,Nx-1+GhostGridX
!$$$$$$                x(i+1+GhostGridX,1:Ny+2*GhostGridY) = REAL(i,long)/REAL((Nx-1),long)*Lx;
!$$$$$$             END DO
!$$$$$$         CASE DEFAULT
!$$$$$$             PRINT *,'Error: Specified grid for x-direction invalid.'
!$$$$$$             STOP
!$$$$$$     END SELECT
!$$$$$$ ENDIF
!$$$$$$ IF (Ny==1) THEN
!$$$$$$     y(1:Nx+2*GhostGridX,1) = Ly;
!$$$$$$ ELSE
!$$$$$$     SELECT CASE (GridY)
!$$$$$$         CASE (0) ! Even Grid
!$$$$$$             DO i = -GhostGridY,Ny-1+GhostGridY
!$$$$$$                y(1:Nx+2*GhostGridX,i+1+GhostGridY) = REAL(i,long)/REAL((Ny-1),long)*Ly;
!$$$$$$             END DO
!$$$$$$         CASE DEFAULT
!$$$$$$             PRINT *,'Error: Specified grid for y-direction invalid.'
!$$$$$$             STOP
!$$$$$$     END SELECT
!$$$$$$ END IF
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
!ALLOCATE(FineGrid%dsigma((Nx+2*GhostGridX)*(Ny+2*GhostGridY)*(Nz+GhostGridZ),5))
!FineGrid%dsigma=zero
END SUBROUTINE SetupCompDomain_test
