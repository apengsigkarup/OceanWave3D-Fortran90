SUBROUTINE MGRestriction3DExplicit(PHIfine,PHIcoarse,FineGrid,CoarseGrid,GhostGridX,GhostGridY,GhostGridZ,Nzf,Nxf,Nyf,Nzc,Nxc,Nyc)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
REAL(KIND=long), DIMENSION(:,:), POINTER :: xf,yf
REAL(KIND=long), DIMENSION(:), POINTER :: zf
TYPE (Level_def) :: FineGrid, CoarseGrid ! Grid information
INTEGER :: xcoarsen, ycoarsen, zcoarsen, GhostGridX,GhostGridY, GhostGridZ
INTEGER :: Nxf, Nyf, Nzf, Nxc, Nyc, Nzc, N, NS
REAL(KIND=long), DIMENSION(Nzf,Nxf,Nyf) :: PHIfine, tmp3D
REAL(KIND=long), DIMENSION(Nzc,Nxc,Nyc) :: PHIcoarse
INTEGER :: kc,jc,ic
INTEGER :: ifine, jfine, kfine
REAL(KIND=long) :: dx_W, dx_E, dy_S, dy_N, dz_D, dz_U, dxc, dyc, dzc

! Initialize output array
PHIcoarse = zero
tmp3D     = PHIfine

! USE POINTERS TO ALREADY ALLOCATED GRID VECTORS
xf => FineGrid%x;   yf => FineGrid%y;   zf => FineGrid%z

! coarsenig strategi is based on the physical grid sizes
Nxf = Nxf-2*GhostGridX
Nxc = Nxc-2*GhostGridX
Nyf = Nyf-2*GhostGridY
Nyc = Nyc-2*GhostGridY
Nzf = Nzf-GhostGridZ
Nzc = Nzc-GhostGridZ

! Determine coarsening strategy for each dimension
xcoarsen = 0; ycoarsen = 0; zcoarsen = 0
IF (Nxf>Nxc) THEN; xcoarsen = 1; ENDIF
IF (Nyf>Nyc) THEN; ycoarsen = 1; ENDIF
IF (Nzf>Nzc) THEN; zcoarsen = 1; ENDIF

! *******************************************
!			 X-COARSEN
! *******************************************

IF (xcoarsen==1) THEN
!  print*,'X-coarsen'
    IF (GhostGridX==1) THEN
        ! We have Ghost Layers present
        ! OMIT Ghost Point Layer's from restriction

	! interior points, full stencil
	DO ifine = 1+GhostGridX , Nxf+GhostGridX, 2
		dx_W = xf(ifine,1+GhostGridY) - xf(ifine-1,1+GhostGridY) ! FIXME: y-index
		dx_E = xf(ifine+1,1+GhostGridY) - xf(ifine,1+GhostGridY) ! FIXME: y-index
!        print*,'dx_W=',dx_W
!        print*,'dx_E=',dx_E
		dxc  = one/(dx_W + dx_E)
		DO kfine = 1+GhostGridZ , Nzf+GhostGridZ-1 ! FREE SURFACE LAYER NOT INCLUDED SINCE THERE WE USE DIRECT INJECTION
			DO jfine = 1+GhostGridY , Nyf+GhostGridY
				tmp3D(kfine,ifine,jfine) = (half*dxc)*(dx_W*tmp3D(kfine,ifine-1,jfine) + (dx_W+dx_E)*tmp3D(kfine,ifine,jfine) + &
				   dx_E*tmp3D(kfine,ifine+1,jfine))
			END DO
		END DO
	END DO
    ELSE
       WRITE (*,'(A)') 'Error: currently multigrid has not been implemented for grid without a single ghost layer.'
       STOP ! FORCED STOP
    ENDIF
ENDIF

! *******************************************
!			 Y-COARSEN
! *******************************************

IF (ycoarsen==1) THEN
!  print*,'Y-coarsen'
    IF (GhostGridY==1) THEN
        ! We have Ghost Layers present
        ! OMIT Ghost Point Layer's from restriction

	! interior points, full stencil
	DO jfine = 1+GhostGridY , Nyf+GhostGridY, 2
        dy_S = yf(1+GhostGridX,jfine)   - yf(1+GhostGridX,jfine-1) ! FIXME: x-index
        dy_N = yf(1+GhostGridX,jfine+1) - yf(1+GhostGridX,jfine) ! FIXME: x-index
        dyc  = one/(dy_S + dy_N)
		DO kfine = 1+GhostGridZ , Nzf+GhostGridZ-1 ! FREE SURFACE LAYER NOT INCLUDED SINCE THERE WE USE DIRECT INJECTION
			DO ifine = 1+GhostGridX , Nxf+GhostGridX, 1+xcoarsen
				tmp3D(kfine,ifine,jfine) = (half*dyc)*(dy_S*tmp3D(kfine,ifine,jfine-1) + (dy_N+dy_S)*tmp3D(kfine,ifine,jfine) + &
				   dy_N*tmp3D(kfine,ifine,jfine+1))
			END DO
		END DO
	END DO
    ELSE
       WRITE (*,'(A)') 'Error: currently multigrid has not been implemented for grid without a single ghost layer.'
       STOP ! FORCED STOP
    ENDIF
ENDIF

! *******************************************
!			 Z-COARSEN
! *******************************************

IF (zcoarsen==1) THEN
!  print*,'Z-coarsen'
    IF (GhostGridZ==1) THEN
        ! We have Ghost Layers present
        ! OMIT Ghost Point Layer's from restriction

	! interior points, full stencil
	DO kfine = 1+GhostGridZ , Nzf+GhostGridZ-1, 2 ! omit free surface since there we use direct injection
        dz_U = zf(kfine+1) - zf(kfine)
        dz_D = zf(kfine)   - zf(kfine-1)
        dzc  = one/(dz_D + dz_U)
	    DO ifine = 1+GhostGridX , Nxf+GhostGridX, 1 + xcoarsen
			DO jfine = 1+GhostGridY , Nyf+GhostGridY, 1 + ycoarsen
				tmp3D(kfine,ifine,jfine) = (half*dzc)*(dz_D*tmp3D(kfine-1,ifine,jfine) + dz_U*tmp3D(kfine+1,ifine,jfine)) + &
				   half*tmp3D(kfine,ifine,jfine)
			END DO
		END DO
	END DO
    ELSE
       WRITE (*,'(A)') 'Error: currently multigrid has not been implemented for grid without a single ghost layer.'
       STOP ! FORCED STOP
    ENDIF

ENDIF

! DIRECT INJECTION
PHIcoarse(1+GhostGridZ:Nzc+GhostGridZ,1+GhostGridX:Nxc+GhostGridX,1+GhostGridY:Nyc+GhostGridY) = &
	tmp3D(1+GhostGridZ:Nzf+GhostGridZ:1+zcoarsen,1+GhostGridX:Nxf+GhostGridX:1+xcoarsen,1+GhostGridY:Nyf+GhostGridY:1+ycoarsen)

if (0==1) then ! if active we assume that the defect is a result of imposing boundary conditions... does not seem to make sense?? it should be assumed zero at boundary equations??
! APPROXIMATE GHOST POINT VALUES ON COARSE GRID
! This is needed for grids containing ghost points. Wall boundary conditions assumed.
DO kc = 1+GhostGridZ, Nzc+GhostGridZ
    DO jc = 1+GhostGridY, Nyc+GhostGridY
       ic = 1 ! WEST
       PHIcoarse(kc,ic,jc) = PHIcoarse(kc,ic+2,jc)
       ic = Nxc+GhostGridX*2 ! EAST
       PHIcoarse(kc,ic,jc) = PHIcoarse(kc,ic-2,jc)
  END DO
END DO

DO kc = 1+GhostGridZ, Nzc+GhostGridZ
  DO ic = 1+GhostGridX, Nxc+GhostGridX
       jc = 1 ! SOUTH
       PHIcoarse(kc,ic,jc) = PHIcoarse(kc,ic,jc+2)
       jc = Nxc+GhostGridY*2 ! NORTH
       PHIcoarse(kc,ic,jc) = PHIcoarse(kc,ic,jc-2)
  END DO
END DO

kc = 1
DO ic = 1+GhostGridX, Nxc+GhostGridX
  DO jc = 1+GhostGridY, Nyc+GhostGridY
     PHIcoarse(kc,ic,jc) = PHIcoarse(kc+2,ic,jc)
  END DO
END DO
end if
END SUBROUTINE MGRestriction3DExplicit
