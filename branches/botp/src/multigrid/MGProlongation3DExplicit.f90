SUBROUTINE MGProlongation3DExplicit(PHIfine,PHIcoarse,FineGrid,CoarseGrid,GhostGridX,GhostGridY,GhostGridZ,Nzf,Nxf,Nyf,Nzc,Nxc,Nyc)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
REAL(KIND=long), DIMENSION(:,:), POINTER :: xf, yf, xc, yc
REAL(KIND=long), DIMENSION(:), POINTER :: zf, zc
TYPE (Level_def) :: FineGrid, CoarseGrid ! Grid information
INTEGER :: xcoarsen, ycoarsen, zcoarsen, GhostGridX, GhostGridY, GhostGridZ
INTEGER :: Nxf, Nyf, Nzf, Nxc, Nyc, Nzc, N, NS
REAL(KIND=long), DIMENSION(Nzf,Nxf,Nyf) :: PHIfine
REAL(KIND=long), DIMENSION(Nzc,Nxc,Nyc) :: PHIcoarse
INTEGER :: kc,jc,ic, XFAC, YFAC, ZFAC
INTEGER :: ifine, jfine, kfine
REAL(KIND=long) :: dx_W, dx_E, dy_S, dy_N, dz_D, dz_U, dxf, dyf, dzf, dzc

! Initialize output array
PHIfine = zero

! USE POINTERS TO ALREADY ALLOCATED GRID VECTORS
xf => FineGrid%x;   yf => FineGrid%y;   zf => FineGrid%z
xc => CoarseGrid%x; yc => CoarseGrid%y; zc => CoarseGrid%z

! coarsenig strategy is based on the physical grid sizes
Nxf = Nxf-2*GhostGridX
Nxc = Nxc-2*GhostGridX
Nyf = Nyf-2*GhostGridY
Nyc = Nyc-2*GhostGridY
Nzf = Nzf-GhostGridZ
Nzc = Nzc-GhostGridZ

!print*,'Nxf=',Nxf
!print*,'Nxc=',Nxc
!print*,'Nyf=',Nyf
!print*,'Nyc=',Nyc
!print*,'Nzf=',Nzf
!print*,'Nzc=',Nzc
!STOP

! Determine coarsening strategy for each dimension
xcoarsen = 0; ycoarsen = 0; zcoarsen = 0;
XFAC = 1; YFAC = 1; ZFAC = 1;
IF (Nxf>Nxc) THEN; xcoarsen = 1; XFAC=2; ENDIF
IF (Nyf>Nyc) THEN; ycoarsen = 1; YFAC=2; ENDIF
IF (Nzf>Nzc) THEN; zcoarsen = 1; ZFAC=2; ENDIF

! BUILD PROLONGATION OPERATOR
! Carry out x-, y-, and z-sweeps independently by first
! carrying out all direct injections as this will allow
! us to only work on the output array when sweeping in
! each direction.

! DIRECT INJECTION (OMIT FREE SURFACE LAYER)
!print*,'GhostGridX=',GhostGridX
!print*,'GhostGridY=',GhostGridY
!print*,'GhostGridZ=',GhostGridZ
!STOP
PHIfine(1+GhostGridZ:Nzf+GhostGridZ-1:ZFAC,1+GhostGridX:Nxf+GhostGridX:XFAC,1+GhostGridY:Nyf+GhostGridY:YFAC) = &
	PHIcoarse(1+GhostGridZ:Nzc+GhostGridZ-1,1+GhostGridX:Nxc+GhostGridX,1+GhostGridY:Nyc+GhostGridY)

! X-SWEEP
! DO LOOPS HAS BEEN ORDERED FOR EFFICIENCY (DUE TO dx's)
IF (XFAC==2) THEN
  DO jfine = 1+GhostGridY, Nyf+GhostGridY, YFAC
    DO ifine = 2+GhostGridX, Nxf+GhostGridX, XFAC
      dx_W = xf(ifine,jfine)   - xf(ifine-1,jfine)
      dx_E = xf(ifine+1,jfine) - xf(ifine,jfine)
	  dxf  = one/(dx_W + dx_E)
	  DO kfine = 1+GhostGridZ, Nzf+GhostGridZ-1, ZFAC ! FIXME: WE DO NOT PROLONGATE AT THE FREE SURFACE OR AT GHOST LAYER
		! linear interpolation
   		PHIfine(kfine,ifine,jfine) = (dx_E*PHIfine(kfine,ifine-1,jfine) + dx_W*PHIfine(kfine,ifine+1,jfine))*dxf
	  END DO
	END DO
  END DO
END IF
XFAC = 1

! Y-SWEEP
IF (YFAC==2) THEN
  DO ifine = 1+GhostGridX, Nxf+GhostGridX, XFAC
    DO jfine = 2+GhostGridY, Nyf+GhostGridY, YFAC
      dy_S = yf(ifine,jfine)   - yf(ifine,jfine-1)
      dy_N = yf(ifine,jfine+1) - yf(ifine,jfine)
	  dyf  = one/(dy_S + dy_N)
	  DO kfine = 1+GhostGridZ, Nzf+GhostGridZ-1, ZFAC ! FIXME: WE PROLONGATE AT THE FREE SURFACE
	    ! linear interpolation
	    PHIfine(kfine,ifine,jfine) = (dy_N*PHIfine(kfine,ifine,jfine-1) + dy_S*PHIfine(kfine,ifine,jfine+1))*dyf
	  END DO
	END DO
  END DO
END IF
YFAC = 1

! Z-SWEEP
! DO LOOPS HAS BEEN ORDERED FOR EFFICIENCY (DUE TO dz's)
IF (ZFAC==2) THEN
	DO kfine = 2+GhostGridZ, Nzf+GhostGridZ-1, ZFAC ! FIXME: WE DO NOT PROLONGATE AT THE FREE SURFACE
	    dz_D = zf(kfine)   - zf(kfine-1)
	    dz_U = zf(kfine+1) - zf(kfine)
	    dzf  = one/(dz_D + dz_U)
		DO jfine = 1+GhostGridY, Nyf+GhostGridY, YFAC
 			DO ifine = 1+GhostGridX, Nxf+GhostGridX, XFAC
			! linear interpolation
			PHIfine(kfine,ifine,jfine) = (dz_U*PHIfine(kfine-1,ifine,jfine) + dz_D*PHIfine(kfine+1,ifine,jfine))*dzf
		END DO
	END DO
  END DO
END IF

!IF (GhostGridZ==1) THEN
!    ! Approximate symmetric reflection about bottom
!	PHIfine(1,:,:) = PHIfine(3,:,:)
!ENDIF

if (0==1) then
! APPROXIMATE GHOST POINT VALUES ON FINE GRID??
! IMPOSE WALL BOUNDARY CONDITIONS THROUGH SYMMETRY (FIXME: needs to be fixed for curvilinear settings /APEK)
DO kfine = 1+GhostGridZ, Nzf+GhostGridZ, ZFAC
    DO jfine = 1+GhostGridY, Nyf+GhostGridY, YFAC
       ifine = 1 ! WEST
       PHIfine(kfine,ifine,jfine) = PHIfine(kfine,ifine+2,jfine)
       ifine = Nxf+GhostGridX*2 ! EAST
       PHIfine(kfine,ifine,jfine) = PHIfine(kfine,ifine-2,jfine)
  END DO
END DO

DO kfine = 1+GhostGridZ, Nzf+GhostGridZ, ZFAC
  DO ifine = 1+GhostGridX, Nxf+GhostGridX, XFAC
       jfine = 1 ! SOUTH
       PHIfine(kfine,ifine,jfine) = PHIfine(kfine,ifine,jfine+2)
       jfine = Nxf+GhostGridY*2 ! NORTH
       PHIfine(kfine,ifine,jfine) = PHIfine(kfine,ifine,jfine-2)
  END DO
END DO

kfine = 1
DO ifine = 1+GhostGridX, Nxf+GhostGridX, XFAC
  DO jfine = 1+GhostGridY, Nyf+GhostGridY, YFAC
     PHIfine(kfine,ifine,jfine) = PHIfine(kfine+2,ifine,jfine)
  END DO
END DO
end if
END SUBROUTINE MGProlongation3DExplicit
