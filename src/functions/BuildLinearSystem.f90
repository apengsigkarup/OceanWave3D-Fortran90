SUBROUTINE BuildLinearSystem(Nx,Ny,Nz,PHI,output,GridStruct,alpha,beta,gamma)
! By Allan P. Engsig-Karup.
!
! Form: output = A PHI
!
! where A is the system matrix and phi the current guess
!
USE Precision
USE Constants
USE DataTypes
USE GlobalVariables, ONLY: GhostGridX, GhostGridY, GhostGridZ, &
	swenseONOFF, Wavefield, IncWaveType ! FIXME: Wavefield is a problem here?? due to Multigrid algorithm
IMPLICIT NONE
TYPE (Level_def), INTENT(IN) :: GridStruct
INTEGER::Gidx, Gidx2, Nx, Ny, Nz, i, j, k, alpha, beta, gamma
REAL(KIND=long), DIMENSION(Nx*Ny*Nz), INTENT(OUT) :: output
REAL(KIND=long), DIMENSION(Nx*Ny*Nz) :: dpdy, ddpdyy, dpdx, ddpdxx, dpds, ddpdsdx, ddpdsdy, ddpdss, tmp !GD: addition for cross derivatives
REAL(KIND=long), DIMENSION(Nx*Ny*Nz), INTENT(IN) :: PHI
!
INTEGER :: Gidx3 !Index to look for boundaries

!print*,'PHI=',PHI
!read*

Gidx3  = 1
Gidx   = 0
output = zero

! Determine derivatives explicitly
IF (Nx>1) THEN
	CALL DiffXEven(PHI,dpdx,1,Nx,Ny,Nz,GridStruct%DiffStencils,alpha)
	CALL DiffXEven(PHI,ddpdxx,2,Nx,Ny,Nz,GridStruct%DiffStencils,alpha)
    ! DISCRETIZATION CHOICE MADE HERE FOR MIXED DERIVATIVE
	!CALL DiffZArbitrary(dpdx,ddpdsdx,1,Nx,Ny,Nz,GridStruct%DiffStencils,gamma)
    !GD: modification for correct treatment of cross derivatives...
    CALL DiffXEven_CD(PHI,tmp,GridStruct%DiffStencils%IndexesX_XZorXY(:,:,:,:,1),&
    	GridStruct%DiffStencils%StencilsX_XZorXY(:,:,:,:,1),Nx,Ny,Nz,alpha)
    CALL DiffZArbitrary_CD(tmp,ddpdsdx,GridStruct%DiffStencils%IndexesZ_XZorYZ(:,:,:,:,1),&
    	GridStruct%DiffStencils%StencilsZ_XZorYZ(:,:,:,:,1),Nx,Ny,Nz,gamma)
ENDIF

IF (Ny>1) THEN
	CALL DiffYEven(PHI,dpdy,1,Nx,Ny,Nz,GridStruct%DiffStencils,beta)
	CALL DiffYEven(PHI,ddpdyy,2,Nx,Ny,Nz,GridStruct%DiffStencils,beta)
    ! DISCRETIZATION CHOICE MADE HERE FOR MIXED DERIVATIVE
	!CALL DiffZArbitrary(dpdy,ddpdsdy,1,Nx,Ny,Nz,GridStruct%DiffStencils,gamma)
    !GD: modification for correct treatment of cross derivatives...
    CALL DiffYEven_CD(PHI,tmp,GridStruct%DiffStencils%IndexesY_YZorXY(:,:,:,:,1),&
    	GridStruct%DiffStencils%StencilsY_YZorXY(:,:,:,:,1),Nx,Ny,Nz,beta)
    CALL DiffZArbitrary_CD(tmp,ddpdsdy,GridStruct%DiffStencils%IndexesZ_XZorYZ(:,:,:,:,2),&
    	GridStruct%DiffStencils%StencilsZ_XZorYZ(:,:,:,:,2),Nx,Ny,Nz,gamma)
ENDIF

CALL DiffZArbitrary(PHI,dpds,1,Nx,Ny,Nz,GridStruct%DiffStencils,gamma)
CALL DiffZArbitrary(PHI,ddpdss,2,Nx,Ny,Nz,GridStruct%DiffStencils,gamma)

! DISTINGUISH BETWEEN 2D AND 3D PROBLEMS BELOW (FOR EFFICIENCY)
! FIXME: REMOVE INDEX GIDX? FROM LOOPS BY USING ARRAYS
IF (Nx==1) THEN
  i=1
  DO j = 1+GhostGridY, Ny-GhostGridY
	DO k = 1, Nz
!		Gidx = Gidx + 1
		Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz
       IF (k==Nz) THEN
			! FREE SURFACE
			output(Gidx) = PHI(Gidx)
		ELSE IF (k==1) THEN
			! BOTTOM, DIRECT KIN. BOTTOM CONDITION
			output(Gidx) = (GridStruct%dsigmanew(k+GhostGridZ,i,j,5) + &
			   GridStruct%hy(i,j)*GridStruct%dsigmanew(k+GhostGridZ,i,j,4) )*dpds(Gidx+GhostGridZ)+&
			   GridStruct%hy(i,j)*dpdy(Gidx+GhostGridZ)
            ! GD: SWENSE addition comes here
            IF (swenseONOFF/=0) THEN
				output(Gidx) = output(Gidx) + (Wavefield%Pz_I_bp(Gidx3)+GridStruct%hy(i,j)*Wavefield%Py_I_bp(Gidx3))
            	Gidx3 = Gidx3+1
            ENDIF
		ELSE
			! INTERIOR POINTS
			output(Gidx) = ddpdyy(Gidx) + GridStruct%dsigmanew(k,i,j,3)*dpds(Gidx) + &
			   two*( GridStruct%dsigmanew(k,i,j,4)*ddpdsdy(Gidx) ) + &
			   (GridStruct%dsigmanew(k,i,j,4)**2+GridStruct%dsigmanew(k,i,j,5)**2)*ddpdss(Gidx)
		END IF
	END DO
  END DO
ELSE IF (Ny==1) THEN
  j=1
  DO i = 1+GhostGridX, Nx-GhostGridX
	DO k = 1, Nz
!		Gidx = Gidx + 1
		Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz
        IF (k==Nz) THEN
			! FREE SURFACE
			output(Gidx) = PHI(Gidx)
		ELSE IF (k==1) THEN
			! BOTTOM, DIRECT KIN. BOTTOM CONDITION
			output(Gidx) = GridStruct%hx(i,j)*dpdx(Gidx+GhostGridZ) + (GridStruct%dsigmanew(k+GhostGridZ,i,j,5) + &
			   GridStruct%hx(i,j)*GridStruct%dsigmanew(k+GhostGridZ,i,j,2) )*dpds(Gidx+GhostGridZ)
            ! GD: SWENSE addition comes here
            IF (swenseONOFF/=0) THEN
				output(Gidx) = output(Gidx) + (Wavefield%Pz_I_bp(Gidx3)+GridStruct%hx(i,j)*Wavefield%Px_I_bp(Gidx3))
            	Gidx3 = Gidx3+1
            ENDIF
		ELSE
			! INTERIOR POINTS
			output(Gidx) = ddpdxx(Gidx) + two*GridStruct%dsigmanew(k,i,j,2)*ddpdsdx(Gidx) + &
			   (GridStruct%dsigmanew(k,i,j,2)**2+GridStruct%dsigmanew(k,i,j,5)**2)*ddpdss(Gidx) + &
               GridStruct%dsigmanew(k,i,j,3)*dpds(Gidx)
		END IF
	END DO
  END DO
ELSE
  DO j = 1+GhostGridY, Ny-GhostGridY
	DO i = 1+GhostGridX, Nx-GhostGridX
		DO k = 1, Nz
!			Gidx = Gidx + 1
			Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz
            IF (k==Nz) THEN
				! FREE SURFACE
				output(Gidx) = PHI(Gidx)
			ELSE IF (k==1) THEN
				! BOTTOM, DIRECT KIN. BOTTOM CONDITION
				output(Gidx) = (GridStruct%dsigmanew(k+GhostGridZ,i,j,5) + &
				   GridStruct%hx(i,j)*GridStruct%dsigmanew(k+GhostGridZ,i,j,2) + &
				   GridStruct%hy(i,j)*GridStruct%dsigmanew(k+GhostGridZ,i,j,4) )*dpds(Gidx+GhostGridZ)+&
				   GridStruct%hx(i,j)*dpdx(Gidx+GhostGridZ)+GridStruct%hy(i,j)*dpdy(Gidx+GhostGridZ)
                ! GD: SWENSE addition comes here
                IF (swenseONOFF/=0) THEN
                    output(Gidx) = output(Gidx) + (Wavefield%Pz_I_bp(Gidx3)+GridStruct%hx(i,j)*Wavefield%Px_I_bp(Gidx3) &
                    	+GridStruct%hy(i,j)*Wavefield%Py_I_bp(Gidx3))
                    Gidx3 = Gidx3+1
                ENDIF
			ELSE
				! INTERIOR POINTS
				output(Gidx) = ddpdxx(Gidx) + ddpdyy(Gidx) + GridStruct%dsigmanew(k,i,j,3)*dpds(Gidx) + &
				   two*(GridStruct%dsigmanew(k,i,j,2)*ddpdsdx(Gidx) + GridStruct%dsigmanew(k,i,j,4)*ddpdsdy(Gidx) ) + &
				   (GridStruct%dsigmanew(k,i,j,2)**2+GridStruct%dsigmanew(k,i,j,4)**2+GridStruct%dsigmanew(k,i,j,5)**2) &
                   *ddpdss(Gidx)
			END IF
		END DO
	END DO
  END DO
ENDIF
  
! NOTE: IF waves are generated by an inhogogenious Neumann boundary condition
! the array Uneumann below is ~=0. The array is updated from
! waveGenerationFromPaddleSignal.f90 called from Runge_Kutta_4. 
 ! 
  ! West  boundary
  IF (GhostGridX==1) THEN
    i = 1
    DO j = 1+GhostGridY, Ny-GhostGridY
  	  DO k = 1, Nz
      !GD: no Neumann condition on bottom points (in horizontal plane) !DO k = 1+GhostGridZ, Nz
      !But keep this loop from 1 to Nz for SWENSE (use of Gidx3...)
		Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz
		Gidx2  = k + (i-1+1)*Nz + (j-1)*Nx*Nz
		output(Gidx) = dpdx(Gidx2) 
        ! GD: SWENSE addition comes here
        IF (swenseONOFF/=0) THEN
            output(Gidx) = output(Gidx) + Wavefield%Px_I_bp(Gidx3)
          	Gidx3 = Gidx3+1
        ENDIF
	  END DO
    END DO
    ! East  boundary
    i = Nx
    DO j = 1+GhostGridY, Ny-GhostGridY
  	  DO k = 1, Nz
      !GD: no Neumann condition on bottom points (in horizontal plane) !DO k = 1+GhostGridZ, Nz
      !But keep this loop from 1 to Nz for SWENSE (use of Gidx3...)
		Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
		Gidx2 = k + (i-1-1)*Nz + (j-1)*Nx*Nz
		output(Gidx) = dpdx(Gidx2)
        ! GD: SWENSE addition comes here
        IF (swenseONOFF/=0) THEN
            output(Gidx) = output(Gidx) + Wavefield%Px_I_bp(Gidx3)
            Gidx3 = Gidx3+1
        ENDIF
	  END DO
    END DO
  ENDIF
  IF (GhostGridY==1) THEN
    ! South boundary
    j = 1
    DO i = 1+GhostGridX, Nx-GhostGridX
	  DO k = 1, Nz
      !GD: no Neumann condition on bottom points (in horizontal plane) !DO k = 1+GhostGridZ, Nz
      !But keep this loop from 1 to Nz for SWENSE (use of Gidx3...)
		Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz
		Gidx2 = k + (i-1)*Nz + (j-1+1)*Nx*Nz
		output(Gidx) = dpdy(Gidx2)
        ! GD: SWENSE addition comes here
        IF (swenseONOFF/=0) THEN
            output(Gidx) = output(Gidx) + Wavefield%Py_I_bp(Gidx3)
            Gidx3 = Gidx3+1
        ENDIF
	  END DO
    END DO
    ! North boundary
    j = Ny
    DO i = 1+GhostGridX, Nx-GhostGridX
	  DO k = 1, Nz
      !GD: no Neumann condition on bottom points (in horizontal plane) !DO k = 1+GhostGridZ, Nz
      !But keep this loop from 1 to Nz for SWENSE (use of Gidx3...)
		Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz
		Gidx2 = k + (i-1)*Nz + (j-1-1)*Nx*Nz
		output(Gidx) = dpdy(Gidx2)
        ! GD: SWENSE addition comes here
        IF (swenseONOFF/=0) THEN
            output(Gidx) = output(Gidx) + Wavefield%Py_I_bp(Gidx3)
            Gidx3 = Gidx3+1
        ENDIF
	  END DO
    END DO
  ENDIF
  IF (GhostGridX+GhostGridY==2) THEN
    ! Four corner type points in domain (ghost points)
    i = 1; j = 1;
    DO k = 1, Nz
        ! FIXME: change to form below: Gidx = Gidx + 1
		Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
		output(Gidx)  = PHI(Gidx)
    END DO
    i = Nx; j = 1;
    DO k = 1, Nz
		Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
		output(Gidx)  = PHI(Gidx)
    END DO
    i = 1; j = Ny;
    DO k = 1, Nz
		Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
		output(Gidx)  = PHI(Gidx)
	END DO
	i = Nx; j = Ny;
	DO k = 1, Nz
		Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
		output(Gidx)  = PHI(Gidx)
	END DO
    ! GD: FIXME the ghost corners on the bottom not treated ?
    ! FIXME: 2D case have to be treated also... (remove the two corner points at the bottom)
    k=1;j=1
    DO i=1,Nx
      Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
      output(Gidx)  = PHI(Gidx)
    ENDDO
    !
    k=1;j=Ny
    DO i=1,Nx
      Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
      output(Gidx)  = PHI(Gidx)
    ENDDO
    k=1;i=1
    DO j=1,Ny
      Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
      output(Gidx)  = PHI(Gidx)
    ENDDO
    !
    k=1;i=Nx
    DO j=1,Ny
      Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
      output(Gidx)  = PHI(Gidx)
    ENDDO
  ENDIF
  !
  ! GD: 2D cases bottom corners have to be treated
  IF ((GhostGridX==1).AND.(Ny==1)) THEN
    k = 1; j = 1; i=1
    Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz
    output(Gidx)  = PHI(Gidx)
	!
    k = 1; j = 1; i=Nx
    Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz
    output(Gidx)  = PHI(Gidx)
  ENDIF
  IF ((GhostGridY==1).AND.(Nx==1)) THEN
    k = 1; i = 1; j=1
    Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz
    output(Gidx)  = PHI(Gidx)
    !
    k = 1; i = 1; j=Ny
    Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz
    output(Gidx)  = PHI(Gidx)
  ENDIF
  
!print*,'output=',output
!read*
  
  
END SUBROUTINE BuildLinearSystem
