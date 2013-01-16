SUBROUTINE BuildLinearSystemTransformedCurvilinear(Nxg,Nyg,Nzg,PHI,out,CompGrid,alpha,beta,gamma)
!
! Build linear system for transformed laplace problem.
! Find solution in computational domain and if GhostGrid
! then define kinematic bottom constraint on those.
!
! By Allan P. Engsig-Karup.
! Note: it is assumed that Nxg, Nyg, Nzg includes ghost points and that alpha=beta=gamma.
! GD: SWENS inhomogeneous Neumann boundary conditions come here...
!
USE Precision
USE Constants
USE DataTypes
USE GlobalVariables, ONLY: GhostGridX, GhostGridY, GhostGridZ, swenseONOFF, Wavefield ! GD: addition for SWENS ! FIXME: global Wavefield a problem for multigrid!
IMPLICIT NONE
TYPE (Level_def), INTENT(IN) :: CompGrid
INTEGER :: Nxg, Nyg, Nzg, alpha, beta, gamma, kappa
REAL(KIND=long), DIMENSION(2*alpha+1,2*alpha+1,2) :: DiffStencilArray
REAL(KIND=long), DIMENSION(Nzg, Nxg, Nyg), INTENT(IN)  :: PHI
REAL(KIND=long), DIMENSION(Nzg, Nxg, Nyg), INTENT(OUT)  :: out
REAL(KIND=long), DIMENSION(Nzg, Nxg, Nyg) :: tmp, &
	dPHIde, dPHIdee, dPHIdn, dPHIdnn, dPHIden, dPHIds, dPHIdss, dPHIdsde, dPHIdsdn
REAL(KIND=long) :: ex, exx, ey, eyy, nx, nxx, ny, nyy, DMe, DMn, DMs, DM2en, DM2e, DM2n, DM2s, &
		DMse, DMsn, dPHIdxx, dPHIdyy, sigmax, sigmaxxyy, sigmay, sigmaz, NormalX, NormalY, NormalZ, &
		dPHIdx, dPHIdy, dPHIdz
INTEGER :: nnz, i, j, k, iconnect, jconnect, kconnect, diffa, diffb, diffg, Gidx
INTEGER, DIMENSION(Nzg, Nxg, Nyg) :: GidxTable, GidxTableBC

IF (alpha==beta .AND. alpha==gamma) THEN
	kappa = alpha
ELSE
	PRINT*,'Error: alpha /= beta /= gamma. (BuildLinearSystemTransformedCurvilinear)'
	STOP
END IF

DiffStencilArray = CompGrid%CurvilinearStuff%DiffStencils%StencilG

! Global numbering table used for lookups
nnz = 0
DO j = 1 , Nyg
	DO i = 1 , Nxg
		DO k = 1 , Nzg
			nnz = nnz + 1
			GidxTable(k,i,j) = nnz
		END DO
	END DO
END DO

! Determine needed derivatives explicitly
IF (Nxg>1) THEN
	CALL DiffXuniform3D(PHI,dPHIde ,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
	CALL DiffXuniform3D(PHI,dPHIdee,2,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
END IF
IF (Nyg>1) THEN
	CALL DiffYuniform3D(PHI,dPHIdn ,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
	CALL DiffYuniform3D(PHI,dPHIdnn,2,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
END IF
CALL DiffZArbitrary(PHI,dPHIds,1,Nxg,Nyg,Nzg,CompGrid%DiffStencils,kappa)
CALL DiffZArbitrary(PHI,dPHIdss,2,Nxg,Nyg,Nzg,CompGrid%DiffStencils,kappa)
!CALL DiffZuniform3D(PHI,dPHIds ,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
!CALL DiffZuniform3D(PHI,dPHIdss,2,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?

! FIXME: choice has been made with order of differentiation for cross-products here:
!        MUST match the choice in BuildLinearSystemMatrix
IF (Nxg>1 .AND. Nyg>1) THEN
	!tmp = dPHIde
	!CALL DiffYuniform3D(dPHIde,dPHIden,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
	!! GD: modif for correct cross derivatives near corners...
    CALL DiffXuniform3D_CD(PHI,tmp,CompGrid%CurvilinearStuff%DiffStencils%IndexesX_XZorXY(:,:,:,:,2),&
    	CompGrid%CurvilinearStuff%DiffStencils%StencilsX_XZorXY(:,:,:,:,2),Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
    CALL DiffYuniform3D_CD(tmp,dPHIden,CompGrid%CurvilinearStuff%DiffStencils%IndexesY_YZorXY(:,:,:,:,2), &
    	CompGrid%CurvilinearStuff%DiffStencils%StencilsY_YZorXY(:,:,:,:,2),Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
END IF

! Mixed derivatives in transformed domain
! FIXME: perhaps more appropriate to determine derivative in x- and y-directions first due to boundary conditions.
!tmp = dPHIds
!IF (Nxg>1) THEN
!	CALL DiffXuniform3D(tmp,dPHIdsde ,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
!END IF
!IF (Nyg>1) THEN
!	CALL DiffYuniform3D(tmp,dPHIdsdn ,1,DiffStencilArray,Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
!END IF
IF (Nxg>1) THEN
!tmp = dPHIde
!CALL DiffZArbitrary(tmp,dPHIdsde,1,Nxg,Nyg,Nzg,CompGrid%DiffStencils,kappa)
!! GD: modif for correct cross derivatives near corners...
CALL DiffXuniform3D_CD(PHI,tmp, CompGrid%CurvilinearStuff%DiffStencils%IndexesX_XZorXY(:,:,:,:,1),&
	CompGrid%CurvilinearStuff%DiffStencils%StencilsX_XZorXY(:,:,:,:,1),Nxg,Nyg,Nzg,kappa) ! FIXME: perhaps diff should only be on interior nodes?
CALL DiffZArbitrary_CD(tmp,dPHIdsde, CompGrid%CurvilinearStuff%DiffStencils%IndexesZ_XZorYZ(:,:,:,:,1),&
	CompGrid%CurvilinearStuff%DiffStencils%StencilsZ_XZorYZ(:,:,:,:,1),Nxg,Nyg,Nzg,kappa)
END IF
IF (Nyg>1) THEN
!tmp = dPHIdn
!!CALL DiffZArbitrary(tmp,dPHIdsdn,1,Nxg,Nyg,Nzg,CompGrid%DiffStencils,kappa)
!! GD: modif for correct cross derivatives near corners...
CALL DiffYuniform3D_CD(PHI,tmp,CompGrid%CurvilinearStuff%DiffStencils%IndexesY_YZorXY(:,:,:,:,1),&
	CompGrid%CurvilinearStuff%DiffStencils%StencilsY_YZorXY(:,:,:,:,1),Nxg,Nyg,Nzg,kappa)
CALL DiffZArbitrary_CD(tmp,dPHIdsdn,CompGrid%CurvilinearStuff%DiffStencils%IndexesZ_XZorYZ(:,:,:,:,2),&
	CompGrid%CurvilinearStuff%DiffStencils%StencilsZ_XZorYZ(:,:,:,:,2),Nxg,Nyg,Nzg,kappa)
END IF

! Free Surface
k = Nzg
DO j = 1, Nyg
	DO i = 1, Nxg
		out(k,i,j) = PHI(k,i,j)
	END DO
END DO

! Interior
DO j = 2 , Nyg-1
	DO i = 2 , Nxg-1
		DO k = 2 , Nzg-1
			ex  = CompGrid%CurvilinearStuff%ex(i,j)
			exx = CompGrid%CurvilinearStuff%exx(i,j)
			ny  = CompGrid%CurvilinearStuff%ny(i,j)
			nyy = CompGrid%CurvilinearStuff%nyy(i,j)
			ey  = CompGrid%CurvilinearStuff%ey(i,j)
			eyy = CompGrid%CurvilinearStuff%eyy(i,j)
			nx  = CompGrid%CurvilinearStuff%nx(i,j)
			nxx = CompGrid%CurvilinearStuff%nxx(i,j)

			DMe   = dPHIde(k,i,j)
			DMn   = dPHIdn(k,i,j)
			DMs   = dPHIds(k,i,j)
			DM2en = dPHIden(k,i,j)
			DM2e  = dPHIdee(k,i,j)
			DM2n  = dPHIdnn(k,i,j)
			DM2s  = dPHIdss(k,i,j)

			DMse  = dPHIdsde(k,i,j)
			DMsn  = dPHIdsdn(k,i,j)

			dPHIdxx = exx*DMe+nxx*DMn+two*ex*nx*DM2en+ex**2*DM2e+nx**2*DM2n
			dPHIdyy = eyy*DMe+nyy*DMn+two*ey*ny*DM2en+ey**2*DM2e+ny**2*DM2n

			sigmax    = CompGrid%dsigmanew(k,i,j,2)
			sigmaxxyy = CompGrid%dsigmanew(k,i,j,3)
			sigmay    = CompGrid%dsigmanew(k,i,j,4)
			sigmaz    = CompGrid%dsigmanew(k,i,j,5)

			out(k,i,j) = dPHIdxx + dPHIdyy + sigmaxxyy*DMs + two*(sigmax*ex+sigmay*ey)*DMse +&
	            two*(sigmax*nx+sigmay*ny)*DMsn + (sigmax**2 + sigmay**2 + sigmaz**2)*DM2s
		END DO
	END DO
END DO

! Kinematic boundary conditions (imposed through ghost layers)
! Let's do a trick here in order to be able to detect where to impose the kinematic boundary
! conditions easily.
GidxTableBC = GidxTable
GidxTableBC(2:Nzg,1+GhostGridX:Nxg-GhostGridX,1+GhostGridY:Nyg-GhostGridY) = 0 ! interior points + free surface points
GidxTableBC(:,1,1)     = 0 ! southwest corner points
GidxTableBC(:,Nxg,1)   = 0 ! southeast corner points
GidxTableBC(:,Nxg,Nyg) = 0 ! northeast corner points
GidxTableBC(:,1,Nyg)   = 0 ! northwest corner points
IF (Nxg>1 .AND. Nyg>1) THEN
	GidxTableBC(1,1+GhostGridX:Nxg-GhostGridX,1)   = 0 ! west  bottom points
	GidxTableBC(1,1+GhostGridX:Nxg-GhostGridX,Nyg) = 0 ! east  bottom points
	GidxTableBC(1,1,1+GhostGridY:Nyg-GhostGridY)   = 0 ! south bottom points
	GidxTableBC(1,Nxg,1+GhostGridY:Nyg-GhostGridY) = 0 ! north bottom points
END IF
!
DO j = 1 , Nyg

	IF (j<kappa+1) THEN
		diffb = kappa + 1 - j
	ELSE IF (j > Nyg - kappa) THEN
		diffb = (Nyg - kappa) - j
	ELSE
		diffb = 0
	ENDIF

	if (j>1 .AND. j<Nyg) then
		jconnect = 0
	else
		jconnect = SGNINT(diffb)
	end if

	DO i = 1, Nxg

		IF (i<kappa+1) THEN
			diffa = kappa + 1 - i
		ELSE IF (i > Nxg - kappa) THEN
			diffa = (Nxg - kappa) - i
		ELSE
			diffa = 0
		ENDIF

		if (i>1 .AND. i<Nxg) then
			iconnect = 0
		else
			iconnect = SGNINT(diffa)
		end if

		DO k = 1 , Nzg

			Gidx = GidxTable(k,i,j)

			! FIXME: Trick is made here... make more efficient
			IF (Gidx==GidxTableBC(k,i,j)) THEN
				IF (k<kappa+1) THEN
					diffg = kappa + 1 - k
				ELSE IF (k > Nzg - kappa) THEN
					diffg = (Nzg - kappa) - k
				ELSE
					diffg = 0
				ENDIF

				if (k>1) then
					kconnect = 0
				else
					kconnect = SGNINT(diffg)
				end if

				ex  = CompGrid%CurvilinearStuff%ex(i+iconnect,j+jconnect)
				ey  = CompGrid%CurvilinearStuff%ey(i+iconnect,j+jconnect)
				nx  = CompGrid%CurvilinearStuff%nx(i+iconnect,j+jconnect)
				ny  = CompGrid%CurvilinearStuff%ny(i+iconnect,j+jconnect)

				! Normal vectors are defined at the ghost points used to impose the kinematic boundary conditions
				NormalX = CompGrid%CurvilinearStuff%NormalX(k,i,j)
				NormalY = CompGrid%CurvilinearStuff%NormalY(k,i,j)
				NormalZ = CompGrid%CurvilinearStuff%NormalZ(k,i,j)

				DMe   = dPHIde(k+kconnect,i+iconnect,j+jconnect)
				DMn   = dPHIdn(k+kconnect,i+iconnect,j+jconnect)
				DMs   = dPHIds(k+kconnect,i+iconnect,j+jconnect)

				dPHIdx = ex*DMe + nx*DMn
				dPHIdy = ey*DMe + ny*DMn

				sigmaz = CompGrid%dsigmanew(k+kconnect,i+iconnect,j+jconnect,5)
				dPHIdz = sigmaz*DMs
	
				out(k,i,j) = NormalX*dPHIdx + NormalY*dPHIdy + NormalZ*dPHIdz
                ! GD: SWENS addition
                IF(swenseONOFF/=0) THEN
                  out(k,i,j) = out(k,i,j) 	+ NormalX*Wavefield%Px_I_bp(Wavefield%GidxTableBP(k,i,j)) &
                  							+ NormalY*Wavefield%Py_I_bp(Wavefield%GidxTableBP(k,i,j)) &
                                            + NormalZ*Wavefield%Pz_I_bp(Wavefield%GidxTableBP(k,i,j))
                ENDIF

			END IF
		END DO
	END DO
END DO
! Ghost points
! Let's do a trick here in order to be able to detect where to make sure that corner ghost points are
! remain unchanged
GidxTableBC = GidxTable
GidxTableBC(2:Nzg,1:Nxg,1+GhostGridY:Nyg-GhostGridY) = 0 ! interior points
GidxTableBC(2:Nzg,1+GhostGridX:Nxg-GhostGridX,1:Nyg) = 0 ! interior points
GidxTableBC(1,1+GhostGridX:Nxg-GhostGridX,1+GhostGridY:Nyg-GhostGridY)   = 0 ! interior points
DO j = 1, Nyg
	DO i = 1, Nxg
		DO k = 1, Nzg
			Gidx = GidxTable(k,i,j)
			! FIXME: Trick is made here... make more efficient
			IF (Gidx==GidxTableBC(k,i,j)) THEN
				out(k,i,j) = PHI(k,i,j)
			END IF
		END DO
	END DO
END DO

CONTAINS

!>
!! Determine sign of integer
!<
FUNCTION SGNINT(I) RESULT(SGN)
IMPLICIT NONE
INTEGER:: I, SGN
IF (I>0) THEN
	SGN = 1
ELSE IF (I<0) THEN
	SGN = -1
ELSE
	SGN = 0
ENDIF
END FUNCTION SGNINT
END SUBROUTINE BuildLinearSystemTransformedCurvilinear

