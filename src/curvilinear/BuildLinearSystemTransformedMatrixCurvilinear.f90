SUBROUTINE BuildLinearSystemTransformedMatrixCurvilinear(CompGrid,kappa, GhostGridX, GhostGridY, GhostGridZ)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def), INTENT(INOUT) :: CompGrid
REAL(KIND=long) :: NormalX, NormalY, NormalZ
REAL(KIND=long) :: ex, exx, nx, nxx, ey, eyy, ny, nyy
INTEGER :: alpha, beta, kappa, iconnect, jconnect, kconnect, i, j, k, n, diffa, diffb, diffg, nnz,rank,ranka,rankb,TotalPoints,&
	Nxg, Nyg, Nzg, Gidx, Lidx, GhostGridX, GhostGridY, GhostGridZ, stencilmaxrank, dummy
INTEGER, DIMENSION(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY) :: GidxTable, GidxTableBC
INTEGER, DIMENSION(2*kappa+1) :: stencilidx, idx, idxa, idxb, idxg
REAL(KIND=long), DIMENSION(2*kappa+1,2*kappa+1,2) :: stencilweightsX, stencilweightsY !, stencilweightsZ
REAL(KIND=long), DIMENSION(2*kappa+1,2) :: stencilweightsZ
REAL(KIND=long), DIMENSION(2*kappa+1)   :: stencilweights
REAL(KIND=long) :: sigmax, sigmay, sigmaxxyy, sigmaz, stencilweights2
REAL(KIND=long), DIMENSION(:), ALLOCATABLE    :: A_val
INTEGER, DIMENSION(:), ALLOCATABLE :: A_colind, A_rowptr

!print*,'BUILDING a matrix using BuildLinearSystemTransformedMatrixCurvilinear subroutine... '

! ALLOCATE TEMPORARY VECTORS FOR SPARSE OPERATOR ASSEMBLY
rank = 2*kappa+1
Nxg = CompGrid%Nx + 2*GhostGridX
Nyg = CompGrid%Ny + 2*GhostGridY
Nzg = CompGrid%Nz +   GhostGridZ

IF (Nxg>1) THEN
	alpha = kappa
	ranka = rank
ELSE
	alpha = 0
	ranka = 1
END IF
IF (Nyg>1) THEN
	beta  = kappa
	rankb = rank
ELSE
	beta  = 0
	rankb = 1
END IF

TotalPoints = Nxg*Nyg*Nzg
stencilmaxrank = rank**2
! factor of 12 below due to number of terms that needs to be added in laplace operator
dummy = 12*TotalPoints*stencilmaxrank
ALLOCATE( A_val(dummy)    )
ALLOCATE( A_colind(dummy) )
ALLOCATE( A_rowptr(dummy) )
A_val =  zero; A_colind = 0; A_rowptr = 0

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

nnz = 0
DO i = -kappa, kappa
	nnz = nnz + 1
	idx(nnz) = i
END DO

nnz = 0
! Free surface
k = Nzg
DO j = 1+GhostGridY , Nyg-GhostGridY
	DO i = 1+GhostGridX , Nxg-GhostGridX
		nnz = nnz + 1
		Gidx = GidxTable(k,i,j)
		A_val(nnz)    = one
		A_colind(nnz) = Gidx
	    A_rowptr(nnz) = Gidx
	END DO
END DO

stencilweightsX = CompGrid%CurvilinearStuff%DiffStencilsPrecond%StencilG
stencilweightsY = CompGrid%CurvilinearStuff%DiffStencilsPrecond%StencilG

!! If problem is only 2D, then make sure derivatives vanish in last dimension
!IF (Nxg==1) THEN
!	stencilweightsX = zero
!END IF
!IF (Nyg==1) THEN
!	stencilweightsY = zero
!END IF

! Interior
DO j = 1+GhostGridY , Nyg-GhostGridY

	IF (j<beta+1) THEN
		diffb = beta + 1 - j
	ELSE IF (j > Nyg - beta) THEN
		diffb = (Nyg - beta) - j
	ELSE
		diffb = 0
	ENDIF
	idxb = j + idx + diffb

	DO i = 1+GhostGridX , Nxg-GhostGridX

		IF (i<alpha+1) THEN
			diffa = alpha + 1 - i
		ELSE IF (i > Nxg - alpha) THEN
			diffa = (Nxg - alpha) - i
		ELSE
			diffa = 0
		ENDIF
		idxa = i + idx + diffa

		DO k = 1+GhostGridZ , Nzg-1

			IF (k<kappa+1) THEN
				diffg = kappa + 1 - k
			ELSE IF (k > Nzg - kappa) THEN
				diffg = (Nzg - kappa) - k
			ELSE
				diffg = 0
			ENDIF
			stencilweightsZ = CompGrid%DiffStencils%StencilZ(k,1:rank,1:2)
			idxg = k + idx + diffg

			Gidx = GidxTable(k,i,j)

			ex  = CompGrid%CurvilinearStuff%ex(i,j)
			exx = CompGrid%CurvilinearStuff%exx(i,j)
			ey  = CompGrid%CurvilinearStuff%ey(i,j)
			eyy = CompGrid%CurvilinearStuff%eyy(i,j)
			nx  = CompGrid%CurvilinearStuff%nx(i,j)
			nxx = CompGrid%CurvilinearStuff%nxx(i,j)
			ny  = CompGrid%CurvilinearStuff%ny(i,j)
			nyy = CompGrid%CurvilinearStuff%nyy(i,j)

			sigmax    = CompGrid%dsigmanew(k,i,j,2)
			sigmaxxyy = CompGrid%dsigmanew(k,i,j,3)
			sigmay    = CompGrid%dsigmanew(k,i,j,4)
			sigmaz    = CompGrid%dsigmanew(k,i,j,5)

			! first terms
			! d/de
			IF (Nxg>1) THEN
			stencilidx     = GidxTable(k,idxa,j)
			stencilweights = stencilweightsX(1:rank,kappa+1-diffa,1)
			A_val(nnz+1:nnz+rank)    = (exx+eyy)*stencilweights
			A_colind(nnz+1:nnz+rank) = stencilidx
			A_rowptr(nnz+1:nnz+rank) = Gidx
			nnz = nnz + rank

			! d/dee
			stencilidx     = GidxTable(k,idxa,j)
			stencilweights = stencilweightsX(1:rank,kappa+1-diffa,2)
			A_val(nnz+1:nnz+rank)    = (ex**2+ey**2)*stencilweights
			A_colind(nnz+1:nnz+rank) = stencilidx
			A_rowptr(nnz+1:nnz+rank) = Gidx
			nnz = nnz + rank
			END IF

			! d/dn
			IF (Nyg>1) THEN
			stencilidx     = GidxTable(k,i,idxb)
			stencilweights = stencilweightsY(1:rank,kappa+1-diffb,1)
			A_val(nnz+1:nnz+rank)    = (nxx+nyy)*stencilweights
			A_colind(nnz+1:nnz+rank) = stencilidx
		    A_rowptr(nnz+1:nnz+rank) = Gidx
			nnz = nnz + rank

			! d/dnn
			stencilidx     = GidxTable(k,i,idxb)
			stencilweights = stencilweightsY(1:rank,kappa+1-diffb,2)
			A_val(nnz+1:nnz+rank)    = (nx**2+ny**2)*stencilweights
			A_colind(nnz+1:nnz+rank) = stencilidx
		    A_rowptr(nnz+1:nnz+rank) = Gidx
			nnz = nnz + rank
			END IF

			! d/ds
			stencilidx     = GidxTable(idxg,i,j)
			stencilweights = stencilweightsZ(1:rank,1)
			A_val(nnz+1:nnz+rank)    = (sigmaxxyy)*stencilweights
			A_colind(nnz+1:nnz+rank) = stencilidx
			A_rowptr(nnz+1:nnz+rank) = Gidx
			nnz = nnz + rank

			! d/dss
			stencilidx     = GidxTable(idxg,i,j)
			stencilweights = stencilweightsZ(1:rank,2)
			A_val(nnz+1:nnz+rank)    = (sigmax**2+sigmay**2+sigmaz**2)*stencilweights
			A_colind(nnz+1:nnz+rank) = stencilidx
		    A_rowptr(nnz+1:nnz+rank) = Gidx
			nnz = nnz + rank

			! cross-differential terms
			!        MUST match the choice of differentiation order in BuildLinearSystem

            ! d/dedn - terms : d/de first, then d/dn
            IF (Nxg>1 .AND. Nyg>1) THEN
            DO n = 1, rankb
                !Lidx = idxb(n) ! d/dn
                !stencilweights2 = stencilweightsY(n,kappa+1-diffb,1) ! d/dn
                !! GD: correct treatment of cross derivatives
                Lidx = CompGrid%CurvilinearStuff%DiffStencilsPrecond%IndexesY_YZorXY(k,i,j,n,2) ! d/dn
                stencilweights2 = CompGrid%CurvilinearStuff%DiffStencilsPrecond%StencilsY_YZorXY(k,i,j,n,2) ! d/dn

                !stencilidx     = GidxTable(k,idxa,Lidx) ! d/de
                !stencilweights = stencilweightsX(1:rank,kappa+1-diffa,1) ! d/de
                !! GD: correct treatment of cross derivatives
                stencilidx = GidxTable(k,CompGrid%CurvilinearStuff%DiffStencilsPrecond%IndexesX_XZorXY(k,i,Lidx,1:ranka,2),Lidx) ! d/de
                stencilweights = CompGrid%CurvilinearStuff%DiffStencilsPrecond%StencilsX_XZorXY(k,i,Lidx,1:ranka,2) ! d/de

                A_val(nnz+1:nnz+rank)    = stencilweights2*(two*(ex*nx+ey*ny))*stencilweights
                A_colind(nnz+1:nnz+rank) = stencilidx
                A_rowptr(nnz+1:nnz+rank) = Gidx
                nnz = nnz + rank
            END DO
            END IF

            ! d/dsde - terms : d/ds first, then d/de
            IF (Nxg>1) THEN
            DO n = 1, rank
                !Lidx = idxg(n)              ! d/ds
                !stencilweights2 = stencilweightsZ(n,1) ! d/ds
                !! GD: correct treatment of cross derivatives
                Lidx = CompGrid%CurvilinearStuff%DiffStencilsPrecond%IndexesZ_XZorYZ(k,i,j,n,1) ! d/ds
                stencilweights2 = CompGrid%CurvilinearStuff%DiffStencilsPrecond%StencilsZ_XZorYZ(k,i,j,n,1) ! d/ds


                !stencilidx     = GidxTable(Lidx,idxa,j) ! d/de
                !stencilweights = stencilweightsX(1:rank,kappa+1-diffa,1) ! d/de
                !! GD: correct treatment of cross derivatives
                stencilidx = GidxTable(Lidx,CompGrid%CurvilinearStuff%DiffStencilsPrecond%IndexesX_XZorXY(Lidx,i,j,1:ranka,1),j) ! d/de
                stencilweights = CompGrid%CurvilinearStuff%DiffStencilsPrecond%StencilsX_XZorXY(Lidx,i,j,1:ranka,1) ! d/de

                A_val(nnz+1:nnz+rank)    = stencilweights2*(two*(sigmax*ex+sigmay*ey))*stencilweights
                A_colind(nnz+1:nnz+rank) = stencilidx
                A_rowptr(nnz+1:nnz+rank) = Gidx
                nnz = nnz + rank
            END DO
            END IF

            ! d/dsdn - terms : d/ds first, then d/dn
            IF (Nyg>1) THEN
            DO n = 1, rank
                !Lidx = idxg(n)              ! d/ds
                !stencilweights2 = stencilweightsZ(n,1) ! d/ds
                !! GD: correct treatment of cross derivatives
                Lidx = CompGrid%CurvilinearStuff%DiffStencilsPrecond%IndexesZ_XZorYZ(k,i,j,n,2) ! d/ds
                stencilweights2 = CompGrid%CurvilinearStuff%DiffStencilsPrecond%StencilsZ_XZorYZ(k,i,j,n,2) ! d/ds

                !stencilidx     = GidxTable(Lidx,i,idxb) ! d/dn
                !stencilweights = stencilweightsY(1:rank,kappa+1-diffb,1) ! d/dn
                !! GD: correct treatment of cross derivatives
                stencilidx = GidxTable(Lidx,i,CompGrid%CurvilinearStuff%DiffStencilsPrecond%IndexesY_YZorXY(Lidx,i,j,1:rankb,1)) ! d/dn
                stencilweights = CompGrid%CurvilinearStuff%DiffStencilsPrecond%StencilsY_YZorXY(Lidx,i,j,1:rankb,1) ! d/dn

                A_val(nnz+1:nnz+rank)    = stencilweights2*(two*(sigmax*nx+sigmay*ny))*stencilweights
                A_colind(nnz+1:nnz+rank) = stencilidx
                A_rowptr(nnz+1:nnz+rank) = Gidx
                nnz = nnz + rank
            END DO
            END IF
		END DO
	END DO
END DO

! Kinematic boundary conditions (imposed through ghost layers)
! Let's do a trick here in order to be able to detect where to impose the kinematic boundary
! conditions easily.
!
! Impose through all ghost points except corner nodes at the bottom and in the vertical. Ghost points next to free surface included.
GidxTableBC = GidxTable
GidxTableBC(2:Nzg,1+GhostGridX:Nxg-GhostGridX,1+GhostGridY:Nyg-GhostGridY) = 0 ! interior points + free surface points
IF (Nxg>1 .AND. Nyg>1) THEN
	! 3D
	GidxTableBC(1:Nzg,1,1)         = 0 ! southwest corner points
	GidxTableBC(1:Nzg,Nxg,1)       = 0 ! southeast corner points
	GidxTableBC(1:Nzg,Nxg,Nyg)     = 0 ! northeast corner points
	GidxTableBC(1:Nzg,1,Nyg)       = 0 ! northwest corner points
	GidxTableBC(1,1+GhostGridX:Nxg-GhostGridX,1)   = 0 ! west  bottom points
	GidxTableBC(1,1+GhostGridX:Nxg-GhostGridX,Nyg) = 0 ! east  bottom points
	GidxTableBC(1,1,1+GhostGridY:Nyg-GhostGridY)   = 0 ! south bottom points
	GidxTableBC(1,Nxg,1+GhostGridY:Nyg-GhostGridY) = 0 ! north bottom points
ELSE IF (Nxg>1) THEN
	! 2D: xz-plane
	GidxTableBC(1,1,1)     = 0 ! southwest corner
	GidxTableBC(1,Nxg,1)   = 0 ! southeast corner
!	GidxTableBC(Nzg,1,1)   = 0 ! northwest corner
!	GidxTableBC(Nzg,Nxg,1) = 0 ! northeast corner
ELSE ! Nyg>1
	! 2D: yz-plane
	GidxTableBC(1,1,1)     = 0 ! southwest corner
	GidxTableBC(1,1,Nyg)   = 0 ! southeast corner
!	GidxTableBC(Nzg,1,1)   = 0 ! northwest corner
!	GidxTableBC(Nzg,1,Nyg) = 0 ! northeast corner
END IF
!
DO j = 1 , Nyg

	IF (j<beta+1) THEN
		diffb = beta + 1 - j
	ELSE IF (j > Nyg - beta) THEN
		diffb = (Nyg - beta) - j
	ELSE
		diffb = 0
	ENDIF
	idxb = j + idx + diffb

	if (j>1 .AND. j<Nyg) then
		jconnect = 0
	else
		jconnect = SGNINT(diffb)
	end if

	DO i = 1, Nxg

		IF (i<alpha+1) THEN
			diffa = alpha + 1 - i
		ELSE IF (i > Nxg - alpha) THEN
			diffa = (Nxg - alpha) - i
		ELSE
			diffa = 0
		ENDIF
		idxa = i + idx + diffa

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
				idxg = k + idx + diffg	

				ex  = CompGrid%CurvilinearStuff%ex(i+iconnect,j+jconnect)
				ey  = CompGrid%CurvilinearStuff%ey(i+iconnect,j+jconnect)
				nx  = CompGrid%CurvilinearStuff%nx(i+iconnect,j+jconnect)
				ny  = CompGrid%CurvilinearStuff%ny(i+iconnect,j+jconnect)

				if (k>1) then
					kconnect = 0
				else
					kconnect = SGNINT(diffg)
				end if
				sigmaz = CompGrid%dsigmanew(k+kconnect,i+iconnect,j+jconnect,5)

				! Normal vectors are defined at the ghost points used to impose the kinematic boundary conditions
				NormalX = CompGrid%CurvilinearStuff%NormalX(k,i,j)
				NormalY = CompGrid%CurvilinearStuff%NormalY(k,i,j)
				NormalZ = CompGrid%CurvilinearStuff%NormalZ(k,i,j)

				! first terms
				IF (Nxg>1) THEN
				stencilidx     = GidxTable(k,idxa,j) ! d/de
				stencilweights = stencilweightsX(1:rank,kappa+1-diffa+iconnect,1) ! d/de
				A_val(nnz+1:nnz+rank)    = (NormalX*ex+NormalY*ey)*stencilweights
				A_colind(nnz+1:nnz+rank) = stencilidx
				A_rowptr(nnz+1:nnz+rank) = Gidx
				nnz = nnz + rank
				END IF

				IF (Nyg>1) THEN
				stencilidx     = GidxTable(k,i,idxb) ! d/dn
				stencilweights = stencilweightsY(1:rank,kappa+1-diffb+jconnect,1) ! d/dn
				A_val(nnz+1:nnz+rank)    = (NormalX*nx+NormalY*ny)*stencilweights
				A_colind(nnz+1:nnz+rank) = stencilidx
				A_rowptr(nnz+1:nnz+rank) = Gidx
				nnz = nnz + rank
				END IF

				! last terms
				stencilidx     = GidxTable(idxg,i,j) ! d/ds
			    !stencilweights = CompGrid%DiffStencils%StencilZ(2,1:rank,1) ! d/ds
                !GD: modification, error before? why was it taken on k=2?
                stencilweights = CompGrid%DiffStencils%StencilZ(k+kconnect,1:rank,1) ! d/ds
				A_val(nnz+1:nnz+rank)    = (NormalZ*sigmaz)*stencilweights
				A_colind(nnz+1:nnz+rank) = stencilidx
				A_rowptr(nnz+1:nnz+rank) = Gidx
				nnz = nnz + rank
			END IF
		END DO
	END DO
END DO

! Ghost points
! Let's do a trick here in order to be able to detect where to make sure that corner ghost points are
! remain unchanged
GidxTableBC = GidxTable
IF (Nyg>1) THEN
	GidxTableBC(2:Nzg,1:Nxg,1+GhostGridY:Nyg-GhostGridY) = 0 ! interior points
END IF
IF (Nxg>1) THEN
	GidxTableBC(2:Nzg,1+GhostGridX:Nxg-GhostGridX,1:Nyg) = 0 ! interior points
END IF
GidxTableBC(1,1+GhostGridX:Nxg-GhostGridX,1+GhostGridY:Nyg-GhostGridY)   = 0 ! bottom points
DO j = 1, Nyg
	DO i = 1, Nxg
		DO k = 1, Nzg
			Gidx = GidxTable(k,i,j)

			! FIXME: Trick is made here... make more efficient
			IF (Gidx==GidxTableBC(k,i,j)) THEN

				nnz           = nnz + 1				
				A_val(nnz)    = one
				A_colind(nnz) = Gidx
				A_rowptr(nnz) = Gidx

			END IF

		END DO
	END DO
END DO

ALLOCATE(CompGrid%PreconditioningMatrix%val(nnz),CompGrid%PreconditioningMatrix%row_ptr(nnz)&
	,CompGrid%PreconditioningMatrix%col_ind(nnz))

CompGrid%PreconditioningMatrix%val     = A_val(1:nnz)
CompGrid%PreconditioningMatrix%col_ind = A_colind(1:nnz)
CompGrid%PreconditioningMatrix%row_ptr = A_rowptr(1:nnz)
CompGrid%PreconditioningMatrix%nnz     = nnz
CompGrid%PreconditioningMatrix%nrow    = Nxg*Nyg*Nzg

DEALLOCATE( A_val    )
DEALLOCATE( A_colind )
DEALLOCATE( A_rowptr )

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
END SUBROUTINE BuildLinearSystemTransformedMatrixCurvilinear
