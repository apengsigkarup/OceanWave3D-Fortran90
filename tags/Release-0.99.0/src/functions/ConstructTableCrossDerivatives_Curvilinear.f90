SUBROUTINE ConstructTableCrossDerivatives_Curvilinear(CompGrid, DiffStencils, kappa, GhostGridX, GhostGridY, GhostGridZ)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def), INTENT(IN) :: CompGrid
TYPE (Diff_def)  :: DiffStencils
INTEGER :: alpha, beta, kappa, i, j, k, diffa, diffb, diffg, rank,ranka,rankb, Nxg, Nyg, Nzg, GhostGridX, GhostGridY, GhostGridZ
!
INTEGER, DIMENSION(2*kappa+1) :: idx, idxa, idxb, idxg
INTEGER :: Gidx
! GD: to exclude corner ghost points, define a global matrix clearer and more efficient?
INTEGER :: idxa_c, idxb_c, idxg_c, idxa_c2, idxb_c2, idxg_c2
! Local variables...
REAL(KIND=long), DIMENSION(2*kappa+1,2) :: stencilweightsZ
INTEGER, DIMENSION(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY) :: L

Nxg = CompGrid%Nx + 2*GhostGridX
Nyg = CompGrid%Ny + 2*GhostGridY
Nzg = CompGrid%Nz + GhostGridZ

rank = 2*kappa+1
! ALLOCATE the corresponding DiffStencils elements... (assumed here that rank=ranka=rankb...)
ALLOCATE(DiffStencils%IndexesX_XZorXY(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,rank,2))
ALLOCATE(DiffStencils%IndexesY_YZorXY(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,rank,2))
ALLOCATE(DiffStencils%IndexesZ_XZorYZ(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,rank,2))
ALLOCATE(DiffStencils%StencilsX_XZorXY(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,rank,2))
ALLOCATE(DiffStencils%StencilsY_YZorXY(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,rank,2))
ALLOCATE(DiffStencils%StencilsZ_XZorYZ(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,rank,2))
!
idxa(1:rank)=0
idxb(1:rank)=0
idxg(1:rank)=0

DO i = -kappa, kappa
	idx(i+kappa+1) = i
END DO
!
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
! Generate array with global indices for every point in the domain
Gidx = 0
DO j = 1,Nyg
	DO i = 1, Nxg
		DO k = 1,Nzg
			Gidx = Gidx + 1
			L(k,i,j) = Gidx
		END DO
	END DO
END DO
!
idxa_c=0
idxb_c=0
idxg_c=0
idxa_c2=0
idxb_c2=0
idxg_c2=0
!
! Construction of the Indexes and Stencils needed
!
DO j = 1, Nyg !1+GhostGridY , Nyg-GhostGridY !
    IF (j<beta+1) THEN
        diffb = beta + 1 - j
    ELSE IF (j > Nyg - beta) THEN
        diffb = (Nyg - beta) - j
    ELSE
        diffb = 0
    ENDIF
    idxb = j + idx + diffb

    DO i = 1, Nxg !   1+GhostGridX , Nxg-GhostGridX !
        idxa_c2 = 0
        idxb_c2 = 0
        !
        IF (i<alpha+1) THEN
            diffa = alpha + 1 - i
        ELSE IF (i > Nxg - alpha) THEN
            diffa = (Nxg - alpha) - i
        ELSE
            diffa = 0
        ENDIF
        idxa = i + idx + diffa
        !
        IF ((j==1).OR.(j==Nyg)) THEN
          IF(idxa(1)==1) THEN
        	idxa_c2 = 1
          ENDIF
          IF(idxa(ranka)==Nxg) THEN
        	idxa_c2 = -1
          ENDIF
       ENDIF
       IF ((i==1).OR.(i==Nxg)) THEN
          IF(idxb(1)==1) THEN
        	idxb_c2 = 1
          ENDIF
          IF(idxb(rankb)==Nyg) THEN
        	idxb_c2 = -1
          ENDIF
       ENDIF
       !
        IF ((i == 1).OR.(i==Nxg)) THEN
            idxa_c2 = 0
        ENDIF
        IF ((j == 1).OR.(j==Nyg)) THEN
            idxb_c2 = 0
        ENDIF

        DO k = 1, Nzg ! 1+GhostGridZ , Nzg-1 !
           idxg_c  = 0
           idxg_c2 = 0
           idxa_c  = 0
           idxb_c  = 0
           !
           IF (k<kappa+1) THEN
                diffg = kappa + 1 - k
           ELSE IF (k > Nzg - kappa) THEN
                diffg = (Nzg - kappa) - k
           ELSE
                diffg = 0
           ENDIF
           idxg = k + idx + diffg
           !
           IF ((i==1).OR.(i==Nxg)) THEN
              IF(idxg(1)==1) THEN
                idxg_c = 1
              ENDIF
              ! We can use the ghost points on the FS (but not on the bottom: corner points not solved during resolution)
              !IF(idxg(rank)==Nzg) THEN
              !  idxg_c = -1
              !ENDIF
           ENDIF
		   IF ((j==1).OR.(j==Nyg)) THEN
              IF(idxg(1)==1) THEN
                idxg_c2 = 1
              ENDIF
              ! We can use the ghost points on the FS (but not on the bottom: corner points not solved during resolution)
              !IF(idxg(rank)==Nzg) THEN
              !  idxg_c2 = -1
              !ENDIF
           ENDIF
           !IF ((k==1).OR.(k==Nzg)) THEN
           IF ((k==1)) THEN  ! We can use the ghost points on the FS (but not on the bottom: corner points not solved during resolution)
              IF(idxa(1)==1) THEN
                idxa_c = 1
              ENDIF
              IF(idxa(ranka)==Nxg) THEN
                idxa_c = -1
              ENDIF
              IF(idxb(1)==1) THEN
                idxb_c = 1
              ENDIF
              IF(idxb(rankb)==Nyg) THEN
                idxb_c = -1
              ENDIF
           ENDIF
           !
           IF ((i == 1).OR.(i==Nxg)) THEN
            idxa_c = 0
           ELSE ! FIXME: This test has to be checked...
			!idxb_c  = 0
            !idxg_c2 = 0
            !IF (k == Nzg) THEN ! We can use the ghost points on the FS (but not on the bottom: corner points not solved during resolution)
            !   idxb_c  = 0
            !   idxg_c2 = 0
            !ENDIF
           ENDIF
           IF ((j == 1).OR.(j==Nyg)) THEN
            idxb_c = 0
           ELSE ! FIXME: This test has to be checked...
             !idxa_c = 0
             !idxg_c = 0
             !IF (k == Nzg) THEN ! We can use the ghost points on the FS (but not on the bottom: corner points not solved during resolution)
             !  idxa_c = 0
             !  idxg_c = 0
             !ENDIF
           ENDIF
           IF ((k == 1).OR.(k==Nzg)) THEN
            idxg_c = 0
            idxg_c2 = 0
           ELSE ! FIXME: This test has to be checked...
             ! Nothing to do because the corner edges are not resolved...
           ENDIF
           !test...
           !idxa_c=0;idxb_c=0;idxg_c=0;idxg_c2=0
           !idxa_c2=0;idxb_c2=0
           !
           ! For correction in XY derivative...
           !
           ! IndexesX in XY derivatives
           DiffStencils%IndexesX_XZorXY(k,i,j,1:ranka,2) = idxa(1:ranka)+idxa_c2
           ! IndexesY in XY derivatives
           DiffStencils%IndexesY_YZorXY(k,i,j,1:rankb,2) = idxb(1:rankb)+idxb_c2
           ! StencilsX in XY derivatives
           DiffStencils%StencilsX_XZorXY(k,i,j,1:ranka,2) = DiffStencils%StencilG(1:ranka,kappa+1-diffa-idxa_c2,1)
           ! StencilsY in XY derivatives
           DiffStencils%StencilsY_YZorXY(k,i,j,1:rankb,2) = DiffStencils%StencilG(1:rankb,kappa+1-diffb-idxb_c2,1)
           !
           ! For correction in XZ derivative...
           !
           ! IndexesX in XZ derivatives
           DiffStencils%IndexesX_XZorXY(k,i,j,1:ranka,1) = idxa(1:ranka)+idxa_c
           ! IndexesZ in XZ derivatives
           DiffStencils%IndexesZ_XZorYZ(k,i,j,1:rank,1) = idxg(1:rank)+idxg_c
           ! StencilsX in XZ derivatives
           DiffStencils%StencilsX_XZorXY(k,i,j,1:ranka,1) = DiffStencils%StencilG(1:ranka,kappa+1-diffa-idxa_c,1)
           ! StencilsZ in XZ derivatives
           IF (idxg_c == 0) THEN
            	!do nothing (Initial)
                stencilweightsZ = CompGrid%DiffStencils%StencilZ(k,1:rank,1:2)
           ELSE
              CALL DiffStencilsZ_modif(k,CompGrid%z,stencilweightsZ,CompGrid%Nz,GhostGridZ,kappa,idxg_c)
           ENDIF
           DiffStencils%StencilsZ_XZorYZ(k,i,j,1:rank,1) = stencilweightsZ(1:rank,1)
           !
           ! For correction in YZ derivative...
           !
           ! IndexesY in YZ derivatives
           DiffStencils%IndexesY_YZorXY(k,i,j,1:rankb,1) = idxb(1:rankb)+idxb_c
           ! IndexesZ in YZ derivatives
           DiffStencils%IndexesZ_XZorYZ(k,i,j,1:rank,2) = idxg(1:rank)+idxg_c2
          ! StencilsY in YZ derivatives
           DiffStencils%StencilsY_YZorXY(k,i,j,1:rankb,1) = DiffStencils%StencilG(1:rankb,kappa+1-diffb-idxb_c,1)
           ! StencilsZ in YZ derivatives
           IF (idxg_c2 == 0) THEN
                 !do nothing (initial)
                 stencilweightsZ = CompGrid%DiffStencils%StencilZ(k,1:rank,1:2)
           ELSE
              CALL DiffStencilsZ_modif(k,CompGrid%z,stencilweightsZ,CompGrid%Nz,GhostGridZ,kappa,idxg_c2)
           ENDIF
           DiffStencils%StencilsZ_XZorYZ(k,i,j,1:rank,2) = stencilweightsZ(1:rank,1)
       ENDDO
   ENDDO
ENDDO

END SUBROUTINE ConstructTableCrossDerivatives_Curvilinear
