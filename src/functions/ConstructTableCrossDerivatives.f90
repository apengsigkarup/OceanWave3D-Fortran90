SUBROUTINE ConstructTableCrossDerivatives(CompGrid, DiffStencils, kappa, GhostGridX, GhostGridY, GhostGridZ, Precond)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def), INTENT(IN) :: CompGrid
TYPE (Diff_def)  :: DiffStencils
INTEGER :: GhostGridX, GhostGridY, GhostGridZ, Precond, fileop
INTEGER :: alpha, beta, kappa, i, j, k, diffa, diffb, diffg, rank,ranka,rankb, Nxg, Nyg, Nzg
INTEGER, DIMENSION((CompGrid%Nz+GhostGridZ),(CompGrid%Nx+2*GhostGridX),(CompGrid%Ny+2*GhostGridY)) :: L
!
INTEGER, DIMENSION(2*kappa+1) :: idx, idxa, idxb, idxg
INTEGER :: ii, jj, kk, Gidx
! GD: to exclude corner ghost points, define a global matrix clearer and more efficient?
INTEGER :: idxa_c, idxb_c, idxg_c, idxa_c2, idxb_c2, idxg_c2
! Local variables...
REAL(KIND=long), DIMENSION(2*kappa+1,2) :: stencilweightsZ

REAL(KIND=long),DIMENSION(2*kappa+1,2*kappa+1,2*kappa+1) :: Stencil3D
INTEGER,DIMENSION(2*kappa+1,2*kappa+1,2*kappa+1) :: Index3D
INTEGER :: IndexX,IndexZ,IndexY,fullrank
REAL(KIND=long) :: StencilX,StencilY,StencilZ

Nxg = CompGrid%Nx + 2*GhostGridX
Nyg = CompGrid%Ny + 2*GhostGridY
Nzg = CompGrid%Nz + GhostGridZ

!print*,'WARNING (APEK): error in SUBROUTINE ConstructTableCrossDerivatives. It is not setup correctly for 2D cases.'
!STOP

rank = 2*kappa+1
fullrank = rank**3
! ALLOCATE the corresponding DiffStencils elements... (assumed here that rank=ranka=rankb...)
! XY derivatives not used in straight boundaries... => size=1
ALLOCATE(DiffStencils%IndexesX_XZorXY(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,rank,1))
ALLOCATE(DiffStencils%IndexesY_YZorXY(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,rank,1))
ALLOCATE(DiffStencils%IndexesZ_XZorYZ(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,rank,2))
ALLOCATE(DiffStencils%StencilsX_XZorXY(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,rank,1))
ALLOCATE(DiffStencils%StencilsY_YZorXY(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,rank,1))
ALLOCATE(DiffStencils%StencilsZ_XZorYZ(CompGrid%Nz+GhostGridZ,CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,rank,2))
! Allocate the full rank stencils if we are in preconditionning step...
IF (Precond == 1) THEN
   ALLOCATE(DiffStencils%FullRankIndexXZ(Nzg*Nxg*Nyg,fullrank))
   ALLOCATE(DiffStencils%FullRankIndexYZ(Nzg*Nxg*Nyg,fullrank))
   ALLOCATE(DiffStencils%FullRankStencilXZ(Nzg*Nxg*Nyg,fullrank))
   ALLOCATE(DiffStencils%FullRankStencilYZ(Nzg*Nxg*Nyg,fullrank))
ENDIF
! Initialize to zero
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
            idxa_c  = 0
        	idxb_c  = 0
            idxg_c  = 0
            idxg_c2 = 0
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
            !  idxb_c = 0
            !  idxg_c2 = 0
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
           ! Test influence corner treatment...
           !idxa_c=0;idxb_c=0;idxg_c=0;idxg_c2=0;idxa_c2=0;idxb_c2=0
           !
           ! For correction in XY derivative... only useful for curvilinear...
           !
           !! IndexesX in XY derivatives
           !DiffStencils%IndexesX_XZorXY(k,i,j,1:ranka,2) = idxa(1:ranka)+idxa_c2
           !! IndexesY in XY derivatives
           !DiffStencils%IndexesY_YZorXY(k,i,j,1:rankb,2) = idxb(1:rankb)+idxb_c2
           !! StencilsX in XY derivatives
           !DiffStencils%StencilsX_XZorXY(k,i,j,1:ranka,2) = CompGrid%DiffStencils%StencilX(i-idxa_c2,1:ranka,1)
           !! StencilsY in XY derivatives
           !DiffStencils%StencilsY_YZorXY(k,i,j,1:ranka,2) = CompGrid%DiffStencils%StencilY(j-idxb_c2,1:rankb,1)
           !
           ! For correction in XZ derivative...
           !
           ! IndexesX in XZ derivatives
           DiffStencils%IndexesX_XZorXY(k,i,j,1:ranka,1) = idxa(1:ranka)+idxa_c
           ! IndexesZ in XZ derivatives
           DiffStencils%IndexesZ_XZorYZ(k,i,j,1:rank,1) = idxg(1:rank)+idxg_c   		
          ! StencilsX in XZ derivatives
           DiffStencils%StencilsX_XZorXY(k,i,j,1:ranka,1) = CompGrid%DiffStencils%StencilX(i-idxa_c,1:ranka,1)
           ! StencilsZ in XZ derivatives
           IF (idxg_c == 0) THEN
             !do nothing (initial)
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
           DiffStencils%StencilsY_YZorXY(k,i,j,1:rankb,1) = CompGrid%DiffStencils%StencilY(j-idxb_c,1:rankb,1)
           ! StencilsZ in YZ derivatives
           IF (idxg_c2 == 0) THEN
                 !do nothing (initial)
                 stencilweightsZ = CompGrid%DiffStencils%StencilZ(k,1:rank,1:2)
           ELSE
              CALL DiffStencilsZ_modif(k,CompGrid%z,stencilweightsZ,CompGrid%Nz,GhostGridZ,kappa,idxg_c2)
           ENDIF
           DiffStencils%StencilsZ_XZorYZ(k,i,j,1:rank,2) = stencilweightsZ(1:rank,1)
          !
       ENDDO
   ENDDO
ENDDO
!
IF (Precond == 1) THEN
   ! Compute the FullRankStencils to use in the preconditionning Matrix
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
           IF (i<alpha+1) THEN
               diffa = alpha + 1 - i
           ELSE IF (i > Nxg - alpha) THEN
               diffa = (Nxg - alpha) - i
           ELSE
               diffa = 0
           ENDIF
           idxa = i + idx + diffa
           !
           DO k = 1, Nzg ! 1+GhostGridZ , Nzg-1 !
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
              ! Treat the possible 2D case...
              IF(Nxg == 1) THEN
                idxa(:)=1
              ENDIF
              IF(Nyg == 1) THEN
                idxb(:)=1
              ENDIF
              !
              IF(Nxg>1) THEN
                 ! DXZ
                 Stencil3D=zero
                 Index3D=0
                 DO jj = 1,rankb !beta+1-Diffb,beta+1-Diffb ! evaluate only on one point?
                   DO kk = 1, rank
                     IndexZ = DiffStencils%IndexesZ_XZorYZ(k,i,j,kk,1)
                     StencilZ = DiffStencils%StencilsZ_XZorYZ(k,i,j,kk,1)
                     DO ii = 1, ranka
                         IndexX = DiffStencils%IndexesX_XZorXY(IndexZ,i,j,ii,1)
                         stencilX = DiffStencils%StencilsX_XZorXY(IndexZ,i,j,ii,1)
                         !
                         Stencil3D(kk,ii,jj) = StencilX*StencilZ
                         Index3D(kk,ii,jj) = L(IndexZ,IndexX,idxb(jj))
                     END DO
                   END DO
                 END DO
                 DiffStencils%FullRankStencilXZ(L(k,i,j),:) = RESHAPE2(fullrank,Stencil3D)
                 DiffStencils%FullRankIndexXZ(L(k,i,j),:) = PACK(Index3D,.TRUE.)
               ENDIF

				IF (Nyg>1) THEN
                  ! DYZ
                  Stencil3D=zero
                  Index3D=0
                  DO ii = 1,ranka! alpha+1-Diffa,alpha+1-Diffa ! evaluate only on one point?
                    DO kk = 1, rank
                      IndexZ = DiffStencils%IndexesZ_XZorYZ(k,i,j,kk,2)
                      StencilZ = DiffStencils%StencilsZ_XZorYZ(k,i,j,kk,2)
                      DO jj = 1,rankb
                          IndexY = DiffStencils%IndexesY_YZorXY(IndexZ,i,j,jj,1)
                          stencilY = DiffStencils%StencilsY_YZorXY(IndexZ,i,j,jj,1)
                          !
                          Stencil3D(kk,ii,jj) = StencilY*StencilZ
                          Index3D(kk,ii,jj) = L(IndexZ,idxa(ii),IndexY)
                      END DO
                    END DO
                  END DO
                  DiffStencils%FullRankStencilYZ(L(k,i,j),:) = RESHAPE2(fullrank,Stencil3D)
                  DiffStencils%FullRankIndexYZ(L(k,i,j),:) = PACK(Index3D,.TRUE.)
                ENDIF
       ENDDO
     ENDDO
   ENDDO
ENDIF
!
CONTAINS

FUNCTION RESHAPE2(nnz,A) RESULT(B)
! Convert array with nnz elements to one-dimensional array of nnz elements
USE Precision
IMPLICIT NONE
INTEGER:: nnz
REAL(KIND=long), DIMENSION(nnz) :: A, B
B=A
END FUNCTION RESHAPE2
END SUBROUTINE ConstructTableCrossDerivatives
