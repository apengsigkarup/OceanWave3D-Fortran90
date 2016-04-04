SUBROUTINE PreparePreconditioner(PreconditioningMatrix,FineGrid,GhostGridX, GhostGridY, GhostGridZ, &
		alpha, beta, gamma, Precond, CurvilinearONOFF,fileop)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: GhostGridX, GhostGridY, GhostGridZ, alpha, beta, gamma, kappa, Precond, CurvilinearONOFF, iERR
INTEGER :: Nxg,Nyg,Nzg,Nx,Ny,Nz,fileop
TYPE (Level_def)       :: FineGrid
TYPE (SparseArray_COO) :: PreconditioningMatrix
TYPE (Wavefield_FS)    :: tmp_wavefield
TYPE (Diff_def)        :: FullRankStencils
!TYPE (Diff_def) :: DiffStencils
!REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: tmpdsigma
Nx  = FineGrid%Nx
Ny  = FineGrid%Ny
Nz  = FineGrid%Nz
Nxg = FineGrid%Nx + 2*GhostGridX
Nyg = FineGrid%Ny + 2*GhostGridY
Nzg = FineGrid%Nz + GhostGridZ
SELECT CASE (Precond)
	CASE (1) ! Linear LU-preconditining
		IF (CurvilinearONOFF==1) THEN
		    ! FDM stencils - vertical stencils needed for precond. matrix generation
	        CALL PreProcessDiffStencilsZ(FineGrid,FineGrid%DiffStencils,GhostGridZ,gamma,fileop)

			! determine curvilinear transformation weights for the 2D plane
			CALL DetermineCurvilinearTransform2D(FineGrid,alpha,beta,gamma,GhostGridX,GhostGridY,GhostGridZ)
			! determine normal vectors at boundary nodes for the 2D plane boundaries
 			CALL ComputeNormalVectors(FineGrid,GhostGridX,GhostGridY,GhostGridZ)

			! Determine linear sigma-coefficients
			ALLOCATE(FineGrid%dsigmanew(Nzg,Nxg,Nyg,5),STAT=iERR)
            CALL CheckError(iERR,11)
			FineGrid%dsigmanew = zero

	        CALL ALLOCATE_Wavefield_Type(tmp_wavefield, Nx, Ny, Nz, GhostGridX, GhostGridy, GhostGridZ, 0)

		    CALL DetermineTransformationConstantsArray(Nxg,Nyg,Nzg,FineGrid,FineGrid%dsigmanew,tmp_wavefield)
	        CALL DEALLOCATE_Wavefield_Type(tmp_wavefield, Nx, Ny, Nz, 0)

			! Construct preconditioning matrix

			! First, deallocate stencils for computational grid and instead determine these to the order of the
			kappa = alpha
			IF (alpha/=beta) THEN
				! FIXME: just picking the largest of alpha and beta here... perhaps check that they are equal in 3D
				kappa = MAX(alpha,beta)
			END IF
			CALL DetermineGenericStencils(FineGrid%CurvilinearStuff%DiffStencilsPrecond,kappa)
            !
            ! GD: Determine the cross derivatives coefficients
            CALL ConstructTableCrossDerivatives_Curvilinear(FineGrid, FineGrid%CurvilinearStuff%DiffStencilsPrecond, kappa, &
            	GhostGridX, GhostGridY, GhostGridZ)
    		!CALL ConstructTableCrossDerivatives_Curvilinear(FineGrid, FineGrid%CurvilinearStuff%DiffStencils, kappa, &
            !	GhostGridX, GhostGridY, GhostGridZ)
			! now we have both transformation weights, geometic information (normal vectors) and stencils...
			! let's construct a linear preconditioning matrix using this information

			CALL BuildLinearSystemTransformedMatrixCurvilinear(FineGrid, kappa, GhostGridX, GhostGridY, GhostGridZ)
			DEALLOCATE(FineGrid%DiffStencils%StencilZ)
            ! GD: deallocate the cross derivative matrices...
            DEALLOCATE(FineGrid%CurvilinearStuff%DiffStencilsPrecond%IndexesX_XZorXY)
            DEALLOCATE(FineGrid%CurvilinearStuff%DiffStencilsPrecond%IndexesY_YZorXY)
            DEALLOCATE(FineGrid%CurvilinearStuff%DiffStencilsPrecond%IndexesZ_XZorYZ)
            DEALLOCATE(FineGrid%CurvilinearStuff%DiffStencilsPrecond%StencilsX_XZorXY)
            DEALLOCATE(FineGrid%CurvilinearStuff%DiffStencilsPrecond%StencilsY_YZorXY)
            DEALLOCATE(FineGrid%CurvilinearStuff%DiffStencilsPrecond%StencilsZ_XZorYZ)
!$$$$$$             DEALLOCATE(FineGrid%CurvilinearStuff%DiffStencils%IndexesX_XZorXY)
!$$$$$$             DEALLOCATE(FineGrid%CurvilinearStuff%DiffStencils%IndexesY_YZorXY)
!$$$$$$             DEALLOCATE(FineGrid%CurvilinearStuff%DiffStencils%IndexesZ_XZorYZ)
!$$$$$$             DEALLOCATE(FineGrid%CurvilinearStuff%DiffStencils%StencilsX_XZorXY)
!$$$$$$             DEALLOCATE(FineGrid%CurvilinearStuff%DiffStencils%StencilsY_YZorXY)
!$$$$$$             DEALLOCATE(FineGrid%CurvilinearStuff%DiffStencils%StencilsZ_XZorYZ)

		ELSE
			ALLOCATE( FineGrid%dsigmanew(Nzg,Nxg,Nyg,5),STAT=iERR )
            CALL CheckError(iERR,10)
            
		    CALL PreProcessDiffStencils(FineGrid,FineGrid%DiffStencils,GhostGridX,GhostGridY,GhostGridZ, alpha,beta,gamma,fileop)
!            print*,'alpha=',alpha
!            print*,'beta=',beta
!            print*,'gamma=',gamma
!            print*,'Nxg=',Nxg
!            print*,'Nyg=',Nyg
!            print*,'Nzg=',Nzg
!            read*            
			CALL PrepareFullOperatorStencils(FineGrid%DiffStencils,FullRankStencils,alpha,beta,gamma,Nxg,Nyg,Nzg) ! table for generating linear system
			CALL ALLOCATE_Wavefield_Type(tmp_wavefield, Nx, Ny, Nz, GhostGridX, GhostGridy, GhostGridZ, 0)
		    CALL DetermineTransformationConstantsArray(Nxg,Nyg,Nzg,FineGrid,FineGrid%dsigmanew,tmp_wavefield)
            ! GD: Determine the cross derivatives coefficients
            CALL ConstructTableCrossDerivatives(FineGrid, FineGrid%DiffStencils, gamma, GhostGridX, GhostGridY, GhostGridZ, 1)
            !
		    CALL BuildLinearSystemMatrix(Nxg,Nyg,Nzg,PreconditioningMatrix,GhostGridX, GhostGridY, GhostGridZ,&
			     FineGrid%dsigmanew,FineGrid,FullRankStencils,alpha,beta,gamma)
			CALL DEALLOCATE_Wavefield_Type(tmp_wavefield, Nx, Ny, Nz, 0)
            ! GD: deallocate the cross derivative matrices...
            DEALLOCATE(FineGrid%DiffStencils%IndexesX_XZorXY)
            DEALLOCATE(FineGrid%DiffStencils%IndexesY_YZorXY)
            DEALLOCATE(FineGrid%DiffStencils%IndexesZ_XZorYZ)
            DEALLOCATE(FineGrid%DiffStencils%StencilsX_XZorXY)
            DEALLOCATE(FineGrid%DiffStencils%StencilsY_YZorXY)
            DEALLOCATE(FineGrid%DiffStencils%StencilsZ_XZorYZ)
			! Deallocate the full rank stencils
			DEALLOCATE(FineGrid%DiffStencils%FullRankIndexXZ)
			DEALLOCATE(FineGrid%DiffStencils%FullRankIndexYZ)
			DEALLOCATE(FineGrid%DiffStencils%FullRankStencilXZ)
			DEALLOCATE(FineGrid%DiffStencils%FullRankStencilYZ)
		ENDIF
		CALL CleanSparseMatrixCSR(FineGrid%PreconditioningMatrix,FineGrid%PreconditioningMatrix_CSR)
!	filename = "SparseMatrix.bin"
!CALL StoreSparseMatrix(FineGrid%PreconditioningMatrix,filename)
!	filename = "dsigma.bin"
!CALL StoreRealArray(FineGrid%dsigma,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),5,filename)
		PRINT*,'  Preconditioning matrix generated.'
		WRITE(FILEOP,*)'  Preconditioning matrix generated.'
	CASE (3) ! Multigrid
	CASE DEFAULT
		PRINT *,'Error: No default precondtioning strategy. (PreparePreconditioner)'
		STOP
END SELECT
END SUBROUTINE PreparePreconditioner
