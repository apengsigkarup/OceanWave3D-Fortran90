SUBROUTINE PreparePreconditionerOLD(PreconditioningMatrix, FullRankStencils, FineGrid, GhostGridX, GhostGridY, GhostGridZ, &
		   alpha, beta, gamma, Precond)
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: GhostGridX, GhostGridY, GhostGridZ, alpha, beta, gamma, Precond
TYPE (Level_def) :: FineGrid
TYPE (Diff_def) :: FullRankStencils
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: tmpdsigma
TYPE (SparseArray_COO) :: PreconditioningMatrix
TYPE (Wavefield_FS) :: tmp_wavefield
ALLOCATE( tmpdsigma((FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),5) )
IF (Precond==1) THEN
	CALL ALLOCATE_Wavefield_Type(tmp_wavefield, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, GhostGridX, GhostGridy, GhostGridZ, 0)
    CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
	     FineGrid,tmpdsigma,tmp_wavefield)
    CALL BuildLinearSystemMatrix(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,PreconditioningMatrix,&
		 GhostGridX, GhostGridY, GhostGridZ,&
	     tmpdsigma,FineGrid,FullRankStencils,alpha,beta,gamma)
	CALL DEALLOCATE_Wavefield_Type(tmp_wavefield, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, 0)
ELSE
	PRINT*,'Error: Chosen preconditioning strategy not implemented.'
	STOP
END IF
DEALLOCATE( tmpdsigma )
PRINT*,'  Preconditioning matrix generated.'
END SUBROUTINE PreparePreconditionerOLD
