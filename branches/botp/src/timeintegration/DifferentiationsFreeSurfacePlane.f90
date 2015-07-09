SUBROUTINE DifferentiationsFreeSurfacePlane(Wavefield,GhostGridX,GhostGridY,FineGrid,alpha,beta)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
USE GlobalVariables, ONLY: curvilinearONOFF,GhostGridZ
IMPLICIT NONE
INTEGER :: GhostGridX, GhostGridY, alpha, beta
TYPE (Level_def) :: FineGrid
TYPE(Wavefield_FS) :: Wavefield

IF (curvilinearONOFF==0) THEN
    ! Differentiation in the Free Surface plane
	IF (FineGrid%Nx>1) THEN
	   CALL UpdateGhostLayerX(Wavefield%E,Wavefield%SourceEx,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
       		FineGrid%DiffStencils,alpha,GhostGridX,GhostGridY)
       CALL DiffXEven(Wavefield%E,Wavefield%Ex, 1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid%DiffStencils,alpha)
       CALL DiffXEven(Wavefield%E,Wavefield%Exx,2,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid%DiffStencils,alpha)
	   CALL UpdateGhostLayerX(Wavefield%P,Wavefield%SourcePx,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
       		FineGrid%DiffStencils,alpha,GhostGridX,GhostGridY)
       CALL DiffXEven(Wavefield%P,Wavefield%Px, 1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid%DiffStencils,alpha)
	END IF

	IF (FineGrid%Ny>1) THEN
	   CALL UpdateGhostLayerY(Wavefield%E,Wavefield%SourceEy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
       		FineGrid%DiffStencils,beta,GhostGridX,GhostGridY)
       CALL DiffYEven(Wavefield%E,Wavefield%Ey, 1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid%DiffStencils,beta)
       CALL DiffYEven(Wavefield%E,Wavefield%Eyy,2,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid%DiffStencils,beta)
	   CALL UpdateGhostLayerY(Wavefield%P,Wavefield%SourcePy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
       		FineGrid%DiffStencils,beta,GhostGridX,GhostGridY)
       CALL DiffYEven(Wavefield%P,Wavefield%Py, 1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid%DiffStencils,beta)
	END IF
ELSE
    ! Differentiation in the Free Surface plane
    !IF (FineGrid%Nx>1) THEN
    !   CALL UpdateGhostLayerXCurvilinear(Wavefield%E,Wavefield%SourceEx,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
    !        FineGrid,alpha,beta,GhostGridX,GhostGridY)
    !   CALL UpdateGhostLayerXCurvilinear(Wavefield%P,Wavefield%SourcePx,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
    !        FineGrid,alpha,beta,GhostGridX,GhostGridY)
    !END IF
	!
    !IF (FineGrid%Ny>1) THEN
    !   CALL UpdateGhostLayerYCurvilinear(Wavefield%E,Wavefield%SourceEy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
    !        FineGrid,alpha,beta,GhostGridX,GhostGridY)
    !   CALL UpdateGhostLayerYCurvilinear(Wavefield%P,Wavefield%SourcePy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
    !        FineGrid,alpha,beta,GhostGridX,GhostGridY)
    !END IF
    ! Probably needed to impose through varn=0... to prevent from problems if rotation of 90 degrees on boundaries
!    IF ((FineGrid%Nx>1).AND.(FineGrid%Ny>1)) THEN
!       CALL UpdateGhostLayerNCurvilinear(Wavefield%E,Wavefield%SourceEx,Wavefield%SourceEy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
!            FineGrid%Nz+GhostGridZ,FineGrid,alpha,beta,GhostGridX,GhostGridY)
!       CALL UpdateGhostLayerNCurvilinear(Wavefield%P,Wavefield%SourcePx,Wavefield%SourcePy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
!            FineGrid%Nz+GhostGridZ,FineGrid,alpha,beta,GhostGridX,GhostGridY)
!       !
!       CALL UpdateGhostLayerECurvilinear(Wavefield%E,Wavefield%SourceEx, Wavefield%SourceEy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
!            FineGrid%Nz+GhostGridZ,FineGrid,alpha,beta,GhostGridX,GhostGridY)
!       CALL UpdateGhostLayerECurvilinear(Wavefield%P,Wavefield%SourcePx, Wavefield%SourcePy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
!            FineGrid%Nz+GhostGridZ,FineGrid,alpha,beta,GhostGridX,GhostGridY)
!    ENDIF

    IF ((FineGrid%Nx>1).AND.(FineGrid%Ny>1)) THEN
		CALL UpdateGhostLayerCurvilinear(WaveField%E,Wavefield%SourceEx,Wavefield%SourceEy,FineGrid%Nx+2*GhostGridX,&
       FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,FineGrid,alpha,beta,GhostGridX,GhostGridY)
		CALL UpdateGhostLayerCurvilinear(WaveField%P,Wavefield%SourcePx,Wavefield%SourcePy,FineGrid%Nx+2*GhostGridX,&
       FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,FineGrid,alpha,beta,GhostGridX,GhostGridY)
	END IF	

	IF (FineGrid%Nx>1) THEN
!	   CALL UpdateGhostLayerXCurvilinear(Wavefield%E,Wavefield%SourceEx,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
!			FineGrid,alpha,beta,GhostGridX,GhostGridY)
	   CALL DiffXEvenCurvilinear(Wavefield%E,Wavefield%Ex,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,alpha)
	   !CALL DiffXXEvenCurvilinear(Wavefield%E,Wavefield%Exx,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,alpha)
       !GD: modif for correct cross derivatives treatment...
       CALL DiffXXEvenCurvilinear_CD(Wavefield%E,Wavefield%Exx,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,alpha,&
       FineGrid%CurvilinearStuff%DiffStencils%IndexesX_XZorXY(1,:,:,:,2), &
       FineGrid%CurvilinearStuff%DiffStencils%StencilsX_XZorXY(1,:,:,:,2), &
       FineGrid%CurvilinearStuff%DiffStencils%IndexesY_YZorXY(1,:,:,:,2), &
       FineGrid%CurvilinearStuff%DiffStencils%StencilsY_YZorXY(1,:,:,:,2))
!	   CALL UpdateGhostLayerXCurvilinear(Wavefield%P,Wavefield%SourcePx,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
!			FineGrid,alpha,beta,GhostGridX,GhostGridY)
	   CALL DiffXEvenCurvilinear(Wavefield%P,Wavefield%Px,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,alpha)
	END IF

	IF (FineGrid%Ny>1) THEN
!	   CALL UpdateGhostLayerYCurvilinear(Wavefield%E,Wavefield%SourceEy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
!			FineGrid,alpha,beta,GhostGridX,GhostGridY)
	   CALL DiffYEvenCurvilinear(Wavefield%E,Wavefield%Ey,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,alpha)
	   !CALL DiffYYEvenCurvilinear(Wavefield%E,Wavefield%Eyy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,alpha)
       !GD: modif for correct cross derivatives treatment...
       CALL DiffYYEvenCurvilinear_CD(Wavefield%E,Wavefield%Eyy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,alpha,&
       FineGrid%CurvilinearStuff%DiffStencils%IndexesX_XZorXY(1,:,:,:,2), &
       FineGrid%CurvilinearStuff%DiffStencils%StencilsX_XZorXY(1,:,:,:,2), &
       FineGrid%CurvilinearStuff%DiffStencils%IndexesY_YZorXY(1,:,:,:,2), &
       FineGrid%CurvilinearStuff%DiffStencils%StencilsY_YZorXY(1,:,:,:,2))
!	   CALL UpdateGhostLayerYCurvilinear(Wavefield%P,Wavefield%SourcePy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
!			FineGrid,alpha,beta,GhostGridX,GhostGridY)
	   CALL DiffYEvenCurvilinear(Wavefield%P,Wavefield%Py,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,alpha)
	END IF
END IF
END SUBROUTINE
