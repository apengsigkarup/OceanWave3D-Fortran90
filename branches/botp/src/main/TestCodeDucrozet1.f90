SUBROUTINE TestCodeDucrozet1
  !
  ! Test code to check the operators.
  !
  USE GlobalVariables
  USE MGLevels
  IMPLICIT NONE
  EXTERNAL BuildLinearSystem, BuildLinearSystemTransformedCurvilinear
  ! GD: to test the cross derivatives...
  REAL(KIND=long), DIMENSION(:,:,:), ALLOCATABLE :: tmpPHI
  TYPE (Diff_def)        :: FullRankStencils
  INTEGER i, j, k
  IF (curvilinearONOFF==1) THEN
     ! determine curvilinear transformation weights for the 2D plane

     IF (FineGrid%Nx>1 .AND. FineGrid%Ny>1) THEN
        CALL DetermineCurvilinearTransform2D(FineGrid,alpha,beta,gamma,GhostGridX,GhostGridY,GhostGridZ)
        ! determine normal vectors at boundary nodes for the 2D plane boundaries
        CALL ComputeNormalVectors(FineGrid,GhostGridX,GhostGridY,GhostGridZ)
     ELSE
        PRINT*,'Error: curvilinear transformation for the plane cannot be invoked for 2D problems.'
        STOP
     END IF
  ENDIF
     PRINT*,'HELLO1'
  ! Determine linear sigma-coefficients
  tmp2D = zero
  ALLOCATE(FineGrid%dsigmanew(FineGrid%Nz+GhostGridZ,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,5))
  FineGrid%dsigmanew = zero
  ! GD: We use the wavefield_FS type...
  ! GD: Allocate the wavefield
  CALL ALLOCATE_Wavefield_Type(Wavefield_tmp, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, GhostGridX, GhostGridy, GhostGridZ, 0)
  CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
       FineGrid,FineGrid%dsigmanew,Wavefield_tmp)
  CALL DEALLOCATE_Wavefield_Type(Wavefield_tmp, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, 0)
  !CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
  !	FineGrid,FineGrid%dsigmanew,tmp2D,tmp2D,tmp2D,tmp2D,tmp2D)
PRINT*,'HELLO2'

  ! Construct preconditioning matrix

  ! First, deallocate stencils for computational grid and instead determine these to the order of the
  kappa = alphaprecond
  IF (alphaprecond/=betaprecond) THEN
     ! FIXME: just picking the largest of alpha and beta here... perhaps check that they are equal in 3D
     kappa = MAX(alphaprecond,betaprecond)
  END IF
  IF (curvilinearONOFF.EQ.1) THEN
     !
     ! FDM stencils - vertical stencils needed for precond. matrix generation
     CALL PreProcessDiffStencilsZ(FineGrid,FineGrid%DiffStencils,GhostGridZ,kappa)

     CALL DetermineGenericStencils(FineGrid%CurvilinearStuff%DiffStencilsPrecond,kappa)
     ! GD: Determine the cross derivatives coefficients
     CALL ConstructTableCrossDerivatives_Curvilinear(FineGrid, FineGrid%CurvilinearStuff%DiffStencilsPrecond, kappa, &
          GhostGridX, GhostGridY, GhostGridZ)
     ! now we have both transformation weights, geometric information (normal vectors) and stencils...
     ! let's construct a linear preconditioning matrix using this information
  ELSE
     !$$$$$$       CALL PreProcessDiffStencils(FineGrid,FineGrid%DiffStencils,GhostGridX,GhostGridY,GhostGridZ, alpha,beta,gamma)
     !$$$$$$       CALL PrepareFullOperatorStencils(FineGrid%DiffStencils,FullRankStencils,alpha,beta,gamma, &
     !$$$$$$             FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ) ! table for generating linear system
     !$$$$$$
     !$$$$$$       ! GD: Determine the cross derivatives coefficients
     !$$$$$$       CALL ConstructTableCrossDerivatives(FineGrid, FineGrid%DiffStencils, gamma, GhostGridX, GhostGridY, GhostGridZ, 1)
  ENDIF
  ! Test system matrix
  IF (.FALSE.) THEN
     IF (curvilinearONOFF.EQ.1) THEN
        kappa = alpha
        IF (alpha/=beta) THEN
           ! FIXME: just picking the largest of alpha and beta here... perhaps check that they are equal in 3D
           kappa = MAX(alpha,beta)
        END IF
        !CALL PreProcessDiffStencils(FineGrid,FineGrid%DiffStencils,GhostGridX,GhostGridY,GhostGridZ, &
        !   	 alphaprecond,betaprecond,gammaprecond)
        CALL DetermineGenericStencils(FineGrid%CurvilinearStuff%DiffStencilsPrecond,kappa)
        ! FDM stencils - vertical stencils needed for precond. matrix generation
        CALL PreProcessDiffStencilsZ(FineGrid,FineGrid%DiffStencils,GhostGridZ,kappa)
        ! GD: Determine the cross derivatives coefficients
        CALL ConstructTableCrossDerivatives_Curvilinear(FineGrid, FineGrid%CurvilinearStuff%DiffStencils, kappa, &
             GhostGridX, GhostGridY, GhostGridZ)
        ALLOCATE(ee((FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)))
        ALLOCATE(tm((FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)))
        ALLOCATE(A((FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)))
        ee = zero
        tm = zero
        A  = zero
        DO i = 1 , (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)
           ee(i) = one
           !   CALL BuildLinearSystemTransformedCurvilinear(FineGrid, ee, tm,GhostGridX,GhostGridY,GhostGridZ,kappa)
           CALL BuildLinearSystemTransformedCurvilinear(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY, &
                FineGrid%Nz+GhostGridZ,ee,tm,FineGrid,alpha,beta,gamma)
           A(:,i) = tm
           ee(i) = zero
        END DO
        filename = "A.bin"
        CALL StoreRealArray(A,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),filename,formattype)

        ! save the vertical derivative for linear stability analysis...
        ee = zero
        tm = zero
        A  = zero
        DO i = 1 , (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)
           ee(i) = one
           CALL DiffZArbitrary(ee,tm,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ, &
                FineGrid%DiffStencils,gamma)
           A(:,i) = tm
           ee(i) = zero
        END DO
        filename = "DMz.bin"
        CALL StoreRealArray(A,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),filename,formattype)

        DEALLOCATE(ee,tm,A)
        !   print*,'stopped here for now...'
        !stop
     ELSE
        kappa = alpha
        IF (alpha/=beta) THEN
           ! FIXME: just picking the largest of alpha and beta here... perhaps check that they are equal in 3D
           kappa = MAX(alpha,beta)
        END IF

        CALL PreProcessDiffStencils(FineGrid,FineGrid%DiffStencils,GhostGridX,GhostGridY,GhostGridZ, alpha,beta,gamma)
        CALL PrepareFullOperatorStencils(FineGrid%DiffStencils,FullRankStencils,alpha,beta,gamma, &
             (FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),(FineGrid%Nz+GhostGridZ)) ! table for generating linear system

        ! GD: Determine the cross derivatives coefficients
        CALL ConstructTableCrossDerivatives(FineGrid, FineGrid%DiffStencils, gamma, GhostGridX, GhostGridY, GhostGridZ, 1)
        !
        ALLOCATE(ee((FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)))
        ALLOCATE(tm((FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)))
        ALLOCATE(A((FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)))
        ee = zero
        tm = zero
        A  = zero

        DO i = 1 , (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)
           ee(i) = one
           !   CALL BuildLinearSystemTransformedCurvilinear(FineGrid, ee, tm,GhostGridX,GhostGridY,GhostGridZ,kappa)
           !CALL BuildLinearSystemTransformedCurvilinear(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY, &
           !     FineGrid%Nz+GhostGridZ,ee,tm,FineGrid,alpha,beta,gamma)
           CALL BuildLinearSystem((FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),(FineGrid%Nz+GhostGridZ), &
                ee,tm,FineGrid,alpha,beta,gamma)
           A(:,i) = tm
           ee(i) = zero
        END DO
        filename = "A.bin"
        CALL StoreRealArray(A,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),filename,formattype)

        ! save the vertical derivative for linear stability analysis...
        ee = zero
        tm = zero
        A  = zero
        DO i = 1 , (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)
           ee(i) = one
           CALL DiffZArbitrary(ee,tm,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ, &
                FineGrid%DiffStencils,gamma)
           A(:,i) = tm
           ee(i) = zero
        END DO
        filename = "DMz.bin"
        CALL StoreRealArray(A,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),filename,formattype)

        DEALLOCATE(ee,tm,A)
        print*,'stopped here for now...'
        stop
     ENDIF
  ENDIF
PRINT*,'HELLO3'

  CALL BuildLinearSystemTransformedMatrixCurvilinear(FineGrid, alphaprecond, GhostGridX, GhostGridY, GhostGridZ)
PRINT*,'HELLO4'

  CALL CleanSparseMatrixCOO(FineGrid%PreconditioningMatrix)
PRINT*,'HELLO5'

  filename = "SparseMatrix.bin"
  CALL StoreSparseMatrix(FineGrid%PreconditioningMatrix,filename,formattype)
PRINT*,'HELLO6'

  CALL FactorPreconditioner(FineGrid%PreconditioningMatrix, &
       (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ))

  print*,'stopped here for now...'
  STOP
END SUBROUTINE TestCodeDucrozet1
