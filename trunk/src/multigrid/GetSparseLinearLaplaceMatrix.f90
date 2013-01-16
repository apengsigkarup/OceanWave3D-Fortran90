SUBROUTINE GetSparseLinearLaplaceMatrix(Gridlevel,nnz,icoo,jcoo,acoo)
! By Allan P. Engsig-Karup.
  USE GlobalVariables
  USE MGLevels
  IMPLICIT NONE
  EXTERNAL BuildLinearSystem, BuildLinearSystemTransformedCurvilinear
  INTEGER, INTENT (IN)  :: Gridlevel
  INTEGER, INTENT (OUT) :: nnz
  INTEGER, DIMENSION(:), POINTER :: icoo, jcoo, acoo

  CALL Initialize
  CALL ReadInputFileParameters
  CALL SetupCompDomain
  CALL InitializeVariables
  time = dt
  CALL SetupInitialConditions
  IF (Precond==3) THEN ! PREPARE FOR MULTIGRID-PRECONDITIONING
    CALL MGPreProcess ( FineGrid, GhostGridX, GhostGridY, GhostGridZ, MGCoarseningStrategy, alphaprecond, betaprecond, &
    	gammaprecond, Precond, MGmaxgrids, CurvilinearONOFF)
  ELSE
	WRITE(*,*)'Error: Wrong Preconditioning strategy chosen.'
	STOP
  ENDIF
  IF (Gridlevel>=1 .AND. Gridlevel<=MG_N_levels) THEN
	  nnz = arrLevels(Gridlevel)%PreconditioningMatrix%nnz
	  ALLOCATE(icoo(nnz), jcoo(nnz), acoo(nnz), STAT=STAT)
	  CALL CheckError(STAT,7)
	  icoo = arrLevels(Gridlevel)%PreconditioningMatrix%row_ptr
	  jcoo = arrLevels(Gridlevel)%PreconditioningMatrix%col_ind
	  acoo = arrLevels(Gridlevel)%PreconditioningMatrix%val
	  WRITE (*,'(A)') 'Sparse system matrix extracted from a grid with dimensions:'
	  WRITE (*,'(A,I3,I3,I3)') '(Nx,Ny,Nz) = ',arrLevels(Gridlevel)%Nx+2*GhostGridX,arrLevels(Gridlevel)%Ny+2*GhostGridY,&
      	arrLevels(Gridlevel)%Nz+GhostGridZ
  ELSE
	  WRITE(*,*)'Error: GridLevel not in interval [1...MaxGridNumber] In GetSparseLinearLaplaceMatrix.'
	  STOP
  END IF

END SUBROUTINE GetSparseLinearLaplaceMatrix
