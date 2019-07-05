SUBROUTINE iterative_solution(rhss,n_rows,A_times_x,sol,CurrentGrid)
! By Allan P. Engsig-Karup.
USE MGlevels
USE Precision
USE Constants
USE SPK
USE GlobalVariables, ONLY: tstep, RKSTAGES, MINITER, MAXITER, TOTALITER, TOTALITERFS, MINITERNOFS, MAXITERNOFS, CYCLET, &
      MAXIT, RELTOL, RINFO, JOB, IS_HSL, MAXIS, SS, MAXS, CNTL, ICNTL, INFOHSL, alpha, beta, gamma, Precond, GhostGridZ, &
	  FineGrid, COLSCA, ROWSCA, KEEP, workspace, ipar, fpar, RESETSOLVER, STAT, abstol, GMRESmaxiterations, solver, GhostGridX, &
	  GhostGridY
  IMPLICIT NONE
  EXTERNAL A_times_x ! external procedures for matrix-vector product
  INTEGER :: n_rows
  REAL(KIND=long), DIMENSION(n_rows), INTENT(IN)    :: rhss
  REAL(KIND=long), DIMENSION(n_rows), INTENT(INOUT) :: sol
  REAL(KIND=long), DIMENSION(n_rows) :: tmp
  TYPE (Level_def), INTENT(IN) :: CurrentGrid
  INTEGER::nnz, iterations, i
  CHARACTER(len=10) :: text

  ipar=0
  fpar=zero;
  ipar(1)  = 0	! always 0 to start an iterative solver
  ipar(2)  = MIN(Precond,1)	! preconditioning strategy (0=none,1=LEFT,2=RIGHT)
  ipar(3)  = 1	! use convergence test scheme 2
  ipar(5)  = GMRESmaxiterations	! use *GMRES(10) (e.g. FGMRES(10))
  ipar(4)  = (n_rows+3)*(ipar(5)+2) + (ipar(5)+1)*ipar(5)/2	! the size of the workspace
IF (solver==0) THEN
  ! high-order defect correction method
!  print*,'   HIGH-ORDER DEFECT CORRECTION METHOD:'
  ipar(4)  = n_rows*2  ! the size of the workspace
ELSE ! default
  ! GMRES
!  print*,'   GMRES METHOD:'
  ipar(4)  = (n_rows+3)*(ipar(5)+2) + (ipar(5)+1)*ipar(5)/2	! the size of the workspace
ENDIF
IF (RESETSOLVER==0) THEN
  ALLOCATE(workspace(ipar(4)),STAT=STAT) ! workspace for residual storage
  PRINT *,'   ITERATIVE SOLVER WORKSPACE SIZE = ',ipar(4)
  CALL CheckError(STAT,6)
  RESETSOLVER = 1 ! MAKE SURE WE ONLY DO THIS ONCE
ENDIF
  workspace=zero
  ipar(6)  = ipar(5)	! use at most iapr(5) matvec's
  fpar(1)  = reltol	    ! relative tolerance
  fpar(2)  = abstol     ! absolute tolerance 100 times smaller than relative tolerance
  fpar(11) = zero	    ! clearing the FLOPS counter

iterations = 0

! GMRES from SPARSKIT2
10 IF (solver==0) THEN
       ! high-order DC method by APEK
       CALL dc(n_rows,rhss,sol,ipar,fpar,workspace)  ! initial guess given as input using "sol"
   ELSE
       CALL gmres(n_rows,rhss,sol,ipar,fpar,workspace)  ! initial guess given as input using "sol"
   ENDIF
!print*,'*** AFTER ITERATIVE SOLVER ITERATION ***'   
   IF (ipar(1).EQ.1) THEN
	  ! DIRECT MATRIX-VECTOR PRODUCT (Determine residual)
      iterations=iterations+1
!      print*,'residual before=',workspace(ipar(8):ipar(8)+n_rows-1)
!      print*,'solution before=',workspace(ipar(9):ipar(9)+n_rows-1)
!      read*
	  CALL A_times_x(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,workspace(ipar(8)), &
      	   workspace(ipar(9)),FineGrid,alpha,beta,gamma)
!      print*,'residual after=',workspace(ipar(8):ipar(8)+n_rows-1)
!      print*,'solution after=',workspace(ipar(9):ipar(9)+n_rows-1)
!      read*
      GOTO 10
   ELSE IF (ipar(1).EQ.2) THEN
	   print*,'Error ipar(1)==2.'
	   STOP
	  GOTO 10
   ELSE IF (ipar(1).EQ.3) THEN
      ! LEFT PRECONDITIONING with the linear low-order matrix.
	  IF (Precond==1) THEN
		  ! Set right hand side for system A_2 u = r
                  CALL LUSOL(n_rows,workspace(ipar(8)),workspace(ipar(9)),alu,jlu,ju)
	  ELSE IF (Precond==2) THEN
		  ! Gauss-Seidel
!		  CALL GaussSeidel(workspace(ipar(9)),workspace(ipar(8)),arrLevels(MG_N_levels)%IterationMatrix)		  		
		PRINT*,'Error: Gauss-Seidel preconditioner not implemented.'
		STOP
	  ELSE IF (Precond==3) THEN
		  ! Multigrid
!      print*,'residual before precond=',workspace(ipar(8):ipar(8)+n_rows-1)
!      print*,'solution before precond=',workspace(ipar(9):ipar(9)+n_rows-1)
!      read*
		CALL MultiGridSolveLinearSystem(workspace(ipar(8)),workspace(ipar(9)),FineGrid%Nx,FineGrid%Ny,&
		        FineGrid%Nz,GhostgridX,GhostGridY,GhostGridZ,cyclet,maxit,reltol,alpha,beta,gamma)
!      print*,'residual after precond=',workspace(ipar(8):ipar(8)+n_rows-1)
!      print*,'solution after precond=',workspace(ipar(9):ipar(9)+n_rows-1)
!      read*
	  END IF
	  GOTO 10
   ELSE IF (ipar(1).EQ.4) THEN
      PRINT*,'Left preconditioning transposed solve not implemented yet.'
	  GOTO 10
   ELSE IF (ipar(1).EQ.5) THEN
      PRINT*,'Right preconditioning not implemented yet.'
	  GOTO 10
   ELSE IF (ipar(1).EQ.6) THEN
      PRINT*,'Right preconditioning transposed solve not implemented yet.'
	  GOTO 10
   ELSE IF (ipar(1).EQ.10) THEN
	  PRINT*,'My own stopping test routine for gmres not implemented yet.'
	  GOTO 10
   ELSE IF (ipar(1).GT.0) THEN
      PRINT*,'ipar(1) returned by gmres solver is an unspecified code.'
   ELSE
!	  PRINT*,'  Iterative solver used ',iterations,' iterations'
!      PRINT*,'the iterative solver terminated with code =', ipar(1)
	  MINITER     = MIN(iterations,MINITER)
	  MAXITER     = MAX(iterations,MAXITER)
	  IF (tstep>1) THEN
		MINITERNOFS = MIN(iterations,MINITERNOFS)
		MAXITERNOFS = MAX(iterations,MAXITERNOFS)
	  ELSE
	    TOTALITERFS = TOTALITERFS + iterations
	  ENDIF
	  TOTALITER = TOTALITER + iterations
	!  IF (3*MAXITER<GMRESmaxiterations) THEN
	!     PRINT*,'WARNING: GMRESmaxiterations (=',GMRESmaxiterations,' is very large compared to the current maximum iterations. This may reduce the computational efficiency significantly.'
!		 PRINT*,'ERROR: Please correct this setting in parameter setup.'
!		 STOP
!	  ENDIF
   ENDIF

END SUBROUTINE iterative_solution
