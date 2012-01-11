!>
!! Solve linear system Ax=b using multigrid method
!<
RECURSIVE SUBROUTINE MultiGridSolver(x,b,k,gam,cyclet,Nxf,Nyf,Nzf,Nxc,Nyc,Nzc,GhostGridX,GhostGridY,GhostGridZ,alpha,beta,gamma)
USE Constants
USE MGLevels
USE HSL_LU
USE GlobalVariables, ONLY: curvilinearONOFF, filename, formattype
IMPLICIT NONE
INTEGER               :: gam, i, k, GhostGridX, GhostGridY, GhostGridZ, nnz, n_rows, alpha, beta, gamma
INTEGER               :: Nxf,Nyf,Nzf,Nxc,Nyc,Nzc
CHARACTER(len=1)      :: cyclet
REAL(KIND=long), DIMENSION((Nxf+2*GhostGridX)*(Nyf+2*GhostGridY)*(Nzf+GhostGridZ)), INTENT(INOUT) :: x
REAL(KIND=long), DIMENSION((Nxf+2*GhostGridX)*(Nyf+2*GhostGridY)*(Nzf+GhostGridZ)) :: dxf, r, Ax, tmp
REAL(KIND=long), DIMENSION((Nxf+2*GhostGridX)*(Nyf+2*GhostGridY)*(Nzf+GhostGridZ)), INTENT(IN) :: b
REAL(KIND=long), DIMENSION((Nxc+2*GhostGridX)*(Nyc+2*GhostGridY)*(Nzc+GhostGridZ)) :: dxc, rc

! If we are on the coarsest grid:
IF (k == 1) THEN

	! FIXME: LET'S NOT DETERMINE nnz AND Nr EVERY TIME WE SOLVE... BUT NEVERMIND FOR NOW
	nnz = SIZE(arrLevels(k)%PreconditioningMatrix%val)
    n_rows = (arrLevels(k)%Nx+2*GhostGridX)*(arrLevels(k)%Ny+2*GhostGridY)*(arrLevels(k)%Nz+GhostGridZ)

    ! Gaussian elimination:
    JOB = 3 ! Solution:
	x = b; ! Use right hand side as input using the vector x which will contain the output (=solution) on return
	CALL MA41AD(JOB, n_rows, nnz, arrLevels(k)%PreconditioningMatrix%row_ptr, &
	    arrLevels(k)%PreconditioningMatrix%col_ind, arrLevels(k)%PreconditioningMatrix%val, x, COLSCA, ROWSCA, KEEP, &
	    IS_HSL, MAXIS, SS, MAXS, CNTL, ICNTL, INFOHSL, RINFO)
	IF (INFOHSL(1) .LT. 0) THEN ! Check for problems
		PRINT *, 'Problems with MA41, JOB = 3 (Solution), direct LU solver on coarsest level of the Multigrid method.'
    END IF 

    ! Choose cycle type:
    IF (cyclet == 'F') gam = 1

ELSE

	! Test operators
!	tmp = one
!	CALL MGRestriction3DExplicit(tmp,rc,arrLevels(k),arrLevels(k-1),GhostGridX,GhostGridY,GhostGridZ,arrLevels(k)%Nz+GhostGridZ,&
 !   	arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k-1)%Nz+GhostGridZ,arrLevels(k-1)%Nx+2*GhostGridX,&
 !       arrLevels(k-1)%Ny+2*GhostGridY )
!	filename = "tmp.bin"
!	CALL StoreRealArray3D(rc,arrLevels(k-1)%Nx+2*GhostGridX,arrLevels(k-1)%Ny+2*GhostGridY,arrLevels(k-1)%Nz+GhostGridZ,filename,&
!	formattype)
!		rc = one
!			CALL MGProlongation3DExplicit(tmp,rc,arrLevels(k),arrLevels(k-1),GhostGridX,GhostGridY,GhostGridZ,arrLevels(k)%Nz+GhostGridZ,&
!	     arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k-1)%Nz+GhostGridZ,arrLevels(k-1)%Nx+2*GhostGridX,&
 !        arrLevels(k-1)%Ny+2*GhostGridY )
!	filename = "tmp.bin"
!	CALL StoreRealArray3D(tmp,arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k)%Nz+GhostGridZ,filename,formattype)
!	STOP
	
    ! Pre-smoothing;
    DO i = 1 , nu(1)
		CALL GaussSeidel(x,b,arrLevels(k)%IterationMatrix)
		CALL GaussSeidelSpecialPoints(x,b,arrLevels(k)%IterationMatrix,arrLevels(k)%iGhost,arrLevels(k)%mapNp)
!		CALL GaussSeidelARRAYback(x,b,arrLevels(k),GhostGridX, GhostGridY, GhostGridZ, alpha, beta, gamma ) ! Initial guess is always zero, so call this one instead
!		CALL GaussSeidelARRAY(x,b,arrLevels(k),GhostGridZ, alpha, beta, gamma )
    END DO

    ! Residual/defect
!	IF (curvilinearONOFF==1) THEN
!		CALL BuildLinearSystemTransformedCurvilinear(arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k)%Nz+GhostGridZ,x,Ax,arrLevels(k),alpha,beta,gamma)
!	ELSE
!		CALL BuildLinearSystem(arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k)%Nz+GhostGridZ,x,Ax,arrLevels(k),alpha,beta,gamma)
!	ENDIF
	! Sparse matrix-vector product (SPARSKIT2 lib routine)
!	Ax = zero
	CALL amux(arrLevels(k)%IterationMatrix%nrow, x, Ax, arrLevels(k)%IterationMatrix%val,arrLevels(k)%IterationMatrix%col_ind,&
    	arrLevels(k)%IterationMatrix%row_ptr) 
	r = b-Ax
    ! Source on coarser grid:
	! We restrict the residual for the computational domain, but what about the ghost layer?
!	rc = zero
	CALL MGRestriction3DExplicit(r,rc,arrLevels(k),arrLevels(k-1),GhostGridX,GhostGridY,GhostGridZ,arrLevels(k)%Nz+GhostGridZ,&
    	arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k-1)%Nz+GhostGridZ,arrLevels(k-1)%Nx+2*GhostGridX,&
        arrLevels(k-1)%Ny+2*GhostGridY )
    dxc = zero
    DO i = 1, gam
        ! Recursive call to MultiGrid solver 
		CALL MultiGridSolver(dxc,rc,k-1,gam,cyclet,Nxc,Nyc,Nzc,arrLevels(k-1)%Nx,arrLevels(k-1)%Ny,arrLevels(k-1)%Nz,&
		     GhostGridX, GhostGridY, GhostGridZ, alpha, beta, gamma)
        IF (k==MG_N_levels .AND. cyclet == 'W') EXIT 
        IF (k==MG_N_levels .AND. cyclet == 'F') EXIT 
    END DO

    ! Correct the guess (u = u* + u'):
!	dxf = zero
	CALL MGProlongation3DExplicit(dxf,dxc,arrLevels(k),arrLevels(k-1),GhostGridX,GhostGridY,GhostGridZ,arrLevels(k)%Nz+GhostGridZ,&
	     arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k-1)%Nz+GhostGridZ,arrLevels(k-1)%Nx+2*GhostGridX,&
         arrLevels(k-1)%Ny+2*GhostGridY )
	x = x + dxf
    ! impose BCs
    CALL GaussSeidelSpecialPoints(x,b,arrLevels(k)%IterationMatrix,arrLevels(k)%iGhost,arrLevels(k)%mapNp)

    ! Post-smooting;
    DO i = 1 , nu(2)
		CALL GaussSeidel(x,b,arrLevels(k)%IterationMatrix)
!	filename = "u.bin"
!	CALL StoreRealArray3D(x,arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k)%Nz+GhostGridZ,filename,formattype)
!	STOP
		CALL GaussSeidelSpecialPoints(x,b,arrLevels(k)%IterationMatrix,arrLevels(k)%iGhost,arrLevels(k)%mapNp)		
!	filename = "u2.bin"
!	CALL StoreRealArray3D(x,arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k)%Nz+GhostGridZ,filename,formattype)
!	STOP
!		CALL GaussSeidelReverse(x,b,arrLevels(k)%IterationMatrix)
!		CALL GaussSeidelARRAY(x,b,arrLevels(k),GhostGridX,GhostGridY,GhostGridZ, alpha, beta, gamma )
    END DO

END IF

END SUBROUTINE MultiGridSolver
