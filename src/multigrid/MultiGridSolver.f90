!>
!! Solve linear system Ax=b using multigrid method
!! By Allan P. Engsig-Karup.
!<
RECURSIVE SUBROUTINE MultiGridSolver(x,b,k,gam,cyclet,Nxf,Nyf,Nzf,Nxc,Nyc,Nzc,GhostGridX,GhostGridY,GhostGridZ,alpha,beta,gamma)
USE Constants
USE MGLevels
USE HSL_LU
USE GlobalVariables, ONLY: curvilinearONOFF, filename, formattype
IMPLICIT NONE
INTEGER               :: gam, i, k, GhostGridX, GhostGridY, GhostGridZ, nnz, n_rows, alpha, beta, gamma
INTEGER               :: Nxf,Nyf,Nzf,Nxc,Nyc,Nzc, ii, jj, kk, count
CHARACTER(len=1)      :: cyclet
REAL(KIND=long), DIMENSION((Nxf+2*GhostGridX)*(Nyf+2*GhostGridY)*(Nzf+GhostGridZ)), INTENT(INOUT) :: x
REAL(KIND=long), DIMENSION((Nxf+2*GhostGridX)*(Nyf+2*GhostGridY)*(Nzf+GhostGridZ)) :: dxf, r, Ax, tmp
REAL(KIND=long), DIMENSION((Nxf+2*GhostGridX)*(Nyf+2*GhostGridY)*(Nzf+GhostGridZ)), INTENT(IN) :: b
REAL(KIND=long), DIMENSION((Nxc+2*GhostGridX)*(Nyc+2*GhostGridY)*(Nzc+GhostGridZ)) :: dxc, rc

!print*,'k=',k
!print*,'b=',b
!print*,'x=',x
!print*,'***************************************'
!read*

! If we are on the coarsest grid:
IF (k == 1) THEN


if (1==1) THEN
	! FIXME: LET'S NOT DETERMINE nnz AND Nr EVERY TIME WE SOLVE... BUT NEVERMIND FOR NOW
	nnz = SIZE(arrLevels(k)%PreconditioningMatrix%val)
    n_rows = (arrLevels(k)%Nx+2*GhostGridX)*(arrLevels(k)%Ny+2*GhostGridY)*(arrLevels(k)%Nz+GhostGridZ)

    ! Gaussian elimination:
    JOB = 3 ! Solution:
	x = b; ! Use right hand side as input using the vector x which will contain the output (=solution) on return
!print*,'MA41AD:'
!print*,'xIN=',x
!print*,'arrLevels(k)%PreconditioningMatrix%val     = ',arrLevels(k)%PreconditioningMatrix%val
!print*,'arrLevels(k)%PreconditioningMatrix%col_ind = ',arrLevels(k)%PreconditioningMatrix%col_ind
!print*,'arrLevels(k)%PreconditioningMatrix%row_ptr = ',arrLevels(k)%PreconditioningMatrix%row_ptr
print *, 'MultiGrid solver not supported without HSL'
stop
	!CALL MA41AD(JOB, n_rows, nnz, arrLevels(k)%PreconditioningMatrix%row_ptr, &
	 !   arrLevels(k)%PreconditioningMatrix%col_ind, arrLevels(k)%PreconditioningMatrix%val, x, COLSCA, ROWSCA, KEEP, &
	 !   IS_HSL, MAXIS, SS, MAXS, CNTL, ICNTL, INFOHSL, RINFO)
	IF (INFOHSL(1) .LT. 0) THEN ! Check for problems
		PRINT *, 'Problems with MA41, JOB = 3 (Solution), direct LU solver on coarsest level of the Multigrid method.'
    END IF 

!print*,'MA41AD:'
!print*,'xOUT=',x
!read*


ELSE
! Check that smoother works
    DO i = 1 , 200
		CALL GaussSeidel(x,b,arrLevels(k)%IterationMatrix)
    END DO
    
END IF

    ! Choose cycle type:
    IF (cyclet == 'F') gam = 1

ELSE

IF (0==1) THEN
	! Test operators
    ! Restrict -> Prolongate -> check fine grid solution
    tmp = zero
    count = 0
    print*,'Nx=',arrLevels(k)%Nx
    print*,'Ny=',arrLevels(k)%Ny    
    print*,'Nz=',arrLevels(k)%Nz
    DO jj = 1+GhostGridY, arrLevels(k)%Ny+2*GhostGridY      
        DO ii = 1+GhostGridX, arrLevels(k)%Nx+2*GhostGridX
            DO kk = 1+GhostGridZ, arrLevels(k)%Nz
            print*,'(i,j,k) = (',ii,',',jj,',',kk,')'
                count = kk + (jj-1)*(arrLevels(k)%Nx+2*GhostGridX)*(arrLevels(k)%Nz+GhostGridZ)&
                      + (ii-1)*(arrLevels(k)%Nz+GhostGridZ)
                tmp(count) = arrLevels(k)%x(ii,jj)
                tmp(count) = arrLevels(k)%z(kk)
!                print*,'x=',tmp(count)
                print*,'z=',tmp(count)
                print*,'count=',count
            END DO
        END DO
    END DO
    print*,'tmp=',tmp
!	tmp = arrLevels(k)%x
	CALL MGRestriction3DExplicit(tmp,rc,arrLevels(k),arrLevels(k-1),GhostGridX,GhostGridY,GhostGridZ,arrLevels(k)%Nz+GhostGridZ,&
    	arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k-1)%Nz+GhostGridZ,arrLevels(k-1)%Nx+2*GhostGridX,&
        arrLevels(k-1)%Ny+2*GhostGridY )
	filename = "R.bin"
	CALL StoreRealArray3D(rc,arrLevels(k-1)%Nx+2*GhostGridX,arrLevels(k-1)%Ny+2*GhostGridY,arrLevels(k-1)%Nz+GhostGridZ,filename,&
	formattype)
    print*,'rc=',rc
!		rc = one
!        print*,'rc=',rc
        print*,'arrLEvels%x=',arrLevels(k-1)%x
!        rc = arrLevels(k-1)%x
			CALL MGProlongation3DExplicit(tmp,rc,arrLevels(k),arrLevels(k-1),GhostGridX,GhostGridY,GhostGridZ,arrLevels(k)%Nz+GhostGridZ,&
	     arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k-1)%Nz+GhostGridZ,arrLevels(k-1)%Nx+2*GhostGridX,&
         arrLevels(k-1)%Ny+2*GhostGridY )
	filename = "P.bin"
	CALL StoreRealArray3D(tmp,arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k)%Nz+GhostGridZ,filename,formattype)
    print*,'tmp=',tmp
	STOP
END IF

! Check that smoother works
!    DO i = 1 , 200
!		CALL GaussSeidel(x,b,arrLevels(k)%IterationMatrix)
!    END DO
    	
if (1==1) THEN
    ! Pre-smoothing;
		CALL GaussSeidelSpecialPoints(x,b,arrLevels(k)%IterationMatrix,arrLevels(k)%iGhost,arrLevels(k)%mapNp)
!        print*,'x before smoothing=',x
    DO i = 1 , nu(1)
		CALL GaussSeidel(x,b,arrLevels(k)%IterationMatrix)
!        print*,'ighostpoints=',arrLevels(k)%iGhost
!        print*,'mapNp=',arrLevels(k)%mapNp
!        stop
		CALL GaussSeidelSpecialPoints(x,b,arrLevels(k)%IterationMatrix,arrLevels(k)%iGhost,arrLevels(k)%mapNp)
!		CALL GaussSeidelARRAYback(x,b,arrLevels(k),GhostGridX, GhostGridY, GhostGridZ, alpha, beta, gamma ) ! Initial guess is always zero, so call this one instead
!		CALL GaussSeidelARRAY(x,b,arrLevels(k),GhostGridZ, alpha, beta, gamma )
    END DO
!        print*,'x after smoothing=',x

    ! Residual/defect
!	IF (curvilinearONOFF==1) THEN
!		CALL BuildLinearSystemTransformedCurvilinear(arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k)%Nz+GhostGridZ,x,Ax,arrLevels(k),alpha,beta,gamma)
!	ELSE
!		CALL BuildLinearSystem(arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY,arrLevels(k)%Nz+GhostGridZ,x,Ax,arrLevels(k),alpha,beta,gamma)
!	ENDIF
	! Sparse matrix-vector product (SPARSKIT2 lib routine)
!	Ax = zero
!print*,'x=',x
!print*,'b=',b
!print*,'matrx=',arrLevels(k)%IterationMatrix%val
	CALL amux(arrLevels(k)%IterationMatrix%nrow, x, Ax, arrLevels(k)%IterationMatrix%val,arrLevels(k)%IterationMatrix%col_ind,&
    	arrLevels(k)%IterationMatrix%row_ptr) 
	r = b-Ax
!    print*,'r=',r
    
!print*,'========================================================================='
!    read*
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
!    CALL GaussSeidelSpecialPoints(x,b,arrLevels(k)%IterationMatrix,arrLevels(k)%iGhost,arrLevels(k)%mapNp)

    ! Post-smooting;
		CALL GaussSeidelSpecialPoints(x,b,arrLevels(k)%IterationMatrix,arrLevels(k)%iGhost,arrLevels(k)%mapNp)		
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

END IF

END SUBROUTINE MultiGridSolver
