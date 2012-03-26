SUBROUTINE LUfactor(SparseMatrix, Nr)
!
! LU factor matrix given in sparse format.
!
USE Precision
USE Constants
USE DataTypes
USE HSL_LU
IMPLICIT NONE
TYPE (SparseArray_COO) :: SparseMatrix
INTEGER :: Nr, Nxz, nnz, STAT

nnz   = SIZE(SparseMatrix%val)
nxz   = Nr ! total number of grid points

! Analysis
JOB = 1

CALL MA41ID(CNTL,ICNTL,KEEP) ! Sets default parameter values

! Initial parameters and workspace for MA41.
MAXIS = 2*(4*nnz + 14*nxz + 1) ! FIXME: BETTER FORMULA HERE?
ALLOCATE (COLSCA(nxz), ROWSCA(nxz), &
         IS_HSL(MAXIS), INFOHSL(20), RINFO(20), &
         irns(nnz), icns(nnz), STAT=STAT )
CALL CheckError(STAT,4)
! GD: Pointer should be allocated during a call in my compiler
! FIXME: More clever way to do this?  Use a temporary array to prevent from allocation-deallocation ?
ALLOCATE (SS(0)) ! GD

!print*,'Fill-in parameter for controlling fill-in in MA1D:'
!print*,'  CNTL(1)=',CNTL(1)
!CNTL(1) = 1.0_long
!print*,'  redefined, CNTL(1)=',CNTL(1)
!print*,''

CALL MA41AD(JOB, Nr, nnz, SparseMatrix%row_ptr, SparseMatrix%col_ind, SparseMatrix%val, one, COLSCA, ROWSCA, KEEP, &
       IS_HSL, MAXIS, SS, MAXS, CNTL, ICNTL, INFOHSL, RINFO)

DEALLOCATE(SS) ! GD: Deallocation of the corresponding pointer

IF (INFOHSL(1) .LT. 0) THEN ! Error
  PRINT *, 'Problems with MA41, JOB = 1 (Analysis)'
  PRINT *, '   INFOHSL(1) = ', INFOHSL(1)
  STOP
END IF

IF (INFOHSL(1) .GT. 0) THEN ! Warning
  PRINT *, 'Warning from MA41, JOB = 1 (Analysis)'
  PRINT *, '   INFOHSL(1) = ', INFOHSL(1)
  PRINT *, '   INFOHSL(2) = ', INFOHSL(2)
END IF

WRITE (*,'(/A)') '  Factorization of preconditioner:'

MAXS = INFOHSL(8) ! Minimum size allowable from analysis
PRINT *, '     Real workspace for LU factors:  MAXS = ', MAXS
ALLOCATE ( SS(MAXS), STAT=STAT)
CALL CheckError(STAT,5)

JOB = 2 ! Factorization

CALL MA41AD(JOB, Nr, nnz, SparseMatrix%row_ptr, SparseMatrix%col_ind, SparseMatrix%val, one, COLSCA, ROWSCA, KEEP, &
   IS_HSL, MAXIS, SS, MAXS, CNTL, ICNTL, INFOHSL, RINFO)

IF (INFOHSL(1) .LT. 0) THEN
   WRITE (*,'(A)') 'Problems with MA41, JOB = 2 (Numerical factorization)'
   STOP
END IF

WRITE (*,'(A,I6)') '     Size of real space used to store the LU factors = ', INFOHSL(9)
WRITE (*,'(A,I6)') '     Size of integer space used to store the LU factors = ', INFOHSL(10)
WRITE (*,'(A,I6/)') '     Number of off-diagonal pivots = ', INFOHSL(12)

JOB = 3 ! Solution

END SUBROUTINE LUfactor
