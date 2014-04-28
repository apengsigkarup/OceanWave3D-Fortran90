SUBROUTINE ConvertCSRtoCOO(tmpSM,SM)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (SparseArray_CSR) :: tmpSM
TYPE (SparseArray_COO) :: SM
INTEGER :: nnz, nrow, ierr, job
!INTEGER, DIMENSION(:), ALLOCATABLE :: iwk, indu

SM%nnz = tmpSM%nnz

ALLOCATE( SM%val(tmpSM%nnz)        )
ALLOCATE( SM%col_ind(tmpSM%nnz)    )
ALLOCATE( SM%row_ptr(tmpSM%nnz)    )

! FIXME: job=3 used here, but can be made more efficient, by only redirecting pointers?
job = 3
CALL csrcoo(tmpSM%nrow,job,tmpSM%row_ptr(tmpSM%nrow+1) - tmpSM%row_ptr(1),&
			tmpSM%val,tmpSM%col_ind,tmpSM%row_ptr,&
			tmpSM%row_ptr(tmpSM%nrow+1) - tmpSM%row_ptr(1),&
			SM%val,SM%row_ptr,SM%col_ind,ierr)

! Release CSR matrix
DEALLOCATE( tmpSM%val, tmpSM%row_ptr, tmpSM%col_ind)

END SUBROUTINE ConvertCSRtoCOO
