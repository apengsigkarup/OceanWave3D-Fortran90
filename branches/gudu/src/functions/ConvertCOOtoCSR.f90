SUBROUTINE ConvertCOOtoCSR(SM,tmpSM)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (SparseArray_CSR) :: tmpSM
TYPE (SparseArray_COO) :: SM
INTEGER :: nnz, nrow !ierr, nrow, job
!INTEGER, DIMENSION(:), ALLOCATABLE :: iwk, indu

! Convert to CSR format and back to COO format to remove duplicate entries
nrow       = SM%nrow
tmpSM%nrow = nrow
nnz        = SM%nnz
tmpSM%nnz  = nnz

ALLOCATE( tmpSM%val(tmpSM%nnz)        )
ALLOCATE( tmpSM%col_ind(tmpSM%nnz)    )
ALLOCATE( tmpSM%row_ptr(tmpSM%nrow+1) )

CALL coocsr(nrow,nnz,SM%val,SM%row_ptr,SM%col_ind,tmpSM%val,tmpSM%col_ind,tmpSM%row_ptr)

! Release COO matrix
DEALLOCATE( SM%val, SM%row_ptr, SM%col_ind)

END SUBROUTINE ConvertCOOtoCSR
