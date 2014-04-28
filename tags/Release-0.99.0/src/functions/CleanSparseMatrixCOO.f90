SUBROUTINE CleanSparseMatrixCOO(SM)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (SparseArray_CSR) :: tmpSM
TYPE (SparseArray_COO) :: SM
INTEGER :: nnz, ierr, nrow, job
INTEGER, DIMENSION(:), ALLOCATABLE :: iwk, indu

! Convert to CSR format and back to COO format to remove duplicate entries
nrow       = SM%nrow
tmpSM%nrow = nrow
nnz        = SM%nnz
tmpSM%nnz  = nnz

ALLOCATE( tmpSM%val(tmpSM%nnz)        )
ALLOCATE( tmpSM%col_ind(tmpSM%nnz)    )
ALLOCATE( tmpSM%row_ptr(tmpSM%nrow+1) )

CALL coocsr(nrow,nnz,SM%val,SM%row_ptr,SM%col_ind,tmpSM%val,tmpSM%col_ind,tmpSM%row_ptr)
!print * ,'precond matrix converted from COO to CSR format.'
!print * ,'nnz = ', tmpsM%row_ptr(tmpSM%nrow+1) - tmpSM%row_ptr(1)

! Clean duplicates

ALLOCATE( iwk(tmpSM%nrow+1) )
ALLOCATE( indu(tmpSM%nrow)  )

CALL clncsr(1,1,tmpSM%nrow,tmpSM%val,tmpSM%col_ind,tmpSM%row_ptr,indu,iwk)

!print * ,'precond matrix cleaned for duplicates.'
!print * ,'nnz = ', tmpSM%row_ptr(tmpSM%nrow+1) - tmpSM%row_ptr(1)

! Convert back to COO format
DEALLOCATE( SM%val, SM%row_ptr, SM%col_ind)
SM%nnz = tmpSM%row_ptr(tmpSM%nrow+1) - tmpSM%row_ptr(1)

ALLOCATE( SM%val(SM%nnz) )
ALLOCATE( SM%col_ind(SM%nnz) )
ALLOCATE( SM%row_ptr(SM%nnz) )

! FIXME: job=3 used here, but can be made more efficient, by only redirecting pointers?
job = 3
CALL csrcoo(tmpSM%nrow,job,tmpSM%row_ptr(tmpSM%nrow+1) - tmpSM%row_ptr(1),&
			tmpSM%val,tmpSM%col_ind,tmpSM%row_ptr,&
			tmpSM%row_ptr(tmpSM%nrow+1) - tmpSM%row_ptr(1),&
			SM%val,SM%row_ptr,SM%col_ind,ierr)

!print * ,'precond matrix converted from CSR to COO format.'

!filename = "SparseMatrix.bin"
!CALL StoreSparseMatrix(FineGrid%PreconditioningMatrix,filename,formattype)
!WRITE(*,'(A)') 'Sparse matrix stored.'
DEALLOCATE(iwk, indu)
DEALLOCATE(tmpSM%val, tmpSM%row_ptr, tmpSM%col_ind)

END SUBROUTINE CleanSparseMatrixCOO
