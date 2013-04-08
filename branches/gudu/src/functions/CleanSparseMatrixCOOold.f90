SUBROUTINE CleanSparseMatrixCOOold(SM)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE

TYPE (SparseArray_COO) :: SM
INTEGER :: nnz, i, nnz_new
INTEGER, DIMENSION(:), ALLOCATABLE :: rowptr, colind
REAL(KIND=long), DIMENSION(:), ALLOCATABLE :: val
!REAL(KIND=long) :: tol

!tol = 1e-14_long
nnz = SIZE(SM%val)

ALLOCATE(val(nnz),rowptr(nnz),colind(nnz))

nnz_new = 0
DO i = 1, nnz
!	IF (ABS(SM%val(i))<tol) THEN
	IF (SM%val(i).EQ.zero) THEN
		! THROW ELEMENT AWAY, SO DO NOTHING
	ELSE
		! KEEP ELEMENT
  	    nnz_new = nnz_new + 1
		val(nnz_new)    = SM%val(i)
		colind(nnz_new) = SM%col_ind(i)
        rowptr(nnz_new) = SM%row_ptr(i)
	END IF
END DO

DEALLOCATE(SM%val,SM%col_ind,SM%row_ptr)
ALLOCATE(SM%val(nnz_new),SM%col_ind(nnz_new),SM%row_ptr(nnz_new))
SM%val     = val(1:nnz_new)
SM%col_ind = colind(1:nnz_new)
SM%row_ptr = rowptr(1:nnz_new)
SM%nnz     = nnz_new;
DEALLOCATE(val,rowptr,colind)
END SUBROUTINE CleanSparseMatrixCOOold
