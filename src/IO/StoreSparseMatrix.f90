SUBROUTINE StoreSparseMatrix(SparseMatrix,filename,formattype)
! By Allan P. Engsig-Karup
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: formattype
TYPE (SparseArray_CSR) :: SparseMatrix
CHARACTER(len=40) :: filename,form
SELECT CASE (formattype)
	CASE (1)
		! Relatively efficient for large data storage
		! information about records stored
        form="unformatted"
	CASE DEFAULT
		! Most efficient for large data storage
		! no information about records stored
        form="binary"
END SELECT
OPEN (unit=22, file=filename,form=form)
WRITE(22) SparseMatrix%nnz, SparseMatrix%row_ptr, SparseMatrix%col_ind, SparseMatrix%val
CLOSE(22)
END SUBROUTINE StoreSparseMatrix
