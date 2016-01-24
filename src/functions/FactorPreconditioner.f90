SUBROUTINE FactorPreconditioner(PreconditioningMatrix, Nfs, fileop)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: Nfs, fileop
TYPE (SparseArray_COO) :: PreconditioningMatrix
CALL LUfactor(PreconditioningMatrix, Nfs, fileop)
END SUBROUTINE FactorPreconditioner
