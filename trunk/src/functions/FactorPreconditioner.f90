SUBROUTINE FactorPreconditioner(PreconditioningMatrix, Nfs)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: Nfs
TYPE (SparseArray_COO) :: PreconditioningMatrix
CALL LUfactor(PreconditioningMatrix, Nfs)
END SUBROUTINE FactorPreconditioner
