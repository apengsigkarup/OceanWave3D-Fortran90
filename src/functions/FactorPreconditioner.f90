SUBROUTINE FactorPreconditioner(PreconditioningMatrix, Nfs)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: Nfs
TYPE (SparseArray_COO) :: PreconditioningMatrix
! Original preconditioner using HSL-MA41
!
!CALL LUfactor(PreconditioningMatrix, Nfs) 

! Alternative preconditioner using SPARSKIT
!
CALL LUfactor_spk(PreconditioningMatrix, Nfs)
END SUBROUTINE FactorPreconditioner
