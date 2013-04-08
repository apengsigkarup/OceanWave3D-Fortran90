SUBROUTINE DetermineGenericStencils(DiffStencils,kappa)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER kappa
TYPE (Diff_def) :: DiffStencils
ALLOCATE(DiffStencils%StencilG(2*kappa+1,2*kappa+1,2))
CALL DetermineGenericStencilsUniform(DiffStencils%StencilG,kappa)
END SUBROUTINE DetermineGenericStencils
