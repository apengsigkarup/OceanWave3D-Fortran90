SUBROUTINE LUfactor(SparseMatrix, Nr, fileop)
!
! LU factor matrix given in sparse format.
!
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
USE SPK
USE GlobalVariables, ONLY: GMRESmaxiterations
IMPLICIT NONE
TYPE (SparseArray_CSR) :: SparseMatrix
INTEGER :: Nr, Nxz, nnz, STAT,i,j, fileop

REAL(8) :: droptol
INTEGER :: lfil, IWK,workspaceSize

nnz   = SIZE(SparseMatrix%val)
nxz   = Nr ! total number of grid points


lfil = 5
IWK  = nnz + lfil + 1 
droptol = 0.005
workspaceSize =(Nr+3)*(GMRESmaxiterations+2) + (GMRESmaxiterations+1) * GMRESmaxiterations/2

ALLOCATE(alu(IWK),jlu(IWK),ju(Nr),w(workspaceSize),jw(3*Nr))
print *, 'workspaceSize = ', workspaceSize


CALL ILUT(Nr,SparseMatrix%val,SparseMatrix%col_ind,SparseMatrix%row_ptr,lfil,droptol,alu,jlu,ju,IWK,w,jw,ierr)

END SUBROUTINE LUfactor
