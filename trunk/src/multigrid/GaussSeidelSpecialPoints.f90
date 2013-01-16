SUBROUTINE GaussSeidelSpecialPoints(x,b,SMCSR,iMap,Np)
! By Allan P. Engsig-Karup.
! Update the global points specified in iMap
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: nnz, nrow, Gidx, i, Np, j
TYPE (SparseArray_CSR) :: SMCSR
REAL(KIND=long), DIMENSION(SMCSR%nnz), INTENT(INOUT) :: x
REAL(KIND=long), DIMENSION(SMCSR%nnz), INTENT(IN) :: b
REAL(KIND=long) :: diagcoef
INTEGER :: IMap(Np)
nrow = SMCSR%nrow
nnz  = SMCSR%nnz
DO j = 1, Np
	Gidx = IMap(j)
	i = Gidx
	diagcoef = SMCSR%val(SMCSR%row_ptr(i))
!	IF (SMCSR%row_ptr(i+1) - SMCSR%row_ptr(i) == 1) THEN
!		x(Gidx) = b(Gidx)/diagcoef
!	ELSE
		x(Gidx) = ( b(Gidx) - DOT_PRODUCT( x( SMCSR%col_ind( SMCSR%row_ptr(i)+1 : SMCSR%row_ptr(i+1)-1 ) ),&
						       SMCSR%val( SMCSR%row_ptr(i)+1 : SMCSR%row_ptr(i+1)-1 ) ) )/diagcoef
!	END IF
END DO
END SUBROUTINE GaussSeidelSpecialPoints
