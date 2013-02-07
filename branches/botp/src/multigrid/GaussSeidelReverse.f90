SUBROUTINE GaussSeidelReverse(x,b,SMCSR)
! By Allan P. Engsig-Karup.
! FIXME: always do free surface first??
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: nnz, nrow, Gidx, i
TYPE (SparseArray_CSR) :: SMCSR
REAL(KIND=long), DIMENSION(SMCSR%nnz), INTENT(INOUT) :: x
REAL(KIND=long), DIMENSION(SMCSR%nnz), INTENT(IN) :: b
REAL(KIND=long) :: diagcoef
nrow = SMCSR%nrow
nnz  = SMCSR%nnz
DO Gidx = nrow,1,-1
	i = Gidx
!	print*,'i=',i
!	print*,'SMCSR%row_ptr(i)=',SMCSR%row_ptr(i)
	diagcoef = SMCSR%val(SMCSR%row_ptr(i))
	IF (SMCSR%row_ptr(i+1) - SMCSR%row_ptr(i) == 1) THEN
		x(Gidx) = b(Gidx)/diagcoef
	ELSE
		x(Gidx) = ( b(Gidx) - DOT_PRODUCT( x( SMCSR%col_ind( SMCSR%row_ptr(i)+1 : SMCSR%row_ptr(i+1)-1 ) ),&
						       SMCSR%val( SMCSR%row_ptr(i)+1 : SMCSR%row_ptr(i+1)-1 ) ) )/diagcoef
	END IF
END DO
END SUBROUTINE GaussSeidelReverse
