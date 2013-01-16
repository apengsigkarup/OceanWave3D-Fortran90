SUBROUTINE CSRdiaREORDER(nrow,nnz,a,ja,ia)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: nrow, nnz, i, k, didx
REAL(KIND=long), DIMENSION(nnz) :: a
REAL(KIND=long)					:: da
INTEGER, DIMENSION(nnz)         :: ja
INTEGER, DIMENSION(nrow+1)      :: ia
DO i = 1, nrow
	da = zero
	didx = 0
!print*,'i=',i
!print*,'row=',a(ia(i):ia(i+1)-1)
!print*,'rowidx=',ja(ia(i):ia(i+1)-1)
	DO k = ia(i), ia(i+1)-1 ! global indexes
		IF (ja(k) == i) THEN
			! found diagonal element
			da   = a(k)  ! diagonal value
			didx = k ! global index
			! we can break the k-loop here because we found the diagonal

		! just interchange first row entry with diagonal entry (not necesarily the most cache efficient way, but the entries can still be spread)
		ja(didx)  = ja(ia(i)) ! ia(i) is global index of first element in row
		a(didx)   = a(ia(i))
		ja(ia(i)) = i
		a(ia(i))  = da
			EXIT ! k-loop
		END IF
	END DO
!	IF (didx /= 0) THEN
!	ELSE
!		print*,'Error: diagonal not found. (CSRdiaREORDER)'
!	END IF
!print*,'row=',a(ia(i):ia(i+1)-1)
!print*,'rowidx=',ja(ia(i):ia(i+1)-1)
!read*
END DO
END SUBROUTINE CSRdiaREORDER
