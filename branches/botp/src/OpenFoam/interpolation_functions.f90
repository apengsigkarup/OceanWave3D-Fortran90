MODULE interpolation_functions
! Module containing functions necessary for 3DInterpolation
CONTAINS
INTEGER FUNCTION fact(n)
      ! Retruns the factorial of n
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      INTEGER p,i
      p = 1   
      do i = 1, n
	    p = p * i
      end do
      fact = p
END FUNCTION fact

! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
function inv(A) result(Ainv)
IMPLICIT NONE
  INTEGER, PARAMETER :: dp = selected_real_kind(15, 307)
  real(KIND=dp), dimension(:,:), intent(in) :: A
  real(KIND=dp), dimension(size(A,1),size(A,2)) :: Ainv

  real(KIND=dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: N,M,info,i
  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A(:,:)

  M = size(A,1)
  N = size(A,2)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(M, N, Ainv, M, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(M, Ainv, M, ipiv, work, M, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv 
END MODULE interpolation_functions
