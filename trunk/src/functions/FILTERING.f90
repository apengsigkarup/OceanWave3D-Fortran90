SUBROUTINE FILTERING(Nx,Ny,U,filterNP,filterALPHA,filtercoefficients,oddeven)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: Nx, Ny, i, j,k, filterALPHA, filterNP, oddeven
REAL(KIND=long) :: U(Nx,Ny), tmp(Nx,Ny), filtercoefficients(filterNP), idx(filterNP)

DO k = 1 , 2
  SELECT CASE (MOD(oddeven+k,2))
	CASE (0)
     IF (Nx>1) THEN
	  tmp = U
!  ! LEFT
!  DO j = 1, Ny
!    DO i = 1, filterALPHA
!       U(i,j) = DOT_PRODUCT(tmp(i+(-filterALPHA:filterALPHA),j),filtercoefficients)
!    END DO
!  END DO
	  ! INTERIOR
	  DO j = 1, Ny
		DO i = 1+filterALPHA, Nx-filterALPHA
			U(i,j) = DOT_PRODUCT(tmp(i-filterALPHA:i+filterALPHA,j),filtercoefficients)
		END DO
	  END DO
	 ENDIF
	CASE (1)
	  IF (Ny>1) THEN
		tmp = U
		! INTERIOR
		DO j = 1+filterALPHA, Ny-filterALPHA
			DO i = 1, Nx
				U(i,j) = DOT_PRODUCT(tmp(i,j-filterALPHA:j+filterALPHA),filtercoefficients)
			END DO
		END DO
	  ENDIF
  END SELECT
END DO
!PRINT*,'  SG-FILTER applied.'
END SUBROUTINE
