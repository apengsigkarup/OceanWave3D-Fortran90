SUBROUTINE GammaFunctions(x,n,dir,ftype,gam,param)
!
! Gamma functions for relaxation
!
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
IMPLICIT NONE
INTEGER::n, ftype
REAL(KIND=long), DIMENSION(n) :: x, gam, tmp
REAL(KIND=long) :: param
INTEGER :: dir
IF (n>1) THEN
	tmp = (x - x(1))/(x(n)-x(1))
ELSE
	tmp = x
ENDIF
SELECT CASE (ABS(ftype))
	CASE (9) ! Allan's sponge filter
		gam = one - tmp**param
	CASE (10) ! third order poly. with clamped BCs
		gam = three*tmp**2-two*tmp**3
        CASE (11)
                gam=(exp(tmp**param)-1)/(exp(one)-one)
END SELECT
IF (dir==-1) THEN ! REVERSE
	gam = gam(n:1:-1)  ! FIXME: not just for x
ENDIF
END SUBROUTINE GammaFunctions
