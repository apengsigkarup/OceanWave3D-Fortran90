SUBROUTINE GammaFunctions(x,n,dir,ftype,gam,param)
!
! Gamma functions for relaxation.
! HBB: added several other options. 3/2019  
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
	CASE (9) ! Allan's sponge filter: 1->0
		gam = one - tmp**param
	CASE (10) ! third order poly. with clamped BCs: 0->1
		gam = three*tmp**2-two*tmp**3
        CASE (11) ! exponential: 0->1
                gam=(exp(tmp**param)-1)/(exp(one)-one)
        CASE (12) ! gamma function, x0=x(n): 0->1
                gam=exp(-((x-x(n))/(x(n)-x(1)))**2/(2*param**2))
        CASE (13) ! 1-gamma function, x0=x(n): 1->0
                gam=1-exp(-((x-x(n))/(x(n)-x(1)))**2/(2*param**2))
        CASE (14) ! exponential: 1->0
                gam=1-(exp(tmp**param)-1)/(exp(one)-one)
END SELECT
IF (dir==-1) THEN ! REVERSE
	gam = gam(n:1:-1)  ! FIXME: not just for x
ENDIF
END SUBROUTINE GammaFunctions
