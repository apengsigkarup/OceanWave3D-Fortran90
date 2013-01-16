FUNCTION disper(g,T,h,n) RESULT(kh)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: n
REAL(KIND=long) :: g, T, h, kh, w
! determine the dimensionless dispersion parameter kh
! for wave equations using the dispersion relation.
w  = 2*pi/T
kh = w**2*h/g
kh = khsolve(g,w,kh,h,n)

CONTAINS

FUNCTION khsolve(g,w,kh,h,n) RESULT(khout)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: i, n
REAL(KIND=long) :: g, T, h, kh, w, khout
! to be used within the function disper for iterative solution of dispersion relation
DO i = 1,n
    kh = SQRT(w**2/g*h*kh*COSH(kh)/SINH(kh)) ! COTH(x) = COSH(x)/SINH(x)
END DO
khout = kh
END FUNCTION

END FUNCTION
