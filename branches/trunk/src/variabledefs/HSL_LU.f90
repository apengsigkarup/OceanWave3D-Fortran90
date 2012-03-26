!>
!! GLOBAL VARIABLES NEEDED FOR THE HARWELL LIBRARY ROUTINES
!<
MODULE HSL_LU
USE Precision
IMPLICIT NONE
INTEGER :: JOB
REAL(KIND=long) :: CNTL(10)
INTEGER :: ICNTL(20), KEEP(50), MAXS, MAXIS, ierr
INTEGER,ALLOCATABLE :: ipivb(:),irns(:),icns(:), repi(:), iw(:)
INTEGER, ALLOCATABLE :: INFOHSL(:), IS_HSL(:)
REAL(KIND=long), ALLOCATABLE :: COLSCA(:), ROWSCA(:), SS(:), RINFO(:)
END MODULE
