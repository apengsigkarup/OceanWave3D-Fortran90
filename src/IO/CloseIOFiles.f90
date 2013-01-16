SUBROUTINE CloseIOFiles
! By Allan P. Engsig-Karup.
USE GlobalVariables
IMPLICIT NONE
INTEGER i
!  DO i=1,SIZE(FILEIP)
!     CLOSE (FILEIP(i))
!  END DO
!  DO I=1,SIZE(FILEOP)
!     CLOSE (FILEOP(I))
!  END DO
! GD: SIZE test not relevant on non-assigned array, here size(FILEIP)=4
! Transform on allocated array and test on allocation or if fileip(i)=i
! and fileop(i)=i+4, following test is working...
  DO i=1,MAXLOC(FILEIP,1)
     CLOSE (FILEIP(i))
  END DO
  DO I=1,MAXLOC(FILEOP,1)
     CLOSE (FILEOP(I))
  END DO
END SUBROUTINE CloseIOFiles
