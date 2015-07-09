!>
!! This module contains definitions to control precision in code.
!! By Allan P. Engsig-Karup.
!<
MODULE Precision
IMPLICIT NONE
! Precision of the internal variables (i.e. not for external libraries)
INTEGER, PARAMETER :: sp = selected_real_kind(6, 37)
INTEGER, PARAMETER :: dp = selected_real_kind(15, 307)
INTEGER, PARAMETER :: qp = selected_real_kind(33, 4931)
INTEGER, PARAMETER :: long=SELECTED_REAL_KIND(15, 307)
END MODULE Precision
