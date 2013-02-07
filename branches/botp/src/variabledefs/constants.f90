!>
!! This module contains definitions for various useful constants.
!! By Allan P. Engsig-Karup.
!<
MODULE Constants
USE Precision
IMPLICIT NONE
! Constant parameters
REAL(KIND=long), PARAMETER :: zero=0.0_long, half=0.5_long,          &
       one=1.0_long,  two=2.0_long, three=3._long,  four=4.0_long,   &
       five=5.0_long, six=6.0_long, seven=7.0_long, eight=8.0_long,  &
       nine=9.0_long, ten=10.0_long
REAL(KIND=long) :: pi, deg2rad
END MODULE
