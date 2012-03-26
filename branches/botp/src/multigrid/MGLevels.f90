MODULE MGLevels
!
! This module defines the multigrid data structure
!
! Allan P. Engsig-Karup, 23 Aug 2007.
!
USE Precision
USE DataTypes

! For multilevels allocate an array of levels
INTEGER :: MG_N_levels, nu(2)
TYPE (Level_def), DIMENSION(:), ALLOCATABLE :: arrLevels
END MODULE MGLevels
