SUBROUTINE waveMaker3D(i0,i1,x,time,time0,Ea,Pa)
!
! Computes surface elevation and velocity potential at time i, for all points in
! the relaxation xone
!
! By Bo Terp Paulsen, botp@mek.dtu.dk
USE GlobalVariables, ONLY: RandomWave3D, FineGrid
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: i0, i1, nx, itime
REAL(KIND=long), DIMENSION(i1-i0+1) :: x, Ea, Pa
REAL(KIND=long) :: time, time0
!

xtract information about computational grid
!

                ! Nearest neighbour interpolation !TODO: implement higher order
                ! scheme.
                itime = (time-time0)/dt0+1

! Get information about wave gauge grid 
!

! Since both grids are structured compute x and y interpolation stencils
!


RETURN

END SUBROUTINE waveMaker3D
