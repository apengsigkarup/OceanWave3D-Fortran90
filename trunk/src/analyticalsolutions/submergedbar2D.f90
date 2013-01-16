SUBROUTINE SubmergedBar_2D(FineGrid,x0,GhostGridX)
!
! Define bottom and bottom gradients for submerged bar test (2D)
!
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: FineGrid
REAL(KIND=long)  :: x,x0
INTEGER :: Nx, i, j, GhostGridX

Nx = FineGrid%Nx+2*GhostGridX
j=1

DO i = 1 , Nx
  	x = FineGrid%x(i,j)-x0
    IF (x<=6.0_long) THEN
        FineGrid%h(i,j)   = 0.4_long
        FineGrid%hx(i,j)  = zero
        FineGrid%hxx(i,j) = zero
    ELSE IF (x<=12.0_long) THEN
        FineGrid%h(i,j)   = 0.4_long - 1.0_long/20.0_long*(x-6.0_long)
        FineGrid%hx(i,j)  = -1.0_long/20.0_long
        FineGrid%hxx(i,j) = zero
    ELSE IF (x<=14.0_long) THEN
        FineGrid%h(i,j)   = 0.1_long
        FineGrid%hx(i,j)  = zero
        FineGrid%hxx(i,j) = zero
    ELSE IF (x<=17.0_long) THEN
        FineGrid%h(i,j)   = 0.1_long + 1.0_long/10.0_long*(x-14.0_long)
        FineGrid%hx(i,j)  = 1.0_long/10.0_long
        FineGrid%hxx(i,j) = zero
    ELSE
        FineGrid%h(i,j)   = 0.4_long
        FineGrid%hx(i,j)  = zero
        FineGrid%hxx(i,j) = zero
    ENDIF
END DO


END SUBROUTINE SubmergedBar_2D
