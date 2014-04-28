SUBROUTINE Zero_Corners(A,B,Nx,Ny)
! By Allan P. Engsig-Karup.
USE PRECISION
USE Constants
IMPLICIT NONE
INTEGER :: Nx,Ny
REAL(KIND=long), DIMENSION(Nx,Ny) :: A,B
!
A(1,1)   = zero
A(Nx,1)  = zero
A(1,Ny)  = zero
A(Nx,Ny) = zero
!
B(1,1)   = zero
B(Nx,1)  = zero
B(1,Ny)  = zero
B(Nx,Ny) = zero

END SUBROUTINE Zero_Corners
