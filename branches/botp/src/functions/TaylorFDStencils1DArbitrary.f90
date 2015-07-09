SUBROUTINE TaylorFDStencils1DArbitrary(alp,bet,order,Stencils,x)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
INTEGER :: r, alp, bet, m, n, order, count
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: Mat
REAL(KIND=long), DIMENSION(alp+bet+1) :: Stencils, x
r = alp + bet + 1; ! rank of points in stencil
ALLOCATE(Mat(r,order+1))
CALL weights(x(alp+1),x,r-1,r-1,order,Mat)
Stencils = Mat(1:r,order+1)
DEALLOCATE(Mat)
END SUBROUTINE TaylorFDStencils1DArbitrary
