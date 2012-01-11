FUNCTION distdot(n,dx,incx,dy,incy) RESULT(out)
! VECTOR DOT PRODUCT ROUTINE. NEEDED BY SPARSKIT GMRES. WRITTEN BY Allan P. Engsig-Karup.
USE Precision
IMPLICIT NONE
INTEGER i,incx,incy,n
REAL(KIND=long) :: dx(n),dy(n),out
out = DOT_PRODUCT(dx,dy)
!print*,'distdot result = ',out
END FUNCTION distdot
