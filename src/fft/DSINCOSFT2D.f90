      SUBROUTINE DSINCOSFT2D(Y,n1_max,N1,n2,wrk)
!
! This routine is the same as DCOSFT2D except that in the first
! direction a SINFT is performed.
!
      IMPLICIT NONE
      INTEGER n1_max, n1, n2
      REAL*8  Y(N1_max,n2), wrk(n2)
      INTEGER j, k

      DO k=1,n2
         CALL dsinft(y(1,k),n1-1)
      END DO
      DO j=1,n1
         DO k=1,n2
            wrk(k)=y(j,k)
         END DO
         CALL dcosft(wrk,n2-1)
         DO k=1,n2
            y(j,k)=wrk(k)
         END DO
      END DO

      RETURN
      END SUBROUTINE DSINCOSFT2D
