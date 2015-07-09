      SUBROUTINE dcosft2d(Y,n1_max,N1,n2,wrk)
!
! This is a fast cosine transform routine for a 2-d function Y(n1,n2).
! n1-1 and n2-1 must both be powers of 2 and >= 2.
! This function is it's own transform, but the inverse must be
! multiplied by 4/((n1-1)*(n2-1)).   Since the data and transform are
! all real, we simply FCT in x along the rows, then in y along the
! columns using the 1-D routine.
!
      IMPLICIT NONE
      INTEGER n1_max, n1, n2
      REAL*8  Y(N1_max,n2), wrk(n2)
      INTEGER j, k

      DO k=1,n2
         CALL dcosft(y(1,k),n1-1)
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
      END SUBROUTINE DCOSFT2D
