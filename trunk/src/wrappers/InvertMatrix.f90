SUBROUTINE InvertMatrix(M,rank)
! By Allan P. Engsig-Karup.
USE GlobalVariables
IMPLICIT NONE
INTEGER :: rank, lwork, info
REAL(KIND=long), DIMENSION(rank,rank) :: M
REAL(KIND=long), DIMENSION(5*rank**2) :: work
INTEGER, DIMENSION(rank) :: ipiv
rank = SIZE(M,1)
CALL dgetrf(rank,rank,M,rank,ipiv,info)
IF(info/=0)THEN
   PRINT *, 'Problems with L-U of the finite-difference matrix',info
   STOP
END IF
lwork=rank
CALL dgetri(rank,M,rank,ipiv,work,lwork,info)
  IF(info/=0)THEN
  PRINT *, 'Problems with inversion of the finite-difference matrix',info
  STOP
END IF
END SUBROUTINE InvertMatrix
