SUBROUTINE FILTERING(Nx,Ny,U,filterNP,filterALPHA,filtercoefficients,oddeven,GhostGridX,GhostGridY,filtercoefficients2)
!
! By Allan P. Engsig-Karup.  Modified by H. Bingham Feb. 2013.
!
! Filter all the non-ghost point grid points first along lines in x and then along lines 
! in y.  
!
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: Nx, Ny, i, j,k, filterALPHA, filterNP, oddeven, GhostGridX, GhostGridY
REAL(KIND=long) :: U(Nx,Ny), tmp(Nx,Ny), filtercoefficients(filterNP), idx(filterNP)
REAL(long) :: filtercoefficients2(filterNP,filterALPHA)
!
! x-direction smoothing. 
!
IF (Nx>1) THEN
   tmp = U
   DO j = 1+GhostGridY, Ny-GhostGridY
      ! Left boundary points
      DO i = 1,filterALPHA
         U(i+GhostGridX,j) = DOT_PRODUCT(tmp(1+GhostGridX:GhostGridX+filterNP,j),&
              filtercoefficients2(1:filterNP,i))
      END DO
      ! INTERIOR
      DO i = 1+GhostGridX+filterALPHA, Nx-GhostGridX-filterALPHA
         U(i,j) = DOT_PRODUCT(tmp(i-filterALPHA:i+filterALPHA,j),filtercoefficients)
      END DO
      ! Right boundary points
      DO i = 1,filterALPHA 
         U(Nx-GhostGridX-filterALPHA+i,j) = DOT_PRODUCT(tmp(Nx-GhostGridX-(filterNP)+1:Nx-GhostGridX,j),&
              filtercoefficients2(filterNP:1:-1,filterALPHA+1-i))
      END DO
   END DO
ENDIF
!
! y-direction smoothing
!
IF (Ny>1) THEN
   tmp = U
   DO i = 1+GhostGridX, Nx-GhostGridX
      ! Bottom boundary points
      DO j = 1,filterALPHA
         U(i,j+GhostGridY) = DOT_PRODUCT(tmp(i,1+GhostGridY:filterNP+GhostGridY),&
              filtercoefficients2(1:filterNP,j))
      END DO
   ! INTERIOR
      DO j = 1+GhostGridY+filterALPHA, Ny-GhostGridY-filterALPHA
         U(i,j) = DOT_PRODUCT(tmp(i,j-filterALPHA:j+filterALPHA),filtercoefficients)
      END DO
      ! Top boundary points
      DO j = 1,filterALPHA
         U(i,Ny-GhostGridY-filterALPHA+j) = DOT_PRODUCT(tmp(i,Ny-GhostGridY-(filterNP)+1:Ny-GhostGridY),&
              filtercoefficients2(filterNP:1:-1,filterALPHA+1-j))
      END DO
   END DO
ENDIF
!PRINT*,'  SG-FILTER applied.'
END SUBROUTINE
