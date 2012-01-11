SUBROUTINE DiscardGhostLayer(Nx,Ny,Nz,PHI)
!
! Zero bottom layer (e.g. ghost layer).
!
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: Nx, Ny, Nz
REAL(KIND=long) :: PHI(Nz,Nx,Ny)
! Discard all elements in ghost layer
PHI(1,:,:) = zero
END SUBROUTINE DiscardGhostLayer
