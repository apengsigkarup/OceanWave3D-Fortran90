SUBROUTINE UpdateGhostLayerX(var,S,Nx,Ny,DiffStencils,alpha,GhostGridX,GhostGridY)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Diff_def) :: DiffStencils
INTEGER :: Nx, Ny, rank, alpha, i,j, GhostGridX, GhostGridY
REAL(KIND=long), DIMENSION(Nx,Ny) :: var, S
REAL(KIND=long), DIMENSION(2*alpha+1) :: Stencil
rank = 2*alpha+1
! two boundaries
i=1
DO j = 1+GhostGridY, Ny-GhostGridY
    Stencil = DiffStencils%StencilX(i+GhostGridX,1:rank,1)
    var(i,j) = (S(i+GhostGridX,j)-DOT_PRODUCT(Stencil(2:rank),var(2:rank,j) )) / Stencil(1)
END DO
i = Nx
DO j = 1+GhostGridY, Ny-GhostGridY
    Stencil = DiffStencils%StencilX(i-GhostGridX,1:rank,1)
    var(i,j) = (S(i-GhostGridX,j)-DOT_PRODUCT(Stencil(1:rank-1),var(i-rank+1:i-1,j) )) / Stencil(rank)
END DO
END SUBROUTINE UpdateGhostLayerX
