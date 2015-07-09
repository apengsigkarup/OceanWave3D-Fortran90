SUBROUTINE UpdateGhostLayerY(var,S,Nx,Ny,DiffStencils,beta,GhostGridX,GhostGridY)
! By Allan P. Engsig-Karup.
USE Precision
USE DataTypes
IMPLICIT NONE
TYPE (Diff_def) :: DiffStencils
INTEGER :: Nx, Ny, rank, beta, i,j, GhostGridX, GhostGridY
REAL(KIND=long), DIMENSION(Nx,Ny) :: var, S
REAL(KIND=long), DIMENSION(2*beta+1) :: Stencil
rank = 2*beta+1
! two boundaries
j=1
DO i = 1+GhostGridX, Nx-GhostGridX
    Stencil = DiffStencils%StencilY(j+GhostGridY,1:rank,1)
    var(i,j) = (S(i,j+GhostGridY)-DOT_PRODUCT(Stencil(2:rank),var(i,2:rank) )) / Stencil(1)
END DO
j = Ny
DO i = 1+GhostGridX, Nx-GhostGridX
    Stencil = DiffStencils%StencilY(j-GhostGridY,1:rank,1)
    var(i,j) = (S(i,j-GhostGridY)-DOT_PRODUCT(Stencil(1:rank-1),var(i,j-rank+1:j-1) )) / Stencil(rank)
END DO
END SUBROUTINE UpdateGhostLayerY
