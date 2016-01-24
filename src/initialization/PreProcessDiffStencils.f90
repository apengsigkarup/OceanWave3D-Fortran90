SUBROUTINE PreProcessDiffStencils(FineGrid,DiffStencils,GhostGridX,GhostGridY,GhostGridZ,alpha,beta,gamma,fileop)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Diff_def), INTENT(OUT) :: DiffStencils
TYPE (Level_def), INTENT(IN) :: FineGrid
INTEGER :: GhostGridX,GhostGridY,GhostGridZ,STAT, fileop
INTEGER :: rank, Nx, Ny, Nz, order, Nzp, Diff, alpha, beta, gamma, i
REAL(KIND=long), DIMENSION(:), ALLOCATABLE :: Stencil
INTEGER, DIMENSION(:), ALLOCATABLE :: idx

! TEMPORARY POINTERS/VARIABLES
Nx = FineGrid%Nx
Ny = FineGrid%Ny
Nz = FineGrid%Nz
!print*,'Nx=',Nx
!print*,'Ny=',Ny
!print*,'Nz=',Nz
! Derivatives in x- and y-directions
!    Neumann conditions imposed by symmetric reflections
!    near walls. Central stencils everywhere.
rank = 2*alpha+1
!print*,'allocation'
ALLOCATE( DiffStencils%StencilX(Nx+2*GhostGridX,rank,2),STAT=STAT)
IF (STAT/=0) THEN
	PRINT*,'Error: Could not allocate memory for differentiation stencil (X).'
	STOP
ENDIF
DiffStencils%StencilX=zero
rank = 2*beta+1
!print*,'allocation'
ALLOCATE( DiffStencils%StencilY(Ny+2*GhostGridY,rank,2),STAT=STAT)
!print*,'hej'
IF (STAT/=0) THEN
	PRINT*,'Error: Could not allocate memory for differentiation stencil (Y).'
	STOP
ENDIF
DiffStencils%StencilY=zero

print*,'Building stencils...'
write(fileop,*)'Building stencils...'
DO order = 1,2
    IF (Nx>1) THEN
!    print*,'FineGrid%x=',FineGrid%x
		CALL BuildStencilsGridX(alpha,order,FineGrid%x,DiffStencils%StencilX(:,:,order),Nx+2*GhostGridX,Ny+2*GhostGridY)
!    print*,'DiffStencils%StencilX(:,:,order) = ',DiffStencils%StencilX(:,:,order)
!    read*
	ENDIF
	IF (Ny>1) THEN
		CALL BuildStencilsGridY(beta,order,FineGrid%y,DiffStencils%StencilY(:,:,order),Nx+2*GhostGridX,Ny+2*GhostGridY)
	ENDIF
END DO
!print*,'DONE!'

! Derivatives in z-direction
!    One-sided stencils near free surface.
!    One-side stencils near bottom, and if ghost points
!    are supported then the directional ghost point
!    is included in stencil.
IF (GhostGridZ==0) THEN; Nzp = Nz; ELSE; Nzp = Nz+1; ENDIF
ALLOCATE( DiffStencils%StencilZ(Nzp,2*gamma+1,2),STAT=STAT)
IF (STAT/=0) THEN
	PRINT*,'Error: Could not allocate memory for differentiation stencil (Z).'
	STOP
ENDIF
ALLOCATE( Stencil(2*gamma+1) )
ALLOCATE( idx(2*gamma+1) )
DO i=-gamma,gamma,1
	idx(i+gamma+1) = i
END DO
DO order=1,2
  DO i = 1, Nzp
	IF (i-1<gamma) THEN
		Diff = gamma-i+1
		CALL TaylorFDStencils1DArbitrary(gamma-Diff,gamma+Diff,order,Stencil,FineGrid%z(i+idx+Diff))
	ELSEIF (i>=Nz-gamma+1) THEN
		Diff = i-(Nzp-gamma)
		CALL TaylorFDStencils1DArbitrary(gamma+Diff,gamma-Diff,order,Stencil,FineGrid%z(i+idx-Diff))
	ELSE
		CALL TaylorFDStencils1DArbitrary(gamma,gamma,order,Stencil,FineGrid%z(i+idx))
	ENDIF
	DiffStencils%StencilZ(i,:,order) = Stencil
  END DO
END DO
DEALLOCATE(Stencil,idx)
END SUBROUTINE PreProcessDiffStencils
