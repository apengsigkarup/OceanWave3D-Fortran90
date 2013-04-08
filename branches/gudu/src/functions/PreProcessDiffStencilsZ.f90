SUBROUTINE PreProcessDiffStencilsZ(FineGrid,DiffStencils,GhostGridZ,gamma)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Diff_def), INTENT(OUT) :: DiffStencils
TYPE (Level_def), INTENT(IN) :: FineGrid
INTEGER :: GhostGridZ,STAT
INTEGER :: rank, Nz, order, Nzp, Diff, gamma, i
REAL(KIND=long), DIMENSION(:), ALLOCATABLE :: Stencil
INTEGER, DIMENSION(:), ALLOCATABLE :: idx

! TEMPORARY POINTERS/VARIABLES
Nz = FineGrid%Nz
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
END SUBROUTINE PreProcessDiffStencilsZ
