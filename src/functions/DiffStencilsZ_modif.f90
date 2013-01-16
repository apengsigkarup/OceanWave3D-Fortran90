SUBROUTINE DiffStencilsZ_modif(k,z,StencilZ,Nz,GhostGridZ,gamma,idxg_c)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
!
INTEGER :: Nz, GhostGridZ, k, gamma
REAL(KIND=long), DIMENSION(1:2*gamma+1,1:2) :: StencilZ
REAL(KIND=long), DIMENSION(1:Nz+GhostGridZ) :: z
!
INTEGER :: Nzp, order, Diff, i
REAL(KIND=long), DIMENSION(:), ALLOCATABLE :: Stencil
INTEGER, DIMENSION(:), ALLOCATABLE :: idx
INTEGER :: idxg_c
!

! TEMPORARY POINTERS/VARIABLES
! Derivatives in z-direction
!    One-sided stencils near free surface.
!    One-side stencils near bottom, and if ghost points
!    are supported then the directional ghost point
!    is included in stencil.
! Correction for correct treatment of corners included through idxg_c
IF (GhostGridZ==0) THEN; Nzp = Nz; ELSE; Nzp = Nz+1; ENDIF
ALLOCATE( Stencil(2*gamma+1) )
ALLOCATE( idx(2*gamma+1) )
DO i=-gamma,gamma,1
	idx(i+gamma+1) = i
END DO
DO order=1,2
  DO i = k,k
	IF (i-1<gamma) THEN
		Diff = gamma-i+1
		CALL TaylorFDStencils1DArbitrary(gamma-Diff-idxg_c,gamma+Diff+idxg_c,order,Stencil,z(i+idx+Diff+idxg_c))
	ELSEIF (i>=Nz-gamma+1) THEN
		Diff = i-(Nzp-gamma)
		CALL TaylorFDStencils1DArbitrary(gamma+Diff-idxg_c,gamma-Diff+idxg_c,order,Stencil,z(i+idx-Diff+idxg_c))
	ELSE
		CALL TaylorFDStencils1DArbitrary(gamma-idxg_c,gamma+idxg_c,order,Stencil,z(i+idx+idxg_c))
	ENDIF
	StencilZ(:,order) = Stencil
  END DO
END DO
DEALLOCATE(Stencil,idx)
END SUBROUTINE DiffStencilsZ_modif
