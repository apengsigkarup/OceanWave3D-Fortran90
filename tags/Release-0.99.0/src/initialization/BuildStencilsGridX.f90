SUBROUTINE BuildStencilsGridX(alpha,order,x,DiffStencils,Nx,Ny)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: rank, alpha, Nx, Ny, order, Diff, i
REAL(KIND=long), DIMENSION(:), ALLOCATABLE   :: Stencil, tmpStencil
!INTEGER, DIMENSION(:), ALLOCATABLE   		 :: idx
INTEGER :: idx(2*alpha+1)
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: c
REAL(KIND=long), DIMENSION(Nx,2*alpha+1)     :: DiffStencils
REAL(KIND=long), DIMENSION(Nx,Ny)            :: x
REAL(KIND=long)                              :: x0
rank = 2*alpha+1
ALLOCATE( Stencil(rank)    )
ALLOCATE( tmpStencil(rank) )
ALLOCATE( c(rank,order+1)  )
IF (Nx==1 .AND. order==0) THEN
    ! Direct interpolation
	Stencil           = zero
	Stencil(alpha+1)  = one
	DiffStencils(1,:) = Stencil
ELSE IF (Nx==1 .AND. order>0) THEN
    ! For a stencil of one point, the derivative is zero
	Stencil           = zero
	Stencil(alpha+1)  = zero
	DiffStencils(1,:) = Stencil
ELSE
!	ALLOCATE( idx(rank) )
	DO i=-alpha,alpha,1
		idx(i+alpha+1) = i
	END DO
	DO i = 1 , alpha
		Diff = alpha-i+1
		! Expansion point where the approximation is to be defined
		x0 = x(i,1) ! FIXME: y-index fixed

		! Determine local finite difference stencil coefficients
		CALL weights(x0,x(i+idx+Diff,1),rank-1,rank-1,order,c) ! FIXME: y-index fixed
		Stencil    = c(1:rank,order+1)
		DiffStencils(i,:) = Stencil
	END DO
	DO i = alpha+1 , Nx-alpha
		! Expansion point where the approximation is to be defined
		x0 = x(i,1) ! FIXME: y-index fixed

		! Determine local finite difference stencil coefficients
		! FIXME: put his outside loop if we do not need flexible grids (e.g. in the case of transformation to an equi-distant grid)
		CALL weights(x0,x(i+idx,1),rank-1,rank-1,order,c) ! FIXME: y-index fixed
		Stencil    = c(1:rank,order+1)
		DiffStencils(i,:) = Stencil
	END DO
	DO i = Nx-alpha+1 , Nx
		Diff = (Nx-alpha)-i
		! Expansion point where the approximation is to be defined
		x0 = x(i,1) ! FIXME: y-index fixed

		! Determine local finite difference stencil coefficients
		CALL weights(x0,x(i+idx+Diff,1),rank-1,rank-1,order,c) ! FIXME: y-index fixed
		Stencil    = c(1:rank,order+1)
		DiffStencils(i,:) = Stencil
	END DO
ENDIF
DEALLOCATE(Stencil,tmpStencil,c)

END SUBROUTINE BuildStencilsGridX
