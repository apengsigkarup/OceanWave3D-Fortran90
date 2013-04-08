SUBROUTINE BuildStencilsGridY(beta,order,y,DiffStencils,Nx,Ny)
! By Allan P. Engsig-Karup.
USE Precision
  USE Constants
  IMPLICIT NONE
  INTEGER :: rank, beta, Nx, Ny, order, Diff, i
  REAL(KIND=long):: Stencil(2*beta+1), tmpStencil(2*beta+1)
  INTEGER :: idx(2*beta+1)
  REAL(KIND=long) :: c(2*beta+1,order+1)
  REAL(KIND=long), DIMENSION(Ny,2*beta+1)      :: DiffStencils
  REAL(KIND=long), DIMENSION(Nx,Ny)            :: y
  REAL(KIND=long)                              :: y0
!
  rank = 2*beta+1
!
  IF (Ny==1 .AND. order==0) THEN
     ! Direct interpolation
     Stencil           = zero
     Stencil(beta+1)   = one
     DiffStencils(1,:) = Stencil
  ELSE IF (Ny==1 .AND. order>0) THEN
     ! For a stencil of one point, the derivative is zero
     Stencil           = zero
     Stencil(beta+1)   = zero
     DiffStencils(1,:) = Stencil
  ELSE
     DO i=-beta,beta,1
        idx(i+beta+1) = i
     END DO
     DO i = 1, Ny
        ! FIXME: Avoid conditional statements below in inner loop. Unroll loops.
        IF (i-1<beta) THEN
           Diff = beta-i+1
        ELSEIF (i>=Ny-beta+1) THEN
           Diff = (Ny-beta)-i
        ELSE
           Diff = 0
        ENDIF
        ! Expansion point where the approximation is to be defined
        y0 = y(1,i) ! FIXME: x-index fixed

        ! Determine local finite difference stencil coefficients
        CALL weights(y0,y(1,i+idx+Diff),rank-1,rank-1,order,c) ! FIXME: y-index fixed
        Stencil    = c(1:rank,order+1)
        DiffStencils(i,:) = Stencil
     END DO
  ENDIF

END SUBROUTINE BuildStencilsGridY
