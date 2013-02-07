SUBROUTINE GaussSeidelARRAY(x,b,GridLevel,GhostGridX,GhostGridY,GhostGridZ, alpha, beta, gamma)
! By Allan P. Engsig-Karup.
!
! Tailored routine to apply the Gauss-Seidel method explicitly
! for solving A x = b
!
! FIXME: IMPLEMENT FOR 2D YZ PROBLEMS, AND 3D
!
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: GhostGridX, GhostGridY, GhostGridZ, i, j, k, Gidx, alpha, beta, gamma, diff, qqq
INTEGER :: Nz
TYPE (Level_def), INTENT(IN) :: GridLevel ! Grid information, stencils, etc.
REAL(KIND=long), DIMENSION(GridLevel%Nx*GridLevel%Ny*(GridLevel%Nz+GhostGridZ)), INTENT(INOUT) :: x
REAL(KIND=long), DIMENSION(GridLevel%Nx*GridLevel%Ny*(GridLevel%Nz+GhostGridZ)), INTENT(IN) :: b
REAL(KIND=long) :: diagcoef, coef, ddpdxx, ddpdyy, ddpdss, dpds, ddpdsdx, ddpdsdy, dpdx, dpdy, rhs
REAL(KIND=long) :: tmpstencil8(8), tmpstencil2(2), STENCILX(2), STENCILY(2)
INTEGER, DIMENSION(3) :: idx, tmpidx

DO i = -1, 1
	idx(i+2) = i
END DO

Nz = GridLevel%Nz+GhostGridZ

IF (GridLevel%Ny==1) THEN
  ! UPDATE FREE SURFACE VALUES
  k = Nz
  j = 1
  DO i = 1, GridLevel%Nx
    x(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)) = b(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1))
  END DO

  j = 1
  ! UPDATE GHOST POINT VALUES
  ! AND SCALE WITH DIAGONAL VALUE TO RETURN GHOST POINT VALUES EXPLICITLY
  k = 1
  tmpstencil2 = GridLevel%DiffStencils%StencilZ(k+GhostGridZ,2:3,1) ! FOR dpds
  DO i = 1, GridLevel%Nx
    ! UPDATE GHOST POINT VALUES
    ! AND SCALE WITH DIAGONAL VALUE TO RETURN GHOST POINT VALUES EXPLICITLY
    coef = GridLevel%dsigmanew(k,i,j,5) + GridLevel%hx(i,j)*GridLevel%dsigmanew(k,i,j,2)
    diagcoef = coef*GridLevel%DiffStencils%StencilZ(k+GhostGridZ,1,1) ! diagonal coefficient for the ghost point in kinematic boundary condition
    dpds = DOT_PRODUCT( tmpstencil2 , x(GridLevel%FullRankStencilsNEW%IndexesZ(k,i,j,2:3)) ) ! ASSUMED 3-POINT STENCIL HERE (SINCE IT IS LOW-ORDER)
    dpdx = DOT_PRODUCT( GridLevel%FullRankStencilsNEW%StencilX(i,2:3,1) , &
	       x(GridLevel%FullRankStencilsNEW%IndexesX(k+GhostGridZ,i,j,2:3)) )
	x(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)) = -( coef*dpds+GridLevel%hx(i,j)*dpdx - &
	       b(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1))) / diagcoef
  END DO

  ! Then update interior points above ghost layer and below free surface
  DO i = 1, GridLevel%Nx
    STENCILX = GridLevel%FullRankStencilsNEW%StencilX(i,2:3,2)
	DO k = 2, Nz-1 ! ONLY INTERIOR POINTS IN DOMAIN ARE UPDATED
	  ! UPDATE TRANSFORMED LAPLACE EQUATIONS - DIAGONAL TERM NOT INCLUDED
	  ! ASSUMED FIRST ORDER STENCIL IN EACH CARTESIAN DIRECTION -> r^2 = 3*3 = 9
      ddpdxx  = DOT_PRODUCT(STENCILX , x(GridLevel%FullRankStencilsNEW%IndexesX(k,i,j,2:3)) )
	  ddpdss  = DOT_PRODUCT(GridLevel%FullRankStencilsNEW%StencilZ(k,2:3,2) , &
	            x(GridLevel%FullRankStencilsNEW%IndexesZ(k,i,j,2:3)) )
	  dpds    = DOT_PRODUCT(GridLevel%FullRankStencilsNEW%StencilZ(k,2:3,1) , &
	            x(GridLevel%FullRankStencilsNEW%IndexesZ(k,i,j,2:3)) )
	  ddpdsdx = DOT_PRODUCT(GridLevel%FullRankStencilsNEW%StencilXZorYZnew(k,i,j,2:9,1) , &
	            x(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,2:9)) )
	  rhs     = ddpdxx + GridLevel%dsigmanew(k,i,j,3)*dpds + two*(GridLevel%dsigmanew(k,i,j,2)*ddpdsdx ) + &
	            (GridLevel%dsigmanew(k,i,j,2)**2+GridLevel%dsigmanew(k,i,j,5)**2)*ddpdss
	  diagcoef = GridLevel%FullRankStencilsNEW%StencilX(i,1,2) + &
	             GridLevel%dsigmanew(k,i,j,3)*GridLevel%FullRankStencilsNEW%StencilZ(k,1,1) + &
				 two*(GridLevel%dsigmanew(k,i,j,2)*GridLevel%FullRankStencilsNEW%StencilXZorYZnew(k,i,j,1,1) ) + &
				 (GridLevel%dsigmanew(k,i,j,2)**2+GridLevel%dsigmanew(k,i,j,5)**2)*GridLevel%FullRankStencilsNEW%StencilZ(k,1,2)
	  x(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)) = (-rhs + &
	            b(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)) )/diagcoef
	END DO
  END DO
ELSE IF (GridLevel%Nx==1) THEN
  ! UPDATE FREE SURFACE VALUES
  k = Nz
  i = 1
  DO j = 1, GridLevel%Ny
    x(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)) = b(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1))
  END DO

  i = 1
  ! UPDATE GHOST POINT VALUES
  ! AND SCALE WITH DIAGONAL VALUE TO RETURN GHOST POINT VALUES EXPLICITLY
  k = 1
  tmpstencil2 = GridLevel%DiffStencils%StencilZ(k+GhostGridZ,2:3,1) ! FOR dpds
  DO j = 1, GridLevel%Ny
    ! UPDATE GHOST POINT VALUES
    ! AND SCALE WITH DIAGONAL VALUE TO RETURN GHOST POINT VALUES EXPLICITLY
    coef = GridLevel%dsigmanew(k,i,j,5) + GridLevel%hy(i,j)*GridLevel%dsigmanew(k,i,j,4)
    diagcoef = coef*GridLevel%DiffStencils%StencilZ(k+GhostGridZ,1,1) ! diagonal coefficient for the ghost point in kinematic boundary condition
    dpds = DOT_PRODUCT( tmpstencil2 , x(GridLevel%FullRankStencilsNEW%IndexesZ(k,i,j,2:3)) ) ! ASSUMED 3-POINT STENCIL HERE (SINCE IT IS LOW-ORDER)
    dpdy = DOT_PRODUCT( GridLevel%FullRankStencilsNEW%StencilY(j,2:3,1) , &
	       x(GridLevel%FullRankStencilsNEW%IndexesY(k+GhostGridZ,i,j,2:3)) )
	x(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)) = -( coef*dpds+GridLevel%hy(i,j)*dpdy - &
	       b(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1))) / diagcoef
  END DO

  ! Then update interior points above ghost layer and below free surface
  DO j = 1, GridLevel%Ny
    STENCILY = GridLevel%FullRankStencilsNEW%StencilY(j,2:3,2)
	DO k = 2, Nz-1 ! ONLY INTERIOR POINTS IN DOMAIN ARE UPDATED
	  ! UPDATE TRANSFORMED LAPLACE EQUATIONS - DIAGONAL TERM NOT INCLUDED
	  ! ASSUMED FIRST ORDER STENCIL IN EACH CARTESIAN DIRECTION -> r^2 = 3*3 = 9
      ddpdyy  = DOT_PRODUCT(STENCILY , x(GridLevel%FullRankStencilsNEW%IndexesY(k,i,j,2:3)) )
	  ddpdss  = DOT_PRODUCT(GridLevel%FullRankStencilsNEW%StencilZ(k,2:3,2) , &
	            x(GridLevel%FullRankStencilsNEW%IndexesZ(k,i,j,2:3)) )
	  dpds    = DOT_PRODUCT(GridLevel%FullRankStencilsNEW%StencilZ(k,2:3,1) , &
	            x(GridLevel%FullRankStencilsNEW%IndexesZ(k,i,j,2:3)) )
	  ddpdsdy = DOT_PRODUCT(GridLevel%FullRankStencilsNEW%StencilXZorYZnew(k,i,j,2:9,2) , &
	            x(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,2:9)) )
	  rhs     = ddpdyy + GridLevel%dsigmanew(k,i,j,3)*dpds + two*(GridLevel%dsigmanew(k,i,j,4)*ddpdsdy ) + &
	            (GridLevel%dsigmanew(k,i,j,4)**2+GridLevel%dsigmanew(k,i,j,5)**2)*ddpdss
	  diagcoef = GridLevel%FullRankStencilsNEW%StencilY(j,1,2) + &
	            GridLevel%dsigmanew(k,i,j,3)*GridLevel%FullRankStencilsNEW%StencilZ(k,1,1) + &
				two*(GridLevel%dsigmanew(k,i,j,4)*GridLevel%FullRankStencilsNEW%StencilXZorYZnew(k,i,j,1,2) ) + &
				(GridLevel%dsigmanew(k,i,j,4)**2+GridLevel%dsigmanew(k,i,j,5)**2)*GridLevel%FullRankStencilsNEW%StencilZ(k,1,2)
	  x(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)) = (-rhs + &
	            b(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)) )/diagcoef

!print*,'Gidx=',GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)
!print*,'ddpdyy=',ddpdyy
!print*,'stencilYY=',
!print*,'ddpdss=',ddpdss
!print*,'dpds=',dpds
!print*,'ddpdsdy=',ddpdsdy
!print*,'rhs=',rhs
!print*,'diagcoef=',diagcoef
!print*,'x(Gidx)=',x(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1))
!read*

	END DO
  END DO
ELSE ! 3D
  ! UPDATE FREE SURFACE VALUES
  k = Nz
  DO j = 1, GridLevel%Ny
    DO i = 1, GridLevel%Nx
      x(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)) = b(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1))
    END DO
  END DO

  ! UPDATE GHOST POINT VALUES
  ! AND SCALE WITH DIAGONAL VALUE TO RETURN GHOST POINT VALUES EXPLICITLY
  k = 1
  tmpstencil2 = GridLevel%DiffStencils%StencilZ(k+GhostGridZ,2:3,1) ! FOR dpds
  DO j = 1, GridLevel%Ny
   STENCILY = GridLevel%FullRankStencilsNEW%StencilY(j,2:3,1)
   DO i = 1, GridLevel%Nx
    ! UPDATE GHOST POINT VALUES
    ! AND SCALE WITH DIAGONAL VALUE TO RETURN GHOST POINT VALUES EXPLICITLY
    coef = GridLevel%dsigmanew(k,i,j,5) + GridLevel%hx(i,j)*GridLevel%dsigmanew(k,i,j,2) + &
	       GridLevel%hy(i,j)*GridLevel%dsigmanew(k,i,j,4)
    diagcoef = coef*GridLevel%DiffStencils%StencilZ(k+GhostGridZ,1,1) ! diagonal coefficient for the ghost point in kinematic boundary condition
    dpds = DOT_PRODUCT( tmpstencil2 , x(GridLevel%FullRankStencilsNEW%IndexesZ(k,i,j,2:3)) ) ! ASSUMED 3-POINT STENCIL HERE (SINCE IT IS LOW-ORDER)
    dpdx = DOT_PRODUCT( GridLevel%FullRankStencilsNEW%StencilX(i,2:3,1) , &
	       x(GridLevel%FullRankStencilsNEW%IndexesX(k+GhostGridZ,i,j,2:3)) )
    dpdy = DOT_PRODUCT( STENCILY , x(GridLevel%FullRankStencilsNEW%IndexesY(k+GhostGridZ,i,j,2:3)) )
	x(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)) = -( coef*dpds+GridLevel%hx(i,j)*dpdx+GridLevel%hy(i,j)*dpdy - &
	       b(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1))) / diagcoef
   END DO
  END DO

  ! Then update interior points above ghost layer and below free surface
  DO j = 1, GridLevel%Ny
   STENCILY = GridLevel%FullRankStencilsNEW%StencilY(j,2:3,2)
   DO i = 1, GridLevel%Nx
    STENCILX = GridLevel%FullRankStencilsNEW%StencilX(i,2:3,2)
	DO k = 2, Nz-1 ! ONLY INTERIOR POINTS IN DOMAIN ARE UPDATED
	  ! UPDATE TRANSFORMED LAPLACE EQUATIONS - DIAGONAL TERM NOT INCLUDED
	  ! ASSUMED FIRST ORDER STENCIL IN EACH CARTESIAN DIRECTION -> r^2 = 3*3 = 9
      ddpdxx  = DOT_PRODUCT(STENCILX , x(GridLevel%FullRankStencilsNEW%IndexesX(k,i,j,2:3)) )
      ddpdyy  = DOT_PRODUCT(STENCILY , x(GridLevel%FullRankStencilsNEW%IndexesY(k,i,j,2:3)) )
	  ddpdss  = DOT_PRODUCT(GridLevel%FullRankStencilsNEW%StencilZ(k,2:3,2) , &
	            x(GridLevel%FullRankStencilsNEW%IndexesZ(k,i,j,2:3)) )
	  dpds    = DOT_PRODUCT(GridLevel%FullRankStencilsNEW%StencilZ(k,2:3,1) , &
	            x(GridLevel%FullRankStencilsNEW%IndexesZ(k,i,j,2:3)) )
	  ddpdsdx = DOT_PRODUCT(GridLevel%FullRankStencilsNEW%StencilXZorYZnew(k,i,j,2:9,1) , &
	            x(GridLevel%FullRankStencilsNEW%IndexesnewXZ(k,i,j,2:9)) )
	  ddpdsdy = DOT_PRODUCT(GridLevel%FullRankStencilsNEW%StencilXZorYZnew(k,i,j,2:9,2) , &
	            x(GridLevel%FullRankStencilsNEW%IndexesnewYZ(k,i,j,2:9)) )
	  rhs     = ddpdxx + ddpdyy + GridLevel%dsigmanew(k,i,j,3)*dpds + two*(GridLevel%dsigmanew(k,i,j,2)*ddpdsdx + &
	            GridLevel%dsigmanew(k,i,j,4)*ddpdsdy ) + (GridLevel%dsigmanew(k,i,j,2)**2+&
				GridLevel%dsigmanew(k,i,j,4)**2+GridLevel%dsigmanew(k,i,j,5)**2)*ddpdss
	  diagcoef = GridLevel%FullRankStencilsNEW%StencilX(i,1,2) + GridLevel%FullRankStencilsNEW%StencilY(j,1,2) + &
	            GridLevel%dsigmanew(k,i,j,3)*GridLevel%FullRankStencilsNEW%StencilZ(k,1,1) + &
				two*(GridLevel%dsigmanew(k,i,j,2)*GridLevel%FullRankStencilsNEW%StencilXZorYZnew(k,i,j,1,1)+&
				GridLevel%dsigmanew(k,i,j,4)*GridLevel%FullRankStencilsNEW%StencilXZorYZnew(k,i,j,1,2) ) + &
				(GridLevel%dsigmanew(k,i,j,2)**2+GridLevel%dsigmanew(k,i,j,4)**2+&
				GridLevel%dsigmanew(k,i,j,5)**2)*GridLevel%FullRankStencilsNEW%StencilZ(k,1,2) ! FIXME: always constant on each grid!
	  x(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)) = (-rhs + &
	            b(GridLevel%FullRankStencilsNEW%Indexesnew(k,i,j,1)) )/diagcoef
	END DO
   END DO
  END DO

ENDIF

END SUBROUTINE GaussSeidelARRAY
