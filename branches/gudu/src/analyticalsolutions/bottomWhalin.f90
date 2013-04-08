!SUBROUTINE BottomWhalin(FineGrid,y0)
!GD: modif
SUBROUTINE BottomWhalin(FineGrid,y0,GhostGridX,GhostGridY)
!
! Define bottom and bottom gradients for whalin's test (3D)
! By Allan P. Engsig-Karup.
!
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: FineGrid
REAL(KIND=long) :: y0, SIGMA, y, x, SIGMAy, SIGMAyy, SIGMAx, SIGMAxx
INTEGER :: Nx, Ny, i, j, GhostGridX, GhostGridY
!Nx = FineGrid%Nx
!Ny = FineGrid%Ny
!GD: modif
Nx = FineGrid%Nx+2*GhostGridX
Ny = FineGrid%Ny+2*GhostGridY

IF (FineGrid%Nx>FineGrid%Ny) THEN
  DO j = 1, Ny
	y = FineGrid%y(1,j) ! FIXME: x-index
	IF (y<y0) y = y0
	SIGMA   = SQRT(6.096_long*y-y*y)
	SIGMAy  = -0.04_long*(6.096_long-two*y)*half/SIGMA
	SIGMAyy = 0.04_long*6.096_long**2/(four*SIGMA**3)
	DO i = 1 , Nx
		x = FineGrid%x(i,1) ! FIXME: y-index
		IF (x<=10.67-SIGMA) THEN
			FineGrid%h(i,j)   = 0.4572_long
			FineGrid%hx(i,j)  = zero
			FineGrid%hy(i,j)  = zero
			FineGrid%hxx(i,j) = zero
			FineGrid%hyy(i,j) = zero
		ELSE IF (x>=18.29-SIGMA) THEN
			FineGrid%h(i,j)   = 0.1524_long
			FineGrid%hx(i,j)  = zero
			FineGrid%hy(i,j)  = zero
			FineGrid%hxx(i,j) = zero
			FineGrid%hyy(i,j) = zero
		ELSE ! 10.67 - SIGMA < x < 18.29 - SIGMA
			FineGrid%h(i,j)   = 0.4572_long + 0.04_long*(10.67_long-SIGMA-x)
			FineGrid%hx(i,j)  = -0.04_long
			FineGrid%hy(i,j)  = SIGMAy
			FineGrid%hxx(i,j) = zero
			FineGrid%hyy(i,j) = SIGMAyy
		ENDIF
	END DO
  END DO
ELSE
  DO i = 1, Nx
	x = FineGrid%x(i,1) ! FIXME: y-index
	IF (x<y0) x = y0
	SIGMA   = SQRT(6.096_long*x-x*x)
	SIGMAx  = -0.04_long*(6.096_long-two*x)*half/SIGMA
	SIGMAxx = 0.04_long*6.096_long**2/(four*SIGMA**3)
	DO j = 1 , Ny
		y = FineGrid%y(1,j) ! FIXME: x-index
		IF (y<=10.67_long-SIGMA) THEN
			FineGrid%h(i,j)   = 0.4572_long
			FineGrid%hx(i,j)  = zero
			FineGrid%hy(i,j)  = zero
			FineGrid%hxx(i,j) = zero
			FineGrid%hyy(i,j) = zero
		ELSE IF (y>=18.29_long-SIGMA) THEN
			FineGrid%h(i,j)   = 0.1524_long
			FineGrid%hx(i,j)  = zero
			FineGrid%hy(i,j)  = zero
			FineGrid%hxx(i,j) = zero
			FineGrid%hyy(i,j) = zero
		ELSE ! 10.67 - SIGMA < y < 18.29 - SIGMA
			FineGrid%h(i,j)   = 0.4572_long + 0.04_long*(10.67_long-SIGMA-y)
			FineGrid%hy(i,j)  = -0.04_long
			FineGrid%hx(i,j)  = SIGMAx
			FineGrid%hyy(i,j) = zero
			FineGrid%hxx(i,j) = SIGMAxx
		ENDIF
	END DO
  END DO
ENDIF
END SUBROUTINE BottomWhalin
