SUBROUTINE DetermineTransformationConstants(Nx,Ny,Nz,FineGrid,dsigma,Wavefield)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: FineGrid
INTEGER :: Nx,Ny,Nz,Gidx
REAL(KIND=long) :: d, z
REAL(KIND=long), DIMENSION(Nx*Ny*Nz,5), INTENT(OUT) :: dsigma
REAL(KIND=long), DIMENSION(Nx,Ny) :: E,Ex,Exx,Ey,Eyy,h,hx,hxx,hy,hyy
REAL(KIND=long), DIMENSION(Nz) :: zgrid
TYPE(Wavefield_FS) :: Wavefield
!Local variables
INTEGER :: i,j,k
h = FineGrid%h
E = Wavefield%E
IF (Nx>1) THEN
	hx  = FineGrid%hx
	hxx = FineGrid%hxx
    Ex  = Wavefield%Ex
    Exx = Wavefield%Exx
ENDIF
IF (Ny>1) THEN
	hy  = FineGrid%hy
	hyy = FineGrid%hyy
    Ey  = Wavefield%Ey
    Eyy = Wavefield%Eyy
ENDIF
! GD: SWENSE addition for nonlinear
! FIXME: during transition period we do not satisfy Ex+Ex_I=0 nor Ey+Ey_I=0...
! The incident terms used here are without time ramp...
! Include this possibility in the NormalZ which will be =/ 0 on horiz. boundaries ?
IF (ASSOCIATED(Wavefield%E_I)) THEN ! This means that swense is ON
  E = E + Wavefield%E_I
  IF (Nx>1) THEN
      Ex  = Ex  + Wavefield%Ex_I
      Exx = Exx + Wavefield%Exx_I
  ENDIF
  IF (Ny>1) THEN
      Ey  = Ey  + Wavefield%Ey_I
      Eyy = Eyy + Wavefield%Eyy_I
  ENDIF
ENDIF
zgrid = FineGrid%z
! Assume that z=0 at below bottom ghost layer.
zgrid(1) = zero
Gidx = 0
IF (Ny==1) THEN
  j = 1
  DO i = 1, Nx
	DO k = 1, Nz
		Gidx = Gidx + 1;
		! Determine z-coordinate (NOTE THAT ACTUALLY sigma = z, see below)
		z = zgrid(k)
		d = E(i,j)+h(i,j)
		! sigma
		dsigma(Gidx,1) = z ! (z+h(i,j))/d
	    ! d(sigma)/dx
		dsigma(Gidx,2) = ((one-dsigma(Gidx,1))*hx(i,j)-dsigma(Gidx,1)*Ex(i,j))/d
!		dsigma(Gidx,2) = (one-dsigma(Gidx,1))/d*hx(i,j)-dsigma(Gidx,1)/d*Ex(i,j)
	    ! d^2(sigma)/dx^2
	    dsigma(Gidx,3) = -two*(one-dsigma(Gidx,1))*(hx(i,j)**2)/d**2+&
		  two*(two*dsigma(Gidx,1)-one)*(Ex(i,j)*hx(i,j))/d**2+&
		  (one-dsigma(Gidx,1))*(hxx(i,j))/d+two*dsigma(Gidx,1)*(Ex(i,j)**2)/d**2-dsigma(Gidx,1)*(Exx(i,j))/d
		! d(sigma)/dz
		dsigma(Gidx,5) = one/d
	END DO
  END DO
ELSE IF (Nx==1) THEN
  i = 1
  DO j = 1, Ny
	DO k = 1, Nz
		Gidx = Gidx + 1;
		! Determine z-coordinate (NOTE THAT ACTUALLY sigma = z, see below)
		z = zgrid(k)
		d = E(i,j)+h(i,j)
		! sigma
		dsigma(Gidx,1) = z ! (z+h(i,j))/d
	    ! d^2(sigma)/dy^2
		dsigma(Gidx,3) = -two*(one-dsigma(Gidx,1))*(hy(i,j)**2)/d**2+&
		  two*(two*dsigma(Gidx,1)-one)*(Ey(i,j)*hy(i,j))/d**2+&
		  (one-dsigma(Gidx,1))*(hyy(i,j))/d+two*dsigma(Gidx,1)*(Ey(i,j)**2)/d**2-dsigma(Gidx,1)*(Eyy(i,j))/d
	    ! d(sigma)/dy
		dsigma(Gidx,4) = ((one-dsigma(Gidx,1))*hy(i,j)-dsigma(Gidx,1)*Ey(i,j))/d
		! d(sigma)/dz
		dsigma(Gidx,5) = one/d
	END DO
  END DO
ELSE ! 3D
  DO j = 1, Ny
	DO i = 1, Nx
		DO k = 1, Nz
			Gidx = Gidx + 1;
			! Determine z-coordinate (NOTE THAT ACTUALLY sigma = z, see below)
			z = zgrid(k)
			d = E(i,j)+h(i,j)
			! sigma
			dsigma(Gidx,1) = z ! (z+h(i,j))/d
		    ! d(sigma)/dx
			dsigma(Gidx,2) = ((one-dsigma(Gidx,1))*hx(i,j)-dsigma(Gidx,1)*Ex(i,j))/d
		    ! d^2(sigma)/dx^2+d^2(sigma)/dy^2
			dsigma(Gidx,3) = -two*(one-dsigma(Gidx,1))*(hx(i,j)**2+hy(i,j)**2)/d**2+&
			  two*(two*dsigma(Gidx,1)-one)*(Ex(i,j)*hx(i,j)+Ey(i,j)*hy(i,j))/d**2+&
			  (one-dsigma(Gidx,1))*(hxx(i,j)+hyy(i,j))/d+two*dsigma(Gidx,1)*(Ex(i,j)**2+&
			  Ey(i,j)**2)/d**2-dsigma(Gidx,1)*(Exx(i,j)+Eyy(i,j))/d
		    ! d(sigma)/dy
			dsigma(Gidx,4) = ((one-dsigma(Gidx,1))*hy(i,j)-dsigma(Gidx,1)*Ey(i,j))/d
			! d(sigma)/dz
			dsigma(Gidx,5) = one/d
		END DO
	END DO
  END DO
ENDIF
END SUBROUTINE DetermineTransformationConstants
