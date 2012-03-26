SUBROUTINE rhsFreeSurface3D_SWENSE(t,Wavefield,g,rhsE,rhsP,Nx,Ny)
!
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx,Ny
REAL(KIND=long) :: g,t
REAL(KIND=long), DIMENSION(Nx,Ny) :: rhsE, rhsP ! E, Ex, Ey, P, Px, Py, W,
!REAL(KIND=long), DIMENSION(Nx,Ny) :: E_I,Ex_I,Ey_I,Et_I,Px_I,Py_I,W_I,Pt_I
TYPE(Wavefield_FS) :: Wavefield
! Local variables
!$$$$$$ REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: PxS_I , PyS_I
REAL(KIND=long), DIMENSION(Nx,Ny) :: PtS_I , PxS_I , PyS_I

IF (Nx==1) THEN
!$$$$$$     ALLOCATE(PyS_I(Nx,Ny))
  	! Surfacic quantities for incident wave field
	PyS_I = Wavefield%Py_I_s + (Wavefield%Ey+Wavefield%Ey_I)*Wavefield%Pz_I_s
    !
	rhsE = - (Wavefield%Ey+Wavefield%Ey_I)*(Wavefield%Py+PyS_I) &
    	   + (Wavefield%W+Wavefield%Pz_I_s)*(one + (Wavefield%Ey+Wavefield%Ey_I)**2) - Wavefield%Et_I                  ! kinematic FS BC
    !
    ! Surfacic quantities for incident wave field
    PtS_I = Wavefield%Pt_I_s + (rhsE+Wavefield%Et_I)*Wavefield%Pz_I_s
    !
	rhsP = -g*(Wavefield%E+Wavefield%E_I) - half*((Wavefield%Py+PyS_I)**2 &
    	   -(Wavefield%W+Wavefield%Pz_I_s)**2*(one + (Wavefield%Ey+Wavefield%Ey_I)**2)) - PtS_I  ! dynamic FS BC!
!$$$$$$     DEALLOCATE(PyS_I)
ELSE IF (Ny==1) THEN
!$$$$$$     ALLOCATE(PxS_I(Nx,Ny))
  	! Surfacic quantities for incident wave field
    PxS_I = Wavefield%Px_I_s + (Wavefield%Ex+Wavefield%Ex_I)*Wavefield%Pz_I_s
    !
	rhsE = -(Wavefield%Ex+Wavefield%Ex_I)*(Wavefield%Px+PxS_I) &
    	   +(Wavefield%W+Wavefield%Pz_I_s)*(one + (Wavefield%Ex+Wavefield%Ex_I)**2) - Wavefield%Et_I                 ! kinematic FS BC
    !
    ! Surfacic quantities for incident wave field
    PtS_I = Wavefield%Pt_I_s + (rhsE+Wavefield%Et_I)*Wavefield%Pz_I_s
    !
	rhsP = -g*(Wavefield%E+Wavefield%E_I) - half*((Wavefield%Px+PxS_I)**2 &
    	   -(Wavefield%W+Wavefield%Pz_I_s)**2*(one + (Wavefield%Ex+Wavefield%Ex_I)**2)) - PtS_I  ! dynamic FS BC!
!$$$$$$     DEALLOCATE(PxS_I)
ELSE
!$$$$$$     ALLOCATE(PxS_I(Nx,Ny),PyS_I(Nx,Ny))
  	! Surfacic quantities for incident wave field
	PxS_I = Wavefield%Px_I_s + (Wavefield%Ex+Wavefield%Ex_I)*Wavefield%Pz_I_s
    PyS_I = Wavefield%Py_I_s + (Wavefield%Ey+Wavefield%Ey_I)*Wavefield%Pz_I_s
    !
	rhsE = -(Wavefield%Ex+Wavefield%Ex_I)*(Wavefield%Px+PxS_I)-(Wavefield%Ey+Wavefield%Ey_I)*(Wavefield%Py+PyS_I) &
    	+ (Wavefield%W+Wavefield%Pz_I_s)*(one + (Wavefield%Ex+Wavefield%Ex_I)**2 + (Wavefield%Ey+Wavefield%Ey_I)**2) &
        - Wavefield%Et_I            ! kinematic FS BC
    !
    ! Surfacic quantities for incident wave field
    PtS_I = Wavefield%Pt_I_s + (rhsE+Wavefield%Et_I)*Wavefield%Pz_I_s
    !
	rhsP = -g*(Wavefield%E+Wavefield%E_I) - half*((Wavefield%Px+PxS_I)**2 + (Wavefield%Py+PyS_I)**2 &
    	- (Wavefield%W+Wavefield%Pz_I_s)**2*(one + (Wavefield%Ex+Wavefield%Ex_I)**2  &
        + (Wavefield%Ey+Wavefield%Ey_I)**2)) - PtS_I       ! dynamic FS BC!
!$$$$$$     DEALLOCATE(PxS_I,PyS_I)
ENDIF

END SUBROUTINE rhsFreeSurface3D_SWENSE
