SUBROUTINE rhsFreeSurface3D(t,Wavefield,g,rhsE,rhsP,Nx,Ny)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx,Ny
REAL(KIND=long) :: g,t
REAL(KIND=long), DIMENSION(Nx,Ny) :: rhsE, rhsP
TYPE(Wavefield_FS) :: Wavefield

IF (Nx==1) THEN
   rhsE = -Wavefield%Ey*Wavefield%Py + Wavefield%W*(one + Wavefield%Ey**2)                   ! kinematic FS BC
   rhsP = -g*Wavefield%E - half*(Wavefield%Py**2 - Wavefield%W**2*(one + Wavefield%Ey**2))  &
        + Wavefield%P0 - Wavefield%Qr_x  ! dynamic FS BC!
ELSE IF (Ny==1) THEN
	rhsE = -Wavefield%Ex*Wavefield%Px + Wavefield%W*(one + Wavefield%Ex**2 ) - Wavefield%Qr_x               ! kinematic FS BC
	rhsP = -g*Wavefield%E - half*(Wavefield%Px**2 - Wavefield%W**2*(one + Wavefield%Ex**2)) + Wavefield%P0 - Wavefield%Mr_t ! dynamic FS BC!
ELSE
	rhsE = -Wavefield%Ex*Wavefield%Px-Wavefield%Ey*Wavefield%Py + Wavefield%W*(one + Wavefield%Ex**2 + Wavefield%Ey**2)         ! kinematic FS BC
	rhsP = -g*Wavefield%E - half*(Wavefield%Px**2 + Wavefield%Py**2 - Wavefield%W**2*(one + Wavefield%Ex**2 + Wavefield%Ey**2))  &
      + Wavefield%P0 ! dynamic FS BC!
ENDIF
!
! Possible pressure damping terms
!
rhsE = rhsE + Wavefield%NuD
rhsP = rhsP + Wavefield%Pd
END SUBROUTINE rhsFreeSurface3D
