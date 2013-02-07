SUBROUTINE ALLOCATE_Wavefield_Type(Wavefield, Nx, Ny, Nz, GhostGridX, GhostGridy, GhostGridZ, swenseONOFF)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
!
INTEGER :: Nx, Ny, Nz, GhostGridX, GhostGridy, GhostGridZ, swenseONOFF
TYPE(WaveField_FS) :: Wavefield
!
!Local variables
INTEGER :: Nx_ext, Ny_ext, Nz_ext, nbp
!
Nx_ext = Nx+2*GhostGridX
Ny_ext = Ny+2*GhostGridY
Nz_ext = Nz+GhostGridZ
! Probably necessary...
NULLIFY(Wavefield%E,Wavefield%Ex,Wavefield%Exx,Wavefield%Ey,Wavefield%Eyy,Wavefield%P0)
NULLIFY(Wavefield%P,Wavefield%Px,Wavefield%Py,Wavefield%W,Wavefield%Qr_x,Wavefield%Mr_t)
NULLIFY(Wavefield%E_I,Wavefield%Ex_I,Wavefield%Exx_I,Wavefield%Ey_I,Wavefield%Eyy_I)
NULLIFY(Wavefield%EtatHist, Wavefield%WHist, Wavefield%NuD, Wavefield%Pd)
NULLIFY(Wavefield%P_I_s,Wavefield%Px_I_s,Wavefield%Py_I_s,Wavefield%Pz_I_s)
NULLIFY(Wavefield%E_I_bp,Wavefield%Ex_I_bp,Wavefield%Ey_I_bp)
NULLIFY(Wavefield%Px_I_bp,Wavefield%Py_I_bp,Wavefield%Pz_I_bp)
NULLIFY(Wavefield%SourceEx,Wavefield%SourcePx)
NULLIFY(Wavefield%SourceEy,Wavefield%SourcePy)
NULLIFY(Wavefield%GidxTableBP)
!
ALLOCATE(Wavefield%E(Nx_ext,Ny_ext), Wavefield%W(Nx_ext,Ny_ext), Wavefield%P(Nx_ext,Ny_ext),  &
     Wavefield%P0(Nx_ext,Ny_ext),Wavefield%Qr_x(Nx_ext,Ny_ext),Wavefield%Mr_t(Nx_ext,Ny_ext), &
     Wavefield%NuD(Nx_ext,Ny_ext), Wavefield%Pd(Nx_ext,Ny_ext))
Wavefield%E=zero; Wavefield%P=zero; Wavefield%W=zero; Wavefield%P0=zero; Wavefield%Qr_x=zero;
Wavefield%Mr_t=zero; Wavefield%NuD=zero; Wavefield%Pd=zero
IF (Nx>1) THEN
	ALLOCATE(Wavefield%Ex(Nx_ext,Ny_ext), Wavefield%Exx(Nx_ext,Ny_ext), Wavefield%Px(Nx_ext,Ny_ext))
	Wavefield%Ex=zero; Wavefield%Exx=zero; Wavefield%Px=zero
ENDIF
IF (Ny>1) THEN
    ALLOCATE(Wavefield%Ey(Nx_ext,Ny_ext), Wavefield%Eyy(Nx_ext,Ny_ext), Wavefield%Py(Nx_ext,Ny_ext))
	Wavefield%Ey=zero; Wavefield%Eyy=zero; Wavefield%Py=zero
ENDIF

ALLOCATE(Wavefield%EtatHist(Nx_ext,Ny_ext,3), Wavefield%WHist(Nx_ext,Ny_ext,3))  ! -hbb Needs extension to 3D.
! ALLOCATE(Wavefield%WHist(Nx_ext,Ny_ext,3))  ! SIGS
! Wavefield%WHist=zero;

IF(swenseONOFF==1) THEN !We allocate the incident wavefields values
    ALLOCATE(Wavefield%E_I(Nx_ext,Ny_ext), Wavefield%Et_I(Nx_ext,Ny_ext), Wavefield%P_I_s(Nx_ext,Ny_ext))
    ALLOCATE(Wavefield%Pz_I_s(Nx_ext,Ny_ext), Wavefield%Pt_I_s(Nx_ext,Ny_ext))
    ! Initialize all variables to zero
    Wavefield%E_I=zero; Wavefield%Et_I=zero; Wavefield%P_I_s=zero; Wavefield%Pz_I_s=zero; Wavefield%Pt_I_s=zero
    IF (Nx>1) THEN
        ALLOCATE(Wavefield%Ex_I(Nx_ext,Ny_ext), Wavefield%Exx_I(Nx_ext,Ny_ext), Wavefield%Px_I_s(Nx_ext,Ny_ext))
        Wavefield%Ex_I=zero; Wavefield%Exx_I=zero; Wavefield%Px_I_s=zero
    ENDIF
    IF (Ny>1) THEN
        ALLOCATE(Wavefield%Ey_I(Nx_ext,Ny_ext), Wavefield%Eyy_I(Nx_ext,Ny_ext), Wavefield%Py_I_s(Nx_ext,Ny_ext))
        Wavefield%Ey_I=zero; Wavefield%Eyy_I=zero; Wavefield%Py_I_s=zero
    ENDIF
    ! Allocate boundary points...
    ! FIXME: treat the 2D cases differently ?
    IF ((Nx>1).AND.(Ny>1)) THEN
        ! Determine number of boundary points
        nbp=2*Nx*Nz_ext+2*Ny*Nz_ext+Nx*Ny
        Wavefield%nbp=nbp
        ALLOCATE(Wavefield%E_I_bp(nbp), Wavefield%Ex_I_bp(nbp), Wavefield%Ey_I_bp(nbp))
        ALLOCATE(Wavefield%Px_I_bp(nbp), Wavefield%Py_I_bp(nbp), Wavefield%Pz_I_bp(nbp))
        Wavefield%E_I_bp=zero; Wavefield%Ex_I_bp=zero; Wavefield%Ey_I_bp=zero
        Wavefield%Px_I_bp=zero; Wavefield%Py_I_bp=zero; Wavefield%Pz_I_bp=zero
    ELSEIF ((Nx>1).AND.(Ny==1)) THEN
        ! Determine number of boundary points
        nbp=2*Nz_ext+Nx
        Wavefield%nbp=nbp
        ALLOCATE(Wavefield%E_I_bp(nbp), Wavefield%Ex_I_bp(nbp))
        ALLOCATE(Wavefield%Px_I_bp(nbp), Wavefield%Pz_I_bp(nbp))
        Wavefield%E_I_bp=zero; Wavefield%Ex_I_bp=zero; Wavefield%Px_I_bp=zero; Wavefield%Pz_I_bp=zero
    ELSEIF ((Nx==1).AND.(Ny>1)) THEN
        ! Determine number of boundary points
        nbp=2*Nz_ext+Ny
        Wavefield%nbp=nbp
        ALLOCATE(Wavefield%E_I_bp(nbp), Wavefield%Ey_I_bp(nbp))
        ALLOCATE(Wavefield%Py_I_bp(nbp), Wavefield%Pz_I_bp(nbp))
        Wavefield%E_I_bp=zero; Wavefield%Ey_I_bp=zero; Wavefield%Py_I_bp=zero; Wavefield%Pz_I_bp=zero
    ELSE
        Wavefield%nbp=0
        STOP
    ENDIF
    ! Index Table of Boundary Points
    ALLOCATE(Wavefield%GidxTableBP(Nz_ext,Nx_ext,Ny_ext))
    Wavefield%GidxTableBP=0
ENDIF
IF (Nx>1) THEN
    ALLOCATE(Wavefield%SourceEx(Nx_ext,Ny_ext), Wavefield%SourcePx(Nx_ext,Ny_ext))
    Wavefield%SourceEx=zero; Wavefield%SourcePx=zero
ENDIF
IF (Ny>1) THEN
    ALLOCATE(Wavefield%SourceEy(Nx_ext,Ny_ext), Wavefield%SourcePy(Nx_ext,Ny_ext))
    Wavefield%SourceEy=zero; Wavefield%SourcePy=zero
ENDIF
END SUBROUTINE ALLOCATE_Wavefield_Type
