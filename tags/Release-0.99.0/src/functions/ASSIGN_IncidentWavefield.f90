SUBROUTINE ASSIGN_IncidentWavefield(Wavefield_1, Wavefield_2, Nx, Ny)

USE DataTypes
IMPLICIT NONE
!
TYPE(WaveField_FS) :: Wavefield_1, Wavefield_2
INTEGER :: Nx, Ny, nbp
!
Wavefield_1%E_I(1:Nx,1:Ny)=Wavefield_2%E_I(1:Nx,1:Ny); Wavefield_1%Et_I(1:Nx,1:Ny)=Wavefield_2%Et_I(1:Nx,1:Ny)
Wavefield_1%P_I_s(1:Nx,1:Ny)=Wavefield_2%P_I_s(1:Nx,1:Ny); Wavefield_1%Pz_I_s(1:Nx,1:Ny)=Wavefield_2%Pz_I_s(1:Nx,1:Ny)
Wavefield_1%Pt_I_s(1:Nx,1:Ny)=Wavefield_2%Pt_I_s(1:Nx,1:Ny);Wavefield_1%nbp=Wavefield_2%nbp; nbp=Wavefield_2%nbp
Wavefield_1%E_I_bp(1:nbp)=Wavefield_2%E_I_bp(1:nbp);Wavefield_1%Pz_I_bp(1:nbp)=Wavefield_2%Pz_I_bp(1:nbp);
Wavefield_1%GidxTableBP=Wavefield_2%GidxTableBP
IF (Nx>1) THEN
    Wavefield_1%Ex_I(1:Nx,1:Ny)=Wavefield_2%Ex_I(1:Nx,1:Ny); Wavefield_1%Exx_I(1:Nx,1:Ny)=Wavefield_2%Exx_I(1:Nx,1:Ny)
    Wavefield_1%Px_I_s(1:Nx,1:Ny)=Wavefield_2%Px_I_s(1:Nx,1:Ny)
    Wavefield_1%Ex_I_bp(1:nbp)=Wavefield_2%Ex_I_bp(1:nbp); Wavefield_1%Px_I_bp(1:nbp)=Wavefield_2%Px_I_bp(1:nbp)
    Wavefield_1%SourceEx(1:Nx,1:Ny)=Wavefield_2%SourceEx(1:Nx,1:Ny); Wavefield_1%SourcePx(1:Nx,1:Ny)=Wavefield_2%SourcePx(1:Nx,1:Ny)
ENDIF
IF (Ny>1) THEN
    Wavefield_1%Ey_I(1:Nx,1:Ny)=Wavefield_2%Ey_I(1:Nx,1:Ny); Wavefield_1%Eyy_I(1:Nx,1:Ny)=Wavefield_2%Eyy_I(1:Nx,1:Ny)
    Wavefield_1%Py_I_s(1:Nx,1:Ny)=Wavefield_2%Py_I_s(1:Nx,1:Ny)
    Wavefield_1%Ey_I_bp(1:nbp)=Wavefield_2%Ey_I_bp(1:nbp); Wavefield_1%Py_I_bp(1:nbp)=Wavefield_2%Py_I_bp(1:nbp)
    Wavefield_1%SourceEy(1:Nx,1:Ny)=Wavefield_2%SourceEy(1:Nx,1:Ny); Wavefield_1%SourcePy(1:Nx,1:Ny)=Wavefield_2%SourcePy(1:Nx,1:Ny)
ENDIF

END SUBROUTINE ASSIGN_IncidentWavefield
