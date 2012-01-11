SUBROUTINE DEALLOCATE_Wavefield_Type(Wavefield, Nx, Ny, Nz, swenseONOFF)

USE DataTypes
IMPLICIT NONE
!
INTEGER :: Nx, Ny, Nz, swenseONOFF
TYPE(WaveField_FS) :: Wavefield
!
!
DEALLOCATE(Wavefield%E, Wavefield%W, Wavefield%P)
IF (Nx>1) THEN
	DEALLOCATE(Wavefield%Ex, Wavefield%Exx, Wavefield%Px)
ENDIF
IF (Ny>1) THEN
    DEALLOCATE(Wavefield%Ey, Wavefield%Eyy, Wavefield%Py)
ENDIF

IF(swenseONOFF==1) THEN !We allocate the incident wavefields values
    DEALLOCATE(Wavefield%E_I, Wavefield%Et_I, Wavefield%P_I_s)
    DEALLOCATE(Wavefield%Pz_I_s, Wavefield%Pt_I_s)
    IF (Nx>1) THEN
        DEALLOCATE(Wavefield%Ex_I, Wavefield%Exx_I, Wavefield%Px_I_s)
    ENDIF
    IF (Ny>1) THEN
        DEALLOCATE(Wavefield%Ey_I, Wavefield%Eyy_I, Wavefield%Py_I_s)
    ENDIF
    ! Allocate boundary points...
    ! FIXME: treat the 2D cases differently ?
    IF ((Nx>1).AND.(Ny>1)) THEN
        DEALLOCATE(Wavefield%E_I_bp, Wavefield%Ex_I_bp, Wavefield%Ey_I_bp)
        DEALLOCATE(Wavefield%Px_I_bp, Wavefield%Py_I_bp, Wavefield%Pz_I_bp)
    ELSEIF ((Nx>1).AND.(Ny==1)) THEN
        DEALLOCATE(Wavefield%E_I_bp, Wavefield%Ex_I_bp)
        DEALLOCATE(Wavefield%Px_I_bp, Wavefield%Pz_I_bp)
    ELSEIF ((Nx==1).AND.(Ny>1)) THEN
        DEALLOCATE(Wavefield%E_I_bp, Wavefield%Ey_I_bp)
        DEALLOCATE(Wavefield%Py_I_bp, Wavefield%Pz_I_bp)
    ENDIF
    DEALLOCATE(Wavefield%GidxTableBP)
ENDIF
IF (Nx>1) THEN
    DEALLOCATE(Wavefield%SourceEx, Wavefield%SourcePx)
ENDIF
IF (Ny>1) THEN
    DEALLOCATE(Wavefield%SourceEy, Wavefield%SourcePy)
ENDIF
!
END SUBROUTINE DEALLOCATE_Wavefield_Type
