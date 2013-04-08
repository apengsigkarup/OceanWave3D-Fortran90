SUBROUTINE DEALLOCATE_Wavefield_Type(Wavefield, Nx, Ny, Nz, swenseONOFF)

USE DataTypes
IMPLICIT NONE
!
INTEGER :: Nx, Ny, Nz, swenseONOFF
INTEGER :: IERR7
TYPE(WaveField_FS) :: Wavefield
!
!
DEALLOCATE(Wavefield%E, Wavefield%W, Wavefield%P, STAT = IERR7)
IF (Nx>1) THEN
	DEALLOCATE(Wavefield%Ex, Wavefield%Exx, Wavefield%Px, STAT = IERR7)
ENDIF
IF (Ny>1) THEN
    DEALLOCATE(Wavefield%Ey, Wavefield%Eyy, Wavefield%Py, STAT = IERR7)
ENDIF

IF(swenseONOFF==1) THEN !We allocate the incident wavefields values
    DEALLOCATE(Wavefield%E_I, Wavefield%Et_I, Wavefield%P_I_s, STAT = IERR7)
    DEALLOCATE(Wavefield%Pz_I_s, Wavefield%Pt_I_s, STAT = IERR7)
    IF (Nx>1) THEN
        DEALLOCATE(Wavefield%Ex_I, Wavefield%Exx_I, Wavefield%Px_I_s, STAT = IERR7)
    ENDIF
    IF (Ny>1) THEN
        DEALLOCATE(Wavefield%Ey_I, Wavefield%Eyy_I, Wavefield%Py_I_s, STAT = IERR7)
    ENDIF
    ! Allocate boundary points...
    ! FIXME: treat the 2D cases differently ?
    IF ((Nx>1).AND.(Ny>1)) THEN
        DEALLOCATE(Wavefield%E_I_bp, Wavefield%Ex_I_bp, Wavefield%Ey_I_bp, STAT = IERR7)
        DEALLOCATE(Wavefield%Px_I_bp, Wavefield%Py_I_bp, Wavefield%Pz_I_bp, STAT = IERR7)
    ELSEIF ((Nx>1).AND.(Ny==1)) THEN
        DEALLOCATE(Wavefield%E_I_bp, Wavefield%Ex_I_bp, STAT = IERR7)
        DEALLOCATE(Wavefield%Px_I_bp, Wavefield%Pz_I_bp, STAT = IERR7)
    ELSEIF ((Nx==1).AND.(Ny>1)) THEN
        DEALLOCATE(Wavefield%E_I_bp, Wavefield%Ey_I_bp, STAT = IERR7)
        DEALLOCATE(Wavefield%Py_I_bp, Wavefield%Pz_I_bp, STAT = IERR7)
    ENDIF
    DEALLOCATE(Wavefield%GidxTableBP, STAT = IERR7)
ENDIF
IF (Nx>1) THEN
    DEALLOCATE(Wavefield%SourceEx, Wavefield%SourcePx, STAT = IERR7)
ENDIF
IF (Ny>1) THEN
    DEALLOCATE(Wavefield%SourceEy, Wavefield%SourcePy, STAT = IERR7)
ENDIF
!
END SUBROUTINE DEALLOCATE_Wavefield_Type
