SUBROUTINE RungeKutta(RKFUNC,RHSFUNC,time,dt,WaveField,RHSE,RHSP)
  
  !-----------------------------------------------------------------------------
  ! Global parameters:
  !-----------------------------------------------------------------------------
  USE GlobalVariables, ONLY: long,FineGrid,GhostGridX,GhostGridY,WaveField_FS

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! External subroutines:
  !-----------------------------------------------------------------------------
  EXTERNAL RKFUNC,RHSFUNC

  !-----------------------------------------------------------------------------
  ! Input/output parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long)       :: time,dt
  TYPE(WaveField_FS)    :: WaveField
  REAL(KIND=long),DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: RHSE,RHSP

  !-----------------------------------------------------------------------------
  ! Call the Runge Kutta function pointed to by RKFUNC with the right hand side
  ! function pointed to by RHSFUNC:
  !-----------------------------------------------------------------------------
  CALL RKFUNC(RHSFUNC,time,dt,WaveField,RHSE,RHSP);

  END SUBROUTINE RungeKutta
