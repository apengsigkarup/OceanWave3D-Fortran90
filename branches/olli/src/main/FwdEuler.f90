SUBROUTINE FwdEuler(RHSFUNC,time,dt,WaveField,RHSE,RHSP)
  
  !-----------------------------------------------------------------------------
  ! Global parameters:
  !-----------------------------------------------------------------------------
  USE GlobalVariables, ONLY: long,FineGrid,ghostgridx,ghostgridy,ghostgridz,&
  WaveField_FS

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! External subroutines:
  !-----------------------------------------------------------------------------
  EXTERNAL RHSFUNC

  !-----------------------------------------------------------------------------
  ! Input/output parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long)       :: time,dt
  TYPE(WaveField_FS)    :: WaveField
  REAL(KIND=long), DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: RHSE,RHSP

  !-----------------------------------------------------------------------------
  ! Summation of right hand sides:
  !-----------------------------------------------------------------------------
  WaveField%E = WaveField%E + dt*RHSE
  WaveField%P = WaveField%P + dt*RHSP

  !-----------------------------------------------------------------------------
  ! Right hand side for the new stage and initial condition for next time step:
  !-----------------------------------------------------------------------------
  print*,"FwdEuler: before RHS"
  CALL RHSFUNC(time+dt,WaveField,RHSE,RHSP);
  print*,"FwdEuler: after RHS"

  END SUBROUTINE FwdEuler
