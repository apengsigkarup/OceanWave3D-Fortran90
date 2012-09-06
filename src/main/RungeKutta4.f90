SUBROUTINE RungeKutta4(RHSFUNC,time,dt,WaveField,RHSE1,RHSP1)
  
  !-----------------------------------------------------------------------------
  ! Global parameters:
  !-----------------------------------------------------------------------------
  USE GlobalVariables, ONLY: long,FineGrid,ghostgridx,ghostgridy,ghostgridz,&
  WaveField_FS,swenseONOFF,WaveField_tmp

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! External subroutines:
  !-----------------------------------------------------------------------------
  EXTERNAL RHSFUNC,ALLOCATE_WaveField_Type

  !-----------------------------------------------------------------------------
  ! Input parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long)       :: time,dt
  TYPE(WaveField_FS)    :: WaveField
  REAL(KIND=long),DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: RHSE1,RHSP1

  !-----------------------------------------------------------------------------
  ! Local parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long)       ::  time_tmp
  REAL(KIND=long),DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: RHSE2,RHSE3,RHSE4,RHSP2,RHSP3,RHSP4

  !-----------------------------------------------------------------------------
  ! Allocation:
  !-----------------------------------------------------------------------------
  IF (.NOT.ASSOCIATED(Wavefield_tmp%E)) THEN
    CALL ALLOCATE_Wavefield_Type(Wavefield_tmp,FineGrid%Nx,FineGrid%Ny,&
    FineGrid%Nz,GhostGridX,GhostGridy,GhostGridZ,swenseONOFF)
  ENDIF
  WaveField_tmp%P0   =   WaveField%P0

  !-----------------------------------------------------------------------------
  ! Call the right hand side function for each stage:
  !-----------------------------------------------------------------------------
  time_tmp          =   time        + dt/2.0
  WaveField_tmp%E   =   WaveField%E + dt/2.0*RHSE1
  WaveField_tmp%P   =   WaveField%P + dt/2.0*RHSP1
  CALL RHSFUNC(time_tmp,WaveField_tmp,RHSE2,RHSP2)
  time_tmp          =   time        + dt/2.0
  WaveField_tmp%E   =   WaveField%E + dt/2.0*RHSE2
  WaveField_tmp%P   =   WaveField%P + dt/2.0*RHSP2
  CALL RHSFUNC(time_tmp,WaveField_tmp,RHSE3,RHSP3)
  time_tmp          =   time        + dt
  WaveField_tmp%E   =   WaveField%E + dt*RHSE3
  WaveField_tmp%P   =   WaveField%P + dt*RHSP3
  CALL RHSFUNC(time_tmp,WaveField_tmp,RHSE4,RHSP4)

  !-----------------------------------------------------------------------------
  ! Summation of right hand sides:
  !-----------------------------------------------------------------------------
  WaveField%E = WaveField%E + dt/6.0*(RHSE1 + 2.0*RHSE2 + 2.0*RHSE3 + RHSE4)
  WaveField%P = WaveField%P + dt/6.0*(RHSP1 + 2.0*RHSP2 + 2.0*RHSP3 + RHSP4)

  !-----------------------------------------------------------------------------
  ! Right hand side for the new stage and initial condition for next time step:
  !-----------------------------------------------------------------------------
  time_tmp = time + dt
  CALL RHSFUNC(time_tmp,WaveField,RHSE1,RHSP1)

  END SUBROUTINE RungeKutta4
