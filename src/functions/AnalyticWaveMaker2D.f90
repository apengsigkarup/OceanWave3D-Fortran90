SUBROUTINE AnalyticWaveMaker2D(nx,x,time,Ea,Pa)
!
! A subroutine to feed in an exact incident wave solution at the positions x(1:nx) and time t.
!
USE GlobalVariables, ONLY: SFsol, g, dt, RandomWave, IncWaveType
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: nx, itime
REAL(KIND=long), DIMENSION(nx) :: x, Ea, Pa
REAL(KIND=long) :: time
!
!
IF (IncWaveType==1)THEN
   CALL stream_func_wave_finite(nx,x,time, &
        SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,Ea,Pa)
ELSEIF(IncWaveType==2)THEN
   itime = time/dt+1
   Ea=RandomWave%eta(:,itime); Pa=RandomWave%Phis(:,itime)
END IF

RETURN

END SUBROUTINE AnalyticWaveMaker2D
