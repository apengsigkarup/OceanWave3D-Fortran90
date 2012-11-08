SUBROUTINE AnalyticWaveMaker2D(i0,i1,x,time,time0,Ea,Pa)
!
! A subroutine to feed in an exact incident wave solution at the positions x(1:nx) and time t.
!
USE GlobalVariables, ONLY: SFsol, g, dt, RandomWave, IncWaveType
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: i0, i1, nx, itime
REAL(KIND=long), DIMENSION(i1-i0+1) :: x, Ea, Pa
REAL(KIND=long) :: time, time0
!
!
nx=i1-i0+1
IF (IncWaveType==1)THEN
   CALL stream_func_wave_finite(nx,x,time, &
        SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,Ea,Pa)
ELSEIF(IncWaveType==2)THEN
   itime = (time-time0)/dt+1
   Ea=RandomWave%eta(i0:i1,itime); Pa=RandomWave%Phis(i0:i1,itime)
END IF

RETURN

END SUBROUTINE AnalyticWaveMaker2D
