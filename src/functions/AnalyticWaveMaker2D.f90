SUBROUTINE AnalyticWaveMaker2D(izone,i0,i1,x,RKtime,time,time0,Ea,Pa)
!
! A subroutine to feed in an exact incident wave solution at the positions x(1:nx) and time t.
!
! By Allan P. Engsig-Karup.
USE GlobalVariables, ONLY: SFsol, g, dt, RandomWave, IncWaveType, half, one
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: izone, i0, i1, nx, it0, itime
REAL(KIND=long), DIMENSION(i1-i0+1) :: x, Ea, Pa
REAL(KIND=long) :: RKtime, time, time0, tol, diff
!
!
nx=i1-i0+1; tol=10E-14
IF (IncWaveType==1)THEN
   CALL stream_func_wave_finite(nx,x,RKtime, &
        SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,Ea,Pa)
ELSEIF(IncWaveType==2)THEN
   diff=(RKtime-time)/dt
   If(diff<tol .or. abs(one-diff)<tol)THEN
      itime = nint(RKtime/dt)+1
      Ea=RandomWave(izone)%eta(1:nx,itime); Pa=RandomWave(izone)%Phis(1:nx,itime)
   ELSE
      itime=nint(time/dt)+1
      Ea=half*(RandomWave(izone)%eta(1:nx,itime) + RandomWave(izone)%eta(1:nx,itime+1)) 
      Pa=half*(RandomWave(izone)%Phis(1:nx,itime) + RandomWave(izone)%Phis(1:nx,itime) )
   END If
END IF

RETURN

END SUBROUTINE AnalyticWaveMaker2D
