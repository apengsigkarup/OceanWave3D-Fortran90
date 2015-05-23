SUBROUTINE AnalyticWaveMaker2D(izone,i0,i1,j,x,RKtime,time,Ea,Pa)
!
! A subroutine to feed in an exact incident wave solution at the positions x(1:nx) and time t.
! time is the beginning of this time-interval while RKtime is the exact time position 
! of this R-K stage.  Mid-interval times are assumed to be at t0+dt/2 and are linearly 
! interpolated between the endpoint values.  j is the y-index of this x-grid line.  
!
! By Allan P. Engsig-Karup and Harry B. Bingham.
!
USE GlobalVariables, ONLY: SFsol, g, dt, RandomWave, IncWaveType, half, one,dt0
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: izone, i0, i1, j, k, nx, it0, itime
REAL(KIND=long), DIMENSION(i1-i0+1) :: x, Ea, Pa
REAL(KIND=long) :: RKtime, time, tol, diff
!
!
nx=i1-i0+1; tol=10E-14
k=(j-1)*nx+1   ! The starting index of this line for 3D generation.  
IF (IncWaveType==1)THEN
   CALL stream_func_wave_finite(nx,x,RKtime, &
        SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,Ea,Pa)
ELSEIF(IncWaveType==2)THEN
   diff=(RKtime-time)/dt0
   If(diff<tol .or. abs(one-diff)<tol)THEN
      itime = nint(RKtime/dt0)+1
      Ea=RandomWave(izone)%eta(k:k+nx-1,itime); Pa=RandomWave(izone)%Phis(k:k+nx-1,itime)
   ELSE
      itime=nint(time/dt0)+1
      Ea=half*( RandomWave(izone)%eta(k:k+nx-1,itime) + RandomWave(izone)%eta(k:k+nx-1,itime+1)  ) 
      Pa=half*( RandomWave(izone)%Phis(k:k+nx-1,itime) + RandomWave(izone)%Phis(k:k+nx-1,itime+1) )
   END If
END IF

RETURN

END SUBROUTINE AnalyticWaveMaker2D
