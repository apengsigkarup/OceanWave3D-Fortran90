SUBROUTINE stream_func_set_up(grav,depth,period,i_wavel_or_per,wave_l,    &
     wavenum,wave_h,i_euler_or_stokes,e_or_s_vel,i_deep_or_finite,        &
     n_h_steps,n_four_modes,nwrk,yy,zz)

! From HBB: utilities.f90

  USE Precision

  IMPLICIT NONE

  INTEGER i_wavel_or_per, i_euler_or_stokes, i_deep_or_finite, n_h_steps, &
       n_four_modes, nwrk
  REAL(KIND=long) :: grav, depth, period, wave_l, wavenum, e_or_s_vel,      &
       wave_h, yy(*), zz(*)

! Workspace for the stream function coefficient calculation
  CHARACTER(len=10) :: deep_or_finite, wavel_or_per, Euler_or_Stokes
  INTEGER :: iwrk(nwrk)
  REAL(KIND=long) :: hoverd, height, zero=0._long
  REAL(KIND=long) :: wrk1(0:nwrk), wrk2(0:nwrk), wrk3(nwrk), wrk4(nwrk,2), &
       wrk5(nwrk), wrk6(nwrk), wrk7(nwrk,nwrk), wrk8(nwrk)

  IF((WAVENUM*DEPTH).GT.9.0D0) THEN
     I_DEEP_OR_FINITE=0
     deep_or_finite='deep'
     hoverd=zero
  ELSE
     I_DEEP_OR_FINITE=1
     deep_or_finite='finite'
     hoverd=wave_H/depth
  ENDIF
  IF(i_Euler_or_Stokes==0)THEN
     Euler_or_Stokes='Euler'
  ELSE
     Euler_or_Stokes='Stokes'
  END IF
  IF(i_wavel_or_per==0)THEN
     wavel_or_per='wavelength'
     height=wave_H/wave_l
  ELSE
     wavel_or_per='period'
     height=wave_H/(grav*period**2)
  END IF
  E_or_S_vel = E_or_S_vel/(grav*wave_H)**0.5
! Compute the Stream Function coefficients.
! This subroutine is basically Fenton's
! original code, with output turned off, turned into a subroutine.
! On return the zz's and the yy's hold the solution.
  call stream_func_coeffs(deep_or_finite,hoverd,wavel_or_per,     &
       height,Euler_or_Stokes,E_or_S_vel,n_four_modes,n_H_steps,  &
       nwrk,wrk1,wrk2,wrk3,wrk4,wrk5,wrk6,wrk7,wrk8,iwrk,zz,yy)
  If(i_wavel_or_per==0)THEN
     period=zz(3)/sqrt(grav*wavenum)
  ELSE IF(i_wavel_or_per==1)THEN
     wavenum=zz(1)/depth
  END If

  RETURN
END SUBROUTINE stream_func_set_up
