SUBROUTINE random_wave_coefficients_3D(i_spec, nt, beta0, dt, dx, Tp, Hs, depth, &
     grav, inc_wave_file, kh_max, seed, seed2, eta0, beta, n_cut)
  !-----------------------------------------------------------------------
  !
  ! Build the coefficients for a pseudo-random 3D wave with direction BETA0 to the 
  ! x-axis, over NT time steps at tim-step DT.  The water depth at the wavemaker 
  ! is DEPTH.  The wave can be based on a long-crested P-M (i_spec=30) or 
  ! JONSWAP (ispec=31) spectrum at angle BETA0 to the x-axis, or it can be 
  ! a JONSWAP with NORMAL directional spreading about the angle BETA0 (i_spec=33).  
  ! The coefficients are returned in COEF ready to be moved around in space and 
  ! FFT'ed by the subroutine random_wave_signal_3D.f90. 
  ! 
  ! When i_spec=-30, a mono-chromatic wave is generated with period Tp, 
  ! height Hs and angle BETA0; so nothing is done here.  
  !
  !----------------------------------------------------------------------
  IMPLICIT none
  CHARACTER(len=30) inc_wave_file, header
  integer, parameter :: long=selected_real_kind(12,99)
  integer i, i_spec, seed, seed2, nt, n_cut
  real(kind=long) :: Tp, Hs, dt, depth, grav, kh_max, dx
  real(kind=long) :: eta0(nt), beta(nt)

  ! Local variables
  real(kind=long) :: factor, domega, beta0, spec, zero=0._long,  one=1._long, &
       two=2._long, pi, twopi
  !
  ! Some parameters 
  !
  pi=acos(-one)
  twopi=two*pi
  !
  ! 
  ! Generate the random wave.
  !
  ! The FT factor. 
  !
  factor = two/nt
  !
  ! The Fourier coefficients are built based on SI dimensional numbers.  Omega is
  ! radian frequency.
  !
  domega = twopi / (nt * dt)
  !
  !
  ! Generate the Fourier components of the wave elevation at the origin.
  ! Note that the same spectral parameters and SEED values always give the same
  ! incident wave.
  !
  IF(i_spec==30)THEN 
     ! A P-M spectrum
     CALL build_coeff(eta0, nt, Hs, Tp, dt, seed, seed2, 0)
  ELSEIf(i_spec==31)THEN
     ! A JONSWAP spectrum
     CALL build_coeff(eta0, nt, Hs, Tp, dt, seed, seed2, 1)
  ELSEIF(i_spec==33)THEN
     ! A JONSWAP spectrum with a normal spreading
     CALL build_coeff_3D(eta0, beta, nt, Hs, Tp, dt, seed, seed2, beta0)
  ELSEIF(i_spec==-30)THEN
     ! Mono-chromatic waves, nothing to do here.  
     RETURN
  END IF
  !
  ! We truncate the spectrum at k_max.
  !
  n_cut=nint(sqrt(grav*kh_max/depth*tanh(kh_max))/domega)+1
  Print *, 'Truncating the spectrum at T=',twopi/((n_cut-1)*domega), &
       ' which corresponds to kh=',kh_max
  print *, 'The deep water resolution of this wave is ',2*pi*depth/(kh_max*dx)+one, &
       ' points per wavelength.'
  Print *, ' '
  eta0(2*n_cut-1:nt)=zero
  !
  !Print out the wave maker elevation.
  !
  CALL drealft (eta0, nt, -1)

  open(21,file='eta0_irregular',status='unknown')
  IF(i_spec==30)THEN
     write(21,32)Hs,Tp,twopi/((n_cut-1)*domega),kh_max,beta0
31   format('% P-M spectrum H_s=',e12.4,' T_p=',e12.4,' truncated at T=',e12.4,&
          ' which is kh=',e12.4,'.  Heading angle is ,',f10.2)
  ELSEIF(i_spec==31)THEN
     write(21,32)Hs,Tp,twopi/((n_cut-1)*domega),kh_max,beta0
32   format('% JONSWAP spectrum H_s=',e12.4,' T_p=',e12.4,' truncated at T=',e12.4,&
          ' which is kh=',e12.4,'.  Heading angle is ,',f10.2)
  ELSEIF(i_spec==33)THEN
     write(21,33)Hs,Tp,twopi/((n_cut-1)*domega),kh_max,beta0
33   format('% JONSWAP with NORMAL spreading H_s=',e12.4,' T_p=',e12.4, &
          ' truncated at T=',e12.4,&
          ' which is kh=',e12.4,'.  Heading angle is ,',f10.2)
  END IF
  do i=1,nt
     write(21,*)(i-1)*dt,eta0(i)
  end do
  close(21)
  !
  ! FFT back to the coefficients and print them out
  !
  CALL drealft (eta0, nt, 1)
  eta0=eta0*factor
  open(21,file='eta0_coeffs',status='unknown')
  write(21,34)
34 format('# f [Hz], mag(c), c_r, c_i')
  write(21,35)zero,eta0(1),eta0(1),eta0(2)
35 format(4e16.5)
  do i=3,nt-1,2
     write(21,35) (i-1)/2 * domega/twopi, sqrt(eta0(i)**2+eta0(i+1)**2),eta0(i),eta0(i+1)
  end do
  close(21)

RETURN
END SUBROUTINE random_wave_coefficients_3D
