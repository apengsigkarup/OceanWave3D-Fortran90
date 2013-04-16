SUBROUTINE random_wave_signal_3D(i_spec, nt, nx, ny, beta0, x0, y0, x, y, dt, Tp, Hs, depth, &
     grav, inc_wave_file, kh_max, seed, seed2, eta, Vs, time0)
  !-----------------------------------------------------------------------
  !
  ! Generate a pseudo-random long crested wave with direction beta to the 
  ! x-axis, over NT time steps at NX by NY points over the computational
  ! free-surface.  The x-positions of the points in the wave generation 
  ! zone are given by x(1:nx,1:ny) and y(1:nx,1:ny), 
  ! and the time-step dt is uniform.  The water depth at the wavemaker (x0,y0) 
  ! is DEPTH.  The wave can be based on a long-crested P-M (i_spec=30) or 
  ! JONSWAP (ispec=31) spectrum at angle BETA0 to the x-axis, or it can be 
  ! a JONSWAP with NORMAL directional spreading about the angle beta (i_spec=33).  
  ! The wave is moved around in space using linear theory in the frequency 
  ! domain.  eta & Vs are returned with the time series
  ! of elevation and the potential at z=0 in the wavemaker zone.  
  ! 
  ! When i_spec=-30, a mono-chromatic wave is generated with period Tp, 
  ! height Hs and angle BETA.  
  !
  !----------------------------------------------------------------------
  IMPLICIT none
  CHARACTER(len=30) inc_wave_file, header
  integer, parameter :: long=selected_real_kind(12,99)
  integer i_spec, seed, seed2, nt, nx, ny
  real(kind=long) :: Tp, Hs, dt, depth, time0, x0, y0, grav
  real(kind=long) :: x(nx,ny), y(nx,ny), Vs(nx,ny,nt,*), eta(nx,ny,nt)

  ! Local variables
  INTEGER :: i, j, k, ns_inc, n_cut, ndat
  real(kind=long) :: cosb, sinb, factor, domega, realpt, beta0, dx,       &
       imagpt, spec, phase, nu, coslnu, sinlnu, magx2, x2(3), x1(2),     &
       kinf, omega, dist, realarg, imagarg, reala, imaga,                &
       zero=0._long,  one=1._long, two=2._long, pi,twopi,udum=1.0_long,  &
       velfact, tanhkhi, dt_inc, dum, kh_max, phifact, kh, t, amp
  !
  ! Local workspace
  !
  real(kind=long),Allocatable :: eta0(:), work_eta(:), work_vel(:), beta(:)
  allocate( eta0(nt), work_eta(nt), work_vel(nt), beta(nt) )
  !
  !
  ! Some parameters 
  !
  pi=acos(-one)
  twopi=two*pi
  !
  IF(i_spec >= 0)THEN
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
     END IF
    
     ! We truncate the spectrum at k_max.
     n_cut=nint(sqrt(grav*kh_max/depth*tanh(kh_max))/domega)+1
     Print *, 'Truncating the spectrum at T=',twopi/((n_cut-1)*domega), &
          ' which corresponds to kh=',kh_max
     dx=max(x(2,1)-x(1,1),y(1,2)-y(1,1))
     print *, 'The deep water resolution of this wave is ',2*pi*depth/(kh_max*dx)+one, &
          ' points per wavelength.'
     Print *, ' '
     eta0(2*n_cut-1:nt)=zero
     !hbb Print out the elevations at the origin to check.
     CALL drealft (eta0, nt, -1)
     !
     open(21,file='eta0_irregular',status='unknown')
     IF(i_spec==30)THEN
        write(21,32)Hs,Tp,twopi/((n_cut-1)*domega),kh_max
31      format('% P-M spectrum H_s=',e12.4,' T_p=',e12.4,' truncated at T=',e12.4,&
             ' which is kh=',e12.4)
     ELSEIF(i_spec==31)THEN
        write(21,32)Hs,Tp,twopi/((n_cut-1)*domega),kh_max
32      format('% JONSWAP spectrum H_s=',e12.4,' T_p=',e12.4,' truncated at T=',e12.4,&
             ' which is kh=',e12.4)
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
34   format('# f [Hz], mag(c), c_r, c_i')
     write(21,35)zero,eta0(1),eta0(1),eta0(2)
35   format(4e16.5)
     do i=3,nt-1,2
        write(21,35) (i-1)/2 * domega/twopi, sqrt(eta0(i)**2+eta0(i+1)**2),eta0(i),eta0(i+1)
     end do
     close(21)
     !
     ! Translate these Fourier amplitudes around in space to generate the
     ! eta and Vs coefficients over some portion of the free-surface.
     !
     DO k=1,ny
        DO j = 1, nx
           !
           ! Take care of omega=0 and  the Nyquist frequency.
           !
           work_eta(1) = zero
           work_eta(2) = zero
           work_vel(1) = zero
           work_vel(2) = zero
           !
           ! The rest of frequency space.
           !
           DO i = 2, n_cut
              omega = (i - 1) * domega
              kinf = omega * omega / grav
              IF (depth.gt.zero) then
                 ! Calling this routine with udum=1. gives the SI dimensional wavenumber.
                 CALL fdwvnum (depth, omega, kinf, udum, grav, nu)
              ELSE
                 nu = kinf
              ENDIF
              !
              ! Compute the distance (dist) from the origin of the coordinates (x0) to
              ! the point where the wave elevation is to be obtained, measured along
              ! the direction of wave propogation.
              !
              x1 (1) = cos (beta(i))
              x1 (2) = sin (beta(i))
              x2 (1) = x(j,k)-x0
              x2 (2) = y(j,k)-y0
              !
              ! Beta is the angle between the (positive sense of the) incident wave
              ! vector and the (positive sense of the) (x,y) point vector.
              !
              dist = x1 (1) * x2 (1) + x1 (2) * x2 (2)
              !
              ! Multiply each Fourier coefficient
              ! by e^-i*nu*l to shift it to the point (x,y).  The coefficients are defined
              ! such that eta(t) = Re B(omega) e^(omega t - k theta).
              !
              if(nu*depth > 15.0)then
                 ! this wave component is in infinite depth (to double-precision).
                 velfact = omega
              else
                 ! this wave component is in finite depth.
                 tanhkhi = cosh(nu*depth)/sinh(nu*depth)
                 velfact = omega * tanhkhi
              end if
              phifact=grav/omega
              !
              coslnu = cos (nu * dist)
              sinlnu = sin (nu * dist)
              reala = eta0(2*i-1)
              imaga = eta0(2*i)
              realarg = (reala * coslnu - imaga * sinlnu)
              imagarg = (coslnu * imaga + sinlnu * reala)
              work_eta(2*i-1) = realarg
              work_eta(2*i) = imagarg
              !
              ! The surface velocity potential
              !
              work_vel(2*i-1) = phifact * imagarg
              work_vel(2*i) =  -phifact * realarg
           enddo
           do i=n_cut+1,nt/2
              work_eta(2*i-1)=zero
              work_eta(2*i) = zero
              work_vel(2*i-1) = zero
              work_vel(2*i) =  zero
           end do
           ! Inverse transform to get the time-series at this point (j,k).
           CALL drealft (work_eta, nt, -1)
           CALL drealft (work_vel, nt, -1)
           eta(j,k,:)=work_eta(:)
           Vs(j,k,:,1)=work_vel(:)
        enddo
     end DO
  ELSE
     !
     ! For the monochromatic case things are much simpler.
     !
     omega=twopi/Tp
     kinf=omega**2/grav; 
     ! Get the wavenumber nu
     CALL fdwvnum(depth,omega,kinf,one,grav,nu)
     amp=Hs/two
     phifact=grav/omega
     !
     do i=1,nt
        t=time0+(i-1)*dt
        eta0(i)=amp*cos(-omega*t)
        do k=1,ny
           do j=1,nx
              x1 (1) = cos (beta0)
              x1 (2) = sin (beta0)
              x2 (1) = x(j,k)-x0
              x2 (2) = y(j,k)-y0
              !
              ! Beta0 is the angle between the (positive sense of the) incident wave
              ! vector and the (positive sense of the) (x,y) point vector.
              !
              dist = x1 (1) * x2 (1) + x1 (2) * x2 (2)
              eta(j,k,i)=amp*cos(nu*dist-omega*t)
              Vs(j,k,i,1)=phifact*amp*sin(nu*dist-omega*t)
           end do
        end do
     end do
     !
     !  Write the wavemaker signal.
     !
     open(21,file='eta0_irregular',status='unknown')
     write(21,36)Hs,Tp,depth
36      format('% Mono-chromatic wave with H=',e12.4,' T=',e12.4,' at x=',e12.4,'.')
     do i=1,nt
        write(21,*)(i-1)*dt,eta0(i)
     end do
     close(21)
  END IF

  DEALLOCATE(eta0, work_eta, work_vel, beta)
444 FORMAT ()
445 FORMAT(3e12.5)
  RETURN
END SUBROUTINE random_wave_signal_3D
