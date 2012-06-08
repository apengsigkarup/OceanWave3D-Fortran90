SUBROUTINE random_wave_signal(i_spec, n1, n2, j0, dx, dt, Tp, Hs, depth, &
     grav, inc_wave_file, kh_max, seed, seed2, eta, Vs, etat0, Vst0, nf)
  !-----------------------------------------------------------------------
  !
  ! Generate a psuedo-random long crested wave over N1 time steps at
  ! N2 points over the 1D computational free-surface which has a total of NF
  ! points.  The grid spacing dx and the time step dt are assumed to be uniform 
  ! in the wave generation zone.  The water depth at the wavemaker (grid point j0) 
  ! is depth.  The wave can be based on a P-M or JONSWAP spectrum, or an 
  ! input time series read from inc_wave_file, depending on i_spec.  The 
  ! wave is moved around in space using linear theory in the frequency 
  ! domain.  eta & Vs are returned with time series
  ! of elevation and potential at z=0 in the wavemaker zone (j=1,n2).
  ! etat0 & Vst0 hold the initial conditions over the entire domain for the P-M or
  ! JONSWAP case but not for an input wave (this is not implemented yet !HBB).  
  ! 
  ! When i_spec=-1, a mono-chromatic wave is generated with period Tp and 
  ! height Hs.  
  !
  !----------------------------------------------------------------------
  IMPLICIT none
  CHARACTER(len=30) inc_wave_file, header
  integer, parameter :: long=selected_real_kind(12,99)
  integer i_spec, seed, seed2, n1, n2, j0, nf
  real(kind=long) :: Tp, Hs, dx, dt, depth
  real(kind=long) :: Vs(n2,n1,*), eta(n2,n1), Vst0(nf,*), etat0(nf)

  ! Local variables
  INTEGER :: i, j, ns_inc, n_cut, ndat
  real(kind=long) :: cosb, sinb, factor, domega, realpt, beta, grav,     &
       imagpt, spec, phase, nu, coslnu, sinlnu, magx2, x2(3), x1(2),     &
       kinf, omega, dist, fn, realarg, imagarg, reala, imaga,            &
       zero=0._long,  one=1._long, two=2._long, pi,twopi,udum=1.0_long,  &
       velfact, tanhkhi, dt_inc, dum, kh_max, phifact, kh, t, x, amp
  !
  ! Local workspace
  !
  real(kind=long),Allocatable :: eta0(:), work_eta(:), work_vel(:)
  allocate( eta0(n1), work_eta(n1), work_vel(n1) )
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
     factor = two/n1
     !
     ! The Fourier coefficients are built based on SI dimensional numbers.  Omega is
     ! radian frequency.
     !
     domega = twopi / (n1 * dt)
     !
     ! Only waves in the x-direction are implemented at this point -hbb
     !
     beta=zero
     cosb = cos (beta )
     sinb = sin (beta )
     !
     ! Generate the Fourier components of the wave elevation at the origin.
     ! Note that the same spectral parameters and SEED values always give the same
     ! incident wave.
     !
     fn = zero
     IF(i_spec==2)THEN
        open(21,file=inc_wave_file,status='old')
        READ(21,'(A)',err=15)header
        READ(21,*)dt_inc
        IF(ABS(dt-dt_inc)>1.e-6)THEN
           print *, 'random_wave_signal.f90: The .inp and .iwf time steps do not agree.'
           print *, header
           print *, dt,dt_inc
           stop
        END IF
        do i=1,n1
           READ(21,*,end=13)eta0(i)
           do j=2,ns_inc
              READ(21,*)dum
           END do
        end do
        go to 14
13      ndat=i-1
        print *, ' Found ',ndat,' data points in the incident wave file.  Padding with zeros up to ',n1
        do i=ndat+1,n1
           eta0(i)=zero
        end do
        go to 14
15      print *, 'random_wave_signal.f90: No header line in the .iwf'
        stop
14      CLOSE(21)
        CALL drealft (eta0, n1, 1)
        !
        ! Scale the FFT
        !
        eta0=eta0*factor
     ELSE
        CALL build_coeff(eta0, n1, Hs, Tp, dt, seed, seed2, i_spec)
     END IF
     ! We truncate the spectrum at k_max.
     n_cut=nint(sqrt(grav*kh_max/depth*tanh(kh_max))/domega)+1
     Print *, 'Truncating the spectrum at T=',twopi/((n_cut-1)*domega), &
          ' which corresponds to kh=',kh_max
     print *, 'The deep water resolution of this wave is ',2*pi*depth/(kh_max*dx)+one, &
          ' points per wavelength.'
     Print *, ' '
     eta0(2*n_cut-1:n1)=zero
     !hbb Print out the elevations at the origin to check.
     CALL drealft (eta0, n1, -1)
     !
     open(21,file='eta0_irregular',status='unknown')
     IF(i_spec==0)THEN
        write(21,32)Hs,Tp,twopi/((n_cut-1)*domega),kh_max
31      format('% P-M spectrum H_s=',e12.4,' T_p=',e12.4,' truncated at T=',e12.4,&
             ' which is kh=',e12.4)
     ELSEIF(i_spec==1)THEN
        write(21,32)Hs,Tp,twopi/((n_cut-1)*domega),kh_max
32      format('% JONSWAP spectrum H_s=',e12.4,' T_p=',e12.4,' truncated at T=',e12.4,&
             ' which is kh=',e12.4)
     ELSEIF(i_spec==2)THEN
        write(21,33)inc_wave_file
33      format('% Incident wave from the file: ',a30)
     END IF
     do i=1,n1
        write(21,*)(i-1)*dt,eta0(i)
     end do
     close(21)
     !
     ! FFT back to the coefficients and print them out
     !
     CALL drealft (eta0, n1, 1)
     eta0=eta0*factor
     open(21,file='eta0_coeffs',status='unknown')
     write(21,34)
34   format('# f [Hz], mag(c), c_r, c_i')
     write(21,35)zero,eta0(1),eta0(1),eta0(2)
35   format(4e16.5)
     do i=3,n1-1,2
        write(21,35) (i-1)/2 * domega/twopi, sqrt(eta0(i)**2+eta0(i+1)**2),eta0(i),eta0(i+1)
     end do
     close(21)
     !
     ! Translate these Fourier amplitudes around in space to generate the
     ! eta and Vs coefficients over some portion of the free-surface.
     !
     DO j = 1, n2
        !
        ! Compute the distance (dist) from the origin of the coordinates (j0) to
        ! the point where the wave elevation is to be obtained, measured along
        ! the direction of wave propogation.
        !
        x1 (1) = cos (beta)
        x1 (2) = sin (beta)
        x2 (1) = (j-j0)*dx
        x2 (2) = zero
        !
        ! Beta is the angle between the (positive sense of the) incident wave
        ! vector and the (positive sense of the) (x,y) point vector.
        !
        dist = x1 (1) * x2 (1) + x1 (2) * x2 (2)
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
           ! The horizontal surface velocity
           !
           !        work_vel(2*i-1) = x1(1) * velfact * realarg
           !        work_vel(2*i) =  x1(1) * velfact * imagarg
           !
           ! The surface velocity potential
           !
           work_vel(2*i-1) = phifact * imagarg
           work_vel(2*i) =  -phifact * realarg
        enddo
        do i=n_cut+1,n1/2
           work_eta(2*i-1)=zero
           work_eta(2*i) = zero
           work_vel(2*i-1) = zero
           work_vel(2*i) =  zero
        end do
        ! Inverse transform to get the time-series at this j-point.
        CALL drealft (work_eta, n1, -1)
        CALL drealft (work_vel, n1, -1)
        eta(j,:)=work_eta(:)
        Vs(j,:,1)=work_vel(:)
     enddo
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
     do i=1,n1
        t=(i-1)*dt
        eta0(i)=amp*cos(-omega*t)
        do j=1,n2
           x=(j-j0)*dx
           eta(j,i)=amp*cos(nu*x-omega*t)
           Vs(j,i,1)=phifact*amp*sin(nu*x-omega*t)
        end do
     end do
     !
     !  Write the wavemaker signal.
     !
     open(21,file='eta0_irregular',status='unknown')
     write(21,36)Hs,Tp
36      format('% Mono-chromatic wave with H=',e12.4,' T=',e12.4,'.')
     do i=1,n1
        write(21,*)(i-1)*dt,eta0(i)
     end do
     close(21)
  END IF

  DEALLOCATE(eta0, work_eta, work_vel)
444 FORMAT ()
445 FORMAT(3e12.5)
  RETURN
END SUBROUTINE random_wave_signal
