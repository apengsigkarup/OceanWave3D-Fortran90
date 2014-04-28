SUBROUTINE random_wave_signal(i_spec, nt, nx, x0, x, dt, Tp, Hs, depth, &
     grav, inc_wave_file, kh_max, seed, seed2, eta, phiS, eta0, n_cut, time0)
  !-----------------------------------------------------------------------
  !
  ! Generate a pseudo-random long crested wave with direction beta to the 
  ! x-axis, over NT time steps at NX points over the 1D computational
  ! free-surface which has a total of NF points.  
  ! The x-positions of the points in the wave generation zone are given by x(1:nx) 
  ! and the time-step dt is uniform.  The water depth at the wavemaker (x0) 
  ! is depth.  The wave can be based on a P-M (i_spec=0) or JONSWAP (ispec=1) 
  ! spectrum, or an input time series read from inc_wave_file (i_spec=2).  
  ! The wave is moved around in space using linear theory in the frequency 
  ! domain.  eta & phiS are returned with the time series
  ! of elevation and potential at z=0 in the wavemaker zone (j=1,nx).
  ! etat0 & phiSt0 hold the initial conditions over the entire domain for the P-M or
  ! JONSWAP case but not for an input wave (this is not implemented yet !HBB).  
  ! 
  ! When i_spec=-1, a mono-chromatic wave is generated with period Tp and 
  ! height Hs.  
  !
  !----------------------------------------------------------------------
  IMPLICIT none
  CHARACTER(len=30) inc_wave_file, header
  integer, parameter :: long=selected_real_kind(12,99)
  integer i_spec, seed, seed2, nt, nx, n_cut
  real(kind=long) :: Tp, Hs, dt, depth, time0, x0, grav
  real(kind=long) :: x(nx), phiS(nx,nt), eta(nx,nt), eta0(nt)

  ! Local variables
  INTEGER :: i, j
  real(kind=long) :: cosb, sinb, factor, domega, realpt, beta, dx,       &
       imagpt, spec, phase, nu, coslnu, sinlnu, magx2, x2(3), x1(2),     &
       kinf, omega, dist, realarg, imagarg, reala, imaga,                &
       zero=0._long,  one=1._long, two=2._long, pi,twopi,udum=1.0_long,  &
       velfact, tanhkhi, dt_inc, dum, kh_max, phifact, kh, t, amp
  !
  ! Local workspace
  !
  real(kind=long),Allocatable :: work_eta(:), work_vel(:)
  allocate( work_eta(nt), work_vel(nt) )
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
     ! Only waves in the x-direction are implemented at this point -hbb
     !
     beta=zero
     cosb = cos (beta )
     sinb = sin (beta )
     !
     !
     ! Translate the Fourier amplitudes around in space to generate the
     ! eta and phiS coefficients over some portion of the free-surface.
     !
     DO j = 1, nx
        !
        ! Compute the distance (dist) from the origin of the coordinates (x0) to
        ! the point where the wave elevation is to be obtained, measured along
        ! the direction of wave propogation.
        !
        x1 (1) = cos (beta)
        x1 (2) = sin (beta)
        x2 (1) = x(j)-x0
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
        do i=n_cut+1,nt/2
           work_eta(2*i-1)=zero
           work_eta(2*i) = zero
           work_vel(2*i-1) = zero
           work_vel(2*i) =  zero
        end do
        ! Inverse transform to get the time-series at this j-point.
        CALL drealft (work_eta, nt, -1)
        CALL drealft (work_vel, nt, -1)
        eta(j,:)=work_eta(:)
        phiS(j,:)=work_vel(:)
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
     do i=1,nt
        t=time0+(i-1)*dt
        eta0(i)=amp*cos(-omega*t)
        do j=1,nx
           dist=x(j)-x0
           eta(j,i)=amp*cos(nu*dist-omega*t)
           phiS(j,i)=phifact*amp*sin(nu*dist-omega*t)
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

  DEALLOCATE(work_eta, work_vel)
444 FORMAT ()
445 FORMAT(3e12.5)
  RETURN
END SUBROUTINE random_wave_signal
