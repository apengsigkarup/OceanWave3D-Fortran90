SUBROUTINE random_wave_signal_3D(i_spec, eta0, n_cut, beta, nt, nx, ny, beta0, x0, y0, x, y,   & 
     dt, Tp, Hs, depth, grav, inc_wave_file, kh_max, seed, seed2, eta, phiS, time0 )
  !-----------------------------------------------------------------------
  !
  ! Generate a pseudo-random long crested wave with direction BETA0 to the 
  ! x-axis, over NT time steps at NX by NY points over the computational
  ! free-surface.  The x-positions of the points in the wave generation 
  ! zone are given by x(1:nx,1:ny) and y(1:nx,1:ny), 
  ! and the time-step DT is uniform.  The water depth at the wavemaker (x0,y0) 
  ! is DEPTH.  The wave can be based on a long-crested P-M (i_spec=30) or 
  ! JONSWAP (ispec=31) spectrum at angle BETA0 to the x-axis, or it can be 
  ! a JONSWAP with NORMAL directional spreading about the angle BETA0 (i_spec=33).  
  ! The wave is moved around in space using linear theory in the frequency 
  ! domain.  eta & phiS are returned with the time series
  ! of elevation and the potential at z=0 in the wavemaker zone.  
  ! 
  ! When i_spec=-30, a mono-chromatic wave is generated with period Tp, 
  ! height Hs and angle BETA0.  
  !
  !
  ! By Harry B. Bingham
  !
  !----------------------------------------------------------------------
  IMPLICIT none
  CHARACTER(len=30) inc_wave_file, header
  integer, parameter :: long=selected_real_kind(12,99)
  integer i_spec, seed, seed2, nt, nx, ny, n_cut
  real(kind=long) :: Tp, Hs, dt, depth, time0, x0, y0, grav
  real(kind=long) :: x(nx,ny), y(nx,ny), phiS(nx,ny,nt), eta(nx,ny,nt), &
       eta0(nt), beta(nt)

  ! Local variables
  INTEGER :: i, j, k, ns_inc, ndat
  real(kind=long) :: factor, domega, realpt, beta0, dx,                  &
       imagpt, spec, phase, nu, coslnu, sinlnu, magx2, x2(3), x1(2),     &
       kinf, omega, dist, realarg, imagarg, reala, imaga,                &
       zero=0._long,  one=1._long, two=2._long, pi,twopi,udum=1.0_long,  &
       velfact, tanhkhi, dt_inc, dum, kh_max, phifact, kh, t, amp, deg2rad
  !
  ! Local workspace
  !
  real(kind=long),Allocatable :: work_eta(:), work_phi(:)
  allocate( work_eta(nt), work_phi(nt))
  !
  !
  ! Some parameters 
  !
  pi=acos(-one)
  twopi=two*pi
  deg2rad=pi/180._long
  !
  IF(i_spec >= 0)THEN
     ! 
     ! Generate the random wave.
     !
     ! The FT factor. 
     !
     factor = two/nt
     !
     ! The Fourier coefficients have been built based on SI dimensional numbers.  Omega is
     ! radian frequency.
     !
     domega = twopi / (nt * dt)
     !
     ! Translate these Fourier amplitudes around in space to generate the
     ! eta and phiS coefficients over this portion of the free-surface.
     !
     DO k=1,ny
        DO j = 1, nx
           !
           ! Take care of omega=0 and  the Nyquist frequency.
           !
           work_eta(1) = zero
           work_eta(2) = zero
           work_phi(1) = zero
           work_phi(2) = zero
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
              x1 (1) = cos (deg2rad*beta(i))
              x1 (2) = sin (deg2rad*beta(i))
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
              work_phi(2*i-1) = phifact * imagarg
              work_phi(2*i) =  -phifact * realarg
           enddo
           do i=n_cut+1,nt/2
              work_eta(2*i-1)=zero
              work_eta(2*i) = zero
              work_phi(2*i-1) = zero
              work_phi(2*i) =  zero
           end do
           ! Inverse transform to get the time-series at this point (j,k).
           CALL drealft (work_eta, nt, -1)
           CALL drealft (work_phi, nt, -1)
           eta(j,k,:)=work_eta(:)
           phiS(j,k,:)=work_phi(:)
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
              x1 (1) = cos (deg2rad*beta0)
              x1 (2) = sin (deg2rad*beta0)
              x2 (1) = x(j,k)-x0
              x2 (2) = y(j,k)-y0
              !
              ! Beta0 is the angle between the (positive sense of the) incident wave
              ! vector and the (positive sense of the) (x,y) point vector.
              !
              dist = x1 (1) * x2 (1) + x1 (2) * x2 (2)
              eta(j,k,i)=amp*cos(nu*dist-omega*t)
              phiS(j,k,i)=phifact*amp*sin(nu*dist-omega*t)
           end do
        end do
     end do
     !
     !  Write the wavemaker signal.
     !
     open(21,file='eta0_irregular',status='unknown')
     write(21,36)Hs,Tp,depth,beta0
36      format('% Mono-chromatic wave with H=',e12.4,' T=',e12.4,' at h=',e12.4,  &
             ' at heading angle ',f10.2,'.')
     do i=1,nt
        write(21,*)(i-1)*dt,eta0(i)
     end do
     close(21)
  END IF

  DEALLOCATE(work_eta, work_phi)
444 FORMAT ()
445 FORMAT(3e12.5)
  RETURN
END SUBROUTINE random_wave_signal_3D
