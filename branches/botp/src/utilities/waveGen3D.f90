SUBROUTINE waveGen3D(nt, n1, j0, dx_fg, depth, grav, inc_wave_file, kh_max)
     !---------------------------------------------------------------------
     ! Generate three-dimensional wave field from input file.
     ! 
     ! Inputs: nt:= Number of Fourier coefficients in time.
     !         n1:= Number of grid points in x-direction of the relaxation zone
     !         j0:= Index indicating location of wave gauge
     !         dx:= Grid size increment 
     !         dy:= Grid size increment 
     !         depth:= Water depth. Waves are generated at constant depth -
     !                 shoaling not included.
     !         grav := Magnitude of gravitational acceleration
     !         inc_wave_file:= Name of incident wave file
     !         kh_max := Water depth at which the signal is truncated
     !
     ! Output: eta:= Free surface elevation at n1*n2 points. 
     !         Vs:= Velocity potential at the free surface
     !         
     !
     !
     ! Written by: Bo Terp Paulsen, botp@mek.dtu.dk
     !----------------------------------------------------------------------
     USE GlobalVariables, ONLY: RandomWave3D, alpha,beta
     IMPLICIT none
     include "/usr/local/include/fftw3.f" ! FIXME: Hardcoded path to FFTW - should not be necessary
     CHARACTER(len=30) :: inc_wave_file, header
     INTEGER, PARAMETER :: long=selected_real_kind(12,99)
     INTEGER :: nt,n1, j0
     REAL(KIND=long) :: depth, grav, dx_fg, dy, kh_max

     ! Local variables
     INTEGER :: i, j,k, ns_inc, n_cut, ndat, M,Mnorm
     INTEGER :: n2_data,counter,l
     DOUBLE COMPLEX :: arg,I1
     REAL(KIND=long) :: domega, nu, coslnu, sinlnu, kinf, omega, dist, &
          realarg, imagarg, reala, imaga, zero=0._long,  one=1._long, two=2._long, & 
          pi, twopi, udum=1.0_long, velfact, tanhkhi, dt_inc, phifact,  df, &
          Lym, Ly,fN
          
     ! Local workspace
     !
     DOUBLE COMPLEX, ALLOCATABLE :: aFFT(:), work_eta(:), work_vel(:)
     INTEGER(KIND=8) plan, planb, planc
     REAL(KIND=long) :: f(nt/2+1), time(nt),ky1,kym,kx,eta1D(nt)
     REAL(KIND=long),ALLOCATABLE :: eta0(:,:),OUTFFT(:)
     REAL(KIND=long),aLLOCATABLE :: y(:)

     
     ! Some parameters 
     !
     pi=acos(-one)
     twopi=two*pi
     I1 = COMPLEX(0,1)
     ndat = nt

     ! Read surface elevation from file "inc_wave_file"
     !
     OPEN(21,FILE=inc_wave_file,STATUS='old')
     READ(21,'(A)',err=15)header
     READ(21,*)dt_inc, n2_data
15   PRINT *, 'random_wave_signal.f90: No header line in the .iwf'

     ! Once the header is read, we need to allocate internal fields
     !
     ALLOCATE( eta0(nt,n2_data),aFFT(nt/2+1),y(n2_data), &
               work_eta(nt/2+1),work_vel(nt/2+1), &
               RandomWave3D%eta(n1,n2_data,nt),OUTFFT(nt), &
               RandomWave3D%Phis(n1,n2_data,nt), &
               RandomWave3D%y(n2_data), RandomWave3D%x(n1))
     eta0 = zero

     ! Now we read the rest of the file
     !
     READ(21,*)(y(j),j=1,n2_data)
     DO i=1,nt
         READ(21,*,end=13) (eta0(i,j),j=1,n2_data)
     ENDDO
     GO TO 14
13   ndat=i-1
14   CLOSE(21)
   
    IF(ndat < nt) THEN 
        PRINT *, 'Time series for incident 3D-waves is to short, ndat = ', &
        ndat , '<nt = ' , nt ,'. Reduce number of timesteps or zeropad signal.'  
    stop
    ENDIF

    ! We save some information in RandomWave3D for global usage
    !
    RandomWave3D%n2_data = n2_data
    RandomWave3D%y = y
    RandomWave3D%dy = RandomWave3D%y(2)-RandomWave3D%y(1)

    ! Prepare for 1D-Fourier transformation. We are using FFTW
    !
    Mnorm = nt!n2_data       ! Number of measurements in timeseries [#]
    fN = 1./(2.*dt_inc)  ! Nyquist frequency [1/s]
    df = 2*fN/nt         ! Frequency increments
    domega = 2*pi*df      
    time = (/((i-1)*dt_inc,i=1,SIZE(time))/)  ! Time vector for plotting purpose only
    f = (/((i-1)*df,i=1,SIZE(f))/)            ! Frequency vector for plotting purpose only


    ! Prepare for truncation
    ! From the dispersion relation we get
    !
    n_cut=nint(sqrt(grav*kh_max/depth*tanh(kh_max))/domega)+1
    PRINT *, 'Truncating the spectrum at T=',twopi/((n_cut-1)*domega), &
         ' which corresponds to kh=',kh_max
    PRINT *, 'The deep water resolution of this wave is ',2*pi*depth/(kh_max*dx_fg)+one, &
         ' points per wavelength.'
    PRINT *, ' '
   

   DO m = 1,n2_data ! We loop over the y-direction
    ! Make a real to complex Fourier transformation
    !
    call dfftw_plan_dft_r2c_1d(plan,nt,eta0(:,m),aFFT,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan, eta0(:,m), aFFT)
    call dfftw_destroy_plan(plan)

    ! Scale output
    !
    aFFT = aFFT/Mnorm

    ! Now we truncate the spectrum at k_max.
    !
    aFFT(n_cut:nt/2+1)=zero 

    ! Translate these Fourier amplitudes around in space to generate the
    ! eta and Vs coefficients over some portion of the free-surface.
    !
    ! To reduce memory requirements we create a coarse local grid for the wave
    ! gauge measurements. During execution we interpolate onto FineGrid.
    DO j = 1, n1  ! Loop over relaxation zone (x-direction)

    RandomWave3D%x(j) = (j0-2)*dx_fg*(j-1)/(n1-1) ! Uniform grid of n1 points from 0 to j0

        ! Compute the distance (dist) from the origin of the local coordinate system (j0) to
        ! the point where the wave elevation is to be obtained, measured in the 
        ! x-direction
        !
        !dist = (j-j0)*dx 
        dist = -(j0-2)*dx_fg + RandomWave3D%x(j) 

        ! Take care of omega=0 and the Nyquist frequency.
        !
        work_eta(1) = zero
        work_vel(1) = zero
        
        ! The rest of frequency space.
        !
        DO i = 2, n_cut 
           omega = (i - 1) * domega

           ! Compute the wave number by solving the linear dispersion 
           ! relation 
           kinf = omega * omega / grav
           IF (depth.gt.zero) then
              CALL fdwvnum (depth, omega, kinf, udum, grav, nu)
           ELSE
              nu = kinf
           ENDIF

           phifact=grav/omega

           work_eta(i) = aFFT(i)*exp(-I1*nu*dist)
           work_vel(i) = I1*phifact*aFFT(i)*exp(-I1*nu*dist)


        ENDDO 

        DO i=n_cut+1,nt/2+1
           work_eta(i) =zero
           work_eta(i) = zero
           work_vel(i) = zero
           work_vel(i) = zero
        ENDDO

        ! We reuse OUTFFT
        OUTFFT = zero

        ! Make a complex to real Fourier transformation of the free surface elevation
        ! 
        call dfftw_plan_dft_c2r_1d(planb,nt,work_eta,OUTFFT,FFTW_ESTIMATE)
        call dfftw_execute_dft_c2r(planb,work_eta,OUTFFT)
        call dfftw_destroy_plan(planb)

        ! Saving translated signal in global struct
        !
        DO i=1,nt
                RandomWave3D%eta(j,m,i) = OUTFFT(i)
        ENDDO

        ! We reuse OUTFFT
        OUTFFT = zero

        ! Make a complex to real Fourier transformation of the velocity potential
        !
        call dfftw_plan_dft_c2r_1d(planc,nt,work_vel,OUTFFT,FFTW_ESTIMATE)
        call dfftw_execute_dft_c2r(planc,work_vel, OUTFFT)
        call dfftw_destroy_plan(planc)

        ! Saving translated signal in global struct
        !
        DO i=1,nt
                RandomWave3D%Phis(j,m,i) = OUTFFT(i)
        ENDDO


     ENDDO
 ENDDO
 RandomWave3D%dx = RandomWave3D%x(2)-RandomWave3D%x(1)

        OPEN(21,FILE='etaTRANSFORMED',STATUS='unknown')
         DO i=1,nt
            write(21,*)time(i),RandomWave3D%eta(n1,:,i)
        ENDDO
         close(21)
  DEALLOCATE(eta0, work_eta, work_vel)
  RETURN
END SUBROUTINE waveGen3D