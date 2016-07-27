SUBROUTINE OceanWave3DT0Setup
  !
  ! This subroutine performs all of the initial set-up work before the time-stepping 
  ! begins.  
  !
  ! By Allan P. Engsig-Karup.
  !
  USE GlobalVariables
  USE MGLevels
  IMPLICIT NONE
  EXTERNAL BuildLinearSystem, BuildLinearSystemTransformedCurvilinear
  ! GD: to test the cross derivatives...
  REAL(KIND=long), DIMENSION(:,:,:), ALLOCATABLE :: tmpPHI
  TYPE (Diff_def)        :: FullRankStencils
  INTEGER i, j, k, n_cut, NxT, NyT

  ! 
  ! Set up the input/output file structure, open the log file and output the 
  ! header to the screen and the log file
  !
  CALL Initialize
  WRITE (6,2010)
  WRITE (fileop(1),2010)
  ! Read in the run parameters 
  CALL ReadInputFileParameters

  print *, 'format type 1:',formattype
  !CALL SetupCompDomain
  !GD: test for the curvilinear domain (definition of initial grid) (to uncomment and to FIX later)
  !
  ! Build the rectangular (x,y) domain and the sigma levels, allocate space for the bathymetry.
  !
  IF (curvilinearONOFF==1) THEN
     CALL SetupCompDomain_test
  ELSE
     CALL SetupCompDomain
  ENDIF

  IF (filteringONOFF>0) THEN
     !
     ! Initialization of filtering
     !
     CALL FilterInit(filtercoefficients,filtercoefficients2)
  ENDIF

  !
  ! Allocate space for the solution variables and wavefield.
  !
  print *,'Initializing variables and arrays.'
  write(fileop(1),*)'Initializing variables and arrays.'
  CALL InitializeVariables
  !
  ! We start at time=0 here but if this is a hot start, time0 will be read in SetupInitialConditions.
  !
  time=zero
  !
  ! Set up for relaxations zones and wave generation.
  !
  IF (relaxONOFF>0) THEN
     CALL PreprocessRelaxationZones
     PRINT *,'  Relaxation zones have been setup.'
     WRITE(FILEOP(1),*)'  Relaxation zones have been setup.'
     IF(IncWaveType==2)THEN
        ! SFsol%T and SFsol%L are used to determine dt based on the Cr number, so they are
        ! re-set here for other wave generation types.  -HBB
        SFsol%T=RandomWave(1)%Tp; SFsol%L=g*RandomWave(1)%Tp**2/(two*pi)
     END IF
     ! Set time step size...
     IF (CFL/=zero) THEN
        IF (FineGrid%Nx>1) THEN
           dxmin = dx
        ELSE
           dxmin = two*dy
        ENDIF
        IF (FineGrid%Ny>1) THEN
           dymin = dy
        ELSE
           dymin = two*dx
        ENDIF
        dsigmamin = FineGrid%z(FineGrid%Nz)-FineGrid%z(FineGrid%Nz-1)
        c = SFsol%L / SFsol%T
        dt = CFL*MIN(dxmin,dymin)/c
        PRINT *,''
        PRINT *,'  Time step size modified based on the incident wave. dt = ',dt
        PRINT *,'  Courant number based on min. (dx,dy) and c=',c,', Cr = ',c*dt/MIN(dxmin,dymin)
        PRINT *,'  Discrete anisotropy, ds_min/dx_min = ',dsigmamin/dxmin
        WRITE(FILEOP(1),*)''
        WRITE(FILEOP(1),*)'  Time step size modified based on the incident wave. dt = ',dt
        WRITE(FILEOP(1),*)'  Courant number based on min. (dx,dy) and c=',c,', Cr = ',c*dt/MIN(dxmin,dymin)
        WRITE(FILEOP(1),*)'  Discrete anisotropy, ds_min/dx_min = ',dsigmamin/dxmin
     END IF
     IF (IncWaveType==1 ) THEN
        CALL stream_func_set_up(g,SFsol%h,SFsol%T,SFsol%i_wavel_or_per,SFsol%L,  &
             SFsol%k,SFsol%HH,SFsol%i_euler_or_stokes,SFsol%e_or_s_vel,SFsol%i_deep_or_finite,  &
             SFsol%n_h_steps,SFsol%n_four_modes,SFsol%nwrk,SFsol%yy,SFsol%zz)
        PRINT *,'     SF solution: k=',SFsol%k,',h=',SFsol%h,',H=',SFsol%HH,',T=',SFsol%T,',L=',SFsol%L
        PRINT *,'                  kh=',SFsol%k*SFsol%h,',c=',SFsol%L/SFsol%T
        WRITE(FILEOP(1),*)'     SF solution: k=',SFsol%k,',h=',SFsol%h,',H=',SFsol%HH,',T=',SFsol%T,',L=',SFsol%L
        WRITE(FILEOP(1),*)'                  kh=',SFsol%k*SFsol%h,',c=',SFsol%L/SFsol%T

        WRITE(6,62)two*pi/SFsol%k,SFsol%T,SFsol%zz(2)/SFsol%k,                  &
             SFsol%zz(5)*sqrt(g/SFsol%k), SFsol%zz(6)*sqrt(g/SFsol%k)
        WRITE(fileop(1),62)two*pi/SFsol%k,SFsol%T,SFsol%zz(2)/SFsol%k,                  &
             SFsol%zz(5)*sqrt(g/SFsol%k), SFsol%zz(6)*sqrt(g/SFsol%k)
62      FORMAT(' The incident wave is a stream function wave with L=',  &
             e12.6,/,'  T=',e12.6,' H=',e12.6,' u_E=',e12.6,' and u_S=',     &
             e12.6,//)
     ELSEIF(IncWaveType==2)THEN
        !
        ! A linear regular or irregular wave
        !
        print *, ' '
        write(fileop(1),*) ' '
        IF(RandomWave(1)%ispec==-1)THEN
           WRITE(6,70)RandomWave(1)%Tp,RandomWave(1)%Hs
           WRITE(fileop(1),70)RandomWave(1)%Tp,RandomWave(1)%Hs
70         FORMAT(' The incident wave is a linear mono-chromatic wave with',/,&
                ' T=', e10.4,' and H=',e10.4,'.',//)
        ELSEIF(RandomWave(1)%ispec==0)THEN
           WRITE(6,71)RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%seed,RandomWave(1)%seed2
           WRITE(fileop(1),71)RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%seed,RandomWave(1)%seed2
71         FORMAT(' The incident wave is a 2D P-M spectrum with',/,&
                ' T_p=', e10.4,' and H_s=',e10.4,', seed values are:',/,2i10,//)
        ELSEIF(RandomWave(1)%ispec==1 .or. RandomWave(1)%ispec==3)THEN
           WRITE(6,72)RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%gamma,RandomWave(1)%seed,RandomWave(1)%seed2
           WRITE(fileop(1),72)RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%gamma,RandomWave(1)%seed,RandomWave(1)%seed2
72         FORMAT(' The incident wave is a 2D JONSWAP spectrum with ',/, &
                'T_p=',  e10.4,', H_s=',e10.4,' and gamma=',e10.4, ' and seed values are:',/,2i10,//)
        ELSEIF(RandomWave(1)%ispec==2)THEN
           WRITE(6,73)RandomWave(1)%inc_wave_file
           WRITE(fileop(1),73)RandomWave(1)%inc_wave_file
73         FORMAT(' The incident wave will be read from file ',a30,/)
        ELSEIF(RandomWave(1)%ispec==-30)THEN
           WRITE(6,74)RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%beta0
           WRITE(fileop(1),74)RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%beta0
74         FORMAT(' The incident wave is a linear, long-crested mono-chromatic wave with',/,&
                ' T=', e10.4,' and H=',e10.4,' at angle ',f10.2,' deg. to the x-axis.',//)
        ELSEIF(RandomWave(1)%ispec==30)THEN
           WRITE(6,75)RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%seed,RandomWave(1)%seed2,&
                RandomWave(1)%beta0
           WRITE(fileop(1),75)RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%seed,RandomWave(1)%seed2,&
                RandomWave(1)%beta0
75         FORMAT(' The incident wave is a 2D P-M spectrum with',/,&
                ' T_p=', e10.4,' and H_s=',e10.4,', seed values are:',/,2i10,/, &
                ' at angle ',f10.2,' deg. to the x-axis.',// )
        ELSEIF(RandomWave(1)%ispec==31)THEN
           WRITE(6,76)RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%gamma,RandomWave(1)%seed, & 
                RandomWave(1)%seed2,RandomWave(1)%beta0
           WRITE(fileop(1),76)RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%gamma,            &
                RandomWave(1)%seed,RandomWave(1)%seed2,RandomWave(1)%beta0
76         FORMAT(' The incident wave is a 2D JONSWAP spectrum with ',/, &
                'T_p=',  e10.4,', H_s=',e10.4,' and gamma=',e10.4,', seed values are:',/,2i10,/,  &
                ' at angle ',f10.2,' deg. to the x-axis.',// )
        ELSEIF(RandomWave(1)%ispec==33)THEN
           WRITE(6,77)RandomWave(1)%beta0,RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%gamma,  &
                RandomWave(1)%seed,RandomWave(1)%seed2
           WRITE(fileop(1),77)RandomWave(1)%beta0,RandomWave(1)%Tp,RandomWave(1)%Hs,              &
                RandomWave(1)%gamma,RandomWave(1)%seed,RandomWave(1)%seed2
77         FORMAT(' The incident wave is a 3D JONSWAP spectrum with Normal spreading at heading angle ',e10.4,/, &
                ' deg. to the x-axis. T_p=',  e10.4,' H_s=',e10.4,' and gamma=',e10.4,  &
                ', seed values are:',/,2i10,//)

         ! 3D random waves with a cos^2s0 spreading

        ELSEIF(RandomWave(1)%ispec==34)THEN
           WRITE(6,78)RandomWave(1)%beta0, RandomWave(1)%S0, RandomWave(1)%Tp, RandomWave(1)%Hs, &
                RandomWave(1)%seed,RandomWave(1)%seed2,RandomWave(1)%gamma
           WRITE(fileop(1),78)RandomWave(1)%beta0, RandomWave(1)%S0, RandomWave(1)%Tp, RandomWave(1)%Hs, &
                RandomWave(1)%seed,RandomWave(1)%seed2,RandomWave(1)%gamma
78         FORMAT(' The incident wave is a 3D JONSWAP spectrum with Cos^(2S) spreading at heading angle ',e10.4,/, &
                ' deg. to the x-axis. S=',e10.4,', T_p=',  e10.4,' H_s=',e10.4,', seed values are:',/,2i10,', gamma=',e10.4//)


        ELSE
           PRINT *, 'ERROR:  RandomWave%ispec must be -1,0,1,2 for 2D waves; or -30,30,31,33,34 for 3D waves.'
           STOP
        END IF
        !
        IF(RandomWave(1)%ispec >=0) THEN
           !
           ! Build the wave using the nearest power of two which is greater than nsteps since we
           ! are using FFT's based on powers of two.
           !
           i=NINT(LOG(real(Nsteps,long))/LOG(two))
           IF(two**i<Nsteps)i=i+1
           n_fft=2**i;
        ELSE
           !
           ! For a monochromatic wave we use Nsteps.
           !
           n_fft=Nsteps
        END IF
        !
        ! Linear wave generation is only implemented for X-directed wave generation,
        ! but angle to the x-axis can be chosen arbitrarily.  -HBB
        !
        ! First compute the Fourier coefficients of the wave.  Find the first 
        ! relaxation zone for generation and base the coefficients on those 
        ! parameters (though they all should be the same.) 
        !
        Do i=1,relaxNo
           If (RelaxZones(i)%XorYgen=='X' .or. RelaxZones(i)%WavegenONOFF==1) exit
        END Do
        If (i>relaxNo) then
           print *, 'Inconsistent relaxation/generation parameters, no x-generation zone '
           print *, 'found even though 3D waves have been asked for.'  
           stop
        end If
        !
        ! Regardless of which parameters we're using, store the coefficients with relaxation 
        ! zone 1.  
        !
        Allocate( RandomWave(1)%eta0(n_fft), RandomWave(1)%beta(n_fft) )
        ! 
        Call random_wave_coefficients( RandomWave(1)%ispec, n_fft, RandomWave(1)%beta0,     &
             dt, dx, RandomWave(1)%Tp, RandomWave(1)%Hs, RandomWave(1)%h0, g,               &
             RandomWave(1)%inc_wave_file, RandomWave(1)%kh_max, RandomWave(1)%seed,         &
             RandomWave(1)%seed2, RandomWave(1)%eta0, RandomWave(1)%beta, RandomWave(1)%S0, &
             n_cut,RandomWave(1)%gamma )
        !
        ! Count the total number of grid points in each generation zone and allocate space
        ! for eta and phiS and compute the elevation and surface potential time-histories.  
        !  
        If (abs(RandomWave(1)%ispec)<30) Then
           ! 2D waves along the x-axis.  
           DO i=1,relaxNo
              If(RelaxZones(i)%XorYgen=='X' .AND. RelaxZones(i)%WavegenONOFF==1) THEN
                 n_wavem=RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1
                 RandomWave(i)%dx=FineGrid%x(RelaxZones(i)%idx(1)+1,1)  &
                      -FineGrid%x(RelaxZones(i)%idx(1),1)
                 !
                 print *, 'Zone ',i,':'
                 print *, 'The generated wave is centered at (x,y)=(',RandomWave(i)%x0,',',&
                      RandomWave(i)%y0,      &
                      ') in a depth of',RandomWave(i)%h0,' This generation zone contains ',&
                      n_wavem, ' grid points.'
                 write(fileop(1),*) 'Zone ',i,':'
                 write(fileop(1),*) 'The generated wave is centered at (x,y)=(',RandomWave(i)%x0,& 
                      ',',RandomWave(i)%y0,      &
                      ') in a depth of',RandomWave(i)%h0,' This generation zone contains ',&
                      n_wavem, ' grid points.'
                 RandomWave(i)%nf=FineGrid%nx+2*GhostGridX
                 !
                 ALLOCATE(RandomWave(i)%eta(n_wavem,n_fft), RandomWave(i)%Phis(n_wavem,n_fft) )

                 CALL random_wave_signal(RandomWave(1)%ispec, n_fft, n_wavem, RandomWave(1)%x0, &
                      FineGrid%x(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),1),                  &
                      dt, RandomWave(1)%Tp, RandomWave(1)%Hs, RandomWave(1)%h0,                 &
                      g, RandomWave(1)%inc_wave_file, RandomWave(1)%kh_max, RandomWave(1)%seed, &
                      RandomWave(1)%seed2, RandomWave(i)%eta, RandomWave(i)%Phis,               &
                      RandomWave(1)%eta0, n_cut, time0 )
              END If
           END DO
        ELSE
           !
           ! 3D waves at angle beta0 to the x-axis.  
           !
           !
           ! Now use the coefficients to get the incident wave time series at all points in 
           ! all wave making relaxation zones.  
           !
           DO i=1,relaxNo
              If(RelaxZones(i)%XorYgen=='X' .AND. RelaxZones(i)%WavegenONOFF==1) THEN
                 RandomWave(i)%nx=RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1
                 RandomWave(i)%ny=RelaxZones(i)%idx(4)-RelaxZones(i)%idx(3)+1
                 n_wavem=RandomWave(i)%nx*RandomWave(i)%ny
                 !
                 print *, 'Zone ',i,':'
                 print *, 'The generated wave is centered at (x,y)=(',RandomWave(i)%x0,',',   &
                      RandomWave(i)%y0,      &
                      ') in a depth of',RandomWave(i)%h0,' This generation zone contains ',   &
                      RandomWave(i)%nx, ' by ',RandomWave(i)%ny,' grid points.'
                 write(fileop(1),*) 'Zone ',i,':'
                 write(fileop(1),*) 'The generated wave is centered at (x,y)=(',RandomWave(i)%x0,&
                      ',', RandomWave(i)%y0,      &
                      ') in a depth of',RandomWave(i)%h0,' This generation zone contains ',   &
                      RandomWave(i)%nx, ' by ',RandomWave(i)%ny,' grid points.'
                 !
                 ALLOCATE(RandomWave(i)%eta(n_wavem,n_fft), RandomWave(i)%Phis(n_wavem,n_fft) )

                 CALL random_wave_signal_3D(RandomWave(1)%ispec, RandomWave(1)%eta0, n_cut,     &
                      RandomWave(1)%beta, n_fft, RandomWave(i)%nx, RandomWave(i)%ny,            &
                      RandomWave(1)%beta0, RandomWave(1)%x0, RandomWave(1)%y0,                  &
                      FineGrid%x(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),                     &
                      RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)),       &
                      FineGrid%y(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),                     &
                      RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)),       &
                      dt, RandomWave(1)%Tp, RandomWave(1)%Hs, RandomWave(1)%h0,                 &
                      g, RandomWave(1)%inc_wave_file, RandomWave(1)%kh_max, RandomWave(1)%seed, &
                      RandomWave(1)%seed2, RandomWave(i)%eta, RandomWave(i)%Phis, time0 )
                 !hbb
                 !                 write(201,*)RandomWave(i)
              END If
           END DO
        END If
     ENDIF
  ENDIF
  IF(IncWaveType==3) THEN
     ! Wave generation by flux condition on western wall, botp
     WRITE(6,79)
     WRITE(fileop(1),79)
79   FORMAT(/,' The incident wave is a flux boundary condition applied at the Western boundary.',/) 
     CALL setupWavePaddle()
  ENDIF
  ! Uneumann is in all cases added to the western boundary RHS. Only if
  ! IncWaveType==3 is it non-zero. 
  ! FIXME: Is there a better solution where this field is only loaded if needed?
  ! botp
  ALLOCATE(Uneumann(FineGrid%Nz+GhostGridZ,FineGrid%Ny+2*GhostGridY))
  Uneumann = zero
  !
  ! Set up the initial conditions and the bathymetry data
  !
  print *,'setup ICs...'
  write(fileop(1),*)'setup ICs...'
  CALL SetupInitialConditions
  time=time0
  dt0 = dt ! botp,Used in AnalyticWaveMaker2D.f90 since OpenFoam changes the timestep
  print *,'done with ICs'
  write(fileop(1),*)'done with ICs'
  !
  ! Set up the Pressure Damping Zones if any.
  !
  If (NDampZones /=0) THEN
     Call PreprocessPDampingZones
     print *, ' '
     print *, 'Pressure damping zones are set up'
     print *, ' '
     write(fileop(1),*) ' '
     write(fileop(1),*) 'Pressure damping zones are set up'
     write(fileop(1),*) ' '
  END If
  IF (.FALSE.) THEN
     !
     ! Test code to validate buildlinearsystem subroutines.
     !
     CALL TestCodeDucrozet1
  END IF
  !
  IF (DetermineBottomGradients==1) THEN
     PRINT *,'Warning: determining bottom gradients numerically, this is only implemented for rectangular grids. (HBB)'
     WRITE(FILEOP(1),*)'Warning: determining bottom gradients numerically, this is only implemented for rectangular grids. (HBB)'
     ! Determine the bottom gradients numerically
     CALL PreProcessDiffStencils(FineGrid,FineGrid%DiffStencils,GhostGridX,GhostGridY,GhostGridZ,alpha,beta,gamma,fileop(1))
     NxT=FineGrid%Nx+2*GhostGridX; NyT=FineGrid%NY+2*GhostGridY
     IF (FineGrid%Nx==1) THEN
        FineGrid%hx = zero; FineGrid%hxx = zero;
     ELSE
        CALL DiffXEven(FineGrid%h,FineGrid%hx, 1,NxT,NyT,1,FineGrid%DiffStencils,alpha)
        CALL DiffXEven(FineGrid%h,FineGrid%hxx,2,NxT,NyT,1,FineGrid%DiffStencils,alpha)
     END IF
     IF (FineGrid%Ny==1) THEN
        FineGrid%hy  = zero; FineGrid%hyy = zero;
     ELSE
        CALL DiffYEven(FineGrid%h,FineGrid%hy, 1,NxT,NyT,1,FineGrid%DiffStencils,beta)
        CALL DiffYEven(FineGrid%h,FineGrid%hyy,2,NxT,NyT,1,FineGrid%DiffStencils,beta)
     ENDIF
     DEALLOCATE(FineGrid%DiffStencils%StencilX,FineGrid%DiffStencils%StencilY,FineGrid%DiffStencils%StencilZ) 
  ENDIF
!
! Now that we are sure that we have all bottom gradients, save the bathymetry data file.
!
  Open(FILEOP(16),file='bathymetry.chk',status='unknown')
  Write(FILEOP(16),81)
81 FORMAT('% Bottom bathymetry: ((h(i,j),h_x,h_xx,h_y,h_yy),j=1,Ny),i=1,Nx)')
  Do i=1,FineGrid%Nx+2*GhostGridX
     Do j=1,FineGrid%Ny+2*GhostGridY
        write(Fileop(16),82)FineGrid%h(i,j),FineGrid%hx(i,j),FineGrid%hxx(i,j), &
             FineGrid%hy(i,j),FineGrid%hyy(i,j)
     end Do
  end Do
  close(FILEOP(16))
82 FORMAT(5e16.6)

  IF (Precond==1) THEN ! PREPARE FOR PRECONDITIONING
     ! DETERMINE LOW-ORDER FINITE DIFFERENCE STENCILS
     ! FIXME: make it possible to choose the order of the preconditioner in the input file
     WRITE(*,FMT='(A,I2,A)') '   Preconditioner: DIRECT LU (',2*alphaprecond,' order, linear)'
     WRITE(fileop(1),FMT='(A,I2,A)') '   Preconditioner: DIRECT LU (',2*alphaprecond,' order, linear)'

     CALL PreparePreconditioner(FineGrid%PreconditioningMatrix,FineGrid,GhostGridX, GhostGridY, GhostGridZ, &
          alphaprecond, betaprecond, gammaprecond, Precond, CurvilinearONOFF,fileop(1))
     !    filename = "SparseMatrix.bin"
     !     CALL StoreSparseMatrix(FineGrid%PreconditioningMatrix,filename,formattype)
     !     print*,'Preconditioningmatrix stored in SparseMatrix.bin.'		
     CALL FactorPreconditioner(FineGrid%PreconditioningMatrix_CSR, &
          (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),fileop(1))

  ELSE IF (Precond==3) THEN

     WRITE(6,3000) maxit,cyclet,nu(1),nu(2)

     ! Prepare for Multigrid
     CALL MGPreProcess ( FineGrid, GhostGridX, GhostGridY, GhostGridZ, MGCoarseningStrategy, alphaprecond, betaprecond, &
          gammaprecond, Precond, MGmaxgrids, CurvilinearONOFF)
  ENDIF
  !
  ! DETERMINE HIGH-ORDER FINITE DIFFERENCE STENCILS
  ! Now, determine fullrank stencils for the x- , y- and z- directions;
  print *,'Determine finite difference stencils for the system matrix...'
  write(fileop(1),*)'Determine finite difference stencils for the system matrix...'
  IF (curvilinearONOFF==0) THEN
     CALL PreProcessDiffStencils(FineGrid,FineGrid%DiffStencils,GhostGridX,GhostGridY,GhostGridZ,alpha,beta,gamma,fileop(1))
     ! GD: Determine the cross derivatives coefficients
     CALL ConstructTableCrossDerivatives(FineGrid, FineGrid%DiffStencils, gamma, GhostGridX, GhostGridY, GhostGridZ, 0)
  ELSE
     ! needed for differentiation stencils in the vertical
     ! CALL PreProcessDiffStencils(FineGrid,FineGrid%DiffStencils,GhostGridX,GhostGridY,GhostGridZ,alpha,beta,gamma)
     !GD: change
     CALL PreProcessDiffStencilsZ(FineGrid,FineGrid%DiffStencils,GhostGridZ,gamma)
     kappa = alpha
     IF (alpha/=beta) THEN
        ! FIXME: just picking the largest of alpha and beta here... perhaps check that they are equal in 3D
        kappa = MAX(alpha,beta)
     END IF
     CALL DetermineGenericStencils(FineGrid%CurvilinearStuff%DiffStencils,kappa)

     CALL DetermineCurvilinearTransform2D(FineGrid,alpha,beta,gamma,GhostGridX,GhostGridY,GhostGridZ)
     ! determine normal vectors at boundary nodes for the 2D plane boundaries
     CALL ComputeNormalVectors(FineGrid,GhostGridX,GhostGridY,GhostGridZ)

     ! Determine linear sigma-coefficients
     ALLOCATE(FineGrid%dsigmanew(FineGrid%Nz+GhostGridZ,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,5))
     FineGrid%dsigmanew = zero
     CALL ALLOCATE_Wavefield_Type(Wavefield_tmp, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, GhostGridX, GhostGridy, GhostGridZ, 0)
     CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
          FineGrid,FineGrid%dsigmanew,Wavefield_tmp)
     CALL DEALLOCATE_Wavefield_Type(Wavefield_tmp, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, 0)

     ! GD: Determine the cross derivatives coefficients
     CALL ConstructTableCrossDerivatives_Curvilinear(FineGrid, FineGrid%CurvilinearStuff%DiffStencils, kappa, &
          GhostGridX, GhostGridY, GhostGridZ)
  END IF
  print*,'...done!'
  write(fileop(1),*)'...done!'
  !
  ! GD: Test to define correct initial spatial derivaties...
  CALL DifferentiationsFreeSurfacePlane(Wavefield,GhostGridX,GhostGridY,FineGrid,alpha,beta)



  !************************************************************************
  !
  ! Debugging setups
  !
  !************************************************************************

  IF (0==1) THEN

     ! Output preconditioner
     IF (Precond==1) THEN
        filename = "SparseMatrix.bin"
        CALL StoreSparseMatrix(FineGrid%PreconditioningMatrix,filename,formattype)
        print*,'Preconditioningmatrix stored in SparseMatrix.bin.'		
     END IF

     ! IF SIMULATION IS LINEAR THEN DETERMINE THE SIGMA-COEFFICIENTS FOR THE TRANSFORMED LAPLACE PROBLEM
     IF (LinearONOFF==0) THEN
        CALL ALLOCATE_Wavefield_Type(Wavefield_tmp, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, GhostGridX, GhostGridy, GhostGridZ, 0)
        CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
             FineGrid%Nz+GhostGridZ,FineGrid,FineGrid%dsigmanew,Wavefield_tmp)
        CALL DEALLOCATE_Wavefield_Type(Wavefield_tmp, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, 0)
     ENDIF

     ! Output linear system matrix 
     ALLOCATE(ee((FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)))
     ALLOCATE(tm((FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)))
     ALLOCATE(A((FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
          (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)))
     ee = zero
     tm = zero
     A  = zero
     DO i = 1 , (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)
        ee(i) = one
        IF (curvilinearONOFF==1) THEN
           !   CALL BuildLinearSystemTransformedCurvilinear(FineGrid, ee, tm,GhostGridX,GhostGridY,GhostGridZ,kappa)
           CALL BuildLinearSystemTransformedCurvilinear(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY, &
                FineGrid%Nz+GhostGridZ,ee,tm,FineGrid,alpha,beta,gamma)
        ELSE
           CALL BuildLinearSystem(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY, &
                FineGrid%Nz+GhostGridZ,ee,tm,FineGrid,alpha,beta,gamma)
        END IF
        A(:,i) = tm
        ee(i) = zero
     END DO
     filename = "A.bin"
     CALL StoreRealArray(A,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
          (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),filename,formattype)
     IF (curvilinearONOFF==1) THEN
        print*,'Linear coefficient matrix A (curvilinear routine) stored in A.bin.'		
     ELSE
        print*,'Linear coefficient matrix A stored in A.bin.'		
     END IF
     print*,'curvilinearONOFF=',curvilinearONOFF

     ! save the vertical derivative for linear stability analysis...
     ee = zero
     tm = zero
     A  = zero
     DO i = 1 , (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ)
        ee(i) = one
        CALL DiffZArbitrary(ee,tm,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ, &
             FineGrid%DiffStencils,gamma)
        A(:,i) = tm
        ee(i) = zero
     END DO
     filename = "DMz.bin"
     CALL StoreRealArray(A,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
          (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),filename,formattype)
     print*,'Matrix DMz stored in DMz.bin.'		

     DEALLOCATE(ee,tm,A)
     !   print*,'stopped here for now...'
     stop

  END IF

  !************************************************************************
  !
  ! Step through time and compute the wave flow.
  !
  !************************************************************************

  ! STORE INITIAL CONDITION
  IF (StoreDataONOFF>0) THEN
     ! GD: SWENSE storage if necessary
     IF(swenseONOFF/=0) THEN
        ! Even numbers will be scattered wavefield
        CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%P,FineGrid,0,formattype)
        ! Odd numbers will be total wavefield
        CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E+Wavefield%E_I,Wavefield%P+Wavefield%P_I_s,&
             FineGrid,1,formattype)
     ELSE
        CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%P,FineGrid,0,formattype)
     ENDIF
     ! GD addition : bottom profile
     CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%h,FineGrid%hx,&
          FineGrid,99999,formattype)
  ELSEIF(StoreDataOnOff<0)THEN
     CALL StoreDataAscii(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%P,FineGrid,0)
     ! Also store bottom profile
     CALL StoreDataAscii(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%h,FineGrid%hx,FineGrid,899)
  ENDIF
  !
  ! Open and initialize the kinematics output file(s) if called for
  !
  If(iKinematics/=0)THEN
        IF (formattype==20) THEN ! Store binary kinematics files
             DO i=1,nOutFiles
                        CALL StoreKinematicData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                     FineGrid%Nz+GhostGridZ,i,0)
             END DO
        ELSEIF (formattype==30) THEN ! Store wave gauges in ASCII format
               CALL StoreWaveGauges(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                     FineGrid%Nz+GhostGridZ,2,0) 
        ENDIF
  END If

  ! IF SIMULATION IS LINEAR THEN DETERMINE THE SIGMA-COEFFICIENTS FOR THE TRANSFORMED LAPLACE PROBLEM
  IF (LinearONOFF==0) THEN
     CALL ALLOCATE_Wavefield_Type(Wavefield_tmp, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, GhostGridX, GhostGridy, GhostGridZ, 0)
     CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
          FineGrid%Nz+GhostGridZ,FineGrid,FineGrid%dsigmanew,Wavefield_tmp)
     CALL DEALLOCATE_Wavefield_Type(Wavefield_tmp, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, 0)
  ENDIF
  !
  ! A number of convergence and consistency checks are coded up in this subroutine.  
  !
  IF(0==1) THEN
     Call TestCodeDucrozet2
  END IF
  !
  ! If kinematics output is requested save the initial conditions
  !
  tstep=0
  IF(iKinematics/=0)THEN
       IF (iKinematics==20) THEN ! Store binary kinematics files
             Do i=1,nOutFiles
                IF (tstep+1 >= Output(i)%tbeg .and. tstep+1 <= Output(i)%tend .and.  &
                     mod(tstep,Output(i)%tstride)==0 )THEN
                           CALL StoreKinematicData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                                FineGrid%Nz+GhostGridZ,i,tstep)
                END IF
             END Do
       ELSEIF (iKinematics==30) THEN ! Store wave gauges in ASCII format
               CALL StoreWaveGauges(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                     FineGrid%Nz+GhostGridZ,2,tstep) 
                       
       ENDIF
  END IF

  ! For coupling OceanWave3D with OpenFOAM we need the velocities and free
  ! surface elevation to be available at all times.
  ! FIXME: Is there a better solution where fields are only allocated if needed?
  ! General problem!!
  !botp
  ALLOCATE( &
       UOF(FineGrid%Nz+GhostGridZ,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY), &
       VOF(FineGrid%Nz+GhostGridZ,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY), &
       WOF(FineGrid%Nz+GhostGridZ,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
  ALLOCATE(dOF(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))


2010 FORMAT(/, '*********************************************************',/,&
       '***                                                   ***',/,&
       '*** OceanWave3D - a coastal engineering tool for      ***',/,&
       '*** simulation of nonlinear free surface waves.       ***',/,&
       '*** Copyright (C) 2009 Allan P. Engsig-Karup.         ***',/,&
       '***                                                   ***',/,&
       '*** This OceanWave3D program comes with ABSOLUTELY NO ***',/,&
       '*** WARRANTY. It is distributed under the conditions  ***',/,&
       '*** of the GNU General Public License version 3.      ***',/,&
       '***                                                   ***',/,&
       '***     Software library developed in 2009 by         ***',/,&
       '***                                                   ***',/,&
       '***     Allan P. Engsig-Karup                         ***',/,&
       '***     Guillaume Ducrozet                            ***',/,&
       '***                                                   ***',/,&
       '*** At DTU Informatics                                ***',/,&
       '***    Scientific Computing Section                   ***',/,&
       '***    Technical University of Denmark                ***',/,&
       '***                                                   ***',/,&
       '***     Original software library written in 2007 by  ***',/,&
       '***                                                   ***',/,&
       '***     Allan P. Engsig-Karup                         ***',/,&
       '***     Harry B. Bingham                              ***',/,&
       '***                                                   ***',/,&
       '*** At Department of Mechanical Engineering           ***',/,&
       '***    Coastal, Maritime and Structural Eng. Section  ***',/,&
       '***    Technical University of Denmark                ***',/,&
       '***                                                   ***',/,&
       '*********************************************************',/)

3000 FORMAT('   Preconditioning: MG-',I1,A1,'(',I1,',',I1,')')

END SUBROUTINE OceanWave3DT0Setup
