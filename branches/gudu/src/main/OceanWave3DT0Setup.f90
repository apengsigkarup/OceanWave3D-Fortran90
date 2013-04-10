SUBROUTINE OceanWave3DT0Setup
  !
  ! This subroutine performs all of the initial set-up work before the time-stepping 
  ! begins.  
  !
  ! By Allan P. Engsig-Karup.
  USE GlobalVariables
  USE MGLevels
  IMPLICIT NONE
  EXTERNAL BuildLinearSystem, BuildLinearSystemTransformedCurvilinear
  ! GD: to test the cross derivatives...
  REAL(KIND=long), DIMENSION(:,:,:), ALLOCATABLE :: tmpPHI
  TYPE (Diff_def)        :: FullRankStencils
  INTEGER i, j, k

  ! OUTPUT HEADER TO SCREEN
  WRITE (6,2010)

  CALL Initialize
  CALL ReadInputFileParameters

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
  print*,'do initialization...'
  CALL InitializeVariables
  !
  ! We start at time=0 here but if this is a hot start, time0 will be read in SetupInitialConditions.
  !
  time=zero
!  !
!  ! Set up the initial conditions and the bathymetry data
!  !
!  print*,'setup ICs...'
!  CALL SetupInitialConditions
!  time=time0
!  print*,'done with ICs'
  !
  ! Set up for relaxations zones and wave generation.
  !
  IF (relaxONOFF>0) THEN
     CALL PreprocessRelaxationZones
     PRINT*,'  Relaxation zones have been setup.'
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
        !print*,'dxmin=',dxmin
        !print*,'dymin=',dymin
        !print*,'min dx,dy =',MIN(dxmin,dymin)
        !print*,'c=',c
        !print*,'CFL=',CFL
        dt = CFL*MIN(dxmin,dymin)/c
        PRINT*,''
        PRINT*,'  Time step size modified based on the incident wave. dt = ',dt
        PRINT*,'  Courant number, Cr = ',c*dt/MIN(dxmin,dymin)
        PRINT*,'  Discrete anisotropy, ds/dx = ',dsigmamin/dxmin
     END IF
     IF (IncWaveType==1 ) THEN
        CALL stream_func_set_up(g,SFsol%h,SFsol%T,SFsol%i_wavel_or_per,SFsol%L,  &
             SFsol%k,SFsol%HH,SFsol%i_euler_or_stokes,SFsol%e_or_s_vel,SFsol%i_deep_or_finite,  &
             SFsol%n_h_steps,SFsol%n_four_modes,SFsol%nwrk,SFsol%yy,SFsol%zz)
        PRINT*,'     SF solution: k=',SFsol%k,',h=',SFsol%h,',H=',SFsol%HH,',T=',SFsol%T,',L=',SFsol%L
        PRINT*,'                  kh=',SFsol%k*SFsol%h,',c=',SFsol%L/SFsol%T

        WRITE(6,62)two*pi/SFsol%k,SFsol%T,SFsol%zz(2)/SFsol%k,                  &
             SFsol%zz(5)*sqrt(g/SFsol%k), SFsol%zz(6)*sqrt(g/SFsol%k)
62      FORMAT(' The incident wave is a stream function wave with L=',  &
             e12.6,/,'  T=',e12.6,' H=',e12.6,' u_E=',e12.6,' and u_S=',     &
             e12.6,//)
     ELSEIF(IncWaveType==2)THEN
        !
        ! A linear regular or irregular wave
        !
        print *, ' '
        IF(RandomWave(1)%ispec==-1)THEN
           WRITE(6,70)RandomWave(1)%Tp,RandomWave(1)%Hs
70         FORMAT(' The incident wave is a linear mono-chromatic wave with',/,&
                ' T=', e10.4,' and H=',e10.4,'.',//)
        ELSEIF(RandomWave(1)%ispec==0)THEN
           WRITE(6,71)RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%seed,RandomWave(1)%seed2
71         FORMAT(' The incident wave is a 2D P-M spectrum with',/,&
                ' T_p=', e10.4,' and H_s=',e10.4,', seed values are:',/,2i10,//)
        ELSEIF(RandomWave(1)%ispec==1)THEN
           WRITE(6,72)RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%seed,RandomWave(1)%seed2
72         FORMAT(' The incident wave is a 2D JONSWAP spectrum with ',/, &
                'T_p=',  e10.4,' and H_s=',e10.4,', seed values are:',/,2i10,//)
        ELSEIF(RandomWave(1)%ispec==2)THEN
           WRITE(6,73)RandomWave(1)%inc_wave_file
73         FORMAT(' The incident wave will be read from file ',a30,/)
        ELSEIF(RandomWave(1)%ispec==3)THEN
           WRITE(6,74)RandomWave(1)%beta,RandomWave(1)%Tp,RandomWave(1)%Hs,RandomWave(1)%seed,RandomWave(1)%seed2
74         FORMAT(' The incident wave is a 3D JONSWAP spectrum with Normal spreading at heading angle ',e10.4,/, &
                'to the x-axis. T_p=',  e10.4,' H_s=',e10.4,', seed values are:',/,2i10,//)
        ELSE
           PRINT *, 'ERROR:  RandomWave%ispec must be -1,0,1,2 or 3.'
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
        ! The grid is also assumed to be uniform in the relaxation/generation zone(s) and 
        ! we expect to generate in X-directed relaxation zones.  Other areas are 
        ! not yet supported.  -HBB
        !
        !
        ! Count the total number of grid points in each generation zone and allocate space
        ! for eta and phiS
        !  
        DO i=1,relaxNo
           n_wavem=0
           If(RelaxZones(i)%XorYgen=='X' .AND. RelaxZones(i)%WavegenONOFF==1) THEN
              n_wavem=n_wavem+RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1
              RandomWave(i)%dx=FineGrid%x(RelaxZones(i)%idx(1)+1,1)  &
                   -FineGrid%x(RelaxZones(i)%idx(1),1)
              j_wavem=nint(RandomWave(1)%x0/RandomWave(1)%dx)+2
              n_wavem=n_wavem+GhostGridX ! Include the ghost point at the left boundary.  
              print *, 'The generated wave is centered at x=',FineGrid%x(j_wavem,1),      &
                   ' in a depth of',RandomWave(i)%h0,', the generation zone contains ',   &
                   n_wavem, ' points.'
              RandomWave(i)%nf=FineGrid%nx
              !
              ALLOCATE(RandomWave(i)%eta(n_wavem,n_fft), RandomWave(i)%Phis(n_wavem,n_fft), &
                   RandomWave(i)%eta0(max(n_fft,RandomWave(i)%nf)),                         &
                   RandomWave(i)%Phis0(max(n_fft,RandomWave(i)%nf)) )
              CALL random_wave_signal(RandomWave(i)%ispec, n_fft, n_wavem, j_wavem-1,       &
                   RandomWave(i)%dx, dt, RandomWave(i)%Tp, RandomWave(i)%Hs, RandomWave(i)%h0, &
                   g, RandomWave(i)%inc_wave_file, RandomWave(i)%kh_max, RandomWave(i)%seed, &
                   RandomWave(i)%seed2, RandomWave(i)%eta, RandomWave(i)%Phis,               &
                   RandomWave(i)%eta0, RandomWave(i)%Phis0, RandomWave(i)%nf, time0)
           END If
        END DO
     ENDIF
  ENDIF
  !
  ! Set up the Pressure Damping Zones if any.
  !
  If (NDampZones /=0) THEN
     Call PreprocessPDampingZones
     print *, ' '
     print *, 'Pressure damping zones are set up'
     print *, ' '
  END If
  !
  ! Set up the initial conditions and the bathymetry data
  ! FIXME: is it needed somewhere to have this setup before
  ! PreprocessRelaxationZones?
  !
  print*,'setup ICs...'
  CALL SetupInitialConditions
  time=time0
  print*,'done with ICs'
  IF (.FALSE.) THEN
     !
     ! Test code to validate buildlinearsystem subroutines.
     !
     CALL TestCodeDucrozet1
  END IF
  !
  IF (DetermineBottomGradients==1) THEN
     PRINT*,'Error: no support yet for determining bottom gradients numerically in current implementation. (APEK)'
     STOP
     !	  ! Determine Bottom Gradients numerically
     !	  CALL PreProcessDiffStencils(FineGrid,FineGrid%DiffStencils,GhostGridX,GhostGridY,GhostGridZ,alpha,beta,gamma)
     !	  IF (FineGrid%Nx==1) THEN
     !	     FineGrid%hx = zero; FineGrid%hxx = zero;
     !  	  ELSE
     !         CALL DiffXEven(FineGrid%h,FineGrid%hx, 1,FineGrid%Nx,FineGrid%Ny,1,FineGrid%DiffStencils,alpha)
     !         CALL DiffXEven(FineGrid%h,FineGrid%hxx,2,FineGrid%Nx,FineGrid%Ny,1,FineGrid%DiffStencils,alpha)
     !	  END IF
     !	  IF (FineGrid%Ny==1) THEN
     !	     FineGrid%hy  = zero; FineGrid%hyy = zero;
     !	  ELSE
     !         CALL DiffYEven(FineGrid%h,FineGrid%hy, 1,FineGrid%Nx,FineGrid%Ny,1,FineGrid%DiffStencils,beta)
     !         CALL DiffYEven(FineGrid%h,FineGrid%hyy,2,FineGrid%Nx,FineGrid%Ny,1,FineGrid%DiffStencils,beta)
     !  	  ENDIF
     !   	  DEALLOCATE(FineGrid%DiffStencils%StencilX,FineGrid%DiffStencils%StencilY,FineGrid%DiffStencils%StencilZ)
  ENDIF

  IF (Precond==1) THEN ! PREPARE FOR PRECONDITIONING
     ! DETERMINE LOW-ORDER FINITE DIFFERENCE STENCILS
     ! FIXME: make it possible to choose the order of the preconditioner in the input file
     WRITE(*,FMT='(A,I2,A)') '   Preconditioner: DIRECT LU (',2*alphaprecond,' order, linear)'

     CALL PreparePreconditioner(FineGrid%PreconditioningMatrix,FineGrid,GhostGridX, GhostGridY, GhostGridZ, &
          alphaprecond, betaprecond, gammaprecond, Precond, CurvilinearONOFF)
!    filename = "SparseMatrix.bin"
!     CALL StoreSparseMatrix(FineGrid%PreconditioningMatrix,filename,formattype)
!     print*,'Preconditioningmatrix stored in SparseMatrix.bin.'		
     CALL FactorPreconditioner(FineGrid%PreconditioningMatrix, &
          (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ))

  ELSE IF (Precond==3) THEN

     WRITE(6,3000) maxit,cyclet,nu(1),nu(2)

     ! Prepare for Multigrid
     CALL MGPreProcess ( FineGrid, GhostGridX, GhostGridY, GhostGridZ, MGCoarseningStrategy, alphaprecond, betaprecond, &
          gammaprecond, Precond, MGmaxgrids, CurvilinearONOFF)
  ENDIF
  
  !
  ! DETERMINE HIGH-ORDER FINITE DIFFERENCE STENCILS
  ! Now, determine fullrank stencils for the x- , y- and z- directions;
  print*,'Determine finite difference stencils for the system matrix...'
  IF (curvilinearONOFF==0) THEN
     CALL PreProcessDiffStencils(FineGrid,FineGrid%DiffStencils,GhostGridX,GhostGridY,GhostGridZ,alpha,beta,gamma)
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
!
  ! GD: Test to define correct initial spatail derivaties...
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
  ENDIF
  !
  ! Open and initialize the kinematics output file(s) if called for
  !
  If(iKinematics/=0)THEN
     Do i=1,nOutFiles
        CALL StoreKinematicData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
             FineGrid%Nz+GhostGridZ,i,0)
     END Do
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
     Do i=1,nOutFiles
        IF (tstep+1 >= Output(i)%tbeg .and. tstep+1 <= Output(i)%tend .and.  &
             mod(tstep,Output(i)%tstride)==0 )THEN
           CALL StoreKinematicData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                FineGrid%Nz+GhostGridZ,i,1)
        END IF
     END Do
  END IF

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
