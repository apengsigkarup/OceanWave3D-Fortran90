SUBROUTINE ReadInputFileParameters
  !>
  !! Read model parameters from ascii input file.
  !!
  !! If this file is modified, make sure to include strict tests and
  !! provide detailed error messages to user in case of errors.
  !!
  !! By Allan P. Engsig-Karup.
  !<
  USE GlobalVariables
  USE MGLevels
  IMPLICIT NONE
  INTEGER ios, i, nxIC, nyIC, iflag_phi, ispec, nGenZones, seed, seed2 
  REAL(kind=long) :: xtankIC, ytankIC, t0IC, Tp, Hs, h0, kh_max, x0, &
       y0, beta0, s0, gamma_jonswap
  CHARACTER(len=30):: inc_wave_file

  READ (FILEIP(1),'(A)',ERR=100,IOSTAT=ios) HEAD(1)
  WRITE (*,FMT='(A,A/)') '   Input file with model parameters : ', filenameINPUT
  WRITE (*,FMT='(A,A/)') '   Header title.................... : ', HEAD(1)
  WRITE (fileop(1),FMT='(A,A/)') '   Input file with model parameters : ', filenameINPUT
  WRITE (fileop(1),FMT='(A,A/)') '   Header title.................... : ', HEAD(1)

  READ (FILEIP(1),*) IC
  BACKSPACE(FILEIP(1))
  READ (FILEIP(1),*,ERR=22) IC, IncWaveType
  BACKSPACE(FILEIP(1))
  READ (FILEIP(1),*,ERR=23) IC, IncWaveType, accel_tol_fact

  WRITE(*,FMT='(A/)') '   PARAMETER INPUT'
  WRITE(*,FMT='(A,I3/)') '   Initial Condition (IC), Predefined : ',IC
  WRITE(*,FMT='(A,F8.4)') '   Local filtering downward vertical acceleration limit is: ', &
       accel_tol_fact
  WRITE(*,FMT='(A//)') ' * g m/s^2.  Theoretically breaking should occur between 0.5 and 1.'
  WRITE(fileop(1),FMT='(A/)') '   PARAMETER INPUT'
  WRITE(fileop(1),FMT='(A,I3/)') '   Initial Condition (IC), Predefined : ',IC
  WRITE(fileop(1),FMT='(A,F8.4)') '   Local filtering downward vertical acceleration limit is: ', &
       accel_tol_fact
  WRITE(fileop(1),FMT='(A//)') ' * g m/s^2.  Theoretically breaking should occur between 0.5 and 1.'

  GO TO 24

  ! If no parameters are specified after IC on the same line choose default values
22 IncWaveType=0;  accel_tol_fact=1000. ! Default, no incident wave, no local smoothing.
  GO TO 24
23 accel_tol_fact=1000.  ! Default, no local filtering to prevent breaking.  
24 Continue

  READ (FILEIP(1),*,ERR=105) Lx, Ly, Lz, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, GridX, GridY, GridZ, GhostGridX, &
       GhostGridY, GhostGridZ
  IF(Lz<=0)THEN
     ! If Lz<=0 then we read the bathymetry definition from a file
     BACKSPACE (FILEIP(1))
     READ (FILEIP(1),*) Lx, Ly, Lz, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, GridX, GridY, GridZ, GhostGridX, &
          GhostGridY, GhostGridZ, fname_bottom
     WRITE(*,FMT='(A,A,A)') '   Bathymetry: ',fname_bottom,' (file)'
     WRITE(fileop(1),FMT='(A,A,A)') '   Bathymetry: ',fname_bottom,' (file)'
  ELSE
     WRITE(*,FMT='(A)') '   Bathymetry                         : Predefined (in setup for IC)'
     WRITE(fileop(1),FMT='(A)') '   Bathymetry                         : Predefined (in setup for IC)'
  END IF
  !
  WRITE (*,900) '   (Nx,Ny,Nz)=(',FineGrid%Nx,',',FineGrid%Ny,',',FineGrid%Nz,')'
  IF (FineGrid%Nx>1) THEN
     dx = Lx/(FineGrid%Nx-one)
     WRITE (*,903) '   Size of dx: ',dx
     WRITE (fileop(1),903) '   Size of dx: ',dx
  ELSE
     WRITE (*,FMT='(/A)') '   The problem is without an x-dimension.'
     WRITE (fileop(1),FMT='(/A)') '   The problem is without an x-dimension.'
     GhostGridX=0
  ENDIF
  IF (FineGrid%Ny>1) THEN
     dy = Ly/(FineGrid%Ny-one)
     WRITE (*,903) '   Size of dy: ',dy
     WRITE (fileop(1),903) '   Size of dy: ',dy
  ELSE
     WRITE (*,FMT='(/A)') '   The problem is without a y-dimension.'
     WRITE (fileop(1),FMT='(/A)') '   The problem is without a y-dimension.'
     GhostGridY=0
  ENDIF
  IF (GridZ==0) THEN
     WRITE (*,FMT='(/A)') '   Even node-distribution in Z'
     WRITE (fileop(1),FMT='(/A)') '   Even node-distribution in Z'
  ELSEIF (GridZ==1) THEN
     WRITE (*,FMT='(/A)') '   Uneven node-distribution in Z'
     WRITE (fileop(1),FMT='(/A)') '   Uneven node-distribution in Z'
  ENDIF
  IF (GhostGridZ==1) THEN
     WRITE (*,FMT='(A/)') '   Ghost point layer included below bottom.'
     WRITE (fileop(1),FMT='(A/)') '   Ghost point layer included below bottom.'
  ELSE IF (GhostGridZ==0) THEN
     WRITE (*,FMT='(A/)') '   Kinematic condition will be imposed directly.'
     WRITE (fileop(1),FMT='(A/)') '   Kinematic condition will be imposed directly.'
  ELSE
     WRITE (*,FMT='(A/)') '   Error: GhostGrid not 0 or 1.'; STOP
     WRITE (fileop(1),FMT='(A/)') '   Error: GhostGrid not 0 or 1.'; STOP
  ENDIF

  IF(IC<0)THEN
     OPEN(unit=fileip(3),file='OceanWave3D.init',status='old')
     READ(fileip(3),'(A)')head(3)
     PRINT *,'INPUT:  Initial conditions will be read from OceanWave3D.init with header'
     PRINT *, head(3)
     print *,' '
     WRITE(FILEOP(1),*)'INPUT:  Initial conditions will be read from OceanWave3D.init with header'
     WRITE(FILEOP(1),*) head(3)
     write(fileop(1),*)' '
     READ(fileip(3),*)xtankIC,ytankIC,nxIC,nyIC,t0IC
     IF ((nxIC-FineGrid%Nx .ne. 0) .or. (nyIC-FineGrid%Ny .ne. 0)) THEN
        PRINT *, 'Your input and initial conditions grids must agree.'
        PRINT *, 'nx=',FineGrid%Nx,' nx IC=',nxIC,' ny=',FineGrid%Ny,' ny IC=',nyIC
        STOP
     END IF
  END IF


  READ (FILEIP(1),*) alpha, beta, gamma, alphaprecond,betaprecond,gammaprecond
  ! Use same stencil size in each direction
  !  beta  = alpha
  !  gamma = alpha
  IF (FineGrid%Nx==1) THEN
     alpha = 0
     alphaprecond = 0
  ENDIF
  IF (FineGrid%Ny==1) THEN
     beta = 0
     betaprecond = 0
  ENDIF
  WRITE (*,901) '   Half-width stencils: (alpha,beta,gamma)=(',alpha,',',beta,',',gamma,')'
  WRITE (*,901) '   Half-width stencils: (alpha,beta,gamma)=(',alphaprecond,',',betaprecond,',',gammaprecond,') (Preconditioner)'
  WRITE (fileop(1),901) '   Half-width stencils: (alpha,beta,gamma)=(',alpha,',',beta,',',gamma,')'
  WRITE (fileop(1),901) '   Half-width stencils: (alpha,beta,gamma)=(',alphaprecond,',',betaprecond,',',gammaprecond,') (Preconditioner)'
  IF (2*alpha+1>FineGrid%Nx .AND. FineGrid%Nx>1) THEN
     GOTO 101
  ENDIF
  IF (2*alpha+1>FineGrid%Nx .AND. FineGrid%Ny==1) THEN
     GOTO 101
  ENDIF
  IF (2*beta+1>FineGrid%Ny .AND. FineGrid%Ny>1) THEN
     GOTO 102
  ENDIF
  IF (2*beta+1>FineGrid%Ny .AND. FineGrid%Nx==1) THEN
     GOTO 102
  ENDIF
  IF (2*gamma+1-GhostGridZ>FineGrid%Nz .AND. FineGrid%Nz>1) THEN
     GOTO 103
  ENDIF
  !
  !
  READ (FILEIP(1),*,err=52) Nsteps, dt, timemethod, CFL, extrapolationONOFF, time0
  GO TO 53
52 Backspace(fileip(1))
  !
  ! If time0 is not input then we either set it to zero, or if an initial condition is 
  ! being read in from a file, then the time specified there is used.  
  !
  READ (FILEIP(1),*,err=54) Nsteps, dt, timemethod, CFL, extrapolationONOFF
  If(IC >= 0) THEN
     time0=zero
  ELSE
     time0=t0IC
  END If
  Go To 53
54 print *, 'Error reading the time set up line: Nsteps, dt, time integration scheme, ...'
  stop
53 print *, ' '
  !
  ! If this is a read initial conditions run and the starting times do not agree, take the 
  ! initial condition time.  
  !
  IF(IC<0 .and. abs(time0-t0IC) > 10.E-10)THEN
     time0=t0IC
     print *, ' '
     print *, '**  Your run starting time from OceanWave3D.inp does not agree with the '
     print *,'   one from OceanWave3D.init.  The init value will be used. **'
     print *, ' '
     write(fileop(1),*) ' '
     write(fileop(1),*) '**  Your run starting time from OceanWave3D.inp does not agree with the '
     write(fileop(1),*)'   one from OceanWave3D.init.  The init value will be used. **'
     write(fileop(1),*) ' '
  END IF

  print *, 'Starting time for this run is ',time0
  WRITE (*,902) '   Number of time steps chosen: ', Nsteps
  write(fileop(1),*), 'Starting time for this run is ',time0
  WRITE (fileop(1),902) '   Number of time steps chosen: ', Nsteps
  IF (CFL/=zero) THEN
     ! GD: FIXME for 2D-3D case with propagation along y...
     IF(FineGrid%Nx==1) THEN
        dt = CFL*dy !/c
        WRITE (*,903) '   CFL constant given by user (dt=CFL*dy)', CFL
        WRITE (fileop(1),903) '   CFL constant given by user (dt=CFL*dy)', CFL
     ELSE
        dt = CFL*dx !/c
        WRITE (*,903) '   CFL constant given by user (dt=CFL*dx)', CFL
        WRITE (fileop(1),903) '   CFL constant given by user (dt=CFL*dx)', CFL
     ENDIF
  ENDIF
  WRITE (*,903) '   Size of time increment: ', dt
  WRITE (fileop(1),903) '   Size of time increment: ', dt

  SELECT CASE (timemethod)
  CASE (1) ! Classical RK4
     RKSTAGES = 4
     WRITE (*,*) '   Time-integration method: Classical Runge-Kutta fourth order'
     WRITE (fileop(1),*) '   Time-integration method: Classical Runge-Kutta fourth order'
     IF (extrapolationONOFF==1) THEN
        WRITE (*,*) '    - Optimization of RK scheme using extrapolation on seperate RK stages will be employed.'
        WRITE (fileop(1),*) '    - Optimization of RK scheme using extrapolation on seperate RK stages will be employed.'
     ENDIF
  CASE (2) ! Carpenter & Kennedy low-storage RK45
     !


     RKSTAGES = 5
     WRITE (*,*) '   Time-integration method: Low-storage five-stage Runge-Kutta fourth order (Carpenter & Kennedy)'
     WRITE (fileop(1),*) '   Time-integration method: Low-storage five-stage Runge-Kutta fourth order (Carpenter & Kennedy)'
  CASE DEFAULT
     WRITE (*,*) 'Error: Chosen time integration method not valid.'
     STOP
  END SELECT

  READ (FILEIP(1),*,err=139) g, rho
  Go To 140
139  BACKSPACE(fileip(1))
  READ (FILEIP(1),*,err=139) g
  Print *, '  ** No input value for rho, using 1000.'
  write(fileop(1),*) '  ** No input value for rho, using 1000.'
  rho=1000.
140 continue

  READ (FILEIP(1),*,err=141) solver, Precond, MGCoarseningStrategy, GMRESmaxiterations, reltol, abstol, maxit, cyclet, &
       nu(1), nu(2), MGmaxgrids
  GO TO 142
141 solver=1
  BACKSPACE(fileip(1))
  READ (FILEIP(1),*,err=141) Precond, MGCoarseningStrategy, GMRESmaxiterations, reltol, abstol, maxit, cyclet, &
       nu(1), nu(2), MGmaxgrids
142  SELECT CASE (solver)
      CASE(0)
         WRITE(*,*) '   Defect correction (DC) method is chosen.'
         WRITE(fileop(1),*) '   Defect correction (DC) method is chosen.'
      CASE DEFAULT
         WRITE(*,*) '   GMRES method is chosen.'
         WRITE(fileop(1),*) '   GMRES method is chosen.'
  END SELECT
  IF (Precond==1) THEN
  SELECT CASE (solver)
      CASE(0)
         WRITE(*,*) '   Strategy: DC + LU (order ',2*alphaprecond,')'
         WRITE(fileop(1),*) '   Strategy: DC + LU (order ',2*alphaprecond,')'
      CASE DEFAULT
         WRITE(*,*) '   Strategy: GMRES + LU (order ',2*alphaprecond,')'
         WRITE(fileop(1),*) '   Strategy: GMRES + LU (order ',2*alphaprecond,')'
  END SELECT
  ELSE IF (Precond==3) THEN
  SELECT CASE (solver)
      CASE(0)
         WRITE(*,*) '   Strategy: DC + MG-RB-',cyclet,'(',nu(1),',',nu(2),')'
         WRITE(fileop(1),*) '   Strategy: DC + MG-RB-',cyclet,'(',nu(1),',',nu(2),')'
      CASE DEFAULT
         WRITE(*,*) '   Strategy: GMRES + MG-RB-',cyclet,'(',nu(1),',',nu(2),')'
         WRITE(fileop(1),*) '   Strategy: GMRES + MG-RB-',cyclet,'(',nu(1),',',nu(2),')'
  END SELECT
  END IF
  WRITE(*,*) '   Tolerance levels user-defined. RelTol = ',reltol,' and AbsTol = ',abstol
  WRITE(fileop(1),*) '   Tolerance levels user-defined. RelTol = ',reltol,' and AbsTol = ',abstol
  IF (Precond==3 .AND. GhostGridZ/=1 ) THEN
     GOTO 104
  ENDIF

  WRITE (*,*) ''
  WRITE (fileop(1),*) ''

  ! STREAM FUNCTION SOLUTION PARAMETERS
  READ (FILEIP(1),*) SFsol%HH, SFsol%h, SFsol%L, SFsol%T, SFsol%i_wavel_or_per, SFsol%e_or_s_vel, &
       SFsol%i_euler_or_stokes, SFsol%n_h_steps, SFsol%n_four_modes
  SFsol%k = two*pi/SFsol%L
  SFsol%nwrk=2*SFsol%n_four_modes+10
  SFsol%g = g
  ALLOCATE ( SFsol%yy(SFsol%nwrk), SFsol%zz(SFsol%nwrk) )
  SFsol%yy = zero; SFsol%zz = zero;

  ! DATA STORAGE
  READ (FILEIP(1),*) StoreDataONOFF, formattype
  IF(formattype==20)THEN
     BACKSPACE(FILEIP(1))
     READ (FILEIP(1),*) StoreDataONOFF, iKinematics, formattype, nOutFiles
     Allocate (Output(nOutFiles))
     IF (nOutFiles>10)THEN
        print *, 'Max. 10 kinematics output files at this point.'
        stop
     END IF
     print *, 'Kinematics output requested in ',nOutFiles,' file(s) named "Kinematics_**.bin".'
     print *, ' '
     write(fileop(1),*) 'Kinematics output requested in ',nOutFiles,' file(s) named "Kinematics_**.bin".'
     write(fileop(1),*) ' '
     Do i=1,nOutFiles
        READ (FILEIP(1),*,err=110)Output(i)%xbeg,Output(i)%xend,Output(i)%xstride,Output(i)%ybeg, &
             Output(i)%yend,Output(i)%ystride,Output(i)%tbeg,Output(i)%tend,Output(i)%tstride
        !
        ! Check that the requested output ranges exist on this grid.
        !
        if ( Output(i)%xbeg<1 .or. Output(i)%xend>FineGrid%Nx .or. Output(i)%xbeg > Output(i)%xend) THEN
           Print *, 'ReadInputFileParameters: Kinematics xrange is invalid'
           stop
        end if
        if(Output(i)%ybeg<1 .or. Output(i)%yend>FineGrid%Ny .or. Output(i)%ybeg > Output(i)%yend ) THEN
           Print *, 'ReadInputFileParameters: Kinematics yrange is invalid'
           stop
        end if
        if(Output(i)%tbeg<1 .or. Output(i)%tend>Nsteps .or. Output(i)%tbeg > Output(i)%tend) THEN
           Print *, 'ReadInputFileParameters: Kinematics trange is invalid'
           stop
        end if
        ! Open the required output files
        OPEN (UNIT=FILEOP(i+1),FILE='Kinematics'//fnt(i)//'.bin',          &
             STATUS='UNKNOWN',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
     END Do
  ELSEIF(formattype==30)THEN
     BACKSPACE(FILEIP(1))
     READ (FILEIP(1),*) StoreDataONOFF, iKinematics, formattype, nOutFiles
     Allocate (Output(nOutFiles))
     IF (nOutFiles>10)THEN
        print *, 'Max. 10 wave gauge files at this point.'
        stop
     END IF
     print *, nOutFiles,' wave gauge(s) requested in file(s) named "waveGauges.dat".'
     print *, ' '
     Do i=1,nOutFiles
        READ (FILEIP(1),*,err=110)Output(i)%xbeg,Output(i)%xend,Output(i)%xstride,Output(i)%ybeg, &
             Output(i)%yend,Output(i)%ystride,Output(i)%tbeg,Output(i)%tend,Output(i)%tstride
        !
        ! Check that the requested output ranges exist on this grid.
        !
        if ( Output(i)%xbeg<1 .or. Output(i)%xend>FineGrid%Nx .or. Output(i)%xbeg > Output(i)%xend) THEN
           Print *, 'ReadInputFileParameters: Kinematics xrange is invalid'
           stop
        end if
        if(Output(i)%ybeg<1 .or. Output(i)%yend>FineGrid%Ny .or. Output(i)%ybeg > Output(i)%yend ) THEN
           Print *, 'ReadInputFileParameters: Kinematics yrange is invalid'
           stop
        end if
        IF (Output(i)%xbeg /= Output(i)%xend) THEN
                PRINT *, 'x_begin not equal x_end: Wave gauges only 1 point per             &
                probe'
                stop
        END IF
        IF (Output(i)%ybeg /= Output(i)%yend) THEN
                PRINT *, 'y_begin not equal y_end: Wave gauges only 1 point per             &
                probe'
                stop
        END IF

     END Do
        ! Open the required output files
        OPEN (UNIT=FILEOP(2),FILE='waveGauges.dat',          &
             STATUS='UNKNOWN',FORM='FORMATTED')
  ELSE
     iKinematics=0
  END IF
  go to 111
110 print *, 'ReadInputFile Parameters:  Error reading kinematics file parameters.'
  print *, 'Failed after reading ',i-1,' output file lines.'
  stop
111 continue

  ! LINEAR/NONLINEAR COMPUTATIONS, applied free-surface pressure
  READ (FILEIP(1),*) LinearONOFF, PressureTermONOFF
  IF (LinearONOFF==0) THEN
     WRITE(*,'(A/)') '   Linear model is employed.'
     WRITE(fileop(1),'(A/)') '   Linear model is employed.'
  ELSE
     WRITE(*,'(A/)') '   Fully nonlinear model is employed.'
     WRITE(fileop(1),'(A/)') '   Fully nonlinear model is employed.'
  ENDIF
  !
  IF(PressureTermOnOff==0)THEN
     WRITE(*,'(A/)') ' No free-surface pressure term is being applied.  '
     WRITE(fileop(1),'(A/)') ' No free-surface pressure term is being applied.  '
  ELSEIF(PressureTermOnOff==1)THEN
     WRITE(*,'(A/)') '   A 2D Gaussian surface pressure is being applied.  (See "funPressureTerm.f90".)'
     WRITE(fileop(1),'(A/)') '   A 2D Gaussian surface pressure is being applied.  (See "funPressureTerm.f90".)'
  ELSEIF(PressureTermOnOff==2)THEN
     WRITE(*,'(A/)') '   A 3D Gaussian surface pressure is being applied.  (See "funPressureTerm.f90".)'
     WRITE(fileop(1),'(A/)') '   A 3D Gaussian surface pressure is being applied.  (See "funPressureTerm.f90".)'
  ELSEIF(PressureTermOnOff==3)THEN
     WRITE(*,'(A/)') '   A 3D tanh surface pressure is being applied.  (See "funPressureTerm.f90".)'
     WRITE(fileop(1),'(A/)') '   A 3D tanh surface pressure is being applied.  (See "funPressureTerm.f90".)'
  ELSE
     Print *, 'No pressure field defined for PressureTermOnOff=',PressureTermOnOff
     stop
  END IF
  !
  IF(TimeMethod /= 1)THEN
     IF(IncWaveType>1 )THEN
        Print *, 'Only RK-4 is implemented with this incident wave-type, TimeMethod -> 1.'
        write(fileop(1),*) 'Only RK-4 is implemented with this incident wave-type, TimeMethod -> 1.'
        TimeMethod=1; RKstages=4;
     END IF
     IF(PressureTermOnOff /= 0)THEN
        Print *, 'Only RK-4 is implemented with an applied free-surface pressure. TimeMethod -> 1.'
        write(fileop(1),*) 'Only RK-4 is implemented with an applied free-surface pressure. TimeMethod -> 1.'
        TimeMethod=1; RKstages=4;
     END IF
  END IF


  ! SG-FILTERING
  READ (FILEIP(1),*) filteringONOFF, filterALPHA, filterORDER, sigma_filt(1), sigma_filt(2), sigma_filt(3)
  IF (filteringONOFF>0) THEN
     WRITE(*,*) '   SG(',2*filterALPHA+1,',',filterORDER,')-filtering will be employed after every ',filteringONOFF,' time step.' 
     WRITE(fileop(1),*) '   SG(',2*filterALPHA+1,',',filterORDER,')-filtering will be employed after every ',filteringONOFF,' time step.' 
     filterNP = filterALPHA*2+1
     ALLOCATE(filtercoefficients(filterNP),tmpfilter(max(filterNP,13)))
     filtercoefficients = zero; tmpfilter = zero
     ! GD: addition filtering on boundaries
     ALLOCATE(filtercoefficients2(max(filterNP,13),max(filterALPHA,6)))
     filtercoefficients2 = zero
     IF (filterALPHA>6)Then
        WRITE(*,*)'** WARNING:  Off-centered filtering coefficients are only implemented for filterALPHA<=6.'
        WRITE(*,*)'**           Some points may be left unfiltered. '
        WRITE(fileop(1),*)'** WARNING:  Off-centered filtering coefficients are only implemented for filterALPHA<=6.'
        WRITE(fileop(1),*)'**           Some points may be left unfiltered. '
     end IF

  ENDIF

  !
  ! BREAKING MODEL
  !
  READ (FILEIP(1),*,err=41) BreakMod%i_breaking, BreakMod%T_half, BreakMod%tan_phi_b, &
       BreakMod%tan_phi_0, BreakMod%del_fac, BreakMod%cel_fac, BreakMod%gamma_break
  If (BreakMod%i_breaking /= 0)Then
     print*, ' '
     Print*, 'Breaking model has been turned on.'
     print *, ' '
     write(fileop(1),*) ' '
     write(fileop(1),*) 'Breaking model has been turned on.'
     write(fileop(1),*) ' '
     IF (BreakMod%i_breaking==2) Then
        BreakMod%i_break_time=2*nsteps
        print *, '  Breaking geometry will be computed, but no model is being applied.'
        write(fileop(1),*) '  Breaking geometry will be computed, but no model is being applied.'
     END IF
  end If
  GoTo 42
41 print*, ' '
  Print*, 'No breaking model line found, the feature is off.'
  print*, ' '
  write(fileop(1),*) ' '
  write(fileop(1),*) 'No breaking model line found, the feature is off.'
  write(fileop(1),*) ' '
  BreakMod%i_breaking=0
  BACKSPACE(fileip(1))
42 continue

  ! RELAXATION ZONES
  READ (FILEIP(1),*) relaxONOFF, relaxTransientTime, relaxNo, relaxXorY, relaxDegrees
  IF (relaxONOFF==1) THEN
     WRITE(*,*) '    Total relaxation zones defined: ',relaxNo
     WRITE(fileop(1),*) '    Total relaxation zones defined: ',relaxNo
     ALLOCATE( RelaxZones(relaxNo) )
     ! 
     DO i=1,relaxNo
        READ (FILEIP(1),*, err=43) RelaxZones(i)%BBox(1), RelaxZones(i)%BBox(2), RelaxZones(i)%BBox(3), &
             RelaxZones(i)%BBox(4), RelaxZones(i)%ftype, RelaxZones(i)%param, RelaxZones(i)%XorY, &
             RelaxZones(i)%WavegenOnOff, RelaxZones(i)%XorYgen, RelaxZones(i)%degrees !, RelaxZones(i)%PhiOnOff
        RelaxZones(i)%PhiOnOff=1
        !        print *, i,'yes',RelaxZones(i)%BBox(1),RelaxZones(i)%PhiOnOff
     END DO
     go to 44
     !hbb
     !hbb  I struggled here to get this read to be backward compatible and finally gave up... 
     !hbb  The feature is implemented, but turned off for now.  Pressure damping on the 
     !hbb  velocity is a much better solution.  
     !hbb
     !443  format(4F10.2,I2,F10.2,A1,I2,A1,F10.2,I2)
     !443     format(4F16.6,I8,F16.6,A1,I8,A1,F16.6,I8)
43   print *, 'Error reading the relaxation zone lines.' !  Note the new format that requires '
     ! print *, 'a value for PhiOnOff=0 (off) or 1 (on) at the end of each zone defn. line.'
     stop
44   continue
  ENDIF
  !
  ! Pressure damping zones
  !
  READ (FILEIP(1),*,err=64) PDampingONOFF,NDampZones 
  IF (NDampZones>1)THEN
     print *, 'Only one pressure damping zone is currently supported.'
     stop
  END IF
  IF (PDampingOnOff /=0) then
     ALLOCATE(PDampZones(NDampZones))
     Do i=1,NDampZones
        READ (FILEIP(1),*,err=63) PDampZones(i)%BBox(1), PDampZones(i)%BBox(2),       &
             PDampZones(i)%BBox(3), PDampZones(i)%BBox(4), PDampZones(i)%g0Phi,         &
             PDampZones(i)%g0Eta, PDampZones(i)%type  
     END Do
     go to 65
  END IF
  go to 65
652 format(2I8)
653 format(7F16.6,I8)
63 print *, 'ReadInputFileParameters:  Error reading pressure damping zones line.'
  stop
64 continue
  print *, ' '
  print *, 'No Pressure damping line found, the feature is off.'
  write(fileop(1),*) ' '
  write(fileop(1),*) 'No Pressure damping line found, the feature is off.'
  backspace(FILEIP(1))
65 continue
  print *, ' '
  print *, 'Found ', NDampZones, ' pressure damping zones.'
  write(fileop(1),*) ' '
  write(fileop(1),*) 'Found ', NDampZones, ' pressure damping zones.'

  ! SWENSE line
  READ(FILEIP(1),*) swenseONOFF, swenseTransientTime, swenseDir, West_refl, East_refl, North_refl, South_refl

  ! CURVILINEAR
  READ(FILEIP(1),*) curvilinearONOFF
  IF (curvilinearONOFF==1 .AND. FineGrid%Nx>1 .AND. FineGrid%Ny>1) THEN
     ! NOTE: It is only possible to run the curvilinear model for 3D cases, i.e. 2D cases excluded
     WRITE(*,'(A/)') '   Curvilinear model employed.'
     WRITE(fileop(1),'(A/)') '   Curvilinear model employed.'
  ELSE
     WRITE(*,'(A/)') '   Standard Cartesian model employed.'
     WRITE(fileop(1),'(A/)') '   Standard Cartesian model employed.'
     curvilinearONOFF = 0    ! Make sure it is the standard model which is employed (also for curvilinear 2D choices)
  END IF
!
! Linear mono-chromatic, random wave or flux boundary generation parameters.  
! 

  IF (IncWaveType==3) THEN 
     ! Wave generation with flux condition on western boundary
     READ(FILEIP(1),*,IOSTAT=ios) wave3DFlux%rampTime, wave3DFlux%order, wave3DFlux%inc_wave_file
     IF (ios>0) THEN
        Print *, 'ReadInputFileParameters:  For IncWaveType==3 we need: Ramp &
             time, Order of interpolation in horizontal direction, file with wave paddle signal.'
        stop
     END IF

  ELSEIF (IncWaveType == 2) THEN ! irregular waves
     ! For irregular waves we have four options: 
     ! 0) PM, 
     ! 1) Normal JONSWAP with gamma = 3.3, 
     ! 2) based on input files and 
     ! 3) JONSWAP with variable gamma value
     !
     ! To have a clean interface without too many goto statement,we first check what type we are trying to read

     READ(FILEIP(1),*,IOSTAT=ios) ispec
     Backspace(FILEIP(1)) 
     IF (ispec==0 .or. ispec==1) THEN !Normal PM or JONSWAP spectrum
        READ(FILEIP(1),*,IOSTAT=ios) ispec,  Tp,  Hs,  h0,   &
             kh_max,  seed,  seed2,  x0,  y0
        gamma_jonswap = 3.3

     ELSEIF (ispec == 2) THEN ! 2D irregular waves with input file
        READ(FILEIP(1),*,IOSTAT=ios) ispec,  Tp,  Hs,  h0,   & 
             kh_max,  seed,  seed2,  x0,  y0, &
             inc_wave_file
        print *, x0, y0, inc_wave_file
     ELSEIF (ispec == 3) THEN ! 2D irregular waves with non-standard gamma value
        READ(FILEIP(1),*,IOSTAT=ios) ispec,  Tp,  Hs,  h0,   &
             kh_max,  seed,  seed2,  x0,  y0, gamma_jonswap

     ELSEIF (ispec>=30) THEN ! multi-directional irregular waves
        READ(FILEIP(1),*,ERR=37,END=37,IOSTAT=ios) ispec,  Tp,  Hs,  h0,   & 
             kh_max,  seed,  seed2,  x0,  y0, &
             inc_wave_file, beta0, s0, gamma_jonswap

     END IF

     IF( abs(ispec)<30 ) THEN
        beta0=0
        s0=1.0
     END IF

37   IF(ios>0)THEN
        IF (abs(ispec)<30) THEN
           Print *, 'ReadInputFileParameters:  For IncWaveType==2 we need irregular wave parameters.'
        ELSE
           Print *, 'ReadInputFileParameters:  For 3D waves, abs(ispec)>30, we need a heading angle, a spreading factor and a JONSWAP gamma factor.'
        END IF
        STOP
     END IF
     !
     ! Check that the seed values are negative and if not make them negative
     !
     IF (seed >=0) THEN
        seed = - seed
        write(6,38)seed
        write(fileop(1),38)seed
38      format(/,'**WARNING**: Seed values must be negative, changed seed to ',i10,/)
     end IF
     IF (seed2 >=0) THEN
        seed2 = - seed2
        write(6,39)seed2
        write(fileop(1),39)seed2
39      format(/,'**WARNING**: Seed values must be negative, changed seed2 to ',i10,/)
     end IF
  END IF

  !
  !
  IF (relaxONOFF==1) THEN ! GD: add this test in case of no relaxation zones defined
     Allocate(RandomWave(relaxNo))
     !
     ! For linear wave generation, allocate a RandomWave structure for each zone, so that 
     ! each zone has a RandomWave associated with it.  
     !
     IF (IncWaveType==2) THEN
        ! Put the scalar wave parameters into each structure and count the ones that will 
        ! actually be used.  
        nGenZones=0
        Do i=1,relaxNo
           RandomWave(i)%ispec=ispec; RandomWave(i)%Tp=Tp; RandomWave(i)%Hs=Hs;
           RandomWave(i)%h0=h0; RandomWave(i)%x0=x0; RandomWave(i)%y0=y0;
           RandomWave(i)%seed=seed; RandomWave(i)%seed2=seed2; RandomWave(i)%kh_max=kh_max;
           RandomWave(i)%inc_wave_file=inc_wave_file; RandomWave(i)%beta0=beta0; 
           RandomWave(i)%S0=s0
           RandomWave(i)%gamma = gamma_jonswap
           If(RelaxZones(i)%XorYgen=='X' .AND. RelaxZones(i)%WavegenONOFF==1) THEN
              nGenZones=nGenZones+1
              If( abs(RandomWave(i)%ispec) >= 30 .and. RelaxZones(i)%degrees /= 0) THEN
                 print *, '3D linear wave generation is only supported for RelaxZones%degrees=0'
                 stop
              END If
           END If
        END Do
        Print *, '  Found ',nGenZones, &
             ' generation zones for the linear wave. (Only used for IncWaveType==2)'
        Print *, ' '
        write(fileop(1),*) '  Found ',nGenZones, &
             ' generation zones for the linear wave. (Only used for IncWaveType==2)'
        write(fileop(1),*) ' '
     ENDIF
  ENDIF

  RETURN

  ! ERROR CONTROL MESSAGES

100 WRITE (*,'(A)') 'Error: Cannot read file header of input file.'; STOP
101 WRITE (*,'(A)') 'Error: Stencil for x-direction too large (rank>Nx) or too small (rank=0).'; STOP
102 WRITE (*,'(A)') 'Error: Stencil for y-direction too large (rank>Ny) or too small (rank=0).'; STOP
103 WRITE (*,'(A)') 'Error: Stencil for z-direction too large (rank>Nz) or too small (rank=0).'; STOP
104 WRITE (*,'(A)') 'Error: Multigrid solver can only implemented when ghost layer is invoked.'; STOP
105 WRITE (*,'(A)') 'Error: Cannot read all required parameters of line with Lx, Ly, Lz, ...'; STOP

  ! REUSABLE OUTPUT FORMATS
900 FORMAT (A,I5,A,I5,A,I5,A)
901 FORMAT (A,I2,A,I2,A,I2,A)
902 FORMAT (A,I8)
903 FORMAT (A,E9.4)

END SUBROUTINE ReadInputFileParameters
