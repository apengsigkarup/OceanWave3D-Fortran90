SUBROUTINE SetupInitialConditions
  ! By Allan P. Engsig-Karup.
  USE GlobalVariables
  IMPLICIT NONE
  INTEGER :: i, j, k, Nx, Ny, Nz
  REAL(KIND=long) :: kh, kw, HH, kh_deep, kh_shallow, h_deep, h_shallow, kx, ky, y0
  REAL(KIND=long) :: tmpx(FineGrid%Nx+2*GhostGridX), tmpy(FineGrid%Ny+2*GhostGridY)
  EXTERNAL BeachGen3Dver2, BuildLinearSystem
  ! TEMPORARY POINTERS/VARIABLES
  Nx = FineGrid%Nx
  Ny = FineGrid%Ny
  Nz = FineGrid%Nz
  DetermineBottomGradients = 0
  print*,'IC chosen is case ',IC
  write(fileop(1),*)'IC chosen is case ',IC
  SELECT CASE (IC)

  CASE (-1, 0) ! Initial condition determined by PressureTermOnOff and funPressureTerm.f90 or read from the init file 'OceanWave3D.init'.  
     IF (IC==-1)THEN
        print *, 'Reading the initial conditions. ** Not implemented curvilinear! **'
        write(fileop(1),*) 'Reading the initial conditions. ** Not implemented curvilinear! **'
        
        DO j=1+GhostGridY,FineGrid%Ny+GhostGridY
           Do i=1+GhostGridX,FineGrid%Nx+GhostGridX
              read(fileip(3),*)WaveField%E(i,j),WaveField%P(i,j)
           END Do
        END DO
        
     ELSE
        IF(PressureTermOnOff==0)THEN
           print*, 'Initial condition is still water.'
           write(fileop(1),*) 'Initial condition is still water.'
        ELSEIF(PressureTermOnOff==1)THEN
           print *, 'Initial condition determined by PressureTermOnOff is a 2D stationary Gaussian hump'
           write(fileop(1),*) 'Initial condition determined by PressureTermOnOff is a 2D stationary Gaussian hump'
        ELSEIF(PressureTermOnOff==2)THEN
           print *, 'Initial condition determined by PressureTermOnOff is a 3D moving Gaussian hump'
           write(fileop(1),*) 'Initial condition determined by PressureTermOnOff is a 3D moving Gaussian hump'
        END IF
        Call funInitialFreeSurfaceElevation(g,FineGrid%Nx+2*GhostGridX,&
             FineGrid%Ny+2*GhostGridY,FineGrid,WaveField)
     END IF
     IF (Lz <= 0)THEN
        ! Read in the bottom contours from a file.
        OPEN(unit=fileip(4),file=fname_bottom,status='old')
        READ(fileip(4),'(A)')head(4)
        PRINT *, 'SetUpInitialConditions:  Reading the bottom contours from file:'
        PRINT *,fname_bottom,' with header:'
        PRINT *, head(4)
        PRINT *, ' '
        WRITE(FILEOP(1),*) 'SetUpInitialConditions:  Reading the bottom contours from file:'
        WRITE(FILEOP(1),*)fname_bottom,' with header:'
        WRITE(FILEOP(1),*) head(4)
        WRITE(FILEOP(1),*) ' '
        READ(fileip(4),*)  DetermineBottomGradients ! Read the gradient  flag

        IF(DetermineBottomGradients==1)THEN
           Print *, ' Reading h and computing derivatives of h numerically.'
           write(fileop(1),*) ' Reading h and computing derivatives of h numerically.'
           ! Read in h(x,y)
           do j=1+GhostGridY,ny+GhostGridY
              DO i=1+GhostGridX,nx+GhostGridX
                 READ(fileip(4),*) FineGrid%h(i,j)
              END DO
           END do
        ELSEIF(DetermineBottomGradients == 0) THEN
           ! Read h(x,y),h_x,h_xx,h_y,h_yy
           Print *, ' Reading h,h_x,h_xx,h_y,h_yy.'
           print *, ' '
           write(fileop(1),*) ' Reading h,h_x,h_xx,h_y,h_yy.'
           write(fileop(1),*) ' '
           !
           Do j=1+GhostGridY,ny+GhostGridY
              DO i=1+GhostGridX,nx+GhostGridX
                 READ(fileip(4),*) FineGrid%h(i,j), FineGrid%hx(i,j), FineGrid%hxx(i,j), &
                      FineGrid%hy(i,j), FineGrid%hyy(i,j)
              END DO
           END Do
        END IF
        !
        ! Extend the bathymetry to the ghost points by copying the domain endpoint values.
        !
        FineGrid%h(1,:)=FineGrid%h(1+GhostGridX,:); FineGrid%h(:,1)=FineGrid%h(:,1+GhostGridY);
        FineGrid%h(nx+2*GhostGridX,:)=FineGrid%h(nx+GhostGridX,:);
        FineGrid%h(:,ny+2*GhostGridY)=FineGrid%h(:,ny+GhostGridY)
        FineGrid%hx(1,:)=FineGrid%hx(1+GhostGridX,:); FineGrid%hx(:,1)=FineGrid%hx(:,1+GhostGridY);
        FineGrid%hx(nx+2*GhostGridX,:)=FineGrid%hx(nx+GhostGridX,:);
        FineGrid%hx(:,ny+2*GhostGridY)=FineGrid%hx(:,ny+GhostGridY)
        FineGrid%hxx(1,:)=FineGrid%hxx(1+GhostGridX,:);
        FineGrid%hxx(:,1)=FineGrid%hxx(:,1+GhostGridY);
        FineGrid%hxx(nx+2*GhostGridX,:)=FineGrid%hxx(nx+GhostGridX,:);
        FineGrid%hxx(:,ny+2*GhostGridY)=FineGrid%hxx(:,ny+GhostGridY)
        FineGrid%hy(1,:)=FineGrid%hy(1+GhostGridX,:); FineGrid%hy(:,1)=FineGrid%hy(:,1+GhostGridY);
        FineGrid%hy(nx+2*GhostGridX,:)=FineGrid%hy(nx+GhostGridX,:);
        FineGrid%hy(:,ny+2*GhostGridY)=FineGrid%hy(:,ny+GhostGridY)
        FineGrid%hyy(1,:)=FineGrid%hyy(1+GhostGridX,:);
        FineGrid%hyy(:,1)=FineGrid%hyy(:,1+GhostGridY);
        FineGrid%hyy(nx+2*GhostGridX,:)=FineGrid%hyy(nx+GhostGridX,:);
        FineGrid%hyy(:,ny+2*GhostGridY)=FineGrid%hyy(:,ny+GhostGridY)
        !
        ! Re-set Lz to be max h(1,1)
        !
        Lz=FineGrid%h(1,1)
     ELSE
        FineGrid%h = Lz; FineGrid%hx=zero; FineGrid%hxx=zero; FineGrid%hy=zero; FineGrid%hyy=zero;
     END IF

  CASE (1) ! Mildly nonlinear standing wave, deep water (Agnon & Glozman (1996))
     print *, 'Mildly nonlinear standing wave (deep water) in a rectangular domain.'
     write(fileop(1),*) 'Mildly nonlinear standing wave (deep water) in a rectangular domain.'
     IF (Nx>1) THEN
        FineGrid%h = two; FineGrid%hx=zero; FineGrid%hxx=zero
        CALL nonlinearstandingwave1D(pi,FineGrid%h(1,1),FineGrid%x,tmp2D,Wavefield%W,Wavefield%E, &
             Wavefield%Ex,Wavefield%Exx,&
             (Nx+2*GhostGridX)*(Ny+2*GhostGridY))
     ELSE IF (Ny>1) THEN
        FineGrid%h = two; FineGrid%hy=zero; FineGrid%hyy=zero
        CALL nonlinearstandingwave1D(pi,FineGrid%h(1,1),FineGrid%y,tmp2D,Wavefield%W,Wavefield%E, &
             Wavefield%Ey,Wavefield%Eyy,&
             (Nx+2*GhostGridX)*(Ny+2*GhostGridY))
     ENDIF
     !GD: test in a curvilinear rotated domain...
     IF (curvilinearONOFF == 1) THEN
        !Rotate the physical grid with a given angle...
        tmp2D = FineGrid%x
        FineGrid%x = tmp2D*COS(20*PI/180)-FineGrid%y*SIN(20*PI/180)
        FineGrid%y = tmp2D*SIN(20*PI/180)+FineGrid%y*COS(20*PI/180)
     ENDIF
     !
     !
  CASE (2) ! Shallow water to Deep water (3D)
     kh_deep    = pi
     kh_shallow = half
     kx         = two*pi/Lx
     ky         = two*pi/Ly
     kw         = SQRT(kx**2 + ky**2)
     h_deep     = kh_deep/kw
     h_shallow  = kh_shallow/kw
     CALL BeachGen3Dver2(FineGrid%h,h_deep,h_shallow,Lx,Ly,FineGrid%x,FineGrid%y,Nx,Ny)
     HH   = 0.4_long*h_shallow
     Wavefield%E   =  HH*COS(kx*FineGrid%x)*COS(ky*FineGrid%y);
     Wavefield%Ex  = -HH*(kx)*SIN(kx*FineGrid%x)*COS(ky*FineGrid%y);
     Wavefield%Exx = -HH*(kx)**2*COS(kx*FineGrid%x)*COS(ky*FineGrid%y);
     Wavefield%Ey  = -HH*(ky)*COS(kx*FineGrid%x)*SIN(ky*FineGrid%y);
     Wavefield%Eyy = -HH*(ky)**2*COS(kx*FineGrid%x)*COS(ky*FineGrid%y);
     DetermineBottomGradients = 1
  CASE (3) ! Whalin (3D)
     IF (FineGrid%Nx>FineGrid%Ny) THEN
        y0 = (FineGrid%y(1,2)-FineGrid%y(1,1))/two ! FIXME: x-index
     ELSE
        y0 = (FineGrid%x(2,1)-FineGrid%x(1,1))/two ! FIXME: y-index
     ENDIF
     ! Define bottom
     CALL BottomWhalin(FineGrid,y0,GhostGridX, GhostGridY)
     !
     IF(swenseONOFF==1)THEN
        IF (LinearONOFF==0) THEN
           print*,'it should be a nonlinear case!'
           STOP
           ! Linear incident wavefield
           CALL incident_linear_wf_finite(swenseDir, Wavefield, & !ramp_type,&
                Nx+2*GhostGridX,FineGrid%x,Ny+2*GhostGridY,FineGrid%y,Nz+GhostGridz,FineGrid%z,FineGrid%h,zero,&
                SFsol%k,g,SFsol%HH,SFsol%h)
        ELSE
           ! Nonlinear incident wavefield
           ! Initialize stream function coefficients
           CALL stream_func_set_up(g,SFsol%h,SFsol%T,SFsol%i_wavel_or_per,SFsol%L,  &
                SFsol%k,SFsol%HH,SFsol%i_euler_or_stokes,SFsol%e_or_s_vel,SFsol%i_deep_or_finite,  &
                SFsol%n_h_steps,SFsol%n_four_modes,SFsol%nwrk,SFsol%yy,SFsol%zz)
           CALL incident_wf_finite(swenseDir, Wavefield, &
                Nx+2*GhostGridX,FineGrid%x,Ny+2*GhostGridY,FineGrid%y,Nz+GhostGridz,FineGrid%z,FineGrid%h,zero, &
                SFsol%k,g,SFsol%n_four_modes,SFsol%zz,SFsol%yy)
        ENDIF
     ENDIF
  CASE (4) ! Deep water, flat bottom, linear standing wave
     IF (Nx>0) THEN
        FineGrid%hx    = zero
        FineGrid%hxx   = zero
     ENDIF
     IF (Ny>0) THEN
        FineGrid%hy    = zero
        FineGrid%hyy   = zero
     ENDIF
     IF (Nx>1) THEN
        kh = two*pi
        kw = two*pi/Lx
        FineGrid%h  = kh/kw
        HH = Lx/ten ! initial wave height
        Wavefield%E   = HH/two*COS(kw*FineGrid%x)
        Wavefield%Ex  = -kw*HH/two*SIN(kw*FineGrid%x)
        Wavefield%Exx = -kw**two*HH/two*COS(kw*FineGrid%x)
     ELSE
        kh = two*pi
        kw = two*pi/Ly
        FineGrid%h  = kh/kw
        HH = Ly/ten ! initial wave height
        Wavefield%E   = HH/two*COS(kw*FineGrid%y)
        Wavefield%Ey  = -kw*HH/two*SIN(kw*FineGrid%y)
        Wavefield%Eyy = -kw**two*HH/two*COS(kw*FineGrid%y)
     ENDIF
  CASE (5) ! Shallow water, flat bottom
     IF (Nx>1) THEN
        FineGrid%hx    = zero
        FineGrid%hxx   = zero
     ENDIF
     IF (Ny>1) THEN
        FineGrid%hy    = zero
        FineGrid%hyy   = zero
     ENDIF
     IF (Nx>1) THEN
        kh = half
        kw = two*pi/Lx
        FineGrid%h  = kh/kw
        HH = FineGrid%h(1,1)*0.75_long ! initial wave height
        Wavefield%E   = HH/two*COS(kw*FineGrid%x)
        Wavefield%Ex  = -kw*HH/two*SIN(kw*FineGrid%x)
        Wavefield%Exx = -kw**two*HH/two*COS(kw*FineGrid%x)
     ELSE
        kh = half
        kw = two*pi/Ly
        FineGrid%h  = kh/kw
        HH = FineGrid%h(1,1)*0.75_long ! initial wave height
        Wavefield%E   = HH/two*COS(kw*FineGrid%y)
        Wavefield%Ey  = -kw*HH/two*SIN(kw*FineGrid%y)
        Wavefield%Eyy = -kw**2*HH/two*COS(kw*FineGrid%y)
     ENDIF
  CASE (6) ! Shallow to Deep Water
     !$$$$$$         kh_deep    = pi
     !$$$$$$         kh_shallow = half
     !$$$$$$         kw         = two*pi/Lx
     !$$$$$$         h_deep     = kh_deep/kw
     !$$$$$$         h_shallow  = kh_shallow/kw
     !$$$$$$         HH   = 0.75_long*h_shallow
     !$$$$$$         IF (Nx>1) THEN
     !$$$$$$                 CALL BeachGen2(h_deep,h_shallow,Lx,FineGrid%x,FineGrid%h,FineGrid%hx,FineGrid%hxx,FineGrid%Nx)
     !$$$$$$                 Wavefield%E   = HH/two*COS(kw*FineGrid%x)
     !$$$$$$                 Wavefield%Ex  = -kw*HH/two*SIN(kw*FineGrid%x)
     !$$$$$$                 Wavefield%Exx = -kw**2*HH/two*COS(kw*FineGrid%x)
     !$$$$$$         ELSE
     !$$$$$$                 CALL BeachGen2(h_deep,h_shallow,Lx,FineGrid%y,FineGrid%h,FineGrid%hy,FineGrid%hyy,FineGrid%Ny)
     !$$$$$$                 Wavefield%E   = HH/two*COS(kw*FineGrid%y)
     !$$$$$$                 Wavefield%Ey  = -kw*HH/two*SIN(kw*FineGrid%y)
     !$$$$$$                 Wavefield%Eyy = -kw**2*HH/two*COS(kw*FineGrid%y)
     !$$$$$$         ENDIF
     ! GD: one has to define wave generation and absorption for this test case ?
     kh_deep    = pi
     !kh_shallow = half
     kh_shallow = pi/20.d0
     !kw         = two*pi/Lx
     kw         = two*pi/SFsol%L
     h_deep     = kh_deep/kw
     h_shallow  = kh_shallow/kw
     IF (SFsol%HH.GT.0.75_long*h_shallow) THEN
        print*, 'reduce amplitude of incident wave...'
        stop
     ENDIF
     HH   = SFsol%HH !0.75_long*h_shallow
     IF (Nx>1) THEN
        ! use relaxation indexes to define correctly ...
        ! First zone: generation on the left
        !RelaxZones(1)%idx(1) !beginnning first zone (will be first point in this case)
        !RelaxZones(1)%idx(2) !end first zone (beginning of the slope then)
        !Second zone:absorption on the right
        !RelaxZones(2)%idx(1) !beginnning first zone (end of the slope)
        !RelaxZones(2)%idx(2) !end first zone (last point)
        !
        CALL BeachGen2_bis(h_deep,h_shallow,RelaxZones(1)%idx(2),RelaxZones(2)%idx(1),FineGrid%x,  &
             FineGrid%h,FineGrid%hx,FineGrid%hxx,FineGrid%Nx+2*GhostGridX)
        ! No intialisation one generates wave
        !Wavefield%E   = HH/two*COS(kw*FineGrid%x)
        !Wavefield%Ex  = -kw*HH/two*SIN(kw*FineGrid%x)
        !Wavefield%Exx = -kw**2*HH/two*COS(kw*FineGrid%x)
     ELSE
        ! use relaxation indexes to define correctly ...
        ! First zone: generation on the left
        !RelaxZones(1)%idx(3) !beginnning first zone (will be first point in this case)
        !RelaxZones(1)%idx(4) !end first zone (beginning of the slope then)
        !Second zone:absorption on the right
        !RelaxZones(2)%idx(3) !beginnning first zone (end of the slope)
        !RelaxZones(2)%idx(4) !end first zone (last point)
        CALL BeachGen2_bis(h_deep,h_shallow,RelaxZones(1)%idx(4),RelaxZones(2)%idx(3),FineGrid%y,  &
             FineGrid%h,FineGrid%hy,FineGrid%hyy,FineGrid%Ny+2*GhostGridZ)
        ! No intialisation one generates wave
        !Wavefield%E   = HH/two*COS(kw*FineGrid%y)
        !Wavefield%Ey  = -kw*HH/two*SIN(kw*FineGrid%y)
        !Wavefield%Eyy = -kw**2*HH/two*COS(kw*FineGrid%y)
     ENDIF
     ! Relaxation
     IF (relaxONOFF==1) THEN
        CALL RelaxationModule(Wavefield%E,Wavefield%P,0.d0)
     ENDIF
     print *, 'initialisation of Linear Shoaling'
     write(fileop(1),*) 'initialisation of Linear Shoaling'

  CASE (7) ! Flat bottom, depth defined from SF-wave
     IF (Nx>1) THEN
        FineGrid%hx    = zero
        FineGrid%hxx   = zero
     ENDIF
     IF (Ny>1) THEN
        FineGrid%hy    = zero
        FineGrid%hyy   = zero
     ENDIF
     FineGrid%h         = SFsol%h
     !
     IF (curvilinearONOFF==1) THEN ! curvilinear coordinates... FIXME for make it work whatever transformation
        ! In the curvilinear case, BBox have to be the indexes of the corresponding relaxation zones
        ! FIXME: use Lx and Ly? Depends how x and y are defined...
        DO j=1,FineGrid%Nx+2*GhostGridX
           tmpx(j) = REAL(j-2,long)*(FineGrid%y(2,1)-FineGrid%y(1,1)) !/REAL((FineGrid%Nx+2*GhostGridX-1),long)*(8+PI)*Lx; !FIXME find a way to define in all situation...
        ENDDO
        DO j=1,FineGrid%Ny+2*GhostGridY
           tmpy(j) = -REAL(j-2,long)*(FineGrid%x(1,2)-FineGrid%x(1,1)) !/REAL((FineGrid%Ny+2*GhostGridY-1),long)*(Ly-Lx); !FIXME find a way to define in all situation...
        ENDDO
     ELSE
        tmpx(:) = FineGrid%x(:,1)
        tmpy(:) = FineGrid%y(1,:)
     ENDIF
     IF (relaxONOFF==1) THEN
        IF (relaxXorY=='X' .OR. relaxXorY=='x') THEN
           IF (relaxDegrees==zero) THEN
              CALL stream_func_wave_finite(FineGrid%Nx+2*GhostGridX,tmpx,time,SFsol%n_four_modes,SFsol%zz,SFsol%yy,&
                   SFsol%k,g,Wavefield%E(:,1),Wavefield%P(:,1))
              DO j=1,FineGrid%Ny+2*GhostGridY
                 Wavefield%E(:,j) = Wavefield%E(:,1)
                 Wavefield%P(:,j) = Wavefield%P(:,1)
              END DO
           ELSE
              DO j=1,FineGrid%Ny+2*GhostGridY
                 CALL stream_func_wave_finite(FineGrid%Nx+2*GhostGridX,tmpx*COS(relaxDegrees/180.0_long*pi)+&
                      tmpy(j)*SIN(relaxDegrees/180.0_long*pi),time,SFsol%n_four_modes,SFsol%zz,&
                      SFsol%yy,SFsol%k,g,Wavefield%E(:,1),Wavefield%P(:,1))
                 Wavefield%E(:,j) = Wavefield%E(:,1)
                 Wavefield%P(:,j) = Wavefield%P(:,1)
              END DO
           ENDIF
        ELSE
           IF (relaxDegrees==zero) THEN
              CALL stream_func_wave_finite(FineGrid%Ny+2*GhostGridY,tmpy,time,SFsol%n_four_modes,SFsol%zz,SFsol%yy,&
                   SFsol%k,g,Wavefield%E(1,:),Wavefield%P(1,:))
              DO j=1,FineGrid%Nx+2*GhostGridX
                 Wavefield%E(j,:) = Wavefield%E(1,:)
                 Wavefield%P(j,:) = Wavefield%P(1,:)
              END DO
           ELSE
              DO j=1,FineGrid%Nx+2*GhostGridX
                 CALL stream_func_wave_finite(FineGrid%Ny+2*GhostGridY,-tmpy*SIN(relaxDegrees/180.0_long*pi)+&
                      tmpx(j)*COS(relaxDegrees/180.0_long*pi),time,SFsol%n_four_modes,SFsol%zz,&
                      SFsol%yy,SFsol%k,g,Wavefield%E(1,:),Wavefield%P(1,:))
                 Wavefield%E(j,:) = Wavefield%E(1,:)
                 Wavefield%P(j,:) = Wavefield%P(1,:)
              END DO
           ENDIF
        ENDIF
        ! Relax initial condition three times for full absorption
        CALL RelaxationModule(Wavefield%E,Wavefield%P,time)
        CALL RelaxationModule(Wavefield%E,Wavefield%P,time)
        CALL RelaxationModule(Wavefield%E,Wavefield%P,time)
     ENDIF
     IF(curvilinearONOFF==1) THEN
        ! FIXME: to add
     ELSE
        CALL PreProcessDiffStencils(FineGrid,FineGrid%DiffStencils,GhostGridX,GhostGridY,GhostGridZ,  &
             alpha,beta,gamma,fileop(1))
        IF (FineGrid%Nx>1) THEN
           CALL DiffXEven(Wavefield%E,Wavefield%Ex,1,FineGrid%Nx,FineGrid%Ny,1,FineGrid%DiffStencils,alpha)
           CALL DiffXEven(Wavefield%E,Wavefield%Exx,2,FineGrid%Nx,FineGrid%Ny,1,FineGrid%DiffStencils,alpha)
           CALL DiffXEven(Wavefield%P,Wavefield%Px,1,FineGrid%Nx,FineGrid%Ny,1,FineGrid%DiffStencils,alpha)
        END IF
        IF (FineGrid%Ny>1) THEN
           CALL DiffYEven(Wavefield%E,Wavefield%Ey,1,FineGrid%Nx,FineGrid%Ny,1,FineGrid%DiffStencils,beta)
           CALL DiffYEven(Wavefield%E,Wavefield%Eyy,2,FineGrid%Nx,FineGrid%Ny,1,FineGrid%DiffStencils, &
                beta)
           CALL DiffYEven(Wavefield%P,Wavefield%Py,1,FineGrid%Nx,FineGrid%Ny,1,FineGrid%DiffStencils, &
                beta)
        END IF
     ENDIF
  CASE (8) ! Linear standing wave (2D) + in curvilinear space rotation with an angle
     FineGrid%h = SFsol%h
     IF (Nx>1) THEN
        FineGrid%hx=zero; FineGrid%hxx=zero
     ENDIF
     IF (Ny>1) THEN
        FineGrid%hy=zero; FineGrid%hyy=zero
     ENDIF
     ! FIXME: Does not work for problems with only one spatial dimension in the horizontal /APEK
     CALL linearstandingwave2D(g,zero,zero,SFsol,SFsol%L,SFsol%L,FineGrid%h(1,1),FineGrid%x,FineGrid%y&
          ,Wavefield%P,tmp2D,tmp2D,Wavefield%W,Wavefield%E,Wavefield%Ex,Wavefield%Exx,Wavefield%Ey, &
          Wavefield%Eyy,(Nx+2*GhostGridX)*(Ny+2*GhostGridY))
     IF (curvilinearONOFF == 1) THEN
        !Rotate the physical grid with a given angle...
        tmp2D = FineGrid%x
        FineGrid%x = tmp2D*COS(20*PI/180)-FineGrid%y*SIN(20*PI/180)
        FineGrid%y = tmp2D*SIN(20*PI/180)+FineGrid%y*COS(20*PI/180)
     ENDIF
     print*,'  Linear standing wave parameters:'
     print*,'      kh = ',SFsol%k*SFsol%h
     print*,'      T  = ',SFsol%T
     print*,'      c  = ',SFsol%c
     print*,'      h  = ',SFsol%h
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
     PRINT *,'      dt = ',dt
     PRINT *,'      Cr = ',SFsol%c*dt/MIN(dxmin,dymin)	
     WRITE(FILEOP(1),*)'      dt = ',dt
     WRITE(FILEOP(1),*)'      Cr = ',SFsol%c*dt/MIN(dxmin,dymin)	
  CASE (9) ! Berkhoff (3D)
     CALL BottomBerkoff(FineGrid,GhostGridX,GhostGridY)
     DetermineBottomGradients = 1
  CASE (10) ! GD: SWENSE, Flat bottom, depth defined from SF-wave
     IF (curvilinearONOFF == 1) THEN
        !Rotate the physical grid with a given angle...
        tmp2D = FineGrid%x
        FineGrid%x = tmp2D*COS(20.d0*PI/180.d0)-FineGrid%y*SIN(20.d0*PI/180.d0)
        FineGrid%y = tmp2D*SIN(20.d0*PI/180.d0)+FineGrid%y*COS(20.d0*PI/180.d0)
        swenseDir = 20.d0
     ENDIF
     IF (Nx>1) THEN
        FineGrid%hx    = zero
        FineGrid%hxx   = zero
     ENDIF
     IF (Ny>1) THEN
        FineGrid%hy    = zero
        FineGrid%hyy   = zero
     ENDIF
     FineGrid%h         = SFsol%h
     ! E,P and their derivatives are zero
     !
     ! Initialization of incident wavefield
     PRINT*,'Beginning of SWENSE Initialization'
     IF(swenseONOFF==0)THEN
        PRINT*,'We should use SWENSE... put swenseONOFF=1 !'
        !STOP
     ELSE
        IF (LinearONOFF==0) THEN
           ! Linear incident wavefield
           CALL incident_linear_wf_finite(swenseDir, Wavefield, & !ramp_type,&
                Nx+2*GhostGridX,FineGrid%x,Ny+2*GhostGridY,FineGrid%y,Nz+GhostGridz,FineGrid%z,  &
                FineGrid%h,zero,SFsol%k,g,SFsol%HH,SFsol%h)
        ELSE
           ! Nonlinear incident wavefield
           print*,' this case is not available already... several checks to make...'
           ! Initialize stream function coefficients
           CALL stream_func_set_up(g,SFsol%h,SFsol%T,SFsol%i_wavel_or_per,SFsol%L,  &
                SFsol%k,SFsol%HH,SFsol%i_euler_or_stokes,SFsol%e_or_s_vel,SFsol%i_deep_or_finite,  &
                SFsol%n_h_steps,SFsol%n_four_modes,SFsol%nwrk,SFsol%yy,SFsol%zz)
           CALL incident_wf_finite(swenseDir, Wavefield, &
                Nx+2*GhostGridX,FineGrid%x,Ny+2*GhostGridY,FineGrid%y,Nz+GhostGridz,FineGrid%z,  &
                FineGrid%h,zero, &
                SFsol%k,g,SFsol%n_four_modes,SFsol%zz,SFsol%yy)
        ENDIF
        PRINT*,'End of SWENSE Initialization'
     ENDIF
  CASE (11) ! Flat bottom, depth defined from SF-wave for semi circular channel...
     IF (Nx>1) THEN
        FineGrid%hx    = zero
        FineGrid%hxx   = zero
     ENDIF
     IF (Ny>1) THEN
        FineGrid%hy    = zero
        FineGrid%hyy   = zero
     ENDIF
     FineGrid%h         = SFsol%h
     !
     IF (curvilinearONOFF==1) THEN ! curvilinear coordinates... FIXME for make it work whatever transformation
        ! In the curvilinear case, BBox have to be the indexes of the corresponding relaxation zones
        ! FIXME: use Lx and Ly? Depends how x and y are defined...  save e and n directly ?
        DO j=1,FineGrid%Nx+2*GhostGridX
           tmpx(j) = REAL(j-2,long)*(FineGrid%y(2,1)-FineGrid%y(1,1)) !/REAL((FineGrid%Nx+2*GhostGridX-1),long)*(8+PI)*Lx; !FIXME find a way to define in all situation...
        ENDDO
        DO j=1,FineGrid%Ny+2*GhostGridY
           tmpy(j) = -REAL(j-2,long)*(FineGrid%x(1,2)-FineGrid%x(1,1)) !/REAL((FineGrid%Ny+2*GhostGridY-1),long)*(Ly-Lx); !FIXME find a way to define in all situation...
        ENDDO
     ELSE
        tmpx(:) = FineGrid%x(:,1)
        tmpy(:) = FineGrid%y(1,:)
     ENDIF
     IF (relaxONOFF==1) THEN
        ! Relax initial condition
        CALL RelaxationModule(Wavefield%E,Wavefield%P,time)
        CALL RelaxationModule(Wavefield%E,Wavefield%P,time)
        CALL RelaxationModule(Wavefield%E,Wavefield%P,time)
     ENDIF
     IF(swenseONOFF==0)THEN
        ! Nothing to do here
     ELSE
        IF (LinearONOFF==0) THEN
           ! Linear incident wavefield
           CALL incident_linear_wf_finite(swenseDir, Wavefield, & !ramp_type,&
                Nx+2*GhostGridX,FineGrid%x,Ny+2*GhostGridY,FineGrid%y,Nz+GhostGridz,FineGrid%z,  &
                FineGrid%h,zero,&
                SFsol%k,g,SFsol%HH,SFsol%h)
        ELSE
           ! Nonlinear incident wavefield
           print*,' this case is not available already... several checks to make...'
           ! Initialize stream function coefficients
           CALL stream_func_set_up(g,SFsol%h,SFsol%T,SFsol%i_wavel_or_per,SFsol%L,  &
                SFsol%k,SFsol%HH,SFsol%i_euler_or_stokes,SFsol%e_or_s_vel,SFsol%i_deep_or_finite,  &
                SFsol%n_h_steps,SFsol%n_four_modes,SFsol%nwrk,SFsol%yy,SFsol%zz)
           CALL incident_wf_finite(swenseDir, Wavefield, &
                Nx+2*GhostGridX,FineGrid%x,Ny+2*GhostGridY,FineGrid%y,Nz+GhostGridz,FineGrid%z,  &
                FineGrid%h,zero, SFsol%k,g,SFsol%n_four_modes,SFsol%zz,SFsol%yy)
        ENDIF
        PRINT*,'End of SWENSE Initialization'
     ENDIF
  CASE (12) ! Flat bottom, depth defined from SF-wave to check convergence...
     IF (Nx>1) THEN
        FineGrid%hx    = zero
        FineGrid%hxx   = zero
     ENDIF
     IF (Ny>1) THEN
        FineGrid%hy    = zero
        FineGrid%hyy   = zero
     ENDIF
     FineGrid%h         = SFsol%h
     kw = SFsol%k
     HH = SFsol%HH ! initial wave height
     Wavefield%E   = HH/two*COS(kw*FineGrid%x)*COS(kw*FineGrid%y)
     Wavefield%Ex  = -kw*HH/two*SIN(kw*FineGrid%x)*COS(kw*FineGrid%y)
     Wavefield%Exx = -kw**two*HH/two*COS(kw*FineGrid%x)*COS(kw*FineGrid%y)
     Wavefield%Ey  = -kw*HH/two*COS(kw*FineGrid%x)*SIN(kw*FineGrid%y)
     Wavefield%Eyy = -kw**two*HH/two*COS(kw*FineGrid%x)*COS(kw*FineGrid%y)
     ! Cross derivative Exy saved in Px
     Wavefield%Px  = kw**two*HH/two*SIN(kw*FineGrid%x)*SIN(kw*FineGrid%y)
     ! Put corners to wrong value
     Wavefield%E(1,1)=1000.0_long
     Wavefield%E(Nx+2*GhostGridX,1)=1000.0_long
     Wavefield%E(1,Ny+2*GhostGridY)=1000.0_long
     Wavefield%E(Nx+2*GhostGridX,Ny+2*GhostGridY)=1000.0_long
     !
     ! Test on thecross derivatives XZ and YZ...
     DO j=1,FineGrid%Nz+GhostGridZ
        PHI(j,:,:) = HH/two*COS(kw*FineGrid%x)*COS(kw*FineGrid%y)  &
             *COSH(kw*((FineGrid%z(j)-one)*FineGrid%h+FineGrid%h))/COSH(kw)!COSH(kw*FineGrid%h)
        !DXZ
        LASTPHI(j,:,:,1) = -kw**2*HH/two*SIN(kw*FineGrid%x)*COS(kw*FineGrid%y)  &
             *SINH(kw*((FineGrid%z(j)-one)&
             *FineGrid%h+FineGrid%h))/COSH(kw)!COSH(kw*FineGrid%h)
        !DYZ
        LASTPHI(j,:,:,2) = -kw**2*HH/two*COS(kw*FineGrid%x)*SIN(kw*FineGrid%y) &
             *SINH(kw*((FineGrid%z(j)-one)*FineGrid%h+FineGrid%h))/COSH(kw)!COSH(kw*FineGrid%h)
     ENDDO
     ! Put corners to zero...
     !PHI(1,1,1) = 1000.0_long
     !PHI(FineGrid%Nz+GhostGridZ,1,1) = 1000.0_long
     !PHI(1,Nx+2*GhostGridX,1) = 1000.0_long
     !PHI(FineGrid%Nz+GhostGridZ,Nx+2*GhostGridX,1) = 1000.0_long
     !PHI(1,1,Ny+2*GhostGridY) = 1000.0_long
     !PHI(FineGrid%Nz+GhostGridZ,1,Ny+2*GhostGridY) = 1000.0_long
     !PHI(FineGrid%Nz+GhostGridZ,Nx+2*GhostGridX,Ny+2*GhostGridY) = 1000.0_long
     !PHI(1,Nx+2*GhostGridX,Ny+2*GhostGridY) = 1000.0_long
     ! This is the whole edges of the domain which are not treated OK ?
     PHI(1:FineGrid%Nz+GhostGridZ,1,1) = 1000.0_long
     PHI(1:FineGrid%Nz+GhostGridZ,Nx+2*GhostGridX,1) = 1000.0_long
     PHI(1:FineGrid%Nz+GhostGridZ,1,Ny+2*GhostGridY) = 1000.0_long
     PHI(1:FineGrid%Nz+GhostGridZ,Nx+2*GhostGridX,Ny+2*GhostGridY) = 1000.0_long
     ! Bottom
     PHI(1,1:Nx+2*GhostGridX,1) = 1000.0_long
     PHI(1,1:Nx+2*GhostGridX,Ny+2*GhostGridY) = 1000.0_long
     PHI(1,1,1:Ny+2*GhostGridY) = 1000.0_long
     PHI(1,Nx+2*GhostGridX,1:Ny+2*GhostGridY) = 1000.0_long
     ! Free surface... We can usr free surface points
     !PHI(FineGrid%Nz+GhostGridZ,1:Nx+2*GhostGridX,1) = 1000.0_long
     !PHI(FineGrid%Nz+GhostGridZ,1:Nx+2*GhostGridX,Ny+2*GhostGridY) = 1000.0_long
     !PHI(FineGrid%Nz+GhostGridZ,1,1:Ny+2*GhostGridY) = 1000.0_long
     !PHI(FineGrid%Nz+GhostGridZ,Nx+2*GhostGridX,1:Ny+2*GhostGridY) = 1000.0_long
  CASE(13) ! Linear standing wave (2D) + in curvilinear space rotation with an angle
     FineGrid%h = SFsol%h
     IF (Nx>1) THEN
        FineGrid%hx=zero; FineGrid%hxx=zero
     ENDIF
     IF (Ny>1) THEN
        FineGrid%hy=zero; FineGrid%hyy=zero
     ENDIF
     IF(swenseONOFF==0)THEN
        CALL linearstandingwave2D(g,zero,zero,SFsol,SFsol%L,SFsol%L,FineGrid%h(1,1),FineGrid%x, &
             FineGrid%y&
             ,Wavefield%P,tmp2D,tmp2D,Wavefield%W,Wavefield%E,Wavefield%Ex,Wavefield%Exx,  &
             Wavefield%Ey,Wavefield%Eyy,&
             (Nx+2*GhostGridX)*(Ny+2*GhostGridY))
        IF (curvilinearONOFF == 1) THEN
           !Rotate the physical grid with a given angle...
           tmp2D = FineGrid%x
           FineGrid%x = tmp2D*COS(20.d0*PI/180.d0)-FineGrid%y*SIN(20.d0*PI/180.d0)
           FineGrid%y = tmp2D*SIN(20.d0*PI/180.d0)+FineGrid%y*COS(20.d0*PI/180.d0)
        ENDIF
     ELSE
        IF (curvilinearONOFF == 1) THEN
           !Rotate the physical grid with a given angle...
           tmp2D = FineGrid%x
           FineGrid%x = tmp2D*COS(20.d0*PI/180.d0)-FineGrid%y*SIN(20.d0*PI/180.d0)
           FineGrid%y = tmp2D*SIN(20.d0*PI/180.d0)+FineGrid%y*COS(20.d0*PI/180.d0)
        ENDIF
        IF (LinearONOFF==0) THEN
           ! Linear incident wavefield
           CALL incident_linear_wf_finite_standing(swenseDir, Wavefield, & !ramp_type,&
                Nx+2*GhostGridX,FineGrid%x,Ny+2*GhostGridY,FineGrid%y,Nz+GhostGridz,FineGrid%z, &
                FineGrid%h,zero,&
                SFsol%k,g,SFsol%HH,SFsol%h)
        ELSE
           ! Nonlinear incident wavefield
           print*,' this case is not available already... several checks to make...'
           stop
        ENDIF
        PRINT*,'End of SWENSE Initialization'
     ENDIF
     print*,'  Linear standing wave parameters:'
     print*,'      kh = ',SFsol%k*SFsol%h
     print*,'      T  = ',SFsol%T
     print*,'      c  = ',SFsol%c
     print*,'      h  = ',SFsol%h
     write(fileop(1),*)'  Linear standing wave parameters:'
     write(fileop(1),*)'      kh = ',SFsol%k*SFsol%h
     write(fileop(1),*)'      T  = ',SFsol%T
     write(fileop(1),*)'      c  = ',SFsol%c
     write(fileop(1),*)'      h  = ',SFsol%h
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
     PRINT*,'      dt = ',dt
     PRINT*,'      Cr = ',SFsol%c*dt/MIN(dxmin,dymin)	
     WRITE(FILEOP(1),*)'      dt = ',dt
     WRITE(FILEOP(1),*)'      Cr = ',SFsol%c*dt/MIN(dxmin,dymin)	
  CASE (14) ! Mildly nonlinear standing wave, deep water (Agnon & Glozman (1996)) with SWENSE
     IF (Nx>1) THEN
        FineGrid%h = two; FineGrid%hx=zero; FineGrid%hxx=zero
     ELSE IF (Ny>1) THEN
        FineGrid%h = two; FineGrid%hy=zero; FineGrid%hyy=zero
     ENDIF
     SFsol%h = FineGrid%h(1,1)
     IF(swenseONOFF==0)THEN
        IF (Nx>1) THEN
           CALL nonlinearstandingwave1D(pi,FineGrid%h(1,1),FineGrid%x,tmp2D,Wavefield%W,Wavefield%E, &
                Wavefield%Ex,Wavefield%Exx,&
                (Nx+2*GhostGridX)*(Ny+2*GhostGridY))
        ELSE IF (Ny>1) THEN
           CALL nonlinearstandingwave1D(pi,FineGrid%h(1,1),FineGrid%y,tmp2D,Wavefield%W,Wavefield%E, &
                Wavefield%Ey,Wavefield%Eyy,&
                (Nx+2*GhostGridX)*(Ny+2*GhostGridY))
        ENDIF
     ENDIF
     !GD: test in a curvilinear rotated domain...
     IF (curvilinearONOFF == 1) THEN
        !Rotate the physical grid with a given angle...
        tmp2D = FineGrid%x
        FineGrid%x = tmp2D*COS(20*PI/180)-FineGrid%y*SIN(20*PI/180)
        FineGrid%y = tmp2D*SIN(20*PI/180)+FineGrid%y*COS(20*PI/180)
     ENDIF
     IF (LinearONOFF==1) THEN
        ! Non-Linear incident wavefield
        CALL incident_nonlinear_wf_finite_standing(swenseDir, Wavefield, & !ramp_type,&
             Nx+2*GhostGridX,FineGrid%x,Ny+2*GhostGridY,FineGrid%y,Nz+GhostGridz,FineGrid%z,  &
             FineGrid%h,zero,&
             SFsol%k,g,SFsol%HH,SFsol%h)
     ELSE
        ! Linear incident wavefield
        print*,' this case is not available... for linear case 13'
        stop
     ENDIF
  CASE(15) !2D submerged bar test
     IF (Ny>1) THEN
        print*,'Only 2D case along x for now'
        STOP
     ENDIF
     ! Determine size of left relaxation zones... (assumed that this is defined in the 1st
     ! relaxation zone)
     y0 = FineGrid%x(RelaxZones(1)%idx(2),1)
     ! Define bottom
     IF (Lz>0) THEN
        CALL SubmergedBar_2D(FineGrid,y0,GhostGridX)
     ELSE
      !-------------------------------------------------------------------------                                                            
      ! Read in the bottom contours from a file:                                                                                            
      !-------------------------------------------------------------------------                                                            
      OPEN(unit=fileip(4),file=fname_bottom,status='old')
      READ(fileip(4),'(A)')head(4)
      PRINT *, 'SetUpInitialConditions:  Reading the bottom contours from file:'
      PRINT *,fname_bottom,' with header:'
      PRINT *, head(4)
      PRINT *, ' '
      WRITE(FILEOP(1),*) 'SetUpInitialConditions:  Reading the bottom contours from file:'
      WRITE(FILEOP(1),*)fname_bottom,' with header:'
      WRITE(FILEOP(1),*) head(4)
      WRITE(FILEOP(1),*) ' '
      READ(fileip(4),*)  DetermineBottomGradients ! Read the gradient flag                                                                  

      !-------------------------------------------------------------------------                                                            
      ! Read h(x,y),h_x,h_xx,h_y,h_yy                                                                                                       
      !-------------------------------------------------------------------------                                                            
      Print *, ' Reading h,h_x,h_xx,h_y,h_yy.'
      print *, ' '
      write(fileop(1),*) ' Reading h,h_x,h_xx,h_y,h_yy.'
      write(fileop(1),*) ' '
      Do j=1+GhostGridY,ny+GhostGridY
        DO i=1+GhostGridX,nx+GhostGridX
          READ(fileip(4),*) FineGrid%h(i,j),    &
                            FineGrid%hx(i,j),   &
                            FineGrid%hxx(i,j),  &
                            FineGrid%hy(i,j),   &
                            FineGrid%hyy(i,j)
        END DO
      END DO

      !-------------------------------------------------------------------------                                                            
      ! Extend the bathymetry to the ghost points by copying the domain                                                                     
      ! endpoint values.                                                                                                                    
      !-------------------------------------------------------------------------                                                            
      FineGrid%h(1,:)                   = FineGrid%h(1+GhostGridX,:);
      FineGrid%h(:,1)                   = FineGrid%h(:,1+GhostGridY);
      FineGrid%h(nx+2*GhostGridX,:)     = FineGrid%h(nx+GhostGridX,:);
      FineGrid%h(:,ny+2*GhostGridY)     = FineGrid%h(:,ny+GhostGridY)
      FineGrid%hx(1,:)                  = FineGrid%hx(1+GhostGridX,:);
      FineGrid%hx(:,1)                  = FineGrid%hx(:,1+GhostGridY);
      FineGrid%hx(nx+2*GhostGridX,:)    = FineGrid%hx(nx+GhostGridX,:);
      FineGrid%hx(:,ny+2*GhostGridY)    = FineGrid%hx(:,ny+GhostGridY)
      FineGrid%hxx(1,:)                 = FineGrid%hxx(1+GhostGridX,:);
      FineGrid%hxx(:,1)                 = FineGrid%hxx(:,1+GhostGridY);
      FineGrid%hxx(nx+2*GhostGridX,:)   = FineGrid%hxx(nx+GhostGridX,:);
      FineGrid%hxx(:,ny+2*GhostGridY)   = FineGrid%hxx(:,ny+GhostGridY)
      FineGrid%hy(1,:)                  = FineGrid%hy(1+GhostGridX,:);
      FineGrid%hy(:,1)                  = FineGrid%hy(:,1+GhostGridY);
      FineGrid%hy(nx+2*GhostGridX,:)    = FineGrid%hy(nx+GhostGridX,:);
      FineGrid%hy(:,ny+2*GhostGridY)    = FineGrid%hy(:,ny+GhostGridY)
      FineGrid%hyy(1,:)                 = FineGrid%hyy(1+GhostGridX,:);
      FineGrid%hyy(:,1)                 = FineGrid%hyy(:,1+GhostGridY);
      FineGrid%hyy(nx+2*GhostGridX,:)   = FineGrid%hyy(nx+GhostGridX,:);
      FineGrid%hyy(:,ny+2*GhostGridY)   = FineGrid%hyy(:,ny+GhostGridY)
        
     END IF
     !
     IF(swenseONOFF==1)THEN
        IF (LinearONOFF==0) THEN
           print*,'it should be a nonlinear case!'
           STOP
           ! Linear incident wavefield
           CALL incident_linear_wf_finite(swenseDir, Wavefield, & !ramp_type,&
                Nx+2*GhostGridX,FineGrid%x,Ny+2*GhostGridY,FineGrid%y,Nz+GhostGridz,FineGrid%z,  &
                FineGrid%h,zero,&
                SFsol%k,g,SFsol%HH,SFsol%h)
        ELSE
           ! Nonlinear incident wavefield
           ! Initialize stream function coefficients
           CALL stream_func_set_up(g,SFsol%h,SFsol%T,SFsol%i_wavel_or_per,SFsol%L,  &
                SFsol%k,SFsol%HH,SFsol%i_euler_or_stokes,SFsol%e_or_s_vel,SFsol%i_deep_or_finite,  &
                SFsol%n_h_steps,SFsol%n_four_modes,SFsol%nwrk,SFsol%yy,SFsol%zz)
           CALL incident_wf_finite(swenseDir, Wavefield, &
                Nx+2*GhostGridX,FineGrid%x,Ny+2*GhostGridY,FineGrid%y,Nz+GhostGridz,FineGrid%z,  &
                FineGrid%h,zero, &
                SFsol%k,g,SFsol%n_four_modes,SFsol%zz,SFsol%yy)
        ENDIF
     ENDIF
   CASE(16)

     !--------------------------------------------------------------------------
     ! Henrik Bredmose 2009 JFM paper setup.
     !--------------------------------------------------------------------------

     IF (Lz <= 0) THEN

      !-------------------------------------------------------------------------
      ! Read in the bottom contours from a file:
      !-------------------------------------------------------------------------
      OPEN(unit=fileip(4),file=fname_bottom,status='old')
      READ(fileip(4),'(A)')head(4)
      PRINT *, 'SetUpInitialConditions:  Reading the bottom contours from file:'
      PRINT *,fname_bottom,' with header:'
      PRINT *, head(4)
      PRINT *, ' '
      WRITE(FILEOP(1),*) 'SetUpInitialConditions:  Reading the bottom contours from file:'
      WRITE(FILEOP(1),*)fname_bottom,' with header:'
      WRITE(FILEOP(1),*) head(4)
      WRITE(FILEOP(1),*) ' '
      READ(fileip(4),*)  DetermineBottomGradients ! Read the gradient flag
      
      !-------------------------------------------------------------------------
      ! Read h(x,y),h_x,h_xx,h_y,h_yy
      !-------------------------------------------------------------------------
      Print *, ' Reading h,h_x,h_xx,h_y,h_yy.'
      print *, ' '
      write(fileop(1),*) ' Reading h,h_x,h_xx,h_y,h_yy.'
      write(fileop(1),*) ' '
      Do j=1+GhostGridY,ny+GhostGridY
        DO i=1+GhostGridX,nx+GhostGridX
          READ(fileip(4),*) FineGrid%h(i,j),    &
                            FineGrid%hx(i,j),   &
                            FineGrid%hxx(i,j),  &
                            FineGrid%hy(i,j),   &
                            FineGrid%hyy(i,j)
        END DO
      END DO
        
      !-------------------------------------------------------------------------
      ! Extend the bathymetry to the ghost points by copying the domain 
      ! endpoint values.
      !-------------------------------------------------------------------------
      FineGrid%h(1,:)                   = FineGrid%h(1+GhostGridX,:); 
      FineGrid%h(:,1)                   = FineGrid%h(:,1+GhostGridY);
      FineGrid%h(nx+2*GhostGridX,:)     = FineGrid%h(nx+GhostGridX,:);
      FineGrid%h(:,ny+2*GhostGridY)     = FineGrid%h(:,ny+GhostGridY)
      FineGrid%hx(1,:)                  = FineGrid%hx(1+GhostGridX,:); 
      FineGrid%hx(:,1)                  = FineGrid%hx(:,1+GhostGridY);
      FineGrid%hx(nx+2*GhostGridX,:)    = FineGrid%hx(nx+GhostGridX,:);
      FineGrid%hx(:,ny+2*GhostGridY)    = FineGrid%hx(:,ny+GhostGridY)
      FineGrid%hxx(1,:)                 = FineGrid%hxx(1+GhostGridX,:);
      FineGrid%hxx(:,1)                 = FineGrid%hxx(:,1+GhostGridY);
      FineGrid%hxx(nx+2*GhostGridX,:)   = FineGrid%hxx(nx+GhostGridX,:);
      FineGrid%hxx(:,ny+2*GhostGridY)   = FineGrid%hxx(:,ny+GhostGridY)
      FineGrid%hy(1,:)                  = FineGrid%hy(1+GhostGridX,:); 
      FineGrid%hy(:,1)                  = FineGrid%hy(:,1+GhostGridY);
      FineGrid%hy(nx+2*GhostGridX,:)    = FineGrid%hy(nx+GhostGridX,:);
      FineGrid%hy(:,ny+2*GhostGridY)    = FineGrid%hy(:,ny+GhostGridY)
      FineGrid%hyy(1,:)                 = FineGrid%hyy(1+GhostGridX,:);
      FineGrid%hyy(:,1)                 = FineGrid%hyy(:,1+GhostGridY);
      FineGrid%hyy(nx+2*GhostGridX,:)   = FineGrid%hyy(nx+GhostGridX,:);
      FineGrid%hyy(:,ny+2*GhostGridY)   = FineGrid%hyy(:,ny+GhostGridY)
      
      !-------------------------------------------------------------------------
      ! Re-set Lz to be max h(1,1)
      !-------------------------------------------------------------------------
      Lz=FineGrid%h(1,1)

    ELSE

      FineGrid%h      =   Lz; 
      FineGrid%hx     =   zero; 
      FineGrid%hxx    =   zero; 
      FineGrid%hy     =   zero; 
      FineGrid%hyy    =   zero;

    END IF

    !---------------------------------------------------------------------------
    ! Set initial surface elevation and free surface potential:
    !---------------------------------------------------------------------------
    CALL stream_func_set_up(g,SFsol%h,SFsol%T,SFsol%i_wavel_or_per,SFsol%L, &
                            SFsol%k,SFsol%HH,SFsol%i_euler_or_stokes,       &
                            SFsol%e_or_s_vel,SFsol%i_deep_or_finite,        &
                            SFsol%n_h_steps,SFsol%n_four_modes,SFsol%nwrk,  &
                            SFsol%yy,SFsol%zz)

    CALL incident_wf_finite(swenseDir, Wavefield,               &
                            Nx+2*GhostGridX,FineGrid%x,         &
                            Ny+2*GhostGridY,FineGrid%y,         &
                            Nz+GhostGridz  ,FineGrid%z,         &
                            FineGrid%h,zero,SFsol%k,g,          &
                            SFsol%n_four_modes,SFsol%zz,SFsol%yy)

    !---------------------------------------------------------------------------
    ! Modify the free surface elevation and the free surface potential by the 
    ! envelope function sech(kx(x-x0)/4) = 1/cosh(kx(x-x0)/4):
    !---------------------------------------------------------------------------
    Wavefield%E = 1.0/COSH(SFsol%k*(FineGrid%x - 200.0)/4.0)*Wavefield%E_I
    Wavefield%P = 1.0/COSH(SFsol%k*(FineGrid%x - 200.0)/4.0)*Wavefield%P_I_s

    !---------------------------------------------------------------------------
    ! Calculate the derivatives of the initial condition numerically:
    !---------------------------------------------------------------------------
    CALL PreProcessDiffStencils(FineGrid,FineGrid%DiffStencils,GhostGridX,  &
    GhostGridY,GhostGridZ,alpha,beta,gamma,fileop(1))
    CALL DiffXEven(Wavefield%E,Wavefield%Ex,1,FineGrid%Nx,FineGrid%Ny,1,    &
    FineGrid%DiffStencils,alpha)
    CALL DiffXEven(Wavefield%E,Wavefield%Exx,2,                             &
    FineGrid%Nx,FineGrid%Ny,1,FineGrid%DiffStencils,alpha)
    CALL DiffXEven(Wavefield%P,Wavefield%Px,1,FineGrid%Nx,FineGrid%Ny,1,    &
    FineGrid%DiffStencils,alpha)
  !
  CASE (17) ! Vincent and Briggs(3D)
     CALL BottomVincent(FineGrid,GhostGridX,GhostGridY)
     DetermineBottomGradients = 1
  CASE DEFAULT
    PRINT * ,'Error: Specified initial conditions are not valid.'
    STOP
  END SELECT

  END SUBROUTINE SetupInitialConditions
