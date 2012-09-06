PROGRAM OceanWave3D
  USE GlobalVariables
  USE MGLevels
  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! External subroutines:
  !-----------------------------------------------------------------------------
  EXTERNAL FwdEuler,RungeKutta4
  
  INTEGER i,j,di,dj,i2,j2
  Real(kind=long) Fx,p1,p2,p3,p4,jac1,jac2,jac3,jac4,dxi1,dxi2,rho,dx1,dx2
  Real(kind=long) dxdxi1,dxdxi2,dydxi1,dydxi2,U1,U2,U3
  Real(kind=long), dimension(3) :: a1xa2 

  rho=1000.0

  !-----------------------------------------------------------------------------
  ! This is the main code for OceanWave3D which is split into preliminary 
  ! set up work and the work involved in moving from one time step to the next.
  ! The idea is that these two routines can be called from other similar 
  ! drivers in order to couple the solver to other codes e.g. OpenFOAM.  
  ! The time-step size dt, can be modified here and everything should work.  
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! Initial set up work:
  !-----------------------------------------------------------------------------
  print*,"OceanWave3D: initial setup start"
  CALL OceanWave3DT0Setup
  print*,"OceanWave3D: initial setup finished"
   
  !-----------------------------------------------------------------------------
  ! Initialise the timer 
  !-----------------------------------------------------------------------------
  CALL SYSTEM_CLOCK(CPUinitial, count_rate, count_max)
  
  !-----------------------------------------------------------------------------
  ! Calculation of ship surface area:
  !-----------------------------------------------------------------------------
  Fx = 0.0
  dxi1 = 1.0_long/(FineGrid%Nx-1)
  dxi2 = 1.0_long/(FineGrid%Ny-1)
  do i=1+GhostGridX,FineGrid%Nx+GhostGridX-1
    do j=1+GhostGridY,FineGrid%Ny+GhostGridY-1
      do di=0,1
        do dj=0,1

          i2 = i + di
          j2 = j + dj

          dxdxi1 = 1.0_long/dxi1*FineGrid%CurvilinearStuff%xe(i2,j2)
          dydxi1 = 1.0_long/dxi1*FineGrid%CurvilinearStuff%ye(i2,j2)
          dxdxi2 = 1.0_long/dxi2*FineGrid%CurvilinearStuff%xn(i2,j2)
          dydxi2 = 1.0_long/dxi2*FineGrid%CurvilinearStuff%yn(i2,j2)

          !-------------------------------------------------------------------
          ! Cross product of covariant basis vectors:
          !-------------------------------------------------------------------
          a1xa2(1)   = dydxi1*dzdxi2(i2,j2) - dydxi2*dzdxi1(i2,j2)
          a1xa2(2)   = dzdxi1(i2,j2)*dxdxi2 - dzdxi2(i2,j2)*dxdxi1
          a1xa2(3)   = dxdxi1*dydxi2        - dxdxi2*dydxi1

          !-------------------------------------------------------------------
          ! Quadrature summation of area:
          !-------------------------------------------------------------------
          Fx = Fx  &
          + 0.25*dxi1*dxi2*wgt_ship(i2,j2) &
          * sqrt(a1xa2(1)**2.0 + a1xa2(2)**2.0 + a1xa2(3)**2.0)

        enddo
      enddo
    enddo
  enddo
  print*,'Ship surface area: ', Fx

  !-----------------------------------------------------------------------------
  ! Initalize the Runge-Kutta timestepping by a call to the right hand side 
  ! function:
  !-----------------------------------------------------------------------------
  CALL RHSLinearShip2D(time,WaveField,RHSE,RHSP)

  open(unit=9,file="forcex")
  write(unit=9,fmt='(A)') '# time, force_x'
  close(unit=9)

  !-----------------------------------------------------------------------------
  ! Step through time and solve the problem:
  !-----------------------------------------------------------------------------
  DO tstep=1,Nsteps-1

    print*,"Time = ",time," step = ",tstep
 
    !---------------------------------------------------------------------------
    ! Runge Kutta time step:
    !---------------------------------------------------------------------------
    CALL RungeKutta(RungeKutta4,RHSLinearShip2D,time,dt,WaveField,RHSE,RHSP)

    !---------------------------------------------------------------------------
    ! Filtering:
    !---------------------------------------------------------------------------
    IF (filteringONOFF>0) THEN
      IF (MOD(tstep,filteringONOFF)==0) THEN
        CALL FILTERING(FineGrid%Nx,FineGrid%Ny,Wavefield%E,filterNP,&
        filterALPHA,filtercoefficients,tstep)
        CALL FILTERING(FineGrid%Nx,FineGrid%Ny,Wavefield%P,filterNP,&
        filterALPHA,filtercoefficients,tstep)
      ENDIF
    ENDIF

    !---------------------------------------------------------------------------
    ! Update time:
    !---------------------------------------------------------------------------
    time = time + dt     

    !---------------------------------------------------------------------------
    ! Calculate wave resistance on ship:
    !---------------------------------------------------------------------------
    Fx = 0.0
    dxi1 = 1.0_long/(FineGrid%Nx-1)
    dxi2 = 1.0_long/(FineGrid%Ny-1)
    do i=1+GhostGridX,FineGrid%Nx+GhostGridX-1
      do j=1+GhostGridY,FineGrid%Ny+GhostGridY-1
        do di=0,1
          do dj=0,1

            i2 = i + di
            j2 = j + dj

            dydxi1  = 1.0_long/dxi1*FineGrid%CurvilinearStuff%ye(i2,j2)
            dydxi2  = 1.0_long/dxi2*FineGrid%CurvilinearStuff%yn(i2,j2)

            jac1    = FineGrid%CurvilinearStuff%xe(i2,j2)/dxi1
            jac2    = FineGrid%CurvilinearStuff%yn(i2,j2)/dxi2

            U1      = (1.0_long/jac1*FPx(i2,j2))/Uship 
            U2      = (1.0_long/jac2*FPy(i2,j2))
            U3      = WaveField%W(i2,j2)

            p1      = rho*(Uship*U1 - Wavefield%P0(i2,j2) &
                                    - 0.5*(U1*U1 + U2*U2 + U3*U3))

            !-------------------------------------------------------------------
            ! Cross product of covariant basis vectors:
            !-------------------------------------------------------------------
            a1xa2(1)   = dydxi1*dzdxi2(i2,j2) - dydxi2*dzdxi1(i2,j2)

            !-------------------------------------------------------------------
            ! Quadrature sumation of force:
            !-------------------------------------------------------------------
            Fx = Fx  &
            + 0.25*dxi1*dxi2*p1*wgt_ship(i2,j2)*a1xa2(1)

          enddo
        enddo
      enddo
    enddo
    print*,time
    print*,Fx

    open(unit=9,file="forcex",position="APPEND")
    write(unit=9,fmt='(F24.16,F24.16)') time,Fx
    close(unit=9)

    !---------------------------------------------------------------------------
    ! Save Data:
    !---------------------------------------------------------------------------
    IF (StoreDataONOFF>0) THEN
      IF (MOD(tstep,StoreDataONOFF)==0) THEN
        CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
        Wavefield%E,Wavefield%P,FineGrid,tstep,formattype)
      ENDIF
      ELSEIF(StoreDataOnOff<0)THEN
      IF (MOD(tstep,-StoreDataONOFF)==0) THEN
        CALL StoreDataAscii(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
        Wavefield%E,Wavefield%P,FineGrid,-tstep/StoreDataOnOff)
      END IF
    ENDIF
    !     CALL OceanWave3DTakeATimeStep
  END DO  

  !-----------------------------------------------------------------------------
  ! Finalize the timer 
  !-----------------------------------------------------------------------------
  CALL SYSTEM_CLOCK(count_1, count_rate, count_max)
  WRITE (6,2030) ((count_1-CPUinitial)*one)/count_rate,Nsteps,&
  FineGrid%Nx*FineGrid%Ny*(FineGrid%Nz+GhostGridZ)

  !-----------------------------------------------------------------------------
  ! Upon completion, write the file OceanWave3D.end, which can be used as 
  ! initial conditions for a continued run
  !-----------------------------------------------------------------------------
  PRINT *, 'Writing OceanWave3D.end file, for possible restart.  Change the file name to OceanWave3D.init and change IC to -1 for a hot start.'
  CLOSE(fileip(3))
  OPEN(fileip(3), file = 'OceanWave3D.end', status = 'unknown')
  WRITE(fileip(3),*) 'Initial conditions outputted from a previous simulation.' ! Header
  WRITE(fileip(3),*) FineGrid%x(FineGrid%Nx,1)-FineGrid%x(1,1),FineGrid%y(FineGrid%Nx,1)-FineGrid%y(1,1), &
       FineGrid%Nx, FineGrid%Ny, time  ! Domain size, number of grid points and ending time.
  DO j=1+GhostGridY,FineGrid%Ny+GhostGridY
     DO i=1+GhostGridX,FineGrid%Nx+GhostGridX
        WRITE(fileip(3),*)WaveField%E(i,j),phi(FineGrid%Nz+GhostGridZ,i,j)
     END DO
  END DO
  ! 
  !
  ! OUTPUT JOB STATUS
  WRITE (6,1999)
  !
  ! CLOSE OPEN FILES
  CALL CloseIOFiles
  !
  ! DEALLOCATE
  CALL CloseVariables

1999 FORMAT(/,' JOB IS COMPLETE',/)
2030 FORMAT(/,' Totals: CPUtime = ',F10.3,' sec., Time steps = ',I5,', DOF=',I6,/)

END PROGRAM
