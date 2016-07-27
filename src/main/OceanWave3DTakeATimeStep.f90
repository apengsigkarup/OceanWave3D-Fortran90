SUBROUTINE OceanWave3DTakeATimeStep
!
! Take one time step to move the solution from tstep to tstep+1
!
  USE GlobalVariables
  USE MGLevels
  IMPLICIT NONE

  INTEGER i, j, k
  REAL(kind=long) :: maxEta, maxh


  TOTALITEROLD = TOTALITER

  CALL SYSTEM_CLOCK(count_0, count_rate, count_max)

  ! Print the simulation time to the screen every 10 time steps.
  IF(MOD(tstep,10)==0)THEN
     WRITE(6,2000) tstep, time
     WRITE(fileop(1),2000) tstep, time
  END IF
  !
  ! Use the free-surface conditions to step forward in time and get new
  ! values for eta and phi.
  !
  CALL timeintegration_fsbc
  !
  ! Check dw/dt to look for spots which may require heavy smoothing.  This feature is turned 
  ! off by setting accel_tol_fact to something larger than 100.  It is off when accel_tol_fact 
  ! does not appear in the input file.  ** Only implemented in 2D -HBB **
  !
  IF(FineGrid%ny==1)THEN
     CALL LocalSmoothing2D( FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY, &
          FineGrid%Nz+GhostGridZ,tstep,fileop(12) )
  ELSE
     CALL LocalSmoothing3D( FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY, &
          FineGrid%Nz+GhostGridZ,tstep,fileop(12) )
  END IF
  !
  ! Check for wave breaking at this step and if found update 
  ! the breaking geometry and roller history.  (Note that the breaking dynamics 
  ! is computed at each stage within the dfdt routine.)   
  !
  IF(BreakMod%i_breaking>0 .and. FineGrid%ny==1) THEN
     ! Call the breaking routine to evolve the rollers and output their geometry.
     ! This is only implemented in 2D.  
     CALL detect_breaking(fileop(14),FineGrid%Nx+2*GhostGridX,Wavefield,1)
  end IF
  !
  ! Check for an unstable solution and abort if found 
  !
  maxEta=maxval(abs(wavefield%E)); maxh=maxval(finegrid%h)
  IF ( maxEta .gt. 10.*maxh) Then
     write(6,*)' ********************************************************.' 
     write(6,*)' The solution looks to be going unstable, aborting here.' 
     write(6,*)' eta_max =',maxEta,' > 10 h_max =',10*maxh
     write(6,*)' ********************************************************.' 
     write(fileop(1),*)' ********************************************************.' 
     write(fileop(1),*)' The solution looks to be going unstable, aborting here.' 
     write(fileop(1),*)' eta_max =',maxEta,' > 10 h_max =',10*maxh
     write(fileop(1),*)' ********************************************************.' 
     ! CLOSE OPEN FILES
     CALL CloseIOFiles
     !
     ! DEALLOCATE
     CALL CloseVariables
     stop
  END IF

  time = time + dt

  CALL SYSTEM_CLOCK(count_1, count_rate, count_max)
  !
  ! Print the simulation time to the screen every 10 time steps.
  !
  If (tstep > 1 .and. tstep < 10) then
     IF (solver==0) THEN
       WRITE (6,2007) TOTALITER
       WRITE (6,2008) TOTALITERFS
       WRITE (6,2009) REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep-1,long))
       WRITE (fileop(1),2007) TOTALITER
       WRITE (fileop(1),2008) TOTALITERFS
       WRITE (fileop(1),2009) REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep-1,long))
     ELSE
       WRITE (6,2001) TOTALITER
       WRITE (6,2002) TOTALITERFS
       WRITE (6,2003) REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep-1,long))
       WRITE (fileop(1),2001) TOTALITER
       WRITE (fileop(1),2002) TOTALITERFS
       WRITE (fileop(1),2003) REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep-1,long))
     ENDIF
     WRITE (6,2005) MINITER, MAXITER
     WRITE (fileop(1),2005) MINITER, MAXITER
     IF (tstep>2) THEN
        IF (solver==0) THEN
          WRITE (6,2010) REAL(TOTALITER-TOTALITERFS,long)/(RKSTAGES*REAL(tstep-2,long))
          WRITE (fileop(1),2010) REAL(TOTALITER-TOTALITERFS,long)/(RKSTAGES*REAL(tstep-2,long))
        ELSE
          WRITE (6,2004) REAL(TOTALITER-TOTALITERFS,long)/(RKSTAGES*REAL(tstep-2,long))
          WRITE (fileop(1),2004) REAL(TOTALITER-TOTALITERFS,long)/(RKSTAGES*REAL(tstep-2,long))
        ENDIF
        WRITE (6,2006) MINITERNOFS, MAXITERNOFS
        WRITE (fileop(1),2006) MINITERNOFS, MAXITERNOFS
     ENDIF
  else
     ! Print the simulation time to the screen every 10 time steps.
     IF(MOD(tstep,10)==0)THEN
        print *, 'time step number ',tstep,' t=',time
        write(fileop(1),*) 'time step number ',tstep,' t=',time
        IF (solver==0) THEN
          WRITE (6,2009) REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep-1,long))
          WRITE (fileop(1),2009) REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep-1,long))
        ELSE
          WRITE (6,2003) REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep-1,long))
          WRITE (fileop(1),2003) REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep-1,long))
        END IF
        WRITE (6,2005) MINITER, MAXITER
        WRITE (fileop(1),2005) MINITER, MAXITER
     endif
  endif
 ! IF(MOD(tstep,10)==0)THEN
 !    WRITE (6,2020) ((count_1-count_0)*one)/count_rate
 !    WRITE (6,333) '    MINITER = ',MINITER,', MAXITER = ',MAXITER,', AVG. ITER=', &
 !         REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep,long))
 !    WRITE (6,334) '    STEP : ',tstep,', Step. iterations = ', &
 !         REAL(TOTALITER,long)-REAL(TOTALITEROLD,long)
 ! endif
  IF (StoreDataONOFF>0) THEN
     IF (MOD(tstep,StoreDataONOFF)==0) THEN		
        ! GD: SWENSE storage if necessary
        IF(swenseONOFF/=0) THEN
           ! GD: New variable definition
           ! Even numbers will be scattered wavefield
           CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%P,FineGrid,2*tstep,formattype)
           !CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%Ex+Wavefield%Ex_I,Wavefield%Ey+Wavefield%Ey_I,FineGrid,2*tstep,formattype)
           ! Odd numbers will be total wavefield
           CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E+Wavefield%E_I,&
                Wavefield%P+Wavefield%P_I_s,FineGrid,2*tstep+1,formattype)
        ELSE
           CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%P,FineGrid,tstep,formattype)
!	   CALL StoreDataVTK(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%P,FineGrid,tstep,formattype)
        !
        ! Write the file OceanWave3D.end, which can be used as initial conditions for a hot start
        WRITE(*,FMT='(A)') 'Writing OceanWave3D.end file, for possible restart.'
        WRITE(*,FMT='(A)') 'Change the file name to OceanWave3D.init and change IC to -1 for a hot start.'
        WRITE(fileop(1),FMT='(A)') 'Writing OceanWave3D.end file, for possible restart.'
        WRITE(fileop(1),FMT='(A)') 'Change the file name to OceanWave3D.init and change IC to -1 for a hot start.'
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
        close(fileip(3))
        ENDIF
     ENDIF
  ELSEIF(StoreDataOnOff<0)THEN
     IF (MOD(tstep,-StoreDataONOFF)==0) THEN		
        CALL StoreDataAscii(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%P, &
             FineGrid,-tstep/StoreDataOnOff)
        !
        ! Write the file OceanWave3D.end, which can be used as initial conditions for a hot start
        WRITE(*,FMT='(A)') 'Writing OceanWave3D.end file, for possible restart.'
        WRITE(*,FMT='(A)') 'Change the file name to OceanWave3D.init and change IC to -1 for a hot start.'
        WRITE(fileop(1),FMT='(A)') 'Writing OceanWave3D.end file, for possible restart.'
        WRITE(fileop(1),FMT='(A)') 'Change the file name to OceanWave3D.init and change IC to -1 for a hot start.'
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
        close(fileip(3))
     END IF
!
  ENDIF
  !
  ! If kinematics output is requested save it
  !
  IF(iKinematics/=0)THEN
          
             IF (iKinematics==20) THEN! Store binary kinematics files
                 Do i=1,nOutFiles
                     IF (tstep+1 >= Output(i)%tbeg .and. tstep+1 <= Output(i)%tend .and.  &
                           mod(tstep,Output(i)%tstride)==0 )THEN
                           CALL StoreKinematicData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                           FineGrid%Nz+GhostGridZ,i,tstep+1)
                     ENDIF 
                ENDDO
             ELSEIF (iKinematics==30) THEN ! Store wave gauges in ASCII format
                     
                        CALL StoreWaveGauges(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                        FineGrid%Nz+GhostGridZ,2,tstep+1)
             ENDIF
        
  END IF
  !
!
333 FORMAT (A,I4,A,I4,A,F5.2)
334 FORMAT (A,I4,A,F5.2)

2000 FORMAT('   Loop no. ',I8,'   Time = ',F9.4,' sec.')
2001 FORMAT(' A total of ',I8,' GMRES iterations done.')
2002 FORMAT(' A total of ',I8,' GMRES iterations done in first time step.')
2003 FORMAT(' An average of ',F8.2,' GMRES iterations done per solve.')
2004 FORMAT(' An average of ',F8.2,' GMRES iterations done per solve excluding first time step iteration counts.')
2005 FORMAT('   Minimum no. ', I8,' and maximum no. ',I8,' iterations per solve')
2006 FORMAT('   Minimum no. ', I8,' and maximum no. ',I8,' iterations per solve excluding first time step solves')
2007 FORMAT(' A total of ',I8,' DC iterations done.')
2008 FORMAT(' A total of ',I8,' DC iterations done in first time step.')
2009 FORMAT(' An average of ',F8.2,' DC iterations done per solve.')
2010 FORMAT(' An average of ',F8.2,' DC iterations done per solve excluding first time step iteration counts.')
2020 FORMAT('   Loop time = ',F12.3,' sec.')

END SUBROUTINE OceanWave3DTakeATimeStep
