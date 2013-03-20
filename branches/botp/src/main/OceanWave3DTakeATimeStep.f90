SUBROUTINE OceanWave3DTakeATimeStep()
  USE GlobalVariables
  USE MGLevels
  IMPLICIT NONE

  INTEGER i, j, k


  TOTALITEROLD = TOTALITER

  CALL SYSTEM_CLOCK(count_0, count_rate, count_max)

  ! Print the simulation time to the screen every 10 time steps.
  !IF(MOD(tstep,10)==0)THEN
     !WRITE(6,2000) tstep, time
  !END IF
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
  time = time + dt
print *, time
  CALL SYSTEM_CLOCK(count_1, count_rate, count_max)
  !
  ! Print the simulation time to the screen every 10 time steps.
  !
  !If (tstep > 1 .and. tstep < 10) then
     !IF (solver==0) THEN
       !WRITE (6,2007) TOTALITER
       !WRITE (6,2008) TOTALITERFS
       !WRITE (6,2009) REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep-1,long))
     !ELSE
       !WRITE (6,2001) TOTALITER
       !WRITE (6,2002) TOTALITERFS
       !WRITE (6,2003) REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep-1,long))
     !ENDIF
     !WRITE (6,2005) MINITER, MAXITER
     !IF (tstep>2) THEN
        !IF (solver==0) THEN
          !WRITE (6,2010) REAL(TOTALITER-TOTALITERFS,long)/(RKSTAGES*REAL(tstep-2,long))
        !ELSE
          !WRITE (6,2004) REAL(TOTALITER-TOTALITERFS,long)/(RKSTAGES*REAL(tstep-2,long))
        !ENDIF
        !WRITE (6,2006) MINITERNOFS, MAXITERNOFS
     !ENDIF
  !else
     !! Print the simulation time to the screen every 10 time steps.
     !IF(MOD(tstep,10)==0)THEN
        !print *, 'time step number ',tstep,' t=',time
        !IF (solver==0) THEN
          !WRITE (6,2009) REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep-1,long))
        !ELSE
          !WRITE (6,2003) REAL(TOTALITER,long)/(RKSTAGES*REAL(tstep-1,long))
        !END IF
        !WRITE (6,2005) MINITER, MAXITER
     !endif
  !endif
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
        ENDIF
     ENDIF
  ELSEIF(StoreDataOnOff<0)THEN
     IF (MOD(tstep,-StoreDataONOFF)==0) THEN		
        CALL StoreDataAscii(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%P, &
             FineGrid,-tstep/StoreDataOnOff)
     END IF
  ENDIF
  !
  ! If kinematics output is requested save it
  !
  IF(iKinematics/=0)THEN
     Do i=1,nOutFiles
        IF (tstep+1 >= Output(i)%tbeg .and. tstep+1 <= Output(i)%tend .and.  &
             mod(tstep,Output(i)%tstride)==0 )THEN
           CALL StoreKinematicData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
                FineGrid%Nz+GhostGridZ,i,tstep+1)
        END IF
     END Do
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
