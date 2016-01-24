PROGRAM OceanWave3D
  USE GlobalVariables
  USE MGLevels
  IMPLICIT NONE
  INTEGER i,j
  ! By Allan P. Engsig-Karup.
  ! 
  ! This is the main code for OceanWave3D which is split into preliminary set up work and 
  ! the work involved in moving from one time step to the next.  The idea is that these 
  ! two routines can be called from other similar drivers in order to couple the solver 
  ! to other codes e.g. OpenFOAM.  The time-step size dt, can be modified here and 
  ! everything should work.  
  !
  ! Initial set up work:  
  ! 
  CALL OceanWave3DT0Setup
  ! Initialise the timer 
  CALL SYSTEM_CLOCK(CPUinitial, count_rate, count_max)
  !
  ! Step through time and solve the problem 
  !
  print *, '  Starting to time step.'
  write(fileop(1),*) '  Starting to time step.'
  DO tstep=1,Nsteps-1 
     CALL OceanWave3DTakeATimeStep 
  END DO  
  !
  CALL SYSTEM_CLOCK(count_1, count_rate, count_max)
  WRITE (6,2030) ((count_1-CPUinitial)*one)/count_rate,Nsteps,FineGrid%Nx*FineGrid%Ny*(FineGrid%Nz+GhostGridZ)
  WRITE (fileop(1),2030) ((count_1-CPUinitial)*one)/count_rate,Nsteps,FineGrid%Nx*FineGrid%Ny*(FineGrid%Nz+GhostGridZ)
  !
  !
  ! Upon completion, write the file OceanWave3D.end, which can be used as initial conditions for a continued run
  WRITE(*,FMT='(A)') 'Writing OceanWave3D.end file, for possible restart.'
  WRITE(*,FMT='(A)') 'Change the file name to OceanWave3D.init and change IC to -1 for a hot start.'
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
2030 FORMAT(/,'   Totals: CPUtime = ',F10.3,' sec., Time steps = ',I5,', DOF=',I6,/)

END PROGRAM
