      SUBROUTINE setupwavegenfrompaddle
      ! Prepare for wave generation from paddle signal
      ! Written by Bo Terp Paulsen, botp@mek.dtu.dk
      !

      ! Inputs / outputs
      !
      USE GlobalVariables
      USE DataTypes
      IMPLICIT NONE
      INTEGER :: Nx, Ny, Nz, ndat, i
      REAL(KIND=long) :: dt_inc 
      CHARACTER(len=30) header
      

      ALLOCATE(waveFlux_inci(Nsteps,FineGrid%Ny))
      !ALLOCATE(waveFlux_tmp(Nsteps))

      Nx = FineGrid%Nx + GhostGridX*2
      Ny = FineGrid%Ny + GhostGridY*2
      Nz = FineGrid%Nz + GhostGridZ

      ! Read paddle signal from file
      
      !IF(i_spec==2)THEN
      open(21,file='paddleSignal',status='old')
      READ(21,'(A)',err=15)header
      READ(21,*)dt_inc
          IF(ABS(dt-dt_inc)>1.e-6)THEN
             print *, 'setupWaveField.f90: The .inp and .iwf time steps &
                 do not agree.'
             print *, header
             print *, dt,dt_inc
             stop
          END IF

          DO i=1,Nsteps-1
             READ(21,*,end=13)waveFlux_inci(i,:)
             !DO j=2,ns_inc
                !READ(21,*)dum
             !END do
          END DO
          go to 14
13    ndat=i-1
      print *, ' Found ',ndat,' data points in the incident wave file.  &
          Padding with zeros up to ',Nsteps
      DO i=ndat+1,Nsteps
         waveFlux_inci(i,:)=zero
      END DO
      go to 14
15    print *, 'wave paddle.f90: No header line in the .iwf'
      stop
14    CLOSE(21)

      ALLOCATE(Uneumann(Nz+GhostGridZ,Ny+GhostGridY))

      END SUBROUTINE
