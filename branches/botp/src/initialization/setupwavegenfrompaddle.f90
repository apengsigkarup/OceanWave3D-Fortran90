      SUBROUTINE setupwavegenfrompaddle
      ! Prepare for wave generation from paddle signal
      ! Written by Bo Terp Paulsen, botp@mek.dtu.dk
      !

      ! Inputs / outputs
      !
      USE GlobalVariables
      USE DataTypes
      IMPLICIT NONE
      INTEGER :: Ny, Nz, ndat, i,n2,j,nt
      REAL(KIND=long) :: dt_inc 
      REAL(KIND=long),ALLOCATABLE :: waveFlux_inci(:,:)

      CHARACTER(len=30) header
      

      Ny = FineGrid%Ny + GhostGridY*2
      Nz = FineGrid%Nz + GhostGridZ

      ! Read header
      !
      open(21,file=wave3DFlux%inc_wave_file,status='old')
      READ(21,'(A)',err=15)header
      READ(21,*)dt_inc, nt, n2
      ! Allocate fields
      !
      ALLOCATE(waveFlux_inci(nt,n2),wave3DFlux%time(nt),&
               wave3DFlux%flux(nt,n2),wave3DFlux%y(n2))

      ! Read rest of file
      !
      READ(21,*)(wave3DFlux%y(j),j=1,n2)
      DO i=1,nt
          READ(21,*,end=16) wave3DFlux%time(i),(waveFlux_inci(i,j),j=1,n2)
      ENDDO
      go to 14
15    print *, 'wave paddle.f90: No header line in the .iwf'
      stop
16    print *, 'error reading wave flux data'
      stop
14    CLOSE(21)

      IF(nt*dt_inc < Nsteps*dt) THEN 
        PRINT *, 'Time series for incident 3D-waves is to short, nt = ', &
        nt , '< Nsteps = ' , Nsteps ,'. Reduce number of timesteps or &
        zeropad signal.'  
      stop
      ENDIF

      ! Save data in global struct
      !
      wave3DFlux%flux = waveFlux_inci
      wave3DFlux%n2 = n2
      wave3DFlux%dt = dt_inc

      END SUBROUTINE
