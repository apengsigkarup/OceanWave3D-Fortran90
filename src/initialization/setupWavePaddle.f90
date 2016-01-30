      SUBROUTINE setupwavepaddle
      ! Prepare for wave generation from paddle signal
      ! Written by Bo Terp Paulsen, botp@mek.dtu.dk
      !
      ! Added output info to screen and LOG file and a check for 
      ! the correct uniform spacing of flux points in the y-direction. - hbb@mek.dtu.dk.
      !

      ! Inputs / outputs
      !
      USE GlobalVariables
      USE DataTypes
      IMPLICIT NONE
      INTEGER :: Ny, Nz, ndat, i,n2,j,nt, flag
      REAL(KIND=long) :: dt_inc, dy_local 
      REAL(KIND=long),ALLOCATABLE :: waveFlux_inci(:,:)

      CHARACTER(len=30) header
      

      Ny = FineGrid%Ny + GhostGridY*2
      Nz = FineGrid%Nz + GhostGridZ

      ! Read header
      !
      open(21,file=wave3DFlux%inc_wave_file,status='old')
      READ(21,'(A)',err=15)header
      READ(21,*)dt_inc, nt, n2
      write(6,23),wave3DFlux%inc_wave_file,header
      write(fileop(1),23),wave3DFlux%inc_wave_file,header
23    FORMAT('Reading the wall boundary flux from file:  ',A,'with header ',/,A)
      !
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
15    print *, 'wave paddle.f90: No header line'
      stop
16    print *, 'error reading wave flux data'
      stop
14    CLOSE(21)
      !
      write(6,24)nt,dt_inc,n2,wave3DFlux%y
      write(fileop(1),24)nt,dt_inc,n2,wave3DFlux%y
24    FORMAT('Found nt=',i10,' time steps at dt = ',e12.4,' for ny= ',i10,' points in the y-direction ',/, &
           'at positions: ',/,1000e12.4)
      !
      ! Check that the points are spaced correctly. 
      !
      dy_local=2*wave3DFlux%y(1); 
      flag=0
      do i=2,n2-1
         if(abs(wave3DFlux%y(i)-(wave3DFlux%y(i-1)+dy_local)).ge.1.e-5) then
            flag=1
         end if
      end do
      if (flag .ne. 0)then
         write(6,*)' ** Your wavemaker flux points must be uniformly spaced from y=dy/2 to ymax-dy/2. **'
         stop
      end if

      IF(nt*dt_inc < Nsteps*dt) THEN 
        PRINT *, 'Time series for incident waves is to short, (nt*dt)_solver &
        < (Nsteps*dt)_waveMaker.'
      stop
      ENDIF

      ! Save data in global struct
      !
      wave3DFlux%flux = waveFlux_inci
      wave3DFlux%n2 = n2
      wave3DFlux%dt = dt_inc

      END SUBROUTINE
