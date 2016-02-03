SUBROUTINE setupwavepaddle
  !
  ! Prepare for wave generation from a paddle signal applied as 
  ! a uniform flux along the Western boundary. 
  !
  ! Written by Bo Terp Paulsen, botp@mek.dtu.dk
  !
  ! Feb. 2016: modified by  Harry B. Bingham - hbb@mek.dtu.dk.
  ! Added output info to screen and LOG file and a check for 
  ! the correct uniform spacing of flux points in the y-direction. 
  ! Also added the free surface slope condition based on the linear 
  ! FSBC eta_x = -1/g u_t to be applied in the ghost point update 
  ! in UpdateGhostLayerX.f90. 
  !
  !
  !
  ! Inputs / outputs
  !
  ! The flux boundary conditions are read from the input file and 
  ! loaded into the structure: wave3DFlux with variables: 
  !                                      %nt
  !                                      %n2 (ny)
  !                                      %dt_inc
  !                                      %flux(nt,n2)
  !                                      %etax(nt,n2)
  !
  USE GlobalVariables
  USE DataTypes
  IMPLICIT NONE
  !
  ! local variables
  INTEGER :: Ny, Nz, ndat, i, n2, j, nt, flag
  REAL(KIND=long) :: dt_inc, dy_local, FluxTime(5), FluxTimeC(5,5), gfac
  REAL(KIND=long),ALLOCATABLE :: waveFlux_inci(:,:)  !hbb should remove this...
  CHARACTER(len=30) header


  Ny = FineGrid%Ny + GhostGridY*2
  Nz = FineGrid%Nz + GhostGridZ
  !
  ! Read the input file header 
  !
  open(21,file=wave3DFlux%inc_wave_file,status='old')
  READ(21,'(A)',err=15)header
  READ(21,*)dt_inc, nt, n2
  ! Informative output
  write(6,23),wave3DFlux%inc_wave_file,header
  write(fileop(1),23),wave3DFlux%inc_wave_file,header
23 FORMAT('Reading the wall boundary flux from file:  ',A,'with header ',/,A)
  !
  ! Allocate fields
  !
  ALLOCATE(waveFlux_inci(nt,n2),wave3DFlux%time(nt),&
       wave3DFlux%flux(nt,n2),wave3DFlux%y(n2), wave3DFlux%etax(nt,n2))

  ! Read the rest of file
  !
  READ(21,*)(wave3DFlux%y(j),j=1,n2)
  DO i=1,nt
     READ(21,*,end=16) wave3DFlux%time(i),(waveFlux_inci(i,j),j=1,n2)
  ENDDO
  go to 14
15 print *, 'wave paddle.f90: No header line'
  stop
16 print *, 'error reading wave flux data'
  stop
14 CLOSE(21)
  !
  ! Informative output
  write(6,24)nt,dt_inc,n2,wave3DFlux%y
  write(fileop(1),24)nt,dt_inc,n2,wave3DFlux%y
24 FORMAT('Found nt=',i10,' time steps at dt = ',e12.4,' for ny= ',  &
        i10,' points in the y-direction ',/, 'at positions: ',/,1000e12.4)
  !
  ! Check that the points are spaced correctly in y. 
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
     PRINT *, 'Time series for the boundary flux incident waves is too short, (nt*dt)_solver &
          < (Nsteps*dt)_waveMaker.'
     stop
  ENDIF
  !
  ! Save the data in the global structure
  !
  write(6,28)wave3DFlux%order-1
28 Format('Interpolation in the y-direction will be done to order ',i10,//)

  wave3DFlux%flux = waveFlux_inci
  deallocate(waveFlux_inci)
  wave3DFlux%n2 = n2
  wave3DFlux%dt = dt_inc
  !
  ! Take a 4th-order derivative of the flux in time and build eta_x=-1/g u_t.
  !
  Do i=1,5
     FluxTime(i)=(i-1)*dt_inc
  END Do
  CALL BuildStencilsGridX(2,1,FluxTime,FluxTimeC,5,1)
  gfac=-1._long/g
  Do j=1,n2
     Do i=1,2
        wave3DFlux%etax(i,j)=gfac*Dot_Product(FluxTimeC(i,:),wave3DFlux%flux(1:5,j))
     END Do
     Do i=3,nt-2
        wave3DFlux%etax(i,j)=gfac*Dot_Product(FluxTimeC(3,:),wave3DFlux%flux(i-2:i+2,j))
     END Do
     Do i=nt-1,nt
        wave3DFlux%etax(i,j)=gfac*Dot_Product(FluxTimeC(5-(nt-i),:),wave3DFlux%flux(nt-4:nt,j))
     END Do
  END Do

END SUBROUTINE setupwavepaddle
