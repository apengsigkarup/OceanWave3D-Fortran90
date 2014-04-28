 SUBROUTINE writeOceanwave3D(counter)
 ! Subroutine printing hotstart file for current timestep. 
 !
 ! Written by Bo Terp Paulsen, botp@mek.dtu.dk
 !
 USE GlobalVariables
 IMPLICIT NONE
 INTEGER :: i, j, counter
 character(80) outputFileName

 print *,"Writing OceanWave3D hotstart to folder OCW3Dhotstart/" 
 ! Create time name
 !
 WRITE (outputFileName, '(a, I0)') 'OCW3Dhotstart/OceanWave3D.',counter
  OPEN(fileip(3), file = outputFileName, status = 'unknown')
  WRITE(fileip(3),*) 'Initial conditions outputted from a previous simulation. Time = ',time ! Header
  WRITE(fileip(3),*) FineGrid%x(FineGrid%Nx,1)-FineGrid%x(1,1),FineGrid%y(FineGrid%Nx,1)-FineGrid%y(1,1), &
       FineGrid%Nx, FineGrid%Ny, time  ! Domain size, number of grid points and ending time.
  DO j=1+GhostGridY,FineGrid%Ny+GhostGridY
     DO i=1+GhostGridX,FineGrid%Nx+GhostGridX
        WRITE(fileip(3),*)WaveField%E(i,j),phi(FineGrid%Nz+GhostGridZ,i,j)
     END DO
  END DO

 CLOSE(fileip(3))
  END SUBROUTINE

