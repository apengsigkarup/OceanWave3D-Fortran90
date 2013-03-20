      SUBROUTINE waveGenerationFromPaddleSignal()
      ! 3-dimensional wave generation.
      ! Implementation through inhomogeneous Neumann boundary conditions
      ! specified in BuildLinearSystem.f90
      !
      ! Inputs:
      !
      !
      ! Outputs:
      !         U:= Velocity array size(NZ,NY) (Global variable)
      !
      !
      ! Written by Bo Terp Paulsen, botp@mek.dtu.dk
      !

      ! Inputs / outputs
      !
      USE Precision
      USE GlobalVariables
      IMPLICIT NONE
      REAL(KIND=long) :: U
      ! Local variables
      !
      INTEGER :: j, k, Nx, Ny, Nz
      REAL(KIND=long) :: omega, dz, dy_local,FAC,arg
      !REAL(KIND=long), ALLOCATABLE :: y(:,:), z(:) 

      Nx = FineGrid%Nx + GhostGridX*2
      Ny = FineGrid%Ny + GhostGridY*2
      Nz = FineGrid%Nz + GhostGridZ

      !ALLOCATE(y(Nx,Ny),z(Nz))
      !dy_local = FineGrid%y(1,3)
      dz = FineGrid%z(3)

      !omega = 2*3.1415/1

        IF (time<4) THEN
          FAC = time/4
        ELSE
          FAC = 1.0
        ENDIF


      DO j =1+GhostGridY,FineGrid%Ny + GhostGridY
        DO k = 1+GhostGridZ,FineGrid%Nz + GhostGridZ
        !arg = 3*3.1415*(j-1)/(Ny-1)
              Uneumann(k,j) = FAC*waveFlux_inci(tstep,j-GhostGridY)!    
              !Uneumann(k,j) = omega*0.01*cos(omega*time + arg)*FAC !waveFlux_inci(tstep,1)!    
        END DO
      END DO
      END SUBROUTINE
