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
      REAL(KIND=long) :: omega, dz, dy_local 
      !REAL(KIND=long), ALLOCATABLE :: y(:,:), z(:) 

      Nx = FineGrid%Nx + GhostGridX*2
      Ny = FineGrid%Ny + GhostGridY*2
      Nz = FineGrid%Nz + GhostGridZ

      !ALLOCATE(y(Nx,Ny),z(Nz))
      !dy_local = FineGrid%y(2,1)
      !dz = FineGrid%z(2)

      !omega = 2*3.1415/1
      DO j =1,Ny
        DO k = 1,Nz
              Uneumann(k,j) = waveFlux_inci(tstep,1)!0.05*cos(omega*time)
        END DO
      END DO
      END SUBROUTINE
