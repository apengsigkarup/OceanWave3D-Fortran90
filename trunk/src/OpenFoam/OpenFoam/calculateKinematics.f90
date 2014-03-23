      SUBROUTINE CalculateKinematics()
      ! Subroutine calculating the kinematics and free surface elevation in the entire
      ! domain. 
      !
      ! Input: 
      !
      ! Output: The global variables UOF, VOF, WOF and dOF are updated
      !
      ! Written by Bo Terp Paulsen (botp), botp@mek.dtu.dk                 
      !
      
      USE GlobalVariables
      IMPLICIT NONE
      ! Local variables
      !
      INTEGER ::  i, j, k, Nx, Ny, Nz
      REAL(KIND=long), DIMENSION(:,:), POINTER :: x, y, h, hx, hy, eta, etax, etay
      REAL(KIND=long), DIMENSION(:), POINTER   :: z
      
      ! Assign local grid range
      !
      Nx = FineGrid%Nx+2*GhostGridX
      Ny = FineGrid%Ny+2*GhostGridY
      Nz = FineGrid%Nz+GhostGridZ
      
      ! Assign the local pointers
      !
      x => FineGrid%x; y => FineGrid%y; z => FineGrid%z; h => FineGrid%h
      hx => FineGrid%hx
      hy => FineGrid%hy; eta => WaveField%E; etax => WaveField%Ex
      etay => WaveField%Ey
      
      ! Ensure that matrices are cleaned
      !
      UOF = zero
      VOF = zero
      WOF = zero
      dOF = zero

      ! Shift the horizontal grid point position to account for the ghost points.  
      !
      !
      ! The fluid thickness d=h+eta
      !
      DO j=1,Ny
         DO i=1,Nx
            dOF(i,j)=h(i,j)+eta(i,j);
         end DO
      End DO
      !
      ! Then the velocity potential at all points in this horizontal slice and at all sigma (vertical) locations.  
      !
      !
      ! Then the velocities at all points in this horizontal slice and at all sigma 
      ! (vertical) locations.  
      !
      !
      ! Compute dphi/dsigma
      !
      CALL DiffZArbitrary(phi,WOF,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
           FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
      IF (FineGrid%Nx>1) THEN
         !
         ! Compute dphi/dx 
         !	
         CALL DiffXEven(phi,UOF,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
              FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,alpha)
         IF ( LinearOnOff /= 0) THEN
            !
            ! Add in the chain rule contribution to get the velocity
            !
            DO j=1,Ny
               DO i=1,Nx
                  DO k=1,Nz
                  UOF(k,i,j) = UOF(k,i,j) + ((1-z(k))/dOF(i,j)*hx(i,j)-z(k)/dOF(i,j)*etax(i,j))*WOF(k,i,j)
                  END DO
               END DO
            END DO
         END IF
      ELSE
         UOF=zero
      END IF
      
       
      IF (FineGrid%Ny>1) THEN
         ! dphi/dy
         CALL DiffYEven(phi,VOF,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
              FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,beta)
         IF ( LinearOnOff /= 0) THEN
            Do j=1,Ny
               Do i=1,Nx
                  Do k=1,Nz
                  VOF(k,i,j) = VOF(k,i,j)+((1-z(k))/dOF(i,j)*hy(i,j)-z(k)/dOF(i,j)*etay(i,j))*WOF(k,i,j)
                  END Do
               END Do
            END Do
         END IF
      ELSE
         VOF=zero
      END IF
     
      IF ( LinearOnOff /= 0) THEN
            Do j=1,Ny
               Do i=1,Nx
                  Do k=1,Nz
                  WOF(k,i,j) = WOF(k,i,j)/dOF(i,j)
                  END Do
               END Do
            END Do
         END IF


END SUBROUTINE CalculateKinematics
