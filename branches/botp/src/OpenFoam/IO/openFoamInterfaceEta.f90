      SUBROUTINE openFoamInterface_eta(x0,etaOF)
      ! Subroutine returning the surface elevation at point x0             
      !                                                                    
      !                                                                    
      ! Input:                                                             
      !       x0: Array with dimension(3) [x y z]                         
      !                                                                    
      ! Output:                                                            
      !        etaOF: Return variable containing the surface     
      !        elevation at point x0                                       
      !                                                                    
      ! Written by Bo Terp Paulsen (botp), botp@mek.dtu.dk                 
      !                                                                    

      USE GlobalVariables 
      IMPLICIT NONE
      
      ! Input parameters
      !
      REAL(KIND=long), DIMENSION(:) :: x0(3)

      ! Return variable
      !
      REAL(KIND=long) :: etaOF

      ! Local variables 
      !
      INTEGER :: i, j, dir
      INTEGER, POINTER :: inOrOut
      INTEGER, DIMENSION(:), POINTER :: NN
      REAL(KIND=long), DIMENSION(:), ALLOCATABLE :: etaTemp
      REAL(KIND=long), DIMENSION(:), POINTER :: xStencil, yStencil, zStencil
      
      ! Assign local pointers
      !

      xStencil => interpolation%stencilX 
      yStencil => interpolation%stencilY 
      zStencil => interpolation%stencilZ
      NN       => interpolation%NN
      inOrOut  => interpolation%inOrOut

      ! Calculate local interpolation stencils
      !

      CALL makeStecils(x0) 
      
      IF (FineGrid%Ny == 1) THEN ! Surface elevation defined in 1D
                      
        
              ! Calculate surface elevation at x0
              !
              !etaOF = FineGrid%h(NN(1),NN(2))
              etaOF = 0
              Do i=1,2*alpha + 1 
                    etaOF = etaOF + xStencil(i) * WaveField%E(NN(1)-alpha-1+i,1)
              END DO
        
      ELSE ! Surface elevation defined in 2D

              ! Allocation of temporary field
              !
              ALLOCATE(etaTemp(2*beta+1))

              ! Interpolate temporary filed in x-direction
              !  
              DO j=1,2*beta +1

                  etaTemp(j) = 0

                  DO i=1,2*alpha + 1 
                        etaTemp(j) = etaTemp(j) + xStencil(i) *& 
                                     WaveField%E(NN(1)-alpha-1+i,NN(2)-beta-1+j)
                  END DO

              END DO


              ! Calculate surface elevation at x0 (y-interpolation)
              !
              !etaOF = FineGrid%h(NN(1),NN(2))
              etaOF = 0
              DO i=1,2*beta + 1 
                    etaOF = etaOF + yStencil(i) * etaTemp(i)
              END DO

      END IF

      END SUBROUTINE openFoamInterface_eta

