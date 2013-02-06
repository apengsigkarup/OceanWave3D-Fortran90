      SUBROUTINE openFoamInterface_U(x0,U,V,W)
      ! Subroutine returning the flow velocity in the x-direction at point    x0             
      !                                                                    
      ! Input:                                                             
      !       x0: Array with dimension(3) [x y z]                         
      !                                                                    
      ! Output:                                                            
      !        U: Return variable containing the u-velocity at x0
      !                                                                    
      ! Written by Bo Terp Paulsen (botp), botp@mek.dtu.dk                 
      !                                                                    

      USE GlobalVariables 
      USE OFmodule
      IMPLICIT NONE
      
      ! Input parameters
      !
      REAL(KIND=long), DIMENSION(:) :: x0(3)

      ! Return variable
      !
      REAL(KIND=long) :: U, V, W

      ! Local variables 
      !
      INTEGER :: i, j, k, dir, Nz, Nzp, left, right, counterX, & 
      counterY, counterZ
      INTEGER, POINTER :: inOrOut
      INTEGER, DIMENSION(:), POINTER :: NN
      REAL(KIND=long), DIMENSION(:), POINTER :: xStencil, yStencil, zStencil
      REAL(KIND=long), DIMENSION(:), ALLOCATABLE :: uTemp1D,uTemp
      REAL(KIND=long), DIMENSION(:), ALLOCATABLE :: vTemp1D,vTemp
      REAL(KIND=long), DIMENSION(:), ALLOCATABLE :: wTemp1D,wTemp
      REAL(KIND=long), DIMENSION(:,:), POINTER :: x, y 
      
      Nz = FineGrid%Nz; 
      ! Assign local pointers
      !
      x => FineGrid%x; y => FineGrid%y;  

      xStencil => interpolation%stencilX 
      yStencil => interpolation%stencilY 
      zStencil => interpolation%stencilZ
      NN       => interpolation%NN
      inOrOut  => interpolation%inOrOut

      ! Calculate local interpolation stencils 
      !
      CALL makeStecils(x0) 

      ! Temporary fields for debugging
      !

            !DO i =1,SIZE(VOF,1)
                !VOF(i,:,:) = FineGrid%y(:,:)
                !UOF(i,:,:) = FineGrid%x(:,:)
            !END DO
            
            !DO j = 1,SIZE(VOF,2)
                !DO k = 1,SIZE(VOF,3)
                    !WOF(:,j,k) = FineGrid%z(:)
                !END DO 
            !END DO
      
      IF (inOrOut == 0) THEN
              ! Limits in z-direction
              !
              IF (GhostGridZ==0) THEN; Nzp = Nz; ELSE; Nzp = Nz+1; ENDIF
              IF(NN(3)<=gamma) THEN
                      left =  NN(3) - 1
                      right = 2*gamma  -left
              ELSE
                      IF (NN(3)>=Nzp-gamma) THEN
                             right = Nzp - NN(3) 
                             left = 2*gamma  - right
                      ELSE
                             left = gamma
                             right = gamma
                      END IF
              END IF

          IF (FineGrid%Ny==1) THEN ! The velocity field is defined in 2D (x/z)
              ALLOCATE(uTemp1D(2*gamma+1))
              ALLOCATE(wTemp1D(2*gamma+1))

              ! Interpolation in x-direction, creating temporary field
              ! in z
              !
              counterZ = 1
              DO k=-left,right

                  uTemp1D(counterZ) = 0
                  wTemp1D(counterZ) = 0

                  DO i=1,2*alpha + 1 
                        uTemp1D(counterZ)=uTemp1D(counterZ)+xStencil(i)&
                        *UOF(NN(3)+k,NN(1)-alpha-1+i,NN(2))

                        wTemp1D(counterZ)=wTemp1D(counterZ)+xStencil(i)&
                        *WOF(NN(3)+k,NN(1)-alpha-1+i,NN(2))
                  END DO
              counterZ = counterZ + 1

              END DO

              ! Calculate velocities at x0 (z-interpolation)
              !
              U =0
              V =0 ! 2D field V==0
              W =0
              DO i=1,2*gamma + 1 
                    U = U + zStencil(i) * uTemp1D(i)
                    W = W + zStencil(i) * wTemp1D(i)
              END DO


          ELSE ! Velocity field in 3D

          ! Allocation of temporary fields
          !

           ALLOCATE(uTemp1D(2*gamma+1))
           ALLOCATE(uTemp(2*beta+1))

           ALLOCATE(vTemp1D(2*gamma+1))
           ALLOCATE(vTemp(2*beta+1))

           ALLOCATE(wTemp1D(2*gamma+1))
           ALLOCATE(wTemp(2*beta+1))

           ! Interpolation in x and z direction. Temporary filed in y
           !
          
           DO j = 1,2*beta+1                      

              counterZ = 1
              DO k=-left,right

                  uTemp1D(counterZ) = 0
                  vTemp1D(counterZ) = 0
                  wTemp1D(counterZ) = 0

                  DO i=1,2*alpha + 1 
                      uTemp1D(counterZ) = uTemp1D(counterZ)+xStencil(i)&
                      * UOF(NN(3)+k,NN(1)-alpha-1+i,NN(2)-beta-1+j)

                      vTemp1D(counterZ) = vTemp1D(counterZ)+xStencil(i)&
                      * VOF(NN(3)+k,NN(1)-alpha-1+i,NN(2)-beta-1+j)

                      wTemp1D(counterZ) = wTemp1D(counterZ)+xStencil(i)&
                      * WOF(NN(3)+k,NN(1)-alpha-1+i,NN(2)-beta-1+j)
                  END DO
              counterZ = counterZ + 1

              END DO


              ! Calculate u-velocity at x0 (z-interpolation)
              !
              uTemp(j) =0
              vTemp(j) =0
              wTemp(j) =0
              DO i=1,2*gamma + 1 
                    uTemp(j) = uTemp(j) + zStencil(i) * uTemp1D(i)
                    vTemp(j) = vTemp(j) + zStencil(i) * vTemp1D(i)
                    wTemp(j) = wTemp(j) + zStencil(i) * wTemp1D(i)
              END DO

           END DO

           ! Calculate u-velocity at x0 (y-interpolation)
           !
           U = 0
           V = 0
           W = 0
           DO i=1,2*beta + 1 
                 U = U + yStencil(i) * uTemp(i)
                 V = V + yStencil(i) * vTemp(i)
                 W = W + yStencil(i) * wTemp(i)
           END DO
            

 
      END IF 
      ELSE ! x0 outside water column
              U = 0;
              V = 0;
              W = 0;
      ENDIF
      END SUBROUTINE openFoamInterface_U

