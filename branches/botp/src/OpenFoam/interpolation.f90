      SUBROUTINE interpolation3D(C,x0,x,r,d)
      ! 3 dimensional interpolation kernel.
      !
      ! Inputs:
      !         x0: array with dimension(3) [x0,y0,z0] 
      !         x : array with dimensions(1:r^3,1:3) 
      !             x=[x1, x2,..;y1,y2,..;z1,z2,..]^T
      !         r : size of stencil (INTEGER)
      !         d: characteristic lengthscale, used for scaling in
      !         order to keep the matrix well conditioned 
      !
      ! Output: 
      !         C: A interpolation kernel with dimension(r^3:r^3)
      !         The stencil is assumed to be logically square so that the grid points 
      !         and derivatives start at the lower left-hand corner and follow 
      !         coordinate 2 (y) first and then coordinate 1 (x).  The coefficients 
      !         are thus defined such that 
      !
      !  c(k,1:r^3), k=(m-1)*r+n; corresponds to d^{m+n}/dx^m/dy^n.
      !   eg.:r=2 [f0 df/dy df/dx d^2f/dxdy df/dz df^2/dzdy df^2]
      !
      !
      ! Written by Bo Terp Paulsen, botp@mek.dtu.dk
      !

      ! Inputs / outputs
      !
      USE Interpolation_functions ! Decleration of matrix inversion,
      ! and factorial.
      IMPLICIT NONE
      INTEGER, PARAMETER :: dp = selected_real_kind(15, 307)
      INTEGER :: r,count1,count2
      REAL(KIND=dp), DIMENSION(:,:) :: x(r**3,3),C(r**3,r**3),mat(r**3,r**3)
      REAL(KIND=dp), DIMENSION(:) :: x0(3) 
      REAL(KIND=dp) :: d 
      

      ! Local variables
      !
      INTEGER :: i,j,k,m,n,b
      REAL(KIND=dp) :: factor

      !
      !
      factor = 1/d
      
      ! Generate Taylor coefficient matrix
      !
      count1 = 0
      count2 = 0
      DO i = 1,r ! x
        DO j = 1,r ! y 
          Do k = 1,r ! z
          count1 = count1 +1 
                DO m = 1,r ! x
                  DO n = 1,r ! y
                    DO b = 1,r ! z
                      count2 = count2 +1 
                mat(count1,count2) = ((factor*(x(count1,1)-x0(1)))**(m-1)/fact(m-1))* &
                                     ((factor*(x(count1,2)-x0(2)))**(n-1)/fact(n-1))* &
                                     ((factor*(x(count1,3)-x0(3)))**(b-1)/fact(b-1))
                    END DO
                  END DO
                END DO
      
           count2 = 0
          END DO
        END DO
      END DO

        C = inv(mat)

      ! Rescaling of coefficients
      !
      count1 = 0

      DO i = 1,r ! x
        DO j = 1,r ! y 
          Do k = 1,r ! z
          count1 = count1 +1
          DO count2 =1,r**3
          C(count1,count2) = factor**(i-1)*factor**(j-1)*factor**(k-1)*C(count1,count2)
          END DO

          END DO
        END DO
      END DO
      END SUBROUTINE interpolation3D

      SUBROUTINE interpolation2D(C,x0,x,r,d)
      ! 2 dimensional interpolation kernel.
      !
      ! Inputs:
      !         x0: array with dimension(3) [x0,y0,z0], only x0 and y0
      !         is used
      !         x=[x1, x2,..;y1,y2,..]^T
      !         r : size of stencil (INTEGER)
      !         dx: characteristic lengthscale, used for scaling in
      !         order to keep the matrix well conditioned 
      !
      ! Output: 
      !         C: A interpolation kernel with dimension(r^2:r^2)
      !         The stencil is assumed to be logically square so that the grid points 
      !         and derivatives start at the lower left-hand corner and follow 
      !         coordinate 2 (y) first and then coordinate 1 (x).  The coefficients !
      !         are 
      !         thus defined such that 
      !
      !  c(k,1:r^3), k=(m-1)*r+n; corresponds to d^{m+n}/dx^m/dy^n.
      !   eg.:r=2 [f0 df/dy df/dx d^2f/dxdy df/dz df^2/dzdy df^2]
      !
      !
      ! Written by Bo Terp Paulsen, botp@mek.dtu.dk
      !

      ! Inputs / outputs
      !
      USE Interpolation_functions
      IMPLICIT NONE
      INTEGER, PARAMETER :: dp = selected_real_kind(15, 307)
      INTEGER :: r
      REAL(KIND=dp), DIMENSION(:,:) :: x(r**2,2),C(r**2,r**2),mat(r**2,r**2)
      REAL(KIND=dp), DIMENSION(:) :: x0(3) 
      REAL(KIND=dp) :: d 
      

      ! Local variables
      !
      INTEGER :: i,j,m,n,count1,count2
      REAL(KIND=dp) :: factor

      ! Scaling factor
      !
      factor = 1/d
      
      ! Generate Taylor coefficient matrix
      !
      count1 = 0
      count2 = 0
      DO i = 1,r 
        DO j = 1,r  
          count1 = count1 +1 
                DO m = 1,r 
                  DO n = 1,r 
                      count2 = count2 +1 
           mat(count1,count2) = ((factor*(x(count1,1)-x0(1)))**(m-1)/fact(m-1))* &
                                ((factor*(x(count1,2)-x0(2)))**(n-1)/fact(n-1))
                  END DO
                END DO
           count2 = 0
        END DO
      END DO
        C = inv(mat)

      ! Rescaling of coefficients
      !
      count1 = 0

      DO i = 1,r 
        DO j = 1,r  
          count1 = count1 +1
          DO count2 =1,r**2
          C(count1,count2) = factor**(i-1)*factor**(j-1)*C(count1,count2)
          END DO

        END DO
      END DO
      END SUBROUTINE interpolation2D

      SUBROUTINE interpolation1D(C,x0,x,r)
      ! 1 dimensional interpolation kernel.
      !
      ! Inputs:
      !         x0: array with dimension(3) [x0,y0,z0], Only x0 is used
      !         x : array with dimensions(1:r,1) 
      !             x=[x1, x2,..]^T
      !         r : size of stencil (INTEGER)
      !
      ! Output: 
      !         C: A interpolation kernel with dimension(r:r)
      !
      !
      ! Written by Bo Terp Paulsen, botp@mek.dtu.dk
      !


      ! Inputs / outputs
      !
      USE Interpolation_functions
      IMPLICIT NONE
      INTEGER, PARAMETER :: dp = selected_real_kind(15, 307)
      INTEGER :: r
      REAL(KIND=dp), DIMENSION(:,:) :: C(r,r),mat(r,r)
      REAL(KIND=dp), DIMENSION(:) :: x0(3),x(r) 
      REAL(KIND=dp) :: d 
      

      ! Local variables
      !
      INTEGER :: i,m,count1, count2
      REAL(KIND=dp) :: factor

      ! Scaling factor
      !
      d = x(r)-x(r-1)
      factor = 1/d
      
      ! Generate Taylor coefficient matrix
      !
      count1 = 0
      count2 = 0
      DO i = 1,r 
          count1 = count1 +1 
                DO m = 1,r 
                      count2 = count2 +1 
          mat(count1,count2) = (factor*(x(count1)-x0(1)))**(m-1)/fact(m-1)
                END DO
          count2 = 0
      END DO
        C = inv(mat)

      ! Rescaling of coefficients
      !
      count1 = 0

      DO i = 1,r 
          count1 = count1 +1
          DO count2 =1,r
          C(count1,count2) = factor**(i-1)*C(count1,count2)
          END DO
      END DO
      END SUBROUTINE interpolation1D

      SUBROUTINE TaylorInterpolation&
      (dir,x0,stencilSize)
      ! Evaluating a Taylor polynomial at point x0, with expansion point
      !at nearest neighbour.
      !
      ! Inputs: 
      !        dir: Integer indication in which
      !        directions interpolation is desired. 
      !        x0 : Point to which the values are interpolated.
      !
      ! Output:
      !        TaylorStencil: Interpolation stencil
      !        inOrOut: Flag indicating however the cell is located in
      !        or outside the water column
      !
      ! Written by Bo Terp Paulsen, botp@mek.dtu.dk
      !
      USE Precision
      USE GlobalVariables
      USE Interpolation_functions !Declaration of factorial func. 
      IMPLICIT NONE

      ! Input variables
      !
      INTEGER :: dir
      REAL(KIND=long), DIMENSION(:) :: x0(3)
        
      ! Local variables
      !
      INTEGER :: stencilSize, i
      INTEGER, DIMENSION(:), POINTER :: NN
      REAL(KIND=long) :: sigma, dl
      REAL(KIND=long), DIMENSION(:), ALLOCATABLE :: dlVector
      REAL(KIND=long), DIMENSION(:,:) :: &
      transMAT(stencilSize,stencilSize)
      REAL(KIND=long), DIMENSION(:), POINTER :: xStencil, yStencil, zStencil

      ! Assign local pointers
      !
      xStencil => interpolation%stencilX 
      yStencil => interpolation%stencilY 
      zStencil => interpolation%stencilZ
      NN       => interpolation%NN


     
        IF(dir == 1) THEN
               ALLOCATE(dlVector(2*alpha +1))
               dl = x0(1)-FineGrid%x(NN(1),1)
               DO i = 0,2*alpha
                        dlVector(i+1) = dl**i / fact(i)
               END DO


               xStencil = &
               MATMUL(Interpolation%dx,dlVector)
       ENDIF
      
       IF(dir == 2) THEN
               ALLOCATE(dlVector(2*beta +1))
               dl = x0(2) - FineGrid%y(1,NN(2))
               DO i = 0,2*beta
                        dlVector(i+1) = dl**i / fact(i)
               END DO

               yStencil = &
               MATMUL(Interpolation%dy(:,:),dlVector)

       ENDIF

       IF(dir == 3) THEN
      ! Transformation of z-physical into sigma
      ! 
      sigma = (x0(3) ) / &
              (FineGrid%h(NN(1),NN(2)) + Wavefield%E(NN(1),NN(2)))

        
               ALLOCATE(dlVector(2*gamma + 1))
               dl = (sigma - FineGrid%z(NN(3))) 
               DO i = 0,2*gamma
                        dlVector(i+1) = dl**i / fact(i)
               END DO

               zStencil = &
               MATMUL(Interpolation%dz(:,:,NN(3)),dlVector)

       ENDIF
        
      END SUBROUTINE TaylorInterpolation
