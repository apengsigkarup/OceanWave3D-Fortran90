      SUBROUTINE OpenFoamNearestNeighbourXY(x0)
      ! Nearest Neighbour algorithm: Runs through a structured Cartesian grid in the $x$, $y$ and $z$ direction respectively. 
      ! 
      ! Input: x0:= vector which contains the coordinate $\{x,y,z\}$ of
      ! which the nearest neighbour is wanted.
      !
      ! Output: 
      !

      USE GlobalVariables
      IMPLICIT NONE
      REAL(KIND=long) :: dist, temp
      REAL(KIND=long), DIMENSION(:) :: x0(3) 
      INTEGER i
      INTEGER, DIMENSION(:), POINTER :: NN
      
      ! Assign local pointers
      !
      NN => interpolation%NN
      
      ! Searches the x-direction (positive direction)
      !
      IF (GridX==1) THEN
          dist = 1*10**5 ! Large number infinite distance
          DO i=1,SIZE(FineGrid%x,1) 
            temp = ABS(x0(1) - FineGrid%x(i,1))
                   IF (temp < dist) THEN
                     dist = temp
                     NN(1) = i
                   ELSE 
                     exit
                   END IF
            
          END DO
      ELSE
          NN(1) = nint(x0(1)/dx) + GhostGridX + 1
      END IF
      
      ! Searches the y-direction (positive direction)
      !
      IF(FineGrid%Ny==1) THEN
      NN(2) = 1
      ELSE IF (GridY==0 .AND. FineGrid%Ny>1) THEN
          NN(2) = nint(x0(2)/dy) + GhostGridY + 1
      ELSE
          dist = 1*10**5
          DO i=1,SIZE(FineGrid%y,2) 
            temp = ABS(x0(2) - FineGrid%y(1,i))
                   IF (temp < dist) THEN
                     dist = temp
                     NN(2) = i
                   ELSE
                           exit
                  END IF
            
          END DO
      END IF

      END SUBROUTINE OpenFoamNearestNeighbourXY


      SUBROUTINE OpenFoamNearestNeighbourZ(x0)
      !Nearest Neighbour algorithm: Runs through a structured Cartesian 
      !grid in the z direction respectively. 
      ! 
      ! Input: x0:= vector which contains the coordinate $\{x,y,z\}$ of
      ! which the nearest neighbour is wanted.
      !
      ! Output: 
      !

      USE GlobalVariables
      IMPLICIT NONE
      REAL(KIND=long) :: dist, temp, sigma
      REAL(KIND=long), DIMENSION(:) :: x0(3) 
      INTEGER i
      INTEGER, POINTER :: inOrOut
      INTEGER, DIMENSION(:), POINTER :: NN


      ! Assign local pointers
      !
      NN => interpolation%NN
      inOrOut => interpolation%inOrOut

      ! Transformation of z-physical into sigma
      ! 
      sigma = (x0(3) ) / &
              (FineGrid%h(NN(1),NN(2)) + Wavefield%E(NN(1),NN(2)))

      ! Identify however the cell is located in or outside the water
      ! column.
      IF(sigma>1.00) THEN 
       
        !Flag indicating the cell is outside the water column 
        inOrOut = 1 
        
      ELSE
        ! Flag indicating that the cell is inside the water column
        inOrOut = 0 
      ! Searches the $z$-direction (positive direction) 
      !
      dist = 1*10**5
      DO i=1,SIZE(FineGrid%z,1)
        temp = ABS(sigma - FineGrid%z(i))
               IF (temp < dist) THEN
                 dist = temp
                 NN(3) = i
               ELSE 
                 exit
               END IF
        
      END DO

      END IF
      END SUBROUTINE OpenFoamNearestNeighbourZ

      !SUBROUTINE NearestNeighbourXY(x0,NN)
      !! Nearest Neighbour algorithm: Runs through a structured Cartesian grid in the $x$, $y$ and $z$ direction respectively. 
      !! 
      !! Input: x0:= vector which contains the coordinate $\{x,y,z\}$ of
      !! which the nearest neighbour is wanted.
      !!
      !! Output: 
      !!

      !USE GlobalVariables
      !IMPLICIT NONE
      !REAL(KIND=long) :: dist, temp
      !REAL(KIND=long), DIMENSION(:) :: x0(3) 
      !INTEGER i
      !INTEGER, DIMENSION(:) :: NN(2)
      
      
      !! Searches the $x$-direction (positive direction)
      !!
      !dist = 1*10**5 ! Large number infinite distance
      !DO i=1,SIZE(FineGrid%x,1) 
        !temp = ABS(x0(1) - FineGrid%x(i,1))
               !IF (temp < dist) THEN
                 !dist = temp
                 !NN(1) = i
               !ELSE 
                 !exit
               !END IF
        
      !END DO
      
      !! Searches the $y$-direction (positive direction)
      !!
      !dist = 1*10**5
      !DO i=1,SIZE(FineGrid%y,2) 
        !temp = ABS(x0(2) - FineGrid%y(1,i))
               !IF (temp < dist) THEN
                 !dist = temp
                 !NN(2) = i
               !ELSE
                       !exit
               !END IF
        
      !END DO

      
      !END SUBROUTINE NearestNeighbourXY

