        SUBROUTINE makeStecils(x0)
        ! Subroutine finding nearest neighbour and calculating local
        ! interpolation stencils

        !Input:
        !       x0: Array with dimension(3) [x y z]                         
        !
        !Output: void
        !
        ! Written by Bo Terp Paulsen (botp), botp@mek.dtu.dk
        !

        USE GlobalVariables
 
        !Input parameters
        !
        REAL(KIND=long), DIMENSION(:) :: x0(3)

        ! Local variables
        !
        INTEGER, POINTER :: inOrOut
        INTEGER, DIMENSION(:), POINTER :: NN
        REAL(KIND=long), DIMENSION(:) :: x0Local(3)
        REAL(KIND=long), DIMENSION(:), POINTER :: xStencil, yStencil, zStencil

        ! Assign local pointers
        !
        xStencil => interpolation%stencilX 
        yStencil => interpolation%stencilY 
        zStencil => interpolation%stencilZ
        NN       => interpolation%NN
        inOrOut  => interpolation%inOrOut
        
        ! Find Nearest Neighbour 
        !
        CALL OpenFoamNearestNeighbourXY(x0)
        ! We make a coordinate transformation from physical domain to
        ! OCW3D domain such that $x0(3) \in [0;eta]
        x0Local(1) = x0(1)
        x0Local(2) = x0(2)
        x0Local(3) = x0(3) + FineGrid%h(NN(1),NN(2)) 
        CALL OpenFoamNearestNeighbourZ(x0Local)

        ! Evaluation of local interpolation stencils
        !
        CALL TaylorInterpolation(1,x0Local,SIZE(xStencil))
        CALL TaylorInterpolation(2,x0Local,SIZE(yStencil))
        IF(inOrOut .EQ. 0) THEN
        CALL TaylorInterpolation(3,x0Local,SIZE(zStencil))
        ELSE
        NN(3) = -1000
        END IF
        END SUBROUTINE makeStecils


