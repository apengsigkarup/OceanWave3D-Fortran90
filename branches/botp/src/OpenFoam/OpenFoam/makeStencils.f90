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
        CALL OpenFoamNearestNeighbourXY(x0,NN)
        CALL OpenFoamNearestNeighbourZ(x0,NN,inOrOut)


        ! Evaluation of local interpolation stencils
        !
        CALL TaylorInterpolation(1,x0,xStencil,NN,SIZE(xStencil))
        CALL TaylorInterpolation(2,x0,yStencil,NN,SIZE(yStencil))
        CALL TaylorInterpolation(3,x0,zStencil,NN,SIZE(zStencil))
        END SUBROUTINE makeStecils


