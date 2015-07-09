SUBROUTINE interpolationInitialize()
! Subroutine initializing interpolation kernels for all cells
! in OpenFOAM relaxation zones. 
! This function should be called from OpenFOAM as a pre-process. 
!
! Input : 
!
! Output : The struct Interpolation defined in GlobalVariables should be
!          initialized. 
!
! Written by Bo Terp Paulsen (botp), botp@mek.dtu.dk

USE Precision
USE GlobalVariables
IMPLICIT NONE

! Local variables
!
!INTEGER :: STAT
INTEGER :: Nx, Ny, Nz, order, Nzp, i, rankx, ranky, rankz, x, y, z, &
left, right
INTEGER, DIMENSION(:), ALLOCATABLE :: idx
REAL(KIND=long), DIMENSION(:) :: x0(3)
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: xStencil, yStencil, zStencil
! Temporary pointers/variables
!
Nx = FineGrid%Nx
Ny = FineGrid%Ny
Nz = FineGrid%Nz

! Rank of stencils
!
rankx = 2*alpha+1
ranky = 2*beta+1
rankz = 2*gamma+1

! Allocation of local stencils
!
ALLOCATE( xStencil(rankx,rankx) )
ALLOCATE( yStencil(ranky,ranky) )
ALLOCATE( zStencil(rankz,rankz) )

! Expansion point
!The equidistant grid imply that derivatives for any interior
!point are similar. Notice that values are updated in if statements below
!
x0 = zero


! Initialization of x-derivative. 
!It is assumed that OceanWave3D mesh always is
!larger than the OpenFOAM grid, so no interpolation near the boundaries is
!necessary. 
!
ALLOCATE(Interpolation%dx(rankx,rankx),STAT = STAT)
IF (STAT/=0) THEN
	PRINT*,'Error: Could not allocate memory for interpolation stencil (X).'
	STOP
ENDIF
IF (Nx>1) THEN
        x = GhostGridX + alpha + 2
        x0(1) = FineGrid%x(x,1)! Expansion point, only x0(1) is used!

        !CALL interpolation1D(xStencil,x0,FineGrid%x(x - alpha :x + alpha,1),rankx)
        
        CALL  weights(x0,FineGrid%x(x-alpha:x+alpha,1),2*alpha,2*alpha,2*alpha,xStencil)
        Interpolation%dx = xStencil

        !CALL weights(x0,FineGrid%x(x - alpha :x+alpha,1), &
                !2*alpha,2*alpha+1,2*alpha,xStencil)
        
ENDIF

! Initialization of y-derivative. 
!It is assumed that OceanWave3D mesh always is
!larger than the OpenFOAM grid, so no interpolation near the boundaries is
!necessary. 
!
ALLOCATE(Interpolation%dy(ranky,ranky),STAT = STAT)
Interpolation%dy = zero
IF (STAT/=0) THEN
	PRINT*,'Error: Could not allocate memory for interpolation stencil (Y).'
	STOP
ENDIF
IF (Ny>1) THEN
        
        y = GhostGridY + beta + 2 
        x0(1) = FineGrid%y(1,y)! Expansion point, only x0(1) is used!

        !CALL interpolation1D(yStencil,x0,FineGrid%y(1,y - beta : y + beta),ranky)
         
        CALL  weights(x0,FineGrid%y(1,y-beta:y+beta),2*beta,2*beta,2*beta,yStencil)
        Interpolation%dy = yStencil
ENDIF

! Initialization of z-derivatives
! This part is copied from PreProcessDiffStencilsZ.f90!!!
!
!One-sided stencils near free surface.
!One-side stencils near bottom, and if ghost points
!are supported then the directional ghost point
!is included in stencil.
!
IF (GhostGridZ==0) THEN; Nzp = Nz; ELSE; Nzp = Nz+1; ENDIF
ALLOCATE( Interpolation%dz(rankz,rankz,Nzp),STAT=STAT)
IF (STAT/=0) THEN
	PRINT*,'Error: Could not allocate memory for differentiation stencil (Z).'
	STOP
ENDIF



DO i = 1, Nzp
! Expansion point 
!
        x0(1) = FineGrid%z(i)

! One sided schemes at boundaries
! 
        IF(i<=gamma) THEN
            left = i - 1 
            right = 2*gamma -left
        ELSE
                IF (i>=Nzp-gamma) THEN
                    right = Nzp - i 
                    left = 2*gamma - right
                ELSE
                    left = gamma
                    right = gamma
                END IF

        END IF
! Calculate derivativatives at x0
!
        !CALL interpolation1D(zStencil,x0,FineGrid%z(i - left : i + right),rankz)

        CALL  weights(x0,FineGrid%z(i-left:i+right),2*gamma,2*gamma,2*gamma,zStencil)

! Save interpolation stencil
!
        Interpolation%dz(:,:,i) = zStencil

END DO

! Allocation of local interpolation stencils
        ALLOCATE(interpolation%stencilX(2*alpha+1))
        ALLOCATE(interpolation%stencilY(2*beta+1))
        ALLOCATE(interpolation%stencilZ(2*gamma+1))

        ALLOCATE(interpolation%NN(3))
        ALLOCATE(interpolation%inOrOut)

END SUBROUTINE interpolationInitialize

