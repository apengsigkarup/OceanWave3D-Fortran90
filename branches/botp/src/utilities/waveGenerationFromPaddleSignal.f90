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
      ! Local variables
      !
      INTEGER :: j, k,NNY,itime
      REAL(KIND=long) :: omega, dz, dy_local,FAC,flux,y0,error
      REAL(KIND=long) :: FluxInterp(wave3DFlux%order)


      ! Nearest neighbour interpolation in time !TODO: implement higher order
      ! scheme.
      itime = (time-time0)/dt0+1

      ! Linear ramping in time
      !
      IF (time<wave3DFlux%rampTime) THEN
        FAC = time/4
      ELSE
        FAC = 1.0
      ENDIF

      ! We assume that the wave gauges are uniformly spaced
      !
      dy_local = wave3DFlux%y(2) - wave3DFlux%y(1)

      ! Then we loop over all cells in the y-direction
      !
      DO j =1+GhostGridY,FineGrid%Ny + GhostGridY

          ! We interpolate from wave gauge grid to fineGrid
          !
          y0 = FineGrid%y(1,j)
          NNY = nint(y0/dy_local)+1-wave3DFlux%order/2

          ! One sided stencils at the boundaries
          !
          IF(NNY<1) THEN
            NNY = 1
          ELSEIF((NNY+wave3DFlux%order-1)>wave3DFlux%n2) THEN
            NNY = wave3DFlux%n2 - wave3DFlux%order + 1
          ENDIF

          FluxInterp = wave3DFlux%flux(tstep,NNy:NNy+wave3DFlux%order-1)

          CALL polint(wave3DFlux%y(NNy),FluxInterp,wave3DFlux%order,y0,flux,error)


          ! Loop over the vertical direction and copy into global array
          DO k = 1+GhostGridZ,FineGrid%Nz + GhostGridZ
                  Uneumann(k,j) = FAC*flux    
          END DO
      END DO
      END SUBROUTINE
