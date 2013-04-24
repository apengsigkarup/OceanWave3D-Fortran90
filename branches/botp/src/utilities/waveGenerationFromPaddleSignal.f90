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
      REAL(KIND=long) :: omega, dz, dy_local,FAC,flux,y0,error, &
      tWeight,tmp
      REAL(KIND=long) :: FluxInterp(wave3DFlux%order)
      REAL(KIND=long), PARAMETER :: small = 1e-1



      ! Nearest neighbour interpolation in time !TODO: implement higher order
      ! scheme.
      itime = floor(time/(wave3DFlux%dt))+1

      !itime = (time-zero)/dt0+1

      IF (MOD(time,wave3DFlux%dt)/wave3DFlux%dt<=1/100) THEN
            tWeight = 0;
      ELSE
            tWeight = time/wave3DFlux%dt-itime+1
      ENDIF

      ! Linear ramping in time
      !
      IF (time<wave3DFlux%rampTime) THEN
        FAC = -two/(wave3DFlux%rampTime**3)*time**3+three/(wave3DFlux%rampTime**2)*time**2
      ELSE
        FAC = one
      ENDIF

      ! We assume that the wave gauges are uniformly spaced
      !
      dy_local = wave3DFlux%y(1)*2

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

          FluxInterp = wave3DFlux%flux(itime,NNy:NNy+wave3DFlux%order-1)
          CALL polint(wave3DFlux%y(NNy),FluxInterp,wave3DFlux%order,y0,flux,error)
        
          IF (tWeight/=0) THEN 
              tmp = flux
              FluxInterp = wave3DFlux%flux(itime+1,NNy:NNy+wave3DFlux%order-1)
              CALL polint(wave3DFlux%y(NNy),FluxInterp,wave3DFlux%order,y0,flux,error)

              flux = (1-tWeight)*flux + tmp*tWeight
          ENDIF
          ! Loop over the vertical direction and copy into global array
          DO k = 1+GhostGridZ,FineGrid%Nz + GhostGridZ
                  Uneumann(k,j) = -FAC*flux
          END DO
      END DO
      END SUBROUTINE
