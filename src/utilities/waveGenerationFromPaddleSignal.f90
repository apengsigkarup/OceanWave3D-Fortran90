SUBROUTINE waveGenerationFromPaddleSignal(RKtime)
  !
  ! 3-dimensional wave generation via a Western boundary flux.
  ! This is implemented via inhomogeneous Neumann boundary conditions
  ! added to the RHS vector in subroutine RungeKutta44.f90. In the 
  ! nonlinear case, source terms are also added to the update of the 
  ! ghost points for eta and phi on the free surface. 
  !
  ! Inputs:
  !         RKtime
  !
  ! Outputs:
  !         Uneumann:= Velocity array size(NZ,NY) (Global variable)
  !         wavefield%SourceEx:= Velocity array size(1,NY) (Global variable)
  !         wavefield%SourcePx:= Velocity array size(1,NY) (Global variable)
  !
  !
  ! Written by Bo Terp Paulsen, botp@mek.dtu.dk
  !
  ! Modified Feb. 2016 by Harry Bingham, hbb@mek.dtu.dk.
  !
  !  Added linear interpolation between the enclosing two time steps 
  !  and placed the inhomogeneous conditions for phi^S_x and eta_x into 
  !  the source terms for the ghost point update in UpdateGhostLayerX.f90. 
  !  Also moved the insertion of Uneumann to the RHS in RungeKutta44.f90 
  !  from the matrix-vector product of BuildLinearSystem.f90. 
  !
  ! Inputs / outputs
  !
  USE Precision
  USE GlobalVariables
  IMPLICIT NONE
  !
  Real(kind=long) :: RKtime
  !
  ! Local variables
  !
  INTEGER :: j, k, NNY, itime
  REAL(KIND=long) :: omega, dz, dy_local, FAC, flux, flux0, flux1, y0, error,  &
       tWeight(2), dt_local(2), etax, etax0, etax1
  REAL(KIND=long) :: FluxInterp(wave3DFlux%order)
  !
  ! Find the nearest boundary flux index less than the current time. 
  ! 
  itime = floor(RKtime/(wave3DFlux%dt))+1
  !
  ! Compute the linear interpolation wieghts for the two boundary flux  
  ! time-levels that bracket the current time level. 
  !
  dt_local(1)=rktime-(itime-1)*wave3DFlux%dt
  dt_local(2)=(itime)*wave3DFlux%dt-rktime
  tWeight(1)=dt_local(2)/(dt_local(1)+dt_local(2))
  tWeight(2)=dt_local(1)/(dt_local(1)+dt_local(2))
  !
  ! Linear ramping in time
  !
  IF (RKtime<wave3DFlux%rampTime) THEN
     FAC = -two/(wave3DFlux%rampTime**3)*RKtime**3+three/(wave3DFlux%rampTime**2)*RKtime**2
  ELSE
     FAC = one
  ENDIF
  !
  ! The wave gauges are uniformly spaced in y (this has been checked in setupwavepaddle.f90)
  !
  dy_local = wave3DFlux%y(1)*2

  ! Then we loop over all the fluid cells in the y-direction
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
     !
     ! Interpolate in y at time step itime
     !
     FluxInterp = wave3DFlux%flux(itime,NNy:NNy+wave3DFlux%order-1)
     CALL polint(wave3DFlux%y(NNy),FluxInterp,wave3DFlux%order,y0,flux0,error)
     FluxInterp = wave3DFlux%etax(itime,NNy:NNy+wave3DFlux%order-1)
     CALL polint(wave3DFlux%y(NNy),FluxInterp,wave3DFlux%order,y0,etax0,error)
     !
     ! Interpolate in y at time step itime+1
     !
     FluxInterp = wave3DFlux%flux(itime+1,NNy:NNy+wave3DFlux%order-1)
     CALL polint(wave3DFlux%y(NNy),FluxInterp,wave3DFlux%order,y0,flux1,error)
     FluxInterp = wave3DFlux%etax(itime+1,NNy:NNy+wave3DFlux%order-1)
     CALL polint(wave3DFlux%y(NNy),FluxInterp,wave3DFlux%order,y0,etax1,error)
     !
     ! Linearly interpolate between the two time values to get the value at RKtime.
     !
     flux = tWeight(1)*flux0 + tWeight(2)*flux1
     etax = tWeight(1)*etax0 + tWeight(2)*etax1

     ! Loop over the vertical direction and copy into the global flux array
     DO k = 1+GhostGridZ,FineGrid%Nz + GhostGridZ
        Uneumann(k,j) = FAC*flux 
     END DO
     ! Add the Neumann conditions on the free surface the source terms for the 
     ! Western boundary ghost points 
     wavefield%sourcePx(1+GhostGridX,j) = Uneumann(FineGrid%Nz+GhostGridZ,j)
     wavefield%sourceEx(1+GhostGridX,j) = etax
  END DO
END SUBROUTINE waveGenerationFromPaddleSignal
