SUBROUTINE Runge_Kutta_45(rhsFreeSurface)
! By Allan P. Engsig-Karup.
  USE Precision
  USE Constants
  USE GlobalVariables, ONLY: FineGrid, time, g, dt, GhostGridZ, PHI, &
      gamma, alpha, beta, RHS, rhsE,rhsP, relaxONOFF, LinearONOFF, filteringONOFF, filterALPHA, filterNP, &
	  filtercoefficients, tstep, GhostGridX, GhostGridY, swenseTransientTime, &
      swenseONOFF, swenseDir, SFsol, Wavefield, curvilinearONOFF, filtercoefficients2

  IMPLICIT NONE
  INTEGER :: INTRK
  REAL(KIND=long) :: RKtime
  REAL(KIND=long), DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY) :: resE, resP
  REAL(KIND=long), DIMENSION(6) :: rk4c
  REAL(KIND=long), DIMENSION(5) :: rk4a, rk4b
  EXTERNAL rhsFreeSurface, BuildLinearSystem, BuildLinearSystemTransformedCurvilinear

  ! Runge-Kutta residual storage
  resE = zero; resP = zero

  ! Low storage Runge-Kutta coefficients
  rk4a = (/ 0.0_long ,  -567301805773.0_long/1357537059087.0_long , &
        -2404267990393.0_long/2016746695238.0_long , -3550918686646.0_long/2091501179385.0_long , &
		-1275806237668.0_long/842570457699.0_long /)
  rk4b = (/ 1432997174477.0_long/9575080441755.0_long , 5161836677717.0_long/13612068292357.0_long , &
         1720146321549.0_long/2090206949498.0_long ,&
            3134564353537.0_long/4481467310338.0_long , 2277821191437.0_long/14882151754819.0_long /)
  rk4c = (/ 0.0_long , 1432997174477.0_long/9575080441755.0_long , 2526269341429.0_long/6820363962896.0_long, &
            2006345519317.0_long/3224310063776.0_long , 2802321613138.0_long/2924317926251.0_long , 1._long /)

  ! Runge-Kutta loop starts here
  DO INTRK = 1 , 5

 	RKtime = time + dt*rk4c(INTRK) ! FIXME: should I use INTRK+1 here due to arraysize of rk4c
    resE = rk4a(INTRK)*resE
    resP = rk4a(INTRK)*resP
  	! GD: in SWENSE case, computation of incident wavefield
  	IF (swenseONOFF/=0) THEN
    	IF (LinearONOFF==0) THEN
    	! First stage: no need to recalculate incident wavefield (linear only)
        CALL incident_linear_wf_finite(swenseDir, Wavefield, &
            FineGrid%Nx+2*GhostGridX,FineGrid%x,FineGrid%Ny+2*GhostGridY,FineGrid%y,&
            FineGrid%Nz+GhostGridz,FineGrid%z,FineGrid%h,RKtime,SFsol%k,g,SFsol%HH,SFsol%h)
    	ELSE
            CALL incident_wf_finite(swenseDir, Wavefield, &
	 		FineGrid%Nx+2*GhostGridX,FineGrid%x,FineGrid%Ny+2*GhostGridY,FineGrid%y,FineGrid%Nz+GhostGridz,FineGrid%z,&
            FineGrid%h,RKtime,SFsol%k,g,SFsol%n_four_modes,SFsol%zz,SFsol%yy)
    	ENDIF
  	ENDIF

	CALL rhsFreeSurface(RKtime,Wavefield,g,rhsE,rhsP,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)
    IF(GhostGridX+GhostGridY==2) THEN !GD: 3D case, no time advance of corners values (which are not useful)...
      CALL Zero_Corners(rhsE,rhsP,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)
    ENDIF

    resE = resE + dt*rhsE
    resP = resP + dt*rhsP

    ! Finish Runge-Kutta stage:
    Wavefield%E = Wavefield%E + rk4b(INTRK)*resE
    Wavefield%P = Wavefield%P + rk4b(INTRK)*resP

    ! Relaxation
    IF (relaxONOFF==1) THEN
        CALL RelaxationModule(Wavefield%E,Wavefield%P,RKtime)
    ENDIF

    ! Filtering
    IF (filteringONOFF>0) THEN
      IF (MOD(tstep,filteringONOFF)==0) THEN
!$$$$$$         !GD change we need to have correct boundary value before the computation of filtered datas...
!$$$$$$         CALL DifferentiationsFreeSurfacePlane(Wavefield,GhostGridX,GhostGridY,FineGrid,alpha,beta)
        ! GD: SWENSE changes, the filtering should be done on scattered+incident...
        ! This subroutine apply filtering on E and P with different specific treatment for SWENSE
        !CALL FILTERING(FineGrid%Nx,FineGrid%Ny,Wavefield%E,filterNP,filterALPHA,filtercoefficients,tstep)
        !CALL FILTERING(FineGrid%Nx,FineGrid%Ny,Wavefield%P,filterNP,filterALPHA,filtercoefficients,tstep)
!$$$$$$         CALL FILTERING_NEW(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield,filterNP,filterALPHA,filtercoefficients,&
!$$$$$$             tstep,swenseONOFF,filtercoefficients2, GhostGridX, GhostGridY)
        ! GD : test the filtering assuming symetry along boundaries...
		!CALL FILTERING_SWENSE(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%E_I,filterNP,filterALPHA, &
        !	filtercoefficients,tstep,GhostGridX,GhostGridY)
		!CALL FILTERING_SWENSE(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%P,Wavefield%P_I_s,filterNP,filterALPHA, &
        !	filtercoefficients,tstep,GhostGridX,GhostGridY)
        ! Change the filtering after time-ramp
        IF((time-dt.GE.swenseTransientTime).AND.(time-3*dt.LE.swenseTransientTime)) THEN
			! If needed change the filtering after the time ramp...
        ENDIF
		! GD : filter only scattered wavefield (use 0 instead of swenseONOFF)
        CALL FILTERING_NEW(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield,filterNP,filterALPHA,filtercoefficients,&
            tstep,0,filtercoefficients2, GhostGridX, GhostGridY)

	  ENDIF
    ENDIF

    CALL DifferentiationsFreeSurfacePlane(Wavefield,GhostGridX,GhostGridY,FineGrid,alpha,beta)

    IF (LinearONOFF==1) THEN
!        CALL DetermineTransformationConstants(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
!		     FineGrid%Nz+GhostGridZ,FineGrid,FineGrid%dsigma,Wavefield)
        ! GD: test !FIXME: check that dsigma is not used in Curvilinear (i.e. just dsigmanew) NO keep both VerticalFreeSurfaceVelocity (to optimize)
 !       IF (curvilinearONOFF==1) THEN
            CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
            	FineGrid,FineGrid%dsigmanew,Wavefield)
  !      ENDIF
	ENDIF

	RHS(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY) = &
	     Wavefield%P(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)

    IF (curvilinearONOFF==1) THEN
		CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
		    BuildLinearSystemTransformedCurvilinear,PHI,FineGrid)
    ELSE
		CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
			BuildLinearSystem,PHI,FineGrid)
    END IF

!    CALL VerticalFreeSurfaceVelocity(Wavefield%W,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
!	     PHI,FineGrid%DiffStencils,FineGrid%dsigma(:,5), gamma)
    CALL VerticalFreeSurfaceVelocity(Wavefield%W,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
	     PHI,FineGrid%DiffStencils,FineGrid%dsigmanew(:,:,:,5), gamma)

  END DO

END SUBROUTINE Runge_Kutta_45
