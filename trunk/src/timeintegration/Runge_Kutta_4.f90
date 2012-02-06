SUBROUTINE Runge_Kutta_4(rhsFreeSurface)
  USE Precision
  USE Constants
  USE DataTypes
  USE GlobalVariables, ONLY: FineGrid, time, g, dt, GhostGridZ, PHI, gamma, alpha, beta, &
       RHS, relaxONOFF, LASTPHI, tstep, extrapolationONOFF, LinearONOFF, filteringONOFF, &
       filterALPHA, filterNP, filtercoefficients, GhostGridX, GhostGridY,                &
       swenseONOFF, swenseDir, SFsol, Wavefield , curvilinearONOFF, Wavefield_tmp,       &
       filtercoefficients2, BreakMod, fileop
  IMPLICIT NONE
  !
  EXTERNAL rhsFreeSurface, BuildLinearSystem, BuildLinearSystemTransformedCurvilinear, DiffXEven, DiffYEven, FILTERING, &
       ALLOCATE_Wavefield_Type,ASSIGN_IncidentWavefield,incident_wf_finite, &
       incident_linear_wf_finite,FILTERING_NEW
  REAL(KIND=long) , DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY) :: k1_E, k1_P, k2_E, k2_P, k3_E, k3_P, k4_E, k4_P
  REAL(KIND=long) :: RKtime
  INTEGER :: idxnew, idxold, idxlast
  !$$$$$$   TYPE(wavefield_FS) :: wavefield_tmp ! Temporary wavefield type...
  ! GD: allocate temporary structure of wavefield_FS type
  IF (.NOT.ASSOCIATED(Wavefield_tmp%E)) THEN
     !IF(time==dt)THEN
     CALL ALLOCATE_Wavefield_Type(Wavefield_tmp,FineGrid%Nx,FineGrid%Ny,FineGrid%Nz,GhostGridX,GhostGridy,GhostGridZ,swenseONOFF)
  ENDIF
  !
  ! Last computed solution for PHI
  idxnew  = MOD(tstep-1,2)+1
  idxlast = MOD(tstep  ,2)+1
  !  idxold  = MOD(tstep+1,3)+1
  LASTPHI(:,:,:,idxnew) = PHI

  ! STAGE 1
  ! GD: choice of rhs for SWENSE...
  IF (swenseONOFF/=0) THEN
     IF (LinearONOFF==0) THEN
        ! First stage: no need to recalculate incident wavefield (linear only), Wavefield and Wavefield_tmp are OK
     ELSE
        CALL incident_wf_finite(swenseDir, Wavefield, &
             FineGrid%Nx+2*GhostGridX,FineGrid%x,FineGrid%Ny+2*GhostGridY,FineGrid%y,FineGrid%Nz+GhostGridz,FineGrid%z,&
             FineGrid%h,time,SFsol%k,g,SFsol%n_four_modes,SFsol%zz,SFsol%yy)
        ! GD: Assign the incident wavefield part...
        CALL ASSIGN_IncidentWavefield(Wavefield_tmp, Wavefield, FineGrid%Nx+2*GhostGridX, FineGrid%Ny+2*GhostGridY)
     ENDIF
  ENDIF
  !
  ! Get a possible applied free-surface pressure
  !
  CALL funPressureTerm(time,g,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
       FineGrid,Wavefield)
  !
  ! Get the breaking terms if called for
  !
  IF (BreakMod%i_breaking>0 .and. FineGrid%ny==1) THEN
     CALL detect_breaking(fileop(14),FineGrid%Nx+2*GhostGridX,Wavefield,0)
  END IF

  CALL rhsFreeSurface(time,Wavefield,g,k1_E,k1_P,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)
  IF(GhostGridX+GhostGridY==2) THEN !GD: 3D case, no time advance of corners values (which are not useful)...
     CALL Zero_Corners(k1_E,k1_P,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)
  ENDIF

  Wavefield_tmp%E = Wavefield%E + half*dt*k1_E
  Wavefield_tmp%P = Wavefield%P + half*dt*k1_P
  RKtime = time+half*dt

  ! Relaxation
  IF (relaxONOFF==1) THEN
     ! CALL RelaxationModule(Wavefield_tmp%E,Wavefield_tmp%P,RKtime)
     CALL RelaxationModule_new(Wavefield_tmp%E,Wavefield_tmp%P,RKtime)
  ENDIF
  CALL DifferentiationsFreeSurfacePlane(Wavefield_tmp,GhostGridX,GhostGridY,FineGrid,alpha,beta)

  IF (LinearONOFF==1) THEN
     !	CALL DetermineTransformationConstants(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
     !	     FineGrid%Nz+GhostGridZ,FineGrid,FineGrid%dsigma,Wavefield_tmp)
     ! GD: test !FIXME: check that dsigma is not used in Curvilinear (i.e. just dsigmanew) NO keep both for VerticalFreeSurfaceVelocity (to optimize)
     !   IF (curvilinearONOFF==1) THEN
     CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
          FineGrid,FineGrid%dsigmanew,Wavefield_tmp)
     !  ENDIF
  ENDIF

  RHS(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY) = &
       Wavefield_tmp%P(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)

  IF (extrapolationONOFF==1) THEN ! Improving initial guess for iterative solver
     PHI = 1.5_long*LASTPHI(:,:,:,idxnew) - LASTPHI(:,:,:,idxlast)*half
     !    PHI = three*LASTPHI(:,:,:,idxnew) - three*LASTPHI(:,:,:,idxlast) + LASTPHI(:,:,:,idxold)
  ENDIF

  IF (curvilinearONOFF==1) THEN
     CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
          BuildLinearSystemTransformedCurvilinear,PHI,FineGrid)
  ELSE
     CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
          BuildLinearSystem,PHI,FineGrid)
  END IF

  CALL VerticalFreeSurfaceVelocity(Wavefield_tmp%W,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,PHI,&
       FineGrid%DiffStencils,FineGrid%dsigmanew(:,:,:,5), gamma)

  ! STAGE 2
  ! GD: choice of rhs for SWENSE...
  IF (swenseONOFF/=0) THEN ! Evaluation of incident wavefield
     IF (LinearONOFF==0) THEN
        CALL incident_linear_wf_finite(swenseDir, Wavefield_tmp, &
             FineGrid%Nx+2*GhostGridX,FineGrid%x,FineGrid%Ny+2*GhostGridY,FineGrid%y,FineGrid%Nz+GhostGridz,FineGrid%z,&
             FineGrid%h,time+half*dt,SFsol%k,g,SFsol%HH,SFsol%h)
     ELSE
        !print*,'This has to be done...'
        !PAUSE
        CALL incident_wf_finite(swenseDir, Wavefield_tmp, &
             FineGrid%Nx+2*GhostGridX,FineGrid%x,FineGrid%Ny+2*GhostGridY,FineGrid%y,FineGrid%Nz+GhostGridz,FineGrid%z,FineGrid%h,&
             time+half*dt,SFsol%k,g,SFsol%n_four_modes,SFsol%zz,SFsol%yy)
     ENDIF
     ! Assign incident wavefield values to Wavefield Needed in the buildlinear subroutines FIXME ! (define Wavefield as argument in interative_solution+buildlinear)
     CALL ASSIGN_IncidentWavefield(Wavefield, Wavefield_tmp, FineGrid%Nx+2*GhostGridX, FineGrid%Ny+2*GhostGridY)
  ENDIF
  !
  ! Get a possible applied free-surface pressure
  !
  CALL funPressureTerm(RKtime,g,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
       FineGrid,Wavefield_tmp)
  !
  ! Get the breaking terms if called for
  !
  IF (BreakMod%i_breaking>0 .and. FineGrid%ny==1) THEN
     CALL detect_breaking(fileop(14),FineGrid%Nx+2*GhostGridX,Wavefield_tmp,0)
  END IF
  !
  CALL rhsFreeSurface(time+half*dt,Wavefield_tmp,g,k2_E,k2_P,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)
  IF(GhostGridX+GhostGridY==2) THEN !GD: 3D case, no time advance of corners values (which are not useful)...
     CALL Zero_Corners(k2_E,k2_P,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)
  ENDIF

  Wavefield_tmp%E = Wavefield%E + half*dt*k2_E
  Wavefield_tmp%P = Wavefield%P + half*dt*k2_P
  RKtime = time+half*dt

  ! Relaxation
  IF (relaxONOFF==1) THEN
     !     CALL RelaxationModule(Wavefield_tmp%E,Wavefield_tmp%P,RKtime)
     CALL RelaxationModule_new(Wavefield_tmp%E,Wavefield_tmp%P,RKtime)
  ENDIF
  CALL DifferentiationsFreeSurfacePlane(Wavefield_tmp,GhostGridX,GhostGridY,FineGrid,alpha,beta)
  IF (LinearONOFF==1) THEN
     CALL DetermineTransformationConstants(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,FineGrid,&
          FineGrid%dsigmanew,Wavefield_tmp)
     ! GD: test !FIXME: check that dsigma is not used in Curvilinear (i.e. just dsigmanew) NO keep both for VerticalFreeSurfaceVelocity (to optimize)
     IF (curvilinearONOFF==1) THEN
        CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
             FineGrid,FineGrid%dsigmanew,Wavefield_tmp)
     ENDIF
  ENDIF

  RHS(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY) = &
       Wavefield_tmp%P(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)

  IF (curvilinearONOFF==1) THEN
     CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
          BuildLinearSystemTransformedCurvilinear,PHI,FineGrid)
  ELSE
     CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
          BuildLinearSystem,PHI,FineGrid)
  END IF

  CALL VerticalFreeSurfaceVelocity(Wavefield_tmp%W,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,PHI,&
       FineGrid%DiffStencils,FineGrid%dsigmanew(:,:,:,5), gamma)

  ! STAGE 3
  ! GD: choice of rhs for SWENSE...
  IF (swenseONOFF/=0) THEN
     IF (LinearONOFF==0) THEN
        ! Third stage: no need to recalculate incident wavefield (linear only)
     ELSE
        CALL incident_wf_finite(swenseDir, Wavefield_tmp, &
             FineGrid%Nx+2*GhostGridX,FineGrid%x,FineGrid%Ny+2*GhostGridY,FineGrid%y,FineGrid%Nz+GhostGridz,FineGrid%z,&
             FineGrid%h,time+half*dt,SFsol%k,g,SFsol%n_four_modes,SFsol%zz,SFsol%yy)
        ! Assign incident wavefield values to Wavefield Needed in the buildlinear subroutines FIXME ! (define Wavefield as argument in interative_solution+buildlinear)
        CALL ASSIGN_IncidentWavefield(Wavefield, Wavefield_tmp, FineGrid%Nx+2*GhostGridX, FineGrid%Ny+2*GhostGridY)
     ENDIF
  ENDIF
  !
  ! Get a possible applied free-surface pressure
  !
  CALL funPressureTerm(RKtime,g,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
       FineGrid,Wavefield_tmp)
  IF (BreakMod%i_breaking>0 .and. FineGrid%ny==1) THEN
     CALL detect_breaking(fileop(14),FineGrid%Nx+2*GhostGridX,Wavefield_tmp,0)
  END IF
  !
  CALL rhsFreeSurface(time+half*dt,Wavefield_tmp,g,k3_E,k3_P,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)	
  IF(GhostGridX+GhostGridY==2) THEN !GD: 3D case, no time advance of corners values (which are not useful)...
     CALL Zero_Corners(k3_E,k3_P,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)
  ENDIF

  Wavefield_tmp%E = Wavefield%E + dt*k3_E
  Wavefield_tmp%P = Wavefield%P + dt*k3_P
  RKtime = time+dt

  ! Relaxation
  IF (relaxONOFF==1) THEN
     !     CALL RelaxationModule(Wavefield_tmp%E,Wavefield_tmp%P,RKtime)
     CALL RelaxationModule_new(Wavefield_tmp%E,Wavefield_tmp%P,RKtime)
  ENDIF

  CALL DifferentiationsFreeSurfacePlane(Wavefield_tmp,GhostGridX,GhostGridY,FineGrid,alpha,beta)
  IF (LinearONOFF==1) THEN
     ! GD: New variable definition
     CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,FineGrid, &
          FineGrid%dsigmanew,Wavefield_tmp)
     ! GD: test !FIXME: check that dsigma is not used in Curvilinear (i.e. just dsigmanew) NO keep both for VerticalFreeSurfaceVelocity (to optimize)
     IF (curvilinearONOFF==1) THEN
        CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
             FineGrid,FineGrid%dsigmanew,Wavefield_tmp)
     ENDIF
  ENDIF

  RHS(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY) = &
       Wavefield_tmp%P(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)

  IF (extrapolationONOFF==1 .AND. tstep>1) THEN ! Improving initial guess for iterative solver
     !  PHI = (two*LASTPHI(:,:,:,idxnew) - LASTPHI(:,:,:,idxold))
     PHI = 2.66667_long*PHI - two*LASTPHI(:,:,:,idxnew) + 0.333333_long*LASTPHI(:,:,:,idxlast) !checked OK Lagrange interp.
     !    ! Work only on interior points ? (seems to be more stable...)
     !    PHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY) = &
     !        2.66667_long*PHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY) &
     !         - two*LASTPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,idxnew) &
     !         + 0.333333_long*LASTPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,idxlast)
     !    PHI = 3.2_long*PHI - three*LASTPHI(:,:,:,idxnew) + LASTPHI(:,:,:,idxlast) -0.2_long*LASTPHI(:,:,:,idxold)
  ENDIF

  IF (curvilinearONOFF==1) THEN
     CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
          BuildLinearSystemTransformedCurvilinear,PHI,FineGrid)
  ELSE
     CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
          BuildLinearSystem,PHI,FineGrid)
  END IF

  CALL VerticalFreeSurfaceVelocity(Wavefield_tmp%W,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ, &
       PHI,FineGrid%DiffStencils,FineGrid%dsigmanew(:,:,:,5), gamma)

  ! STAGE 4
  ! GD: choice of rhs for SWENSE...
  IF (swenseONOFF/=0) THEN ! Evaluation of incident wavefield
     IF (LinearONOFF==0) THEN
        CALL incident_linear_wf_finite(swenseDir, Wavefield_tmp, &
             FineGrid%Nx+2*GhostGridX,FineGrid%x,FineGrid%Ny+2*GhostGridY,FineGrid%y,FineGrid%Nz+GhostGridz,FineGrid%z,&
             FineGrid%h,time+dt,SFsol%k,g,SFsol%HH,SFsol%h)
     ELSE
        CALL incident_wf_finite(swenseDir, Wavefield_tmp, &
             FineGrid%Nx+2*GhostGridX,FineGrid%x,FineGrid%Ny+2*GhostGridY,FineGrid%y,FineGrid%Nz+GhostGridz,FineGrid%z,&
             FineGrid%h,time+dt,SFsol%k,g,SFsol%n_four_modes,SFsol%zz,SFsol%yy)
     ENDIF
     ! Assign incident wavefield values to Wavefield (to prevent from 1st RKstep computation + following computation)
     CALL ASSIGN_IncidentWavefield(Wavefield, Wavefield_tmp, FineGrid%Nx+2*GhostGridX, FineGrid%Ny+2*GhostGridY)
  ENDIF
  !
  ! Get a possible applied free-surface pressure
  !
  CALL funPressureTerm(RKtime,g,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
       FineGrid,Wavefield_tmp)
  IF (BreakMod%i_breaking>0 .and. FineGrid%ny==1) THEN
     CALL detect_breaking(fileop(14),FineGrid%Nx+2*GhostGridX,Wavefield_tmp,0)
  END IF
  !
  CALL rhsFreeSurface(time+dt,Wavefield_tmp,g,k4_E,k4_P,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)
  IF(GhostGridX+GhostGridY==2) THEN !GD: 3D case, no time advance of corners values (which are not useful)...
     CALL Zero_Corners(k4_E,k4_P,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)
  ENDIF
  !
  ! Shift the deta/dt history and put in the new value.  
  !
  Wavefield%EtatHist(:,:,3)=Wavefield%EtatHist(:,:,2);  Wavefield%EtatHist(:,:,2)=Wavefield%EtatHist(:,:,1);
  Wavefield%EtatHist(:,:,1)= 1/six*(k1_E+two*k2_E+two*k3_E+k4_E)
  ! Compute the new surface quantities
  Wavefield%E    = Wavefield%E + dt/six*(k1_E+two*k2_E+two*k3_E+k4_E)
  Wavefield%P    = Wavefield%P + dt/six*(k1_P+two*k2_P+two*k3_P+k4_P)
  RKtime = time+dt

  ! Relaxation
  IF (relaxONOFF==1) THEN
     !     CALL RelaxationModule(Wavefield%E,Wavefield%P,RKtime)
     CALL RelaxationModule_new(Wavefield%E,Wavefield%P,RKtime)
  ENDIF
  ! Filtering
  IF (filteringONOFF>0) THEN
     IF (MOD(tstep,filteringONOFF)==0) THEN
        IF (SWENSEOnOff==0) THEN
           !hbb I've turned off the boundary filtering unless the incident-scattered decomposition is applied.  
           CALL FILTERING(FineGrid%Nx,FineGrid%Ny,Wavefield%E,filterNP,filterALPHA,filtercoefficients, &
                tstep)
           CALL FILTERING(FineGrid%Nx,FineGrid%Ny,Wavefield%P,filterNP,filterALPHA,filtercoefficients, &
                tstep)
        ELSE
           ! GD: SWENSE changes, the filtering should be done on scattered+incident...
           ! This subroutine apply filtering on E and P with different specific treatment for SWENSE
           ! GD test to filter only scattered wavefield
           ! HBB: Note that the boundary filters here are hard coded to a 13-point stencil and 
           ! will surely not work correctly for other choices of filterNP.  
           CALL FILTERING_NEW(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield,filterNP, &
                filterALPHA,filtercoefficients,tstep,0,filtercoefficients2, GhostGridX, GhostGridY)
        END IF
     ENDIF
  ENDIF
  CALL DifferentiationsFreeSurfacePlane(Wavefield,GhostGridX,GhostGridY,FineGrid,alpha,beta)

  IF (LinearONOFF==1) THEN
     CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,FineGrid, &
          FineGrid%dsigmanew,Wavefield)
     ! GD: test !FIXME: check that dsigma is not used in Curvilinear (i.e. just dsigmanew) NO keep both for VerticalFreeSurfaceVelocity (to optimize)
     IF (curvilinearONOFF==1) THEN
        CALL DetermineTransformationConstantsArray(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
             FineGrid,FineGrid%dsigmanew,Wavefield)
     ENDIF
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

  CALL VerticalFreeSurfaceVelocity(Wavefield%W,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,PHI, &
       FineGrid%DiffStencils,FineGrid%dsigmanew(:,:,:,5),gamma)
  
  Wavefield%WHist(:,:,3)=Wavefield%WHist(:,:,2);  Wavefield%WHist(:,:,2)=Wavefield%WHist(:,:,1);
  Wavefield%WHist(:,:,1)=Wavefield%W(:,:); 

END SUBROUTINE Runge_Kutta_4
