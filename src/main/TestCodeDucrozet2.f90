SUBROUTINE TestCodeDucrozet2

  USE GlobalVariables
  USE MGLevels
  IMPLICIT NONE
  EXTERNAL BuildLinearSystem, BuildLinearSystemTransformedCurvilinear
  ! GD: to test the cross derivatives...
  REAL(KIND=long), DIMENSION(:,:,:), ALLOCATABLE :: tmpPHI
  TYPE (Diff_def)        :: FullRankStencils
  INTEGER i, j, k

  !
  ! A number of convergence and consistency checks
  !
  IF ((IC==11).AND.(0==1)) THEN
     ALLOCATE(EA(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
     ALLOCATE(WA(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
     ! DO Convergence test
     time = zero*SFsol%T/two

     DO k=1,FineGrid%Nz+GhostGridZ
        DO j=1,FineGrid%Ny+2*GhostGridY
           DO i=1,FineGrid%Nx+2*GhostGridX
              CALL McCamyFuchs54LocalXY(FineGrid%x(i,j),FineGrid%y(i,j),FineGrid%z(k)*SFsol%h-SFsol%h,Lx,&
                   SFsol%HH/two,two*pi/SFsol%L,SFsol%h,g,time,0.0000001_long,100,Wavefield%E(i,j),PHI(k,i,j),WA(i,j))
           ENDDO
        ENDDO
     ENDDO
     Wavefield%P(1:FineGrid%Nx+2*GhostGridX,1:FineGrid%Ny+2*GhostGridY) = &
          PHI(FineGrid%Nz+GhostGridZ,1:FineGrid%Nx+2*GhostGridX,1:FineGrid%Ny+2*GhostGridY)
     print*,'finished initialisation'
     !
     CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E, &
          PHI(FineGrid%Nz+GhostGridZ,:,:),FineGrid,0,formattype)
     !
     CALL DifferentiationsFreeSurfacePlane(Wavefield,GhostGridX,GhostGridY,FineGrid,alpha,beta)
     !RHS(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+2*GhostGridX,1+GhostGridY:FineGrid%Ny+2*GhostGridY) = &
     !   Wavefield%P(1+GhostGridX:FineGrid%Nx+2*GhostGridX,1+GhostGridY:FineGrid%Ny+2*GhostGridY)
     !GD: do not take ghost points...
     RHS(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY) = &
          Wavefield%P(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)

     IF (curvilinearONOFF==1) THEN
        CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             BuildLinearSystemTransformedCurvilinear,PHI,FineGrid)
     ELSE
        CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             BuildLinearSystem,PHI,FineGrid)
     END IF
     print*,'iterations finished'
     CALL VerticalFreeSurfaceVelocity(Wavefield%W,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
          PHI,FineGrid%DiffStencils,FineGrid%dsigmanew(:,:,:,5),gamma)
     print*,'Vertical velocity computed'
     filename = "P.bin"
     CALL StoreRealArray(Wavefield%P,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "WA.bin"
     CALL StoreRealArray(WA,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "W.bin"
     CALL StoreRealArray(Wavefield%W,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "X.bin"
     CALL StoreRealArray(FineGrid%x,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "Y.bin"
     CALL StoreRealArray(FineGrid%y,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     print*,'Arrays stored'

     !CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny),Wavefield%W(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),&
     !		WA(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     ! the resolution of the system involve reflexion on outer wall... let's take only half point to exclude the non-correct ones
     !CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny)/2,Wavefield%W(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny/2+GhostGridY),&
     !		WA(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny/2+GhostGridY),dum)
     CALL maxnorm((FineGrid%Nx),Wavefield%W(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY),&
          WA(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY),dum)
     !dum = dum / (half*SFsol%k*SFsol%HH*SFsol%c)
     WRITE(*,*) ' Max.-norm error in free surface w (on cylinder) = ',dum !,maxval(abs(WA)),maxval(abs(Wavefield%W))
     ! the resolution of the system involve reflexion on outer wall... let's take only half point to exclude the non-correct ones
     CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny)/2,&
          Wavefield%W(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny/2+GhostGridY),&
          WA(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny/2+GhostGridY),dum)
     !
     WRITE(*,*) ' Max.-norm error in free surface w (inside domain) = ',dum !,maxval(abs(WA)),maxval(abs(Wavefield%W))
     !
     DEALLOCATE(EA,WA)
     PRINT*,''
     PRINT*,'Analysis done. Terminating program.'
     stop

  ENDIF

  !$$$$$$ IF (.TRUE.) THEN
  !$$$$$$   ! Test filtering
  !$$$$$$   Wavefield%E_I = zero
  !$$$$$$   Wavefield%E = COS(two*pi/(0.049)*FineGrid%x)+two/ten*COS(10000*two*pi/(0.049)*FineGrid%x)
  !$$$$$$   !
  !$$$$$$   Wavefield%P_I_s = zero
  !$$$$$$   Wavefield%P = zero
  !$$$$$$   CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%P,FineGrid,0,formattype)
  !$$$$$$   CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%P,FineGrid,1,formattype)
  !$$$$$$   !
  !$$$$$$   CALL FILTERING_NEW(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield,filterNP,filterALPHA,filtercoefficients,&
  !$$$$$$             tstep,swenseONOFF,filtercoefficients2, GhostGridX, GhostGridY)
  !$$$$$$   !
  !$$$$$$   CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%P,FineGrid,2,formattype)
  !$$$$$$   CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%E,Wavefield%P,FineGrid,3,formattype)
  !$$$$$$   STOP
  !$$$$$$ ENDIF
  IF ((IC==12).AND.(curvilinearONOFF==0)) THEN ! Convergence analysis, straight boundaries
     ALLOCATE(tmpPHI((FineGrid%Nz+GhostGridZ),(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY)))
     !
     WRITE (*,'(A)') '  Convergence analysis'
     WRITE (*,'(A/)') '========================'
     !
     ALLOCATE(EA(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
     ! Store reference values on the free surface DXZ and DYZ
     CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,LASTPHI(FineGrid%Nz+GhostGridZ,:,:,1), &
          LASTPHI(FineGrid%Nz+GhostGridZ,:,:,2),FineGrid,0,formattype)
     ! Compute derivatives
     ! DXZ
     EA=LASTPHI(FineGrid%Nz+GhostGridZ,:,:,1)
     !CALL DiffXEven(PHI,PHI2,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,alpha)
     !CALL DiffZArbitrary(PHI2,tmpPHI,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
     CALL DiffXEven_CD(PHI,PHI2,FineGrid%DiffStencils%IndexesX_XZorXY(:,:,:,:,1),FineGrid%DiffStencils%StencilsX_XZorXY(:,:,:,:,1),&
          FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,alpha)
     CALL DiffZArbitrary_CD(PHI2,tmpPHI,FineGrid%DiffStencils%IndexesZ_XZorYZ(:,:,:,:,1),&
          FineGrid%DiffStencils%StencilsZ_XZorYZ(:,:,:,:,1),&
          FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,gamma)
     !
     DO k=1,FineGrid%Nz+GhostGridZ
        print*,'apek:array do not conform.'
        STOP
        !  	tmpPHI(k,:,:) = tmpPHI(k,:,:)*FineGrid%dsigma(:,5)
     ENDDO
     !
     !CALL maxnorm((FineGrid%Nz+GhostGridZ)*(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY),LASTPHI(:,:,:,1),tmpPHI,dum)
     CALL maxnorm((FineGrid%Nz)*(FineGrid%Nx)*(FineGrid%Ny),&
          LASTPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,1),&
          tmpPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     PRINT*, ' Max.-norm error in whole domain Dxz(PHI) = ', dum
     PRINT*,''
     !PRINT*,'location maximum...',MAXLOC(LASTPHI(:,:,:,1)-tmpPHI)
     !PRINT*,'location maximum...',MAXLOC(LASTPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,1)-&
     !  tmpPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY))
     CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny),EA(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),&
          tmpPHI(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     PRINT*, ' Max.-norm error on free surface Dxz(PHI) = ', dum
     PRINT*,''

     LASTPHI(:,:,:,1)=tmpPHI
     !
     ! Compute derivatives
     !DYZ
     EA=LASTPHI(FineGrid%Nz+GhostGridZ,:,:,2)
     !CALL DiffYEven(PHI,PHI2,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,alpha)
     !CALL DiffZArbitrary(PHI2,tmpPHI,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
     CALL DiffYEven_CD(PHI,PHI2,FineGrid%DiffStencils%IndexesY_YZorXY(:,:,:,:,1),FineGrid%DiffStencils%StencilsY_YZorXY(:,:,:,:,1),&
          FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,alpha)
     CALL DiffZArbitrary_CD(PHI2,tmpPHI,FineGrid%DiffStencils%IndexesZ_XZorYZ(:,:,:,:,2),&
          FineGrid%DiffStencils%StencilsZ_XZorYZ(:,:,:,:,2),&
          FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,gamma)
     !
     DO k=1,FineGrid%Nz+GhostGridZ
        print*,'apek:array do not conform.'
        STOP
        !  	tmpPHI(k,:,:) = tmpPHI(k,:,:)*FineGrid%dsigma(:,5)
     ENDDO
     !
     !CALL maxnorm((FineGrid%Nz+GhostGridZ)*(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY),LASTPHI(:,:,:,2),tmpPHI,dum)
     CALL maxnorm((FineGrid%Nz)*(FineGrid%Nx)*(FineGrid%Ny),&
          LASTPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,2),&
          tmpPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     PRINT*, ' Max.-norm error in whole domain Dyz(PHI) = ', dum
     PRINT*,''
     !PRINT*,'location maximum...',MAXLOC(LASTPHI(:,:,:,2)-tmpPHI)
     !PRINT*,'location maximum...',MAXLOC(LASTPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,2)-&
     !  tmpPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY))
     CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny),EA(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),&
          tmpPHI(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     PRINT*, ' Max.-norm error on free surface Dyz(PHI) = ', dum
     PRINT*,''
     LASTPHI(:,:,:,2)=tmpPHI

     ! Store computed values
     CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,LASTPHI(FineGrid%Nz+GhostGridZ,:,:,1),&
          LASTPHI(FineGrid%Nz+GhostGridZ,:,:,2),FineGrid,1,formattype)
     !
     PRINT*,'Analysis done. Terminating program.'
     filename = "SparseMatrix.bin"
     CALL StoreSparseMatrix(FineGrid%PreconditioningMatrix,filename,formattype)
     print*,'Sparse matrix stored. stopped here.'
     STOP
  ENDIF

  IF ((IC==12).AND.(curvilinearONOFF==1)) THEN ! Convergence analysis, curvilinear
     !
     WRITE (*,'(A)') '  Convergence analysis'
     WRITE (*,'(A/)') '========================'
     !
     ALLOCATE(EA(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
     ! Store reference values
     CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%Ex,Wavefield%Exx,FineGrid,0,formattype)
     CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%Ey,Wavefield%Eyy,FineGrid,1,formattype)
     ! Compute derivatives
     EA=Wavefield%Ex
     CALL DiffXEvenCurvilinear(Wavefield%E,Wavefield%Ex, FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,&
          alpha)
     CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny),Wavefield%Ex(1+GhostGridX:FineGrid%Nx+GhostGridX,&
          1+GhostGridY:FineGrid%Ny+GhostGridY),EA(1+GhostGridX:FineGrid%Nx+GhostGridX,&
          1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     !dum = dum / (half*SFsol%k*SFsol%HH*SFsol%c)
     PRINT*, ' Max.-norm error in free surface Ex = ', dum
     PRINT*,''
     !
     EA=Wavefield%Exx
     !CALL DiffXXEvenCurvilinear(Wavefield%E,Wavefield%Exx,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,&
     !	alpha)
     CALL DiffXXEvenCurvilinear_CD(Wavefield%E,Wavefield%Exx,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,alpha, &
          FineGrid%CurvilinearStuff%DiffStencils%IndexesX_XZorXY(1,:,:,:,2),&
          FineGrid%CurvilinearStuff%DiffStencils%StencilsX_XZorXY(1,:,:,:,2), &
          FineGrid%CurvilinearStuff%DiffStencils%IndexesY_YZorXY(1,:,:,:,2),&
          FineGrid%CurvilinearStuff%DiffStencils%StencilsY_YZorXY(1,:,:,:,2))
     CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny),Wavefield%Exx(1+GhostGridX:FineGrid%Nx+GhostGridX,&
          1+GhostGridY:FineGrid%Ny+GhostGridY),EA(1+GhostGridX:FineGrid%Nx+GhostGridX,&
          1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     !dum = dum / (half*SFsol%k*SFsol%HH*SFsol%c)
     PRINT*, ' Max.-norm error in free surface Exx = ', dum
     PRINT*, ' Location Max. error in free surface Exx = ', &
          MAXLOC(ABS(Wavefield%Exx(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)-&
          EA(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)))
     PRINT*,''
     !
     EA=Wavefield%Ey
     CALL DiffYEvenCurvilinear(Wavefield%E,Wavefield%Ey, FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,&
          alpha)
     !CALL DiffYEvenCurvilinear_CD(Wavefield%E,Wavefield%Ey, FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid%DiffStencils,alpha, GidxTableCD_alpha, GidxTableCD_beta)
     CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny),Wavefield%Ey(1+GhostGridX:FineGrid%Nx+GhostGridX,&
          1+GhostGridY:FineGrid%Ny+GhostGridY),EA(1+GhostGridX:FineGrid%Nx+GhostGridX,&
          1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     !dum = dum / (half*SFsol%k*SFsol%HH*SFsol%c)
     PRINT*, ' Max.-norm error in free surface Ey = ', dum
     PRINT*,''
     !
     EA=Wavefield%Eyy
     !CALL DiffYYEvenCurvilinear(Wavefield%E,Wavefield%Eyy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,1,FineGrid,&
     !	alpha)
     CALL DiffYYEvenCurvilinear_CD(Wavefield%E,Wavefield%Eyy,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,&
          1,FineGrid,alpha, &
          FineGrid%CurvilinearStuff%DiffStencils%IndexesX_XZorXY(1,:,:,:,2),&
          FineGrid%CurvilinearStuff%DiffStencils%StencilsX_XZorXY(1,:,:,:,2), &
          FineGrid%CurvilinearStuff%DiffStencils%IndexesY_YZorXY(1,:,:,:,2),&
          FineGrid%CurvilinearStuff%DiffStencils%StencilsY_YZorXY(1,:,:,:,2))
     CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny),Wavefield%Eyy(1+GhostGridX:FineGrid%Nx+GhostGridX,&
          1+GhostGridY:FineGrid%Ny+GhostGridY),EA(1+GhostGridX:FineGrid%Nx+GhostGridX,&
          1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     !dum = dum / (half*SFsol%k*SFsol%HH*SFsol%c)
     PRINT*, ' Max.-norm error in free surface Eyy = ', dum
     PRINT*, ' Location Max. error in free surface Eyy = ', &
          MAXLOC(ABS(Wavefield%Eyy(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)-&
          EA(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)))
     PRINT*,''
     ! Store computed values
     CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%Ex,Wavefield%Exx,FineGrid,2,formattype)
     CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,Wavefield%Ey,Wavefield%Eyy,FineGrid,3,formattype)
     !
     ! Check convergence of other cross derivatives...
     !
     ALLOCATE(tmpPHI((FineGrid%Nz+GhostGridZ),(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY)))
     tmpPHI = zero
     !dsde
     !CALL DiffXuniform3D(PHI,PHI2 ,1,FineGrid%CurvilinearStuff%DiffStencils%StencilG,&
     !    FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,alpha)
     !CALL DiffZArbitrary(PHI2,tmpPHI,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,alpha)
     !
     CALL DiffXuniform3D_CD(PHI,PHI2, FineGrid%CurvilinearStuff%DiffStencils%IndexesX_XZorXY(:,:,:,:,1),&
          FineGrid%CurvilinearStuff%DiffStencils%StencilsX_XZorXY(:,:,:,:,1),&
          FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,alpha)
     CALL DiffZArbitrary_CD(PHI2,tmpPHI, FineGrid%CurvilinearStuff%DiffStencils%IndexesZ_XZorYZ(:,:,:,:,1),&
          FineGrid%CurvilinearStuff%DiffStencils%StencilsZ_XZorYZ(:,:,:,:,1),&
          FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,alpha)

     !LASTPHI(:,:,:,2)=tmpPHI
     LASTPHI(:,:,:,1)=tmpPHI
     !dsdn
     !CALL DiffYuniform3D(PHI,PHI2 ,1,FineGrid%CurvilinearStuff%DiffStencils%StencilG,&
     !    FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,alpha)
     !CALL DiffZArbitrary(PHI2,tmpPHI,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,alpha)
     !
     CALL DiffYuniform3D_CD(PHI,PHI2,FineGrid%CurvilinearStuff%DiffStencils%IndexesY_YZorXY(:,:,:,:,1),&
          FineGrid%CurvilinearStuff%DiffStencils%StencilsY_YZorXY(:,:,:,:,1),&
          FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,alpha)
     CALL DiffZArbitrary_CD(PHI2,tmpPHI,FineGrid%CurvilinearStuff%DiffStencils%IndexesZ_XZorYZ(:,:,:,:,2),&
          FineGrid%CurvilinearStuff%DiffStencils%StencilsZ_XZorYZ(:,:,:,:,2),&
          FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,alpha)
     PHI2=zero
     ! Interior
     DO j = 2 , FineGrid%Ny+GhostGridY
        DO i = 2 , FineGrid%Nx+GhostGridX
           DO k = 2 , FineGrid%Nz+GhostGridZ
              ! Dxz
              !PHI2(k,i,j) = FineGrid%dsigmanew(k,i,j,5)*(FineGrid%CurvilinearStuff%ex(i,j)*LASTPHI(k,i,j,2)+FineGrid%CurvilinearStuff%nx(i,j)*tmpPHI(k,i,j))
              ! Dyz
              PHI2(k,i,j) = FineGrid%dsigmanew(k,i,j,5)*(FineGrid%CurvilinearStuff%ey(i,j)*LASTPHI(k,i,j,1)&
                   +FineGrid%CurvilinearStuff%ny(i,j)*tmpPHI(k,i,j))
           END DO
        END DO
     END DO

     tmpPHI = PHI2

     PHI2=zero
     PHI2(2:FineGrid%Nz+GhostGridZ,2:FineGrid%Nx+GhostGridX,2:FineGrid%Ny+GhostGridY) = &
          (LASTPHI(2:FineGrid%Nz+GhostGridZ,2:FineGrid%Nx+GhostGridX,2:FineGrid%Ny+GhostGridY,1)&
          -tmpPHI(2:FineGrid%Nz+GhostGridZ,2:FineGrid%Nx+GhostGridX,2:FineGrid%Ny+GhostGridY))

     !$$$$$$    ! Store computed and reference values
     !$$$$$$   CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,tmpPHI(FineGrid%Nz+GhostGridZ,:,:),LASTPHI(FineGrid%Nz+GhostGridZ,:,:,1),FineGrid,4,formattype)
     !$$$$$$
     !$$$$$$   CALL maxnorm((FineGrid%Nz)*(FineGrid%Nx)*(FineGrid%Ny),&
     !$$$$$$     LASTPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,1),&
     !$$$$$$     tmpPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     !$$$$$$   PRINT*, ' Max.-norm error in whole domain Dxz(PHI) = ', dum
     !$$$$$$   PRINT*, 'Max original quantity...', MAXVAL(ABS(LASTPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,1)))
     !$$$$$$   PRINT*,''
     !$$$$$$   !PRINT*,'location maximum...',MAXLOC(LASTPHI(:,:,:,1)-tmpPHI)
     !$$$$$$   PRINT*,'location maximum...',MAXLOC(LASTPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,1)-&
     !$$$$$$     tmpPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY))
     !$$$$$$   CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny),LASTPHI(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,1),&
     !$$$$$$     tmpPHI(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     !$$$$$$   PRINT*, ' Max.-norm error on free surface Dxz(PHI) = ', dum
     !$$$$$$   PRINT*, 'Max original quantity...',MAXVAL(ABS(LASTPHI(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,1)))
     !$$$$$$   PRINT*,''

     ! Store computed and reference values
     CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,tmpPHI(FineGrid%Nz+GhostGridZ,:,:),&
          LASTPHI(FineGrid%Nz+GhostGridZ,:,:,2),FineGrid,4,formattype)

     CALL maxnorm((FineGrid%Nz)*(FineGrid%Nx)*(FineGrid%Ny),&
          LASTPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,2),&
          tmpPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     PRINT*, ' Max.-norm error in whole domain Dyz(PHI) = ', dum
     !PRINT*, 'Max original quantity...', MAXVAL(ABS(LASTPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,2)))
     !PRINT*,'location maximum error...',MAXLOC(ABS(LASTPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,2)-&
     !  tmpPHI(1+GhostGridZ:FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)))

     PRINT*,''
     CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny),&
          LASTPHI(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,2),&
          tmpPHI(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     PRINT*, ' Max.-norm error on free surface Dyz(PHI) = ', dum
     !PRINT*, 'Max original quantity...',MAXVAL(ABS(LASTPHI(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY,2)))
     PRINT*,''
     ! Test
     CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%CurvilinearStuff%ey,FineGrid%CurvilinearStuff%ny,&
          FineGrid,5,formattype)
     CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%CurvilinearStuff%ex,FineGrid%CurvilinearStuff%nx,&
          FineGrid,6,formattype)

     PRINT*,'Analysis done. Terminating program.'
     !filename = "SparseMatrix.bin"
     !CALL StoreSparseMatrix(FineGrid%PreconditioningMatrix,filename,formattype)
     !print*,'Sparse matrix stored. stopped here.'
     STOP
  ENDIF

  !IF (IC==8 .AND. 1==0) THEN
  IF (IC==8) THEN
     ! GD: Initialization of tstep needed in iterative_solution
     tstep=0
     WRITE (*,'(A)') '  Convergence analysis'
     WRITE (*,'(A/)') '========================'
     ALLOCATE(EA(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
     ALLOCATE(WA(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
     ! DO Convergence test
     time = SFsol%T/four

     IF (CurvilinearONOFF == 1) THEN
        !Rotate back the physical grid with a given angle to determine the reference values...
        tmp2D = FineGrid%x
        FineGrid%x = tmp2D*COS(-20*PI/180)-FineGrid%y*SIN(-20*PI/180)
        FineGrid%y = tmp2D*SIN(-20*PI/180)+FineGrid%y*COS(-20*PI/180)
     ENDIF

     CALL linearstandingwave2D(g,zero,time,SFsol,SFsol%L,SFsol%L,FineGrid%h(1,1),FineGrid%x,FineGrid%y&
          ,Wavefield%P,tmp2D,tmp2D,WA,EA,Wavefield%Ex,Wavefield%Exx,Wavefield%Ey,Wavefield%Eyy, &
          (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY))

     IF (CurvilinearONOFF == 1) THEN
        !Rotate the physical grid with a given angle...
        tmp2D = FineGrid%x
        FineGrid%x = tmp2D*COS(20*PI/180)-FineGrid%y*SIN(20*PI/180)
        FineGrid%y = tmp2D*SIN(20*PI/180)+FineGrid%y*COS(20*PI/180)
     ENDIF

     filename = "P.bin"
     CALL StoreRealArray(Wavefield%P,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)

     !RHS(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+2*GhostGridX,1+GhostGridY:FineGrid%Ny+2*GhostGridY) = &
     !   Wavefield%P(1+GhostGridX:FineGrid%Nx+2*GhostGridX,1+GhostGridY:FineGrid%Ny+2*GhostGridY)
     !GD: do not take ghost points...
     RHS(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY) = &
          Wavefield%P(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)

     IF (curvilinearONOFF==1) THEN
        CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             BuildLinearSystemTransformedCurvilinear,PHI,FineGrid)
     ELSE
        CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             BuildLinearSystem,PHI,FineGrid)
     END IF

     !filename = "SparseMatrix.bin"
     !CALL StoreSparseMatrix(FineGrid%PreconditioningMatrix,filename,formattype)
     !print*,'Sparse matrix stored. stopped here.'

     !	! Convert COO matrix to CSR matrix
     !	CALL ConvertCOOtoCSR(FineGrid%PreconditioningMatrix,FineGrid%IterationMatrix)
     !	! reorder CSR matrix to prepare for gauss-seidel operations
     !	CALL CSRdiaREORDER(FineGrid%IterationMatrix%nrow,FineGrid%IterationMatrix%nnz,&
     !		FineGrid%IterationMatrix%val,FineGrid%IterationMatrix%col_ind,FineGrid%IterationMatrix%row_ptr)
     !	CALL ConvertCSRtoCOO(FineGrid%IterationMatrix,FineGrid%PreconditioningMatrix)
     !filename = "SparseMatrix.bin"
     !CALL StoreSparseMatrix(FineGrid%PreconditioningMatrix,filename,formattype)
     !print*,'Sparse matrix stored. stopped here.'
     !STOP

     !	! Test Gauss-Seidel
     !	PHI = RHS
     !	DO i = 1, 500
     !		CALL GaussSeidel(PHI,RHS,FineGrid%IterationMatrix)
     !	END DO

     CALL VerticalFreeSurfaceVelocity(Wavefield%W,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
          PHI,FineGrid%DiffStencils,FineGrid%dsigmanew(:,:,:,5), gamma)
     filename = "WA.bin"
     CALL StoreRealArray(WA,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "W.bin"
     CALL StoreRealArray(Wavefield%W,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "X.bin"
     CALL StoreRealArray(FineGrid%x,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "Y.bin"
     CALL StoreRealArray(FineGrid%y,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     !	filename = "Z.bin"
     !CALL StoreRealArray(FineGrid%z,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)

     !  CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,W,WA,FineGrid,9999,formattype)
     CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny),Wavefield%W(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),&
          WA(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),dum)
     dum = dum / (half*SFsol%k*SFsol%HH*SFsol%c)
     WRITE(6,651) dum
     CALL GatherData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,dum,GridZ,2*gamma+1,SFsol%k*SFsol%h,&
          formattype)
651  FORMAT(' Max.-norm error in free surface vertical velocity w = ',E10.4)
     DEALLOCATE(EA,WA)
     PRINT*,''
     PRINT*,'Analysis done. Terminating program.'
     STOP
  ENDIF
  !
  IF (IC==13 .AND. 1==0) THEN
     !IF (IC==13) THEN
     ! GD: Initialization of tstep needed in iterative_solution
     tstep=0
     WRITE (*,'(A)') '  Convergence analysis'
     WRITE (*,'(A/)') '========================'
     ALLOCATE(EA(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
     ALLOCATE(WA(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
     ! DO Convergence test
     time = SFsol%T/four

     IF (CurvilinearONOFF == 1) THEN
        !Rotate back the physical grid with a given angle to determine the reference values...
        tmp2D = FineGrid%x
        FineGrid%x = tmp2D*COS(-20.d0*PI/180.d0)-FineGrid%y*SIN(-20.d0*PI/180.d0)
        FineGrid%y = tmp2D*SIN(-20.d0*PI/180.d0)+FineGrid%y*COS(-20.d0*PI/180.d0)
     ENDIF
     ! we just need WA and EA here
     CALL linearstandingwave2D(g,zero,time,SFsol,SFsol%L,SFsol%L,FineGrid%h(1,1),FineGrid%x,FineGrid%y&
          ,tmp2D,tmp2D,tmp2D,WA,EA,tmp2D,tmp2D,tmp2D,tmp2D, &
          (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY))

     IF (CurvilinearONOFF == 1) THEN
        !Rotate the physical grid with a given angle...
        tmp2D = FineGrid%x
        FineGrid%x = tmp2D*COS(20.d0*PI/180.d0)-FineGrid%y*SIN(20.d0*PI/180.d0)
        FineGrid%y = tmp2D*SIN(20.d0*PI/180.d0)+FineGrid%y*COS(20.d0*PI/180.d0)
     ENDIF

     IF (LinearONOFF==0) THEN
        ! Linear incident wavefield
        CALL incident_linear_wf_finite_standing(1, Wavefield, & !ramp_type,&
             FineGrid%Nx+2*GhostGridX,FineGrid%x,FineGrid%Ny+2*GhostGridY,FineGrid%y,FineGrid%Nz+GhostGridz,FineGrid%z,FineGrid%h,&
             time,SFsol%k,g,SFsol%HH,SFsol%h)
     ELSE
        ! Nonlinear incident wavefield
        print*,' this case is not available already... several checks to make...'
        stop
     ENDIF
     !
     CALL DifferentiationsFreeSurfacePlane(Wavefield,GhostGridX,GhostGridY,FineGrid,alpha,beta)

     filename = "P.bin"
     CALL StoreRealArray(Wavefield%P,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)

     !RHS(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+2*GhostGridX,1+GhostGridY:FineGrid%Ny+2*GhostGridY) = &
     !   Wavefield%P(1+GhostGridX:FineGrid%Nx+2*GhostGridX,1+GhostGridY:FineGrid%Ny+2*GhostGridY)
     !GD: do not take ghost points...
     RHS(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY) = &
          Wavefield%P(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)

     IF (curvilinearONOFF==1) THEN
        CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             BuildLinearSystemTransformedCurvilinear,PHI,FineGrid)
     ELSE
        CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             BuildLinearSystem,PHI,FineGrid)
     END IF

     CALL VerticalFreeSurfaceVelocity(Wavefield%W,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
          PHI,FineGrid%DiffStencils,FineGrid%dsigmanew(:,:,:,5), gamma)
     filename = "WA.bin"
     CALL StoreRealArray(WA,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "W.bin"
     !CALL StoreRealArray(Wavefield%W+Wavefield%Pz_I_s,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     CALL StoreRealArray(Wavefield%W,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     !CALL StoreRealArray(Wavefield%Pz_I_s,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "X.bin"
     CALL StoreRealArray(FineGrid%x,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "Y.bin"
     CALL StoreRealArray(FineGrid%y,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     !	filename = "Z.bin"
     !CALL StoreRealArray(FineGrid%z,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)

     !  CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,W,WA,FineGrid,9999,formattype)
     CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny),Wavefield%W(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY) &
          +Wavefield%Pz_I_s(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),&
          WA(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),dum)

     dum = dum / (half*SFsol%k*SFsol%HH*SFsol%c)
     WRITE(6,651) dum
     CALL GatherData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,dum,GridZ,2*gamma+1,SFsol%k*SFsol%h,&
          formattype)
     stop
     DEALLOCATE(EA,WA)
     PRINT*,' '
     PRINT*,'Analysis done. Terminating program.'
     STOP
  ENDIF
  IF (IC==14 .AND. 1==0) THEN
     !IF (IC==13) THEN
     ! GD: Initialization of tstep needed in iterative_solution
     tstep=0
     WRITE (*,'(A)') '  Convergence analysis'
     WRITE (*,'(A/)') '========================'
     ALLOCATE(EA(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
     ALLOCATE(WA(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
     !
     time = zero

     IF (CurvilinearONOFF == 1) THEN
        !Rotate back the physical grid with a given angle to determine the reference values...
        tmp2D = FineGrid%x
        FineGrid%x = tmp2D*COS(-20.d0*PI/180.d0)-FineGrid%y*SIN(-20.d0*PI/180.d0)
        FineGrid%y = tmp2D*SIN(-20.d0*PI/180.d0)+FineGrid%y*COS(-20.d0*PI/180.d0)
     ENDIF
     ! we just need WA and EA here
     CALL nonlinearstandingwave1D(pi,FineGrid%h(1,1),FineGrid%x,tmp2D,WA,EA,tmp2D,tmp2D,&
          (FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY))

     IF (CurvilinearONOFF == 1) THEN
        !Rotate the physical grid with a given angle...
        tmp2D = FineGrid%x
        FineGrid%x = tmp2D*COS(20.d0*PI/180.d0)-FineGrid%y*SIN(20.d0*PI/180.d0)
        FineGrid%y = tmp2D*SIN(20.d0*PI/180.d0)+FineGrid%y*COS(20.d0*PI/180.d0)
     ENDIF

     ! NonLinear incident wavefield
     CALL incident_nonlinear_wf_finite_standing(1, Wavefield, & !ramp_type,&
          FineGrid%Nx+2*GhostGridX,FineGrid%x,FineGrid%Ny+2*GhostGridY,FineGrid%y,FineGrid%Nz+GhostGridz,FineGrid%z,FineGrid%h,&
          time,SFsol%k,g,SFsol%HH,SFsol%h)
     !
     CALL DifferentiationsFreeSurfacePlane(Wavefield,GhostGridX,GhostGridY,FineGrid,alpha,beta)

     filename = "P.bin"
     CALL StoreRealArray(Wavefield%P,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)

     !GD: do not take ghost points...
     RHS(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY) = &
          Wavefield%P(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY)

     IF (curvilinearONOFF==1) THEN
        CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             BuildLinearSystemTransformedCurvilinear,PHI,FineGrid)
     ELSE
        CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
             BuildLinearSystem,PHI,FineGrid)
     END IF

     CALL VerticalFreeSurfaceVelocity(Wavefield%W,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,&
          PHI,FineGrid%DiffStencils,FineGrid%dsigmanew(:,:,:,5), gamma)
     filename = "WA.bin"
     CALL StoreRealArray(WA,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "W.bin"
     !CALL StoreRealArray(Wavefield%W+Wavefield%Pz_I_s,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     CALL StoreRealArray(Wavefield%W,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     !CALL StoreRealArray(Wavefield%Pz_I_s,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "X.bin"
     CALL StoreRealArray(FineGrid%x,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     filename = "Y.bin"
     CALL StoreRealArray(FineGrid%y,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)
     !	filename = "Z.bin"
     !CALL StoreRealArray(FineGrid%z,(FineGrid%Nx+2*GhostGridX),(FineGrid%Ny+2*GhostGridY),filename,formattype)

     !  CALL StoreData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,W,WA,FineGrid,9999,formattype)
     CALL maxnorm((FineGrid%Nx)*(FineGrid%Ny),Wavefield%W(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY) &
          +Wavefield%Pz_I_s(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),&
          WA(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),dum)

     dum = dum / (half*SFsol%k*SFsol%HH*SFsol%c)
     WRITE(6,651) dum
     CALL GatherData(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,dum,GridZ,2*gamma+1,SFsol%k*SFsol%h,&
          formattype)
     stop
     DEALLOCATE(EA,WA)
     PRINT*,' '
     PRINT*,'Analysis done. Terminating program.'
     STOP
  END IF
END SUBROUTINE TestCodeDucrozet2
