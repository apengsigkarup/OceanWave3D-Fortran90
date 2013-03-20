SUBROUTINE RelaxationModule_new(E,P,time,time0)
!
! Apply relaxation zones
!
! By Allan P. Engsig-Karup
USE GlobalVariables, ONLY: RelaxZones, FineGrid, SFsol, g, relaxTransientTime, relaxNo, &
	GhostGridX, GhostGridY, curvilinearONOFF, dt0, RandomWave3D, IncWaveType
USE Precision
USE Constants
IMPLICIT NONE
!REAL(KIND=long), DIMENSION(FineGrid%Nx,FineGrid%Ny) :: E, P
! GD: change
REAL(KIND=long), DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY) :: E, P
REAL(KIND=long) :: time, time0, FAC, etaLocal,x0,y0,yint,dy, etaInterp(RandomWave3d%order,RandomWave3d%order)
INTEGER :: i, j, k, itime,nnx,nny,order,countX,countY
REAL(KIND=long), ALLOCATABLE ::eta3d(:,:), phi3d(:,:)
REAL(KIND=long) :: tmpx(FineGrid%Nx+2*GhostGridX), tmpy(FineGrid%Ny+2*GhostGridY)
CHARACTER*3 CHOUT3

IF (time<relaxTransientTime) THEN
  FAC = time/relaxTransientTime
ELSE
  FAC = one
ENDIF
IF (curvilinearONOFF==1) THEN ! curvilinear coordinates... FIXME for make it work whatever transformation
  ! FIXME: use Lx and Ly? Depends how x and y are defined... save e and n directly ?
  DO j=1,FineGrid%Nx+2*GhostGridX
    !tmpx(j) = REAL(j-2,long)*(FineGrid%y(2,1)-FineGrid%y(1,1)) !/REAL((FineGrid%Nx+2*GhostGridX-1),long)*(8+PI)*Lx; !FIXME find a way to define in all situation...
    tmpx(j) = FineGrid%y(j,1)
    tmpx(j) = REAL(j,long)
  ENDDO
  DO j=1,FineGrid%Ny+2*GhostGridY
    !tmpy(j) = -REAL(j-2,long)*(FineGrid%x(1,2)-FineGrid%x(1,1)) !/REAL((FineGrid%Ny+2*GhostGridY-1),long)*(Ly-Lx); !FIXME find a way to define in all situation...
    tmpy(j) = FineGrid%x(1,j)
    tmpy(j) = REAL(j,long)
  ENDDO
ELSE
  tmpx(:) = FineGrid%x(:,1)
  tmpy(:) = FineGrid%y(1,:)
ENDIF
! FOR EACH RELAXATION ZONE
DO i = 1, relaxNo
   IF (RelaxZones(i)%XorY=='X' .AND. RelaxZones(i)%XorYgen=='X') THEN
      IF (RelaxZones(i)%WavegenONOFF==0) THEN
         ! No incident wave
         IF (RelaxZones(i)%PhiOnOff==1) THEN
            ! Relax both phi and eta
            DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
               E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) =  &
                    E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam)
               P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) =  &
                    P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam)
            END DO
         ELSE
            ! Relax only eta
            DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
               E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) =  &
                    E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam)
            END DO
         END IF
      ELSE
         IF (RelaxZones(i)%degrees==zero) THEN
             IF(IncWaveType==3) THEN !botp

                ! Nearest neighbour interpolation !TODO: implement higher order
                ! scheme.
                itime = (time-time0)/dt0+1

                ALLOCATE(eta3d(RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1,RelaxZones(i)%idx(4)-RelaxZones(i)%idx(3)+1), &
                         phi3d(RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1,RelaxZones(i)%idx(4)-RelaxZones(i)%idx(3)+1))

                ! We loop over all points in the relaxation zone and interpolate
                ! from the coarse wave gauge mesh
                countY = 0
                DO k =RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
                    countY = countY + 1
                    countX = 0
                    DO  j=RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
                    countX = countX +1

                        x0 = FineGrid%x(j,k)
                        y0 = FineGrid%y(j,k)

                        NNx = nint(x0/(RandomWave3d%dx))+1-RandomWave3d%order/2
                        NNY = nint(y0/(RandomWave3d%dy))+1-RandomWave3d%order/2

                        IF(NNX<1) THEN
                            NNX = 1
                        ELSEIF((NNx+RandomWave3d%order-1)>RandomWave3d%n1) THEN
                            NNX = RandomWave3d%n1 - RandomWave3d%order+1
                        ENDIF
                        
                        IF(NNY<1) THEN
                            NNy = 1
                        ELSEIF((NNY+RandomWave3d%order-1)>RandomWave3d%n2_data) THEN
                            NNY = RandomWave3d%n2_data - RandomWave3d%order + 1
                        ENDIF

                        etaInterp = RandomWave3d%eta(NNx:NNx+RandomWave3d%order-1,NNy:NNy+RandomWave3d%order-1,itime)

                        CALL polin2(RandomWave3d%x(NNx),RandomWave3d%y(NNy),etaInterp,RandomWave3d%order&
                                    ,RandomWave3d%order,x0,y0,yint,dy)
                        eta3d(countX,countY) = yint

                        etaInterp = RandomWave3d%phis(NNx:NNx+RandomWave3d%order-1,NNy:NNy+RandomWave3d%order-1,itime)
                        CALL polin2(RandomWave3d%x(NNx),RandomWave3d%y(NNy),etaInterp,RandomWave3d%order&
                                    ,RandomWave3d%order,x0,y0,yint,dy)
                        phi3d(countX,countY) = yint
                    ENDDO
                ENDDO
!open(unit=itime, file='eta3d.'//CHOUT3(INT(itime)),status='unknown')
   !DO j=1,countY
      !WRITE (itime,*) phi3d(:,j)
   !ENDDO
!CLOSE(itime)

                countY = 0
                DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
                   countY = countY + 1

                   E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = &
                        E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam) &
                        + FAC*eta3d(:,countY)*(one-RelaxZones(i)%gam)
                   P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = &
                        P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam) &
                        + FAC*phi3d(:,countY)*(one-RelaxZones(i)%gam)
                 ENDDO
             ELSE
                CALL AnalyticWaveMaker2D(RelaxZones(i)%idx(1),RelaxZones(i)%idx(2),&
                     FineGrid%x(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),1),time, &
                     time0,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
!      	CALL     stream_func_wave_finite(RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1,&
!!               FineGrid%x(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),1),time &
!                ,SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
                 DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
                   E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = &
                        E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam) &
                        + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
                   P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = &
                        P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam) &
                        + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
                END DO
            ENDIF
         ELSE
            DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
               CALL AnalyticWaveMaker2D(RelaxZones(i)%idx(1),RelaxZones(i)%idx(2),&
                    tmpx(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2))*COS(RelaxZones(i)%degrees/180.0_long*pi) &
                    +tmpy(j)*SIN(RelaxZones(i)%degrees/180.0_long*pi),time,time0, &
                    RelaxZones(i)%Ea,RelaxZones(i)%Pa)
!          CALL stream_func_wave_finite(RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1,&
!		       tmpx(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2))&
!			   *COS(RelaxZones(i)%degrees/180.0_long*pi)+tmpy(j)*SIN(RelaxZones(i)%degrees/180.0_long*pi),&
!			   time,SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
               E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)

               P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
            END DO
         ENDIF
      ENDIF
   ELSE IF (RelaxZones(i)%XorY=='Y' .AND. RelaxZones(i)%XorYgen=='Y') THEN ! Y-direction
      IF (RelaxZones(i)%WavegenONOFF==0) THEN
         ! No incident wave
         IF (RelaxZones(i)%PhiOnOff==1) THEN
            ! Relax both phi and eta
            DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
               E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = &
                    E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam)
               P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = &
                    P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam)
            END DO
         ELSE
            ! Relax only eta
            DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
               E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = &
                    E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam)
            END DO
         END IF
      ELSE
         IF (RelaxZones(i)%degrees==zero) THEN
            CALL AnalyticWaveMaker2D(RelaxZones(i)%idx(3),RelaxZones(i)%idx(4),   &
                 tmpy(RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)),time,time0,  &
                 RelaxZones(i)%Ea,RelaxZones(i)%Pa)
           !         CALL stream_func_wave_finite(RelaxZones(i)%idx(4)-RelaxZones(i)%idx(3)+1,&
!		      tmpy(RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)),time,&
!			  SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
            DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
               E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
               P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
            END DO
         ELSE
            DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
               CALL AnalyticWaveMaker2D(RelaxZones(i)%idx(3),RelaxZones(i)%idx(4), &
                    -tmpy(RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))                 &
                    *SIN(RelaxZones(i)%degrees/180.0_long*pi)+tmpx(j)*COS(RelaxZones(i)%degrees/180.0_long*pi), &
                    time,time0,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
!           CALL stream_func_wave_finite(RelaxZones(i)%idx(4)-RelaxZones(i)%idx(3)+1,&
!		        -tmpy(RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))&
!				*SIN(RelaxZones(i)%degrees/180.0_long*pi)+tmpx(j)*COS(RelaxZones(i)%degrees/180.0_long*pi),time,&
!				SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
               E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
               P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
            END DO
         ENDIF
      ENDIF
   ELSE IF (RelaxZones(i)%XorY=='X' .AND. RelaxZones(i)%XorYgen=='Y') THEN
      IF (RelaxZones(i)%WavegenONOFF==0) THEN
         IF (RelaxZones(i)%PhiOnOff==1) THEN
           ! Relax both phi and eta
            DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
               E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) &
                    *(RelaxZones(i)%gam)
               P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) &
                    *(RelaxZones(i)%gam)
            END DO
         ELSE
           ! Relax only eta
            DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
               E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) &
                    *(RelaxZones(i)%gam)
            END DO
         END IF
      ELSE
         IF (RelaxZones(i)%degrees==zero) THEN
            CALL AnalyticWaveMaker2D(RelaxZones(i)%idx(3),RelaxZones(i)%idx(4),   &
                 tmpy(RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)),time,time0, &
                 RelaxZones(i)%Ea,RelaxZones(i)%Pa)
!         CALL stream_func_wave_finite(RelaxZones(i)%idx(4)-RelaxZones(i)%idx(3)+1,&
!			tmpy(RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)),time,&
!			SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
            DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
               k = j - RelaxZones(i)%idx(1) + 1
               E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
               P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
            END DO
         ELSE
            DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
               k = j - RelaxZones(i)%idx(1) + 1
               CALL AnalyticWaveMaker2D(RelaxZones(i)%idx(3),RelaxZones(i)%idx(4),            &
                    -tmpy(RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))                            &
                    *SIN(RelaxZones(i)%degrees/180.0_long*pi)+tmpx(j)*COS(RelaxZones(i)%degrees/180.0_long*pi), &
                    time,time0,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
!           CALL stream_func_wave_finite(RelaxZones(i)%idx(4)-RelaxZones(i)%idx(3)+1,&
!			    -tmpy(RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))&
!				*SIN(RelaxZones(i)%degrees/180.0_long*pi)+tmpx(j)*COS(RelaxZones(i)%degrees/180.0_long*pi),time,&
!				SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
               E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) &
                    *(RelaxZones(i)%gam(k)) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam(k))
               P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) &
                    *(RelaxZones(i)%gam(k)) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam(k))
            END DO
         ENDIF
      ENDIF
   ELSE IF (RelaxZones(i)%XorY=='Y' .AND. RelaxZones(i)%XorYgen=='X') THEN ! Y-direction
      IF (RelaxZones(i)%WavegenONOFF==0) THEN
         IF (RelaxZones(i)%PhiOnOff==1) THEN
            ! Relax both phi and eta
            DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
               E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) &
                    *(RelaxZones(i)%gam)
               P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) &
                    *(RelaxZones(i)%gam)
            END DO
         ELSE
            ! Relax only eta
            DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
               E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) &
                    *(RelaxZones(i)%gam)
            END DO
         END IF
      ELSE
         IF (RelaxZones(i)%degrees==zero) THEN
            CALL AnalyticWaveMaker2D(RelaxZones(i)%idx(1),RelaxZones(i)%idx(2), &
                 tmpx(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2)),time,time0, &
                 RelaxZones(i)%Ea,RelaxZones(i)%Pa)
	    DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
               k = j - RelaxZones(i)%idx(3) + 1
               E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)&
                    *(RelaxZones(i)%gam(k)) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam(k))
               P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)&
                    *(RelaxZones(i)%gam(k)) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam(k))
            END DO
         ELSE
	    DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
               k = j - RelaxZones(i)%idx(3) + 1
               CALL AnalyticWaveMaker2D(RelaxZones(i)%idx(1),RelaxZones(i)%idx(2),            &
                    tmpx(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2))                            &
                    *COS(RelaxZones(i)%degrees/180.0_long*pi)+tmpy(j)*SIN(RelaxZones(i)%degrees/180.0_long*pi),&
                    time,time0,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
               E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)  &
                    *(RelaxZones(i)%gam(k)) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam(k))
               P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)  &
                    *(RelaxZones(i)%gam(k)) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam(k))
 	    END DO
         ENDIF
      ENDIF
   ENDIF
END DO

END SUBROUTINE RelaxationModule_new

