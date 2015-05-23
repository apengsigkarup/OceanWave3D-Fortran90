SUBROUTINE RelaxationModule_new(E,P,RKtime,time)
!
! Apply relaxation in all relaxation zones
!
! By Allan P. Engsig-Karup and Harry B. Bingham
!
USE GlobalVariables, ONLY: RelaxZones, FineGrid, SFsol, g, relaxTransientTime, relaxNo, &
	GhostGridX, GhostGridY, curvilinearONOFF, RandomWave, IncWaveType
USE Precision
USE Constants
IMPLICIT NONE
!REAL(KIND=long), DIMENSION(FineGrid%Nx,FineGrid%Ny) :: E, P
! GD: change
REAL(KIND=long), DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY) :: E, P
REAL(KIND=long) :: RKtime, time, FAC
INTEGER :: i, j, k, j0, j1, k0, k1, krel
REAL(KIND=long) :: tmpx(FineGrid%Nx+2*GhostGridX), tmpy(FineGrid%Ny+2*GhostGridY)

!
! The ramp factor in time for smooth initialization of the wave generation.
!
IF (RKtime<relaxTransienttime) THEN
  FAC = RKtime/relaxTransienttime
ELSE
  FAC = one
ENDIF
!
! Local x- and y-coordinate arrays.
!
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
!
! FOR EACH RELAXATION ZONE
!
DO i = 1, relaxNo
   ! Local index numbers to make things more readable 
   j0=RelaxZones(i)%idx(1); j1=RelaxZones(i)%idx(2)
   k0=RelaxZones(i)%idx(3); k1=RelaxZones(i)%idx(4)
   !
   ! The relaxation is split into four blocks depending on the values of 
   ! XorY (relaxation direction) and XorYgen (wave generation direction)
   !
   IF (RelaxZones(i)%XorY=='X' .AND. RelaxZones(i)%XorYgen=='X') THEN
      IF (RelaxZones(i)%WavegenONOFF==0) THEN
         ! No incident wave -> damping zone
         IF (RelaxZones(i)%PhiOnOff==1) THEN
            ! Relax both phi and eta
            DO j = k0,k1
               E(j0:j1,j) =  &
                    E(j0:j1,j)*(RelaxZones(i)%gam)
               P(j0:j1,j) =  &
                    P(j0:j1,j)*(RelaxZones(i)%gam)
            END DO
         ELSE
            ! Relax only eta
            DO j = k0,k1
               E(j0:j1,j) =  &
                    E(j0:j1,j)*(RelaxZones(i)%gam)
            END DO
         END IF
      ELSE
         !
         ! Get the target eta and phiS and apply to the relaxation zone.  
         !
         IF (RelaxZones(i)%degrees==zero) THEN
            If ( IncWaveType==2 .and. abs(RandomWave(i)%ispec) >= 30 ) THEN
               !
               ! 3D wave generation using linear regular and irregular waves.  
               !
               Do k=k0,k1
                  !
                  ! Load this x-directed line of solutions and relax.  
                  !
                  krel = k - k0 + 1
                  CALL AnalyticWaveMaker2D(i,j0,j1,krel,FineGrid%x(j0:j1,k),RKtime,time, &
                       RelaxZones(i)%Ea,RelaxZones(i)%Pa)
                  E(j0:j1,k) = &
                       E(j0:j1,k)*(RelaxZones(i)%gam) &
                       + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
                  P(j0:j1,k) = &
                       P(j0:j1,k)*(RelaxZones(i)%gam) &
                       + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
               End Do
            ELSE
               ! Pure 2D wave generation
               CALL AnalyticWaveMaker2D(i,j0,j1,1,FineGrid%x(j0:j1,1),RKtime,time, &
                    RelaxZones(i)%Ea,RelaxZones(i)%Pa)
               DO k = k0,k1
                  E(j0:j1,k) = &
                       E(j0:j1,k)*(RelaxZones(i)%gam) &
                       + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
                  P(j0:j1,k) = &
                       P(j0:j1,k)*(RelaxZones(i)%gam) &
                       + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
               END DO
            END If
         ELSE
            DO j = k0,k1
               CALL AnalyticWaveMaker2D(i,j0,j1,1,tmpx(j0:j1)*COS(RelaxZones(i)%degrees/180.0_long*pi) &
                    +tmpy(j)*SIN(RelaxZones(i)%degrees/180.0_long*pi),RKtime,time,  &
                    RelaxZones(i)%Ea,RelaxZones(i)%Pa)
               E(j0:j1,j) = E(j0:j1,j) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)

               P(j0:j1,j) = P(j0:j1,j) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
            END DO
         ENDIF
      ENDIF
   ELSE IF (RelaxZones(i)%XorY=='Y' .AND. RelaxZones(i)%XorYgen=='Y') THEN ! Y-direction
      IF (RelaxZones(i)%WavegenONOFF==0) THEN
         ! No incident wave
         IF (RelaxZones(i)%PhiOnOff==1) THEN
            ! Relax both phi and eta
            DO j = j0,j1
               E(j,k0:k1) = &
                    E(j,k0:k1)*(RelaxZones(i)%gam)
               P(j,k0:k1) = &
                    P(j,k0:k1)*(RelaxZones(i)%gam)
            END DO
         ELSE
            ! Relax only eta
            DO j = j0,j1
               E(j,k0:k1) = &
                    E(j,k0:k1)*(RelaxZones(i)%gam)
            END DO
         END IF
      ELSE
         IF (RelaxZones(i)%degrees==zero) THEN
            CALL AnalyticWaveMaker2D(i,k0,k1,1,tmpy(k0:k1),RKtime,time,   &
                 RelaxZones(i)%Ea,RelaxZones(i)%Pa)
            DO j = j0,j1
               E(j,k0:k1) = E(j,k0:k1) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
               P(j,k0:k1) = P(j,k0:k1) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
            END DO
         ELSE
            DO j = j0,j1
               CALL AnalyticWaveMaker2D(i,k0,k1,1, &
                    -tmpy(k0:k1)                 &
                    *SIN(RelaxZones(i)%degrees/180.0_long*pi)+tmpx(j)*COS(RelaxZones(i)%degrees/180.0_long*pi), &
                    RKtime,time, RelaxZones(i)%Ea,RelaxZones(i)%Pa)
               E(j,k0:k1) = E(j,k0:k1) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
               P(j,k0:k1) = P(j,k0:k1) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
            END DO
         ENDIF
      ENDIF
   ELSE IF (RelaxZones(i)%XorY=='X' .AND. RelaxZones(i)%XorYgen=='Y') THEN
      IF (RelaxZones(i)%WavegenONOFF==0) THEN
         IF (RelaxZones(i)%PhiOnOff==1) THEN
           ! Relax both phi and eta
            DO j = k0,k1
               E(j0:j1,j) = E(j0:j1,j) &
                    *(RelaxZones(i)%gam)
               P(j0:j1,j) = P(j0:j1,j) &
                    *(RelaxZones(i)%gam)
            END DO
         ELSE
           ! Relax only eta
            DO j = k0,k1
               E(j0:j1,j) = E(j0:j1,j) &
                    *(RelaxZones(i)%gam)
            END DO
         END IF
      ELSE
         IF (RelaxZones(i)%degrees==zero) THEN
            CALL AnalyticWaveMaker2D(i,k0,k1,1,tmpy(k0:k1),RKtime,time,  &
                 RelaxZones(i)%Ea,RelaxZones(i)%Pa)
            DO j = j0,j1
               k = j - j0 + 1
               E(j,k0:k1) = E(j,k0:k1) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
               P(j,k0:k1) = P(j,k0:k1) &
                    *(RelaxZones(i)%gam) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
            END DO
         ELSE
            DO j = j0,j1
               k = j - j0 + 1
               CALL AnalyticWaveMaker2D(i,k0,k1,1,            &
                    -tmpy(k0:k1)                              &
                    *SIN(RelaxZones(i)%degrees/180.0_long*pi)+tmpx(j)*COS(RelaxZones(i)%degrees/180.0_long*pi), &
                    RKtime,time, RelaxZones(i)%Ea,RelaxZones(i)%Pa)
               E(j,k0:k1) = E(j,k0:k1) &
                    *(RelaxZones(i)%gam(k)) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam(k))
               P(j,k0:k1) = P(j,k0:k1) &
                    *(RelaxZones(i)%gam(k)) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam(k))
            END DO
         ENDIF
      ENDIF
   ELSE IF (RelaxZones(i)%XorY=='Y' .AND. RelaxZones(i)%XorYgen=='X') THEN ! Y-direction
      IF (RelaxZones(i)%WavegenONOFF==0) THEN
         IF (RelaxZones(i)%PhiOnOff==1) THEN
            ! Relax both phi and eta
            DO j = j0,j1
               E(j,k0:k1) = E(j,k0:k1) &
                    *(RelaxZones(i)%gam)
               P(j,k0:k1) = P(j,k0:k1) &
                    *(RelaxZones(i)%gam)
            END DO
         ELSE
            ! Relax only eta
            DO j = j0,j1
               E(j,k0:k1) = E(j,k0:k1) &
                    *(RelaxZones(i)%gam)
            END DO
         END IF
      ELSE
         IF (RelaxZones(i)%degrees==zero) THEN
            If ( IncWaveType==2 .and. abs(RandomWave(i)%ispec) >= 30 ) THEN
               !
               ! 3D wave generation using linear regular and irregular waves.  
               !
               Do k=k0,k1
                  !
                  ! Load this x-directed line of solutions and relax.  
                  !
                  krel = k - k0 + 1
                  CALL AnalyticWaveMaker2D(i,j0,j1,krel,tmpx(j0:j1),RKtime,time,  &
                       RelaxZones(i)%Ea,RelaxZones(i)%Pa)
                  E(j0:j1,k) = E(j0:j1,k)&
                       *(RelaxZones(i)%gam(krel)) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam(krel))
                  P(j0:j1,k) = P(j0:j1,k)&
                       *(RelaxZones(i)%gam(krel)) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam(krel))
               END DO
            ELSE

               CALL AnalyticWaveMaker2D(i,j0,j1,1,tmpx(j0:j1),RKtime,time,  &
                    RelaxZones(i)%Ea,RelaxZones(i)%Pa)
               DO j = k0,k1
                  k = j - k0 + 1
                  E(j0:j1,j) = E(j0:j1,j)&
                       *(RelaxZones(i)%gam(k)) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam(k))
                  P(j0:j1,j) = P(j0:j1,j)&
                       *(RelaxZones(i)%gam(k)) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam(k))
               END DO
            END IF
         ELSE
	    DO j = k0,k1
               k = j - k0 + 1
               CALL AnalyticWaveMaker2D(i,j0,j1,1,            &
                    tmpx(j0:j1)                               &
                    *COS(RelaxZones(i)%degrees/180.0_long*pi)+tmpy(j)*SIN(RelaxZones(i)%degrees/180.0_long*pi),&
                    RKtime,time, RelaxZones(i)%Ea,RelaxZones(i)%Pa)
               E(j0:j1,j) = E(j0:j1,j)  &
                    *(RelaxZones(i)%gam(k)) + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam(k))
               P(j0:j1,j) = P(j0:j1,j)  &
                    *(RelaxZones(i)%gam(k)) + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam(k))
 	    END DO
         ENDIF
      ENDIF
   ENDIF
END DO

END SUBROUTINE RelaxationModule_new
