SUBROUTINE RelaxationModule(E,P,time)
!
! Apply relaxation zone
!
! By Allan P. Engsig-Karup.
USE GlobalVariables, ONLY: RelaxZones, FineGrid, SFsol, g, relaxTransientTime, relaxNo, &
	GhostGridX, GhostGridY, curvilinearONOFF
USE Precision
USE Constants
IMPLICIT NONE
!REAL(KIND=long), DIMENSION(FineGrid%Nx,FineGrid%Ny) :: E, P
! GD: change
REAL(KIND=long), DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY) :: E, P
REAL(KIND=long) :: time, FAC
INTEGER :: i, j, k
REAL(KIND=long) :: tmpx(FineGrid%Nx+2*GhostGridX), tmpy(FineGrid%Ny+2*GhostGridY)

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
	   DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
         E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam)
         P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam)
 	   END DO
	 ELSE
	  IF (RelaxZones(i)%degrees==zero) THEN
      	CALL stream_func_wave_finite(RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1,&
		     FineGrid%x(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),1),time &
			 ,SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
        !CALL stream_func_wave_finite(RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1,&
		!     tmpx(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2)),time &
		!	 ,SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
	    DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
          E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam) &
		      + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
          P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam) &
		      + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
 	    END DO
	  ELSE
	    DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
          CALL stream_func_wave_finite(RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1,&
		       tmpx(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2))&
			   *COS(RelaxZones(i)%degrees/180.0_long*pi)+tmpy(j)*SIN(RelaxZones(i)%degrees/180.0_long*pi),&
			   time,SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
          E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam) &
		       + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
          P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam) &
		       + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
 	    END DO
	  ENDIF
	 ENDIF
  ELSE IF (RelaxZones(i)%XorY=='Y' .AND. RelaxZones(i)%XorYgen=='Y') THEN ! Y-direction
	 IF (RelaxZones(i)%WavegenONOFF==0) THEN
	   DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
         E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam)
         P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam)
 	   END DO
	 ELSE
	   IF (RelaxZones(i)%degrees==zero) THEN
         CALL stream_func_wave_finite(RelaxZones(i)%idx(4)-RelaxZones(i)%idx(3)+1,&
		      tmpy(RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)),time,&
			  SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
	     DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
           E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam) &
		        + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
           P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam) &
		        + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
 	     END DO
	   ELSE
	     DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
           CALL stream_func_wave_finite(RelaxZones(i)%idx(4)-RelaxZones(i)%idx(3)+1,&
		        -tmpy(RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))&
				*SIN(RelaxZones(i)%degrees/180.0_long*pi)+tmpx(j)*COS(RelaxZones(i)%degrees/180.0_long*pi),time,&
				SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
           E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam) &
		        + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
           P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam) &
		        + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
 	     END DO
	   ENDIF
	 ENDIF
  ELSE IF (RelaxZones(i)%XorY=='X' .AND. RelaxZones(i)%XorYgen=='Y') THEN
	 IF (RelaxZones(i)%WavegenONOFF==0) THEN
	   DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
         E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam)
         P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam)
 	   END DO
	 ELSE
	  IF (RelaxZones(i)%degrees==zero) THEN
         CALL stream_func_wave_finite(RelaxZones(i)%idx(4)-RelaxZones(i)%idx(3)+1,&
			tmpy(RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)),time,&
			SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
	     DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
		   k = j - RelaxZones(i)%idx(1) + 1
           E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam) &
		       + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam)
           P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam) &
		       + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam)
 	     END DO
	   ELSE
	     DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
		   k = j - RelaxZones(i)%idx(1) + 1
           CALL stream_func_wave_finite(RelaxZones(i)%idx(4)-RelaxZones(i)%idx(3)+1,&
			    -tmpy(RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))&
				*SIN(RelaxZones(i)%degrees/180.0_long*pi)+tmpx(j)*COS(RelaxZones(i)%degrees/180.0_long*pi),time,&
				SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
           E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam(k)) &
		       + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam(k))
           P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam(k)) &
		       + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam(k))
 	     END DO
	   ENDIF
	 ENDIF
  ELSE IF (RelaxZones(i)%XorY=='Y' .AND. RelaxZones(i)%XorYgen=='X') THEN ! Y-direction
	 IF (RelaxZones(i)%WavegenONOFF==0) THEN
	   DO j = RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
         E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = E(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam)
         P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4)) = P(j,RelaxZones(i)%idx(3):RelaxZones(i)%idx(4))*(RelaxZones(i)%gam)
 	   END DO
	 ELSE
	   IF (RelaxZones(i)%degrees==zero) THEN
        CALL stream_func_wave_finite(RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1,&
		   tmpx(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2)),time,&
		   SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
	    DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
		  k = j - RelaxZones(i)%idx(3) + 1
          E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam(k)) &
		     + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam(k))
          P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam(k)) &
		     + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam(k))
 	    END DO
	  ELSE
	    DO j = RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
		  k = j - RelaxZones(i)%idx(3) + 1
          CALL stream_func_wave_finite(RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1,&
		     tmpx(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2))&
			 *COS(RelaxZones(i)%degrees/180.0_long*pi)+tmpy(j)*SIN(RelaxZones(i)%degrees/180.0_long*pi),time,&
			 SFsol%n_four_modes,SFsol%zz,SFsol%yy,SFsol%k,g,RelaxZones(i)%Ea,RelaxZones(i)%Pa)
          E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = E(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam(k)) &
		     + FAC*RelaxZones(i)%Ea*(one-RelaxZones(i)%gam(k))
          P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j) = P(RelaxZones(i)%idx(1):RelaxZones(i)%idx(2),j)*(RelaxZones(i)%gam(k)) &
		     + FAC*RelaxZones(i)%Pa*(one-RelaxZones(i)%gam(k))
 	    END DO
	  ENDIF
	 ENDIF
   ENDIF
END DO

END SUBROUTINE RelaxationModule
