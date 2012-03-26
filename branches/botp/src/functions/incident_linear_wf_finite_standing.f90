!>
!! Evaluate the linear incident wave field  from
!! linear theory (eta, phi and derivatives)
!!
!!    Evaluation on the free surface
!!    Evaluation on boundaries
!<
!
subroutine incident_linear_wf_finite_standing(wf_type, Wavefield, & !ramp_type,&
	 nfx,xp,nfy,yp,nfz,zp,bottom,tp,wavenum,grav,H, d)
!
USE PRECISION
USE GlobalVariables, ONLY: GhostGridX,GhostGridY,GhostGridZ, swenseTransientTime, &
     West_refl, East_refl, North_refl, South_refl, curvilinearONOFF, SFSOL
USE DataTypes
USE Constants
!
IMPLICIT NONE
INTERFACE
	SUBROUTINE phi_linear_3D_standing(H,angle,Lx,Ly,x,y,z,d,omega,time,grav,phix,phiy,phiz,phi,phit)
    USE PRECISION
   	REAL(KIND=long), INTENT(IN)  :: H,Lx,Ly,x,y,z,d,omega,time,grav,angle
   	REAL(KIND=long), INTENT(OUT) :: phix,phiy,phiz
   	REAL(KIND=long), INTENT(OUT), OPTIONAL :: phi,phit
   	END SUBROUTINE phi_linear_3D_standing
END INTERFACE
!
INTEGER nfx, nfy, nfz, nbp, wf_type
REAL(KIND=long) tp, wavenum, grav, H, d
! For curvilinear
REAL(KIND=long), DIMENSION(nfx,nfy) :: xp,yp
REAL(KIND=long), DIMENSION(nfz) :: zp
REAL(KIND=long), DIMENSION(nfx,nfy) :: bottom

!
TYPE(Wavefield_FS) :: Wavefield
! Local variables
INTEGER i, j, k, tmpIdx, Gidx3
REAL(KIND=long) temp1,omega,angle
REAL(KIND=long) zbp
REAL(KIND=long) FAC
!
!
! Time ramp and its derivatives
!time_ramp(:)=ramp_value(ramp_type,tp)
IF (swenseTransientTime/=zero) THEN
  IF(tp==zero) THEN
    FAC = zero
  ELSEIF(tp<=swenseTransientTime) THEN
    !FAC = one/swenseTransientTime*tp
    FAC = -two/(swenseTransientTime**3)*tp**3+three/(swenseTransientTime**2)*tp**2
  ELSE
    FAC = one
  ENDIF
ELSE
  FAC=one
ENDIF
!
omega = sqrt(grav*wavenum*tanh(wavenum*d))
omega=wf_type/abs(wf_type)*omega
!
SFSOL%k=SQRT((two*pi/SFSOL%L)**2+(zero*two*pi/SFSOL%L)**2)
SFSOL%T=SQRT((two*pi)**2/(grav*SFSOL%k*TANH(SFSOL%k*d)) ); 
SFSOL%c=SQRT(grav/SFSOL%k*TANH(SFSOL%k*d))
      !
      IF (curvilinearONOFF.EQ.0) THEN
        IF((nfx.EQ.1).OR.(nfy.EQ.1)) THEN
          write(*,*) 'This is a 3d case !!'
          stop
        ENDIF
      ELSE
        ! FIXME: 3D only for now
		IF((nfx.EQ.1).OR.(nfy.EQ.1)) THEN
          write(*,*) 'The 2D case in curvilinear SEENSE has to be done...'
          stop
        ENDIF
      ENDIF
      !
      IF (curvilinearONOFF.EQ.1) THEN
      	angle = 20.0_long*PI/180.0_long
      ELSE
        angle= zero
      ENDIF
	  !   
      DO j=1,nfx
        DO k=1,nfy
            CALL eta_linear_3D_standing(H,angle,SFsol%L,SFsol%L,xp(j,k),yp(j,k),omega,tp,Wavefield%E_I(j,k),Wavefield%Ex_I(j,k), &
            	 Wavefield%Ey_I(j,k),Wavefield%Exx_I(j,k),Wavefield%Eyy_I(j,k),Wavefield%Et_I(j,k))
        	!
        	! Here phi, is computed on the free-surface z=0
        	!    
        	CALL phi_linear_3D_standing(H,angle,SFsol%L,SFsol%L,xp(j,k),yp(j,k),zero,d,omega,tp,grav, &
            	Wavefield%Px_I_s(j,k),Wavefield%Py_I_s(j,k),Wavefield%Pz_I_s(j,k),Wavefield%P_I_s(j,k),Wavefield%Pt_I_s(j,k))
        ENDDO
      END DO
      !
	!
    ! Determine the incident wavefield on the boundary points...
    !
	!
    Wavefield%SourceEx = -FAC*Wavefield%Ex_I
    Wavefield%SourcePx = -FAC*Wavefield%Px_I_s
    Wavefield%SourceEy = -FAC*Wavefield%Ey_I
    Wavefield%SourcePy = -FAC*Wavefield%Py_I_s
    !
  	! FIXME: Assumption here is that the plane is aligned with the x- and y-axes. Extend to
  	!        general curvilinear coordinates.
    ! FIXME: Optimization of loops
    Gidx3 = 1
    !
    ! Bottom boundary
    !
    DO j = 1+GhostGridY, nfy-GhostGridY
		DO i = 1+GhostGridX, nfx-GhostGridX
            ! Wave elevation + derivative (Dimensional)
       		Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(i,j)
        	Wavefield%Ex_I_bp(Gidx3) = FAC*Wavefield%Ex_I(i,j)
            Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(i,j)
            ! transform on the real z-space... Dimensional
            zbp=bottom(i,j)*zp(1+GhostGridZ)-bottom(i,j)
            !zbp=d*zp(2)-d
            CALL phi_linear_3D_standing(FAC*H,angle,SFsol%L,SFsol%L,xp(i,j),yp(i,j),zbp,d,omega,tp,grav, &
        			Wavefield%Px_I_bp(Gidx3),Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(1,i,j) = Gidx3
            !Wavefield%GidxTableBP(1+GhostGridZ,i,j) = Gidx3
			! Increment global index 
            Gidx3 = Gidx3+1
		END DO
	END DO	
    !
   	! West  boundary
	!
    i = 1
	!
    DO j = 1+GhostGridY, nfy-GhostGridY
        DO k = 1, nfz
            ! Wave elevation + derivative
            Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(i+GhostGridX,j)
            Wavefield%Ex_I_bp(Gidx3) = FAC*Wavefield%Ex_I(i+GhostGridX,j)
            Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(i+GhostGridX,j)
            IF (West_refl==0) THEN
				Wavefield%Px_I_bp(Gidx3)=zero; Wavefield%Py_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
            ELSE
                ! transform on the real z-space... Dimensional
                zbp=bottom(i+GhostGridX,j)*zp(k)-bottom(i+GhostGridX,j)
                !zbp=d*zp(k)-d
                CALL phi_linear_3D_standing(FAC*H,angle,SFsol%L,SFsol%L,xp(i+GhostGridX,j),yp(i+GhostGridX,j),zbp,d,omega,tp,grav, &
                	Wavefield%Px_I_bp(Gidx3),Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
            ENDIF
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(k,i,j) = Gidx3
            !Wavefield%GidxTableBP(k,i+GhostGridX,j) = Gidx3
            ! Increment global index    
            Gidx3 = Gidx3+1     
        END DO
    END DO
	! East  boundary
   	i = nfx
    ! index of the first element (=Nz*(Ny-2*GhostY))
    tmpIdx = Gidx3
	!
    DO j = 1+GhostGridY, nfy-GhostGridY
        DO k = 1, nfz
            ! Wave elevation + derivative
            Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(i-GhostGridX,j)
            Wavefield%Ex_I_bp(Gidx3) = FAC*Wavefield%Ex_I(i-GhostGridX,j)
            Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(i-GhostGridX,j)
            IF (East_refl==0) THEN
				Wavefield%Px_I_bp(Gidx3)=zero; Wavefield%Py_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
           	ELSE
                ! transform on the real z-space... Dimensional
                zbp=bottom(i-GhostGridX,j)*zp(k)-bottom(i-GhostGridX,j)
                !zbp=d*zp(k)-d
                CALL phi_linear_3D_standing(FAC*H,angle,SFsol%L,SFsol%L,xp(i-GhostGridX,j),yp(i-GhostGridX,j),zbp,d,omega,tp,grav, &
                    Wavefield%Px_I_bp(Gidx3),Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
            ENDIF
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(k,i,j) = Gidx3
            !Wavefield%GidxTableBP(k,i-GhostGridX,j) = Gidx3
            ! Increment global index        
            Gidx3 = Gidx3+1     
        END DO
    END DO
    !
    ! South boundary
    !
    j = 1 
    DO i = 1+GhostGridX, nfx-GhostGridX
        ! index of the first element (=2*Nz*(Ny-2*GhostY)+i)
        tmpIdx = Gidx3
        ! Wave elevation + derivative (Dimensional)
        Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(i,j+GhostGridY)
        Wavefield%Ex_I_bp(Gidx3) = FAC*Wavefield%Ex_I(i,j+GhostGridY)
        Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(i,j+GhostGridY)
        DO k = 1, nfz
            ! Wave elevation + derivative
            Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
            Wavefield%Ex_I_bp(Gidx3) = Wavefield%Ex_I_bp(tmpIdx)
            Wavefield%Ey_I_bp(Gidx3) = Wavefield%Ey_I_bp(tmpIdx)
            IF (South_refl==0) THEN
                Wavefield%Px_I_bp(Gidx3)=zero; Wavefield%Py_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
            ELSE
                ! transform on the real z-space... Dimensional
                zbp=bottom(i,j+GhostGridY)*zp(k)-bottom(i,j+GhostGridY)
                !zbp=d*zp(k)-d
                CALL phi_linear_3D_standing(FAC*H,angle,SFsol%L,SFsol%L,xp(i,j+GhostGridY),yp(i,j+GhostGridY),zbp,d,omega,tp,grav, &
                    Wavefield%Px_I_bp(Gidx3),Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
            ENDIF
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(k,i,j) = Gidx3
            !Wavefield%GidxTableBP(k,i,j+GhostGridY) = Gidx3
            ! Increment global index 
            Gidx3 = Gidx3+1
        END DO
    END DO
    !
    ! North boundary
    j = nfy 
    DO i = 1+GhostGridX, nfx-GhostGridX
        ! index of the first element (=2*Nz*(Ny-2*GhostY)+...+i)
        tmpIdx = Gidx3
        ! Wave elevation + derivative (Dimensional)
        Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(i,j-GhostGridY)
        Wavefield%Ex_I_bp(Gidx3) = FAC*Wavefield%Ex_I(i,j-GhostGridY)
        Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(i,j-GhostGridY)
        DO k = 1, nfz
            ! Wave elevation + derivative
            Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
            Wavefield%Ex_I_bp(Gidx3) = Wavefield%Ex_I_bp(tmpIdx)
            Wavefield%Ey_I_bp(Gidx3) = Wavefield%Ey_I_bp(tmpIdx)
            IF (North_refl==0) THEN
                Wavefield%Px_I_bp(Gidx3)=zero; Wavefield%Py_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
            ELSE
              ! transform on the real z-space... Dimensional
              zbp=bottom(i,j-GhostGridY)*zp(k)-bottom(i,j-GhostGridY)
              !zbp=d*zp(k)-d
              CALL phi_linear_3D_standing(FAC*H,angle,SFsol%L,SFsol%L,xp(i,j-GhostGridY),yp(i,j-GhostGridY),zbp,d,omega,tp,grav, &
                  Wavefield%Px_I_bp(Gidx3),Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
            ENDIF
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(k,i,j) = Gidx3
            !Wavefield%GidxTableBP(k,i,j-GhostGridY) = Gidx3
            ! Increment global index 
            Gidx3 = Gidx3+1
        END DO
    END DO
    ! FIXME: How to treat the corners ?
    IF (West_refl==0) THEN
      Wavefield%SourceEx(1+GhostGridX,2+GhostGridY:nfy-1-GhostGridY) = zero
      Wavefield%SourcePx(1+GhostGridX,2+GhostGridY:nfy-1-GhostGridY) = zero
      Wavefield%SourceEy(1+GhostGridX,2+GhostGridY:nfy-1-GhostGridY) = zero
      Wavefield%SourcePy(1+GhostGridX,2+GhostGridY:nfy-1-GhostGridY) = zero
    ENDIF
    IF (East_refl==0) THEN
      Wavefield%SourceEx(nfx-GhostGridX,2+GhostGridY:nfy-1-GhostGridY) = zero
      Wavefield%SourcePx(nfx-GhostGridX,2+GhostGridY:nfy-1-GhostGridY) = zero
      Wavefield%SourceEy(nfx-GhostGridX,2+GhostGridY:nfy-1-GhostGridY) = zero
      Wavefield%SourcePy(nfx-GhostGridX,2+GhostGridY:nfy-1-GhostGridY) = zero
    ENDIF
    IF (South_refl==0) THEN
      Wavefield%SourceEx(2+GhostGridX:nfx-1-GhostGridX,1+GhostGridY) = zero
      Wavefield%SourcePx(2+GhostGridX:nfx-1-GhostGridX,1+GhostGridY) = zero
      Wavefield%SourceEy(2+GhostGridX:nfx-1-GhostGridX,1+GhostGridY) = zero
      Wavefield%SourcePy(2+GhostGridX:nfx-1-GhostGridX,1+GhostGridY) = zero
    ENDIF
    IF (North_refl==0) THEN
      Wavefield%SourceEx(2+GhostGridX:nfx-1-GhostGridX,nfy-GhostGridY) = zero
      Wavefield%SourcePx(2+GhostGridX:nfx-1-GhostGridX,nfy-GhostGridY) = zero
      Wavefield%SourceEy(2+GhostGridX:nfx-1-GhostGridX,nfy-GhostGridY) = zero
      Wavefield%SourcePy(2+GhostGridX:nfx-1-GhostGridX,nfy-GhostGridY) = zero
    ENDIF 
    ! Corners are set to zero if both reflexions are zero...
    IF (West_refl==0) THEN
      IF (South_refl==0) THEN
      	Wavefield%SourceEx(1+GhostGridX,1+GhostGridY) = zero
      	Wavefield%SourcePx(1+GhostGridX,1+GhostGridY) = zero
      	Wavefield%SourceEy(1+GhostGridX,1+GhostGridY) = zero
      	Wavefield%SourcePy(1+GhostGridX,1+GhostGridY) = zero
      ENDIF
      IF (North_refl==0) THEN
      	Wavefield%SourceEx(1+GhostGridX,nfy-GhostGridY) = zero
      	Wavefield%SourcePx(1+GhostGridX,nfy-GhostGridY) = zero
      	Wavefield%SourceEy(1+GhostGridX,nfy-GhostGridY) = zero
      	Wavefield%SourcePy(1+GhostGridX,nfy-GhostGridY) = zero
      ENDIF
    ENDIF
    IF (East_refl==0) THEN
      IF (South_refl==0) THEN
      	Wavefield%SourceEx(nfx-GhostGridX,1+GhostGridY) = zero
      	Wavefield%SourcePx(nfx-GhostGridX,1+GhostGridY) = zero
      	Wavefield%SourceEy(nfx-GhostGridX,1+GhostGridY) = zero
      	Wavefield%SourcePy(nfx-GhostGridX,1+GhostGridY) = zero
      ENDIF
      IF (North_refl==0) THEN
      	Wavefield%SourceEx(nfx-GhostGridX,nfy-GhostGridY) = zero
      	Wavefield%SourcePx(nfx-GhostGridX,nfy-GhostGridY) = zero
      	Wavefield%SourceEy(nfx-GhostGridX,nfy-GhostGridY) = zero
      	Wavefield%SourcePy(nfx-GhostGridX,nfy-GhostGridY) = zero
      ENDIF
    ENDIF

RETURN
END SUBROUTINE incident_linear_wf_finite_standing
