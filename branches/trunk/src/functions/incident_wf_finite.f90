!>
!! Evaluate the nonlinear incident wave field  from
!! stream function theory (eta, phi and derivatives)
!!
!!    Evaluation on the free surface
!!    Evaluation on boundaries
!<
!
SUBROUTINE incident_wf_finite(wf_type, Wavefield, &
	 nfx,xp,nfy,yp,nfz,zp,bottom,tp,wavenum,grav,n_four_modes,zz,yy)
!
USE PRECISION
USE GlobalVariables, ONLY: GhostGridX,GhostGridY,GhostGridZ, swenseTransientTime, &
     West_refl, East_refl, North_refl, South_refl, curvilinearONOFF
USE DataTypes
USE Constants
IMPLICIT NONE
!
INTERFACE
	SUBROUTINE phi_nonlinear(k,x,z,omega,time,n_four_modes,zz,phix,phiz,phi,phit)
	USE PRECISION
	REAL(KIND=long), INTENT(IN)  :: k,x,z,omega,time, zz(2*n_four_modes+10) !GD: change zz(*)
	INTEGER :: n_four_modes
	REAL(KIND=long), INTENT(OUT) :: phix,phiz
	REAL(KIND=long), INTENT(OUT), OPTIONAL :: phi,phit
   	END SUBROUTINE phi_nonlinear
	!
    SUBROUTINE phi_nonlinear_3D(k,angle,x,y,z,omega,time,n_four_modes,zz,phix,phiy,phiz,phi,phit)
    USE PRECISION
    INTEGER :: n_four_modes
    REAL(KIND=long), INTENT(IN)  :: k,angle,x,y,z,omega,time,zz(2*n_four_modes+10)
    REAL(KIND=long), INTENT(OUT) :: phix,phiy,phiz
    REAL(KIND=long), INTENT(OUT), OPTIONAL :: phi,phit
    END SUBROUTINE phi_nonlinear_3D
END INTERFACE
!
INTEGER nfx, nfy, nfz, wf_type, n_four_modes
REAL(KIND=long) tp, wavenum, grav, H, d, zz(2*n_four_modes+10), yy(n_four_modes) ! GD: change zz(*), yy(*)
!REAL(KIND=long), DIMENSION(nfx) :: xp
!REAL(KIND=long), DIMENSION(nfy) :: yp
! For curvilinear
REAL(KIND=long), DIMENSION(nfx,nfy) :: xp, yp
REAL(KIND=long), DIMENSION(nfz) :: zp
REAL(KIND=long), DIMENSION(nfx,nfy) :: bottom
!
TYPE(Wavefield_FS) :: Wavefield
! Local variables
INTEGER i, j, k, tmpIdx, Gidx3
REAL(KIND=long) temp1,omega,angle,omega_t,kinv
REAL(KIND=long) zbp
REAL(KIND=long) FAC, tp_adim
REAL(KIND=long) :: ufactor, pfactor
REAL(KIND=long), DIMENSION(nfx,nfy) :: eta_S ! Temporary array (needed ?)
REAL(KIND=long), DIMENSION(2*n_four_modes+10) :: zz2
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
kinv=one/wavenum
ufactor=sqrt(grav*kinv)
pfactor=ufactor*kinv
omega_t=wavenum*ufactor*zz(4)*tp
omega=zz(4)
! Work in nondimensional space
tp_adim = tp*wavenum*ufactor
!      wf_type=1  => Wavefield propagating along x>0,rectangular domain
!      wf_type=-1 => Wavefield propagating along x<0,rectangular domain
!      wf_type=2  => Wavefield propagating along y>0,rectangular domain
!      wf_type=-2 => Wavefield propagating along y<0,rectangular domain
!	   wf_type=3  => boundary-fitted ?
!	   wf_type>3  => 3D regular with a wf_type angle (in degrees)
! FIXME: Is it necessary? (GD)
!
! FIXME: Think about the boundary fitted coordinate transformation...
! How to make the incident wave field in a clever way... Compute on all
! free surface points ? Going through an interpolation between 2 grids ?
omega=wf_type/abs(wf_type)*omega
omega_t=wf_type/abs(wf_type)*omega_t
!

IF((ABS(wf_type).EQ.1).AND.(curvilinearONOFF.EQ.0))THEN ! propagation along x
	IF(nfx==1) THEN
		PRINT*, 'SWENSE Propagation along x with only one point 1 x-point'
        STOP
    ENDIF
	!
    ! Be careful call eta_nonlinear and phi_nonlinear with non-dimensional time !
    DO j=1,nfx
      	CALL eta_nonlinear(wavenum,xp(j,1),omega,tp_adim,n_four_modes,yy,&
        	Wavefield%E_I(j,1),Wavefield%Ex_I(j,1),Wavefield%Exx_I(j,1),Wavefield%Et_I(j,1))
        ! Put back in dimensional space...
        Wavefield%E_I(j,1)   = kinv*Wavefield%E_I(j,1)
        Wavefield%Exx_I(j,1) = wavenum*Wavefield%Exx_I(j,1)
        Wavefield%Et_I(j,1)  = ufactor*Wavefield%Et_I(j,1)
    ENDDO
    ! Expand along y-direction if needed
	IF(nfy>1) THEN
        DO j=1,nfx
           Wavefield%E_I(j,:)   = Wavefield%E_I(j,1)
           Wavefield%Ex_I(j,:)  = Wavefield%Ex_I(j,1)
           Wavefield%Exx_I(j,:) = Wavefield%Exx_I(j,1)
           Wavefield%Et_I(j,:)  = Wavefield%Et_I(j,1)
           !
           Wavefield%Ey_I(j,:)   = zero
           Wavefield%Eyy_I(j,:)  = zero
           Wavefield%Py_I_s(j,:) = zero
        ENDDO
    ENDIF
    !
    DO i=1,nfx
      DO j=1,nfy 
        ! Careful, eta^S can be 3D...
        ! Here phi, is computed on the free-surface z=eta^S+eta^I
        ! eta^I and eta^S are dimensional
        !
        zbp = wavenum*(Wavefield%E(i,j)+Wavefield%E_I(i,j))
        CALL phi_nonlinear(wavenum,xp(i,1),zbp,omega,tp_adim,n_four_modes,zz, &
            Wavefield%Px_I_s(i,j),Wavefield%Pz_I_s(i,j),Wavefield%P_I_s(i,j),Wavefield%Pt_I_s(i,j))
        !
        Wavefield%Px_I_s(i,j)= ufactor*Wavefield%Px_I_s(i,j)
        Wavefield%Pz_I_s(i,j)= ufactor*Wavefield%Pz_I_s(i,j)
        Wavefield%P_I_s(i,j) = pfactor*Wavefield%P_I_s(i,j)
        Wavefield%Pt_I_s(i,j)= ufactor**2*Wavefield%Pt_I_s(i,j)
      END DO
    END DO
    !
    ! Determine the incident wavefield on the boundary points...
    !
    IF(nfy>1) THEN
    	Wavefield%Ey_I_bp = zero
    	Wavefield%Py_I_bp = zero
        !
        Wavefield%SourceEy = zero
        Wavefield%SourcePy = zero
    ENDIF
    Wavefield%SourceEx = -FAC*Wavefield%Ex_I
    Wavefield%SourcePx = -FAC*Wavefield%Px_I_s
    ! Take into account the eventual ramp...
    zz2(1:2*n_four_modes+10) = zz(1:2*n_four_modes+10)
    DO j=1,n_four_modes
    	zz2(n_four_modes+j+10) = FAC*zz(n_four_modes+j+10)
    ENDDO
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
            ! transform on the real z-space... non-dimensional
            ! FIXME: Here, use E_I or FAC*E_I ?? To check !
            zbp=wavenum*((bottom(i,j)+Wavefield%E(i,j)+FAC*Wavefield%E_I(i,j))*zp(1+GhostGridZ)-bottom(i,j))
            !zbp=d*zp(2)-d
            CALL phi_nonlinear(wavenum,xp(i,1),zbp,omega,tp_adim,n_four_modes,zz2, &
        			Wavefield%Px_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
            ! Put back in dimensional space...
    		Wavefield%Px_I_bp(Gidx3)= ufactor*Wavefield%Px_I_bp(Gidx3)
    		Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(1,i,j) = Gidx3
            !Wavefield%GidxTableBP(1-GhostGridZ,i,j) = Gidx3
			! Increment global index 
            Gidx3 = Gidx3+1
		END DO
	END DO
    !
   	! West  boundary
	!
    i = 1
    tmpIdx = Gidx3
    ! Wave elevation + derivative (Dimensional)
    Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(i+GhostGridX,1)
    Wavefield%Ex_I_bp(Gidx3) = FAC*Wavefield%Ex_I(i+GhostGridX,1)
	!
    DO j = 1+GhostGridY, nfy-GhostGridY
        DO k = 1, nfz
            ! Wave elevation + derivative
            Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
            Wavefield%Ex_I_bp(Gidx3) = Wavefield%Ex_I_bp(tmpIdx)
            IF (West_refl==0) THEN
				Wavefield%Px_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
            ELSE
                ! transform on the real z-space... non-dimensional
            	! FIXME: Here, use E_I or FAC*E_I ?? To check !
                zbp=wavenum*((bottom(i+GhostGridX,j)+Wavefield%E(i+GhostGridX,j)+FAC*Wavefield%E_I(i+GhostGridX,j))*zp(k) &
                	-bottom(i+GhostGridX,j))
                !zbp=bottom(i+GhostGridX,j)*zp(k)-bottom(i+GhostGridX,j)
                !zbp=d*zp(k)-d
                CALL phi_nonlinear(wavenum,xp(i+GhostGridX,1),zbp,omega,tp_adim,n_four_modes,zz2, &
        			Wavefield%Px_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                ! Put back in dimensional space...
    			Wavefield%Px_I_bp(Gidx3)= ufactor*Wavefield%Px_I_bp(Gidx3)
    			Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
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
    ! Wave elevation + derivative (Dimensional)
    Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(i-GhostGridX,1)
    Wavefield%Ex_I_bp(Gidx3) = FAC*Wavefield%Ex_I(i-GhostGridX,1)
	!
    DO j = 1+GhostGridY, nfy-GhostGridY
        DO k = 1, nfz
            ! Wave elevation + derivative
            Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
            Wavefield%Ex_I_bp(Gidx3) = Wavefield%Ex_I_bp(tmpIdx)
            IF (East_refl==0) THEN
				Wavefield%Px_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
           	ELSE
              	! transform on the real z-space... non-dimensional
            	! FIXME: Here, use E_I or FAC*E_I ?? To check !
                zbp=wavenum*((bottom(i-GhostGridX,j)+Wavefield%E(i-GhostGridX,j)+FAC*Wavefield%E_I(i-GhostGridX,j))*zp(k) &
                	-bottom(i-GhostGridX,j))
                !zbp=bottom(i-GhostGridX,j)*zp(k)-bottom(i-GhostGridX,j)
                !zbp=d*zp(k)-d
                CALL phi_nonlinear(wavenum,xp(i-GhostGridX,1),zbp,omega,tp_adim,n_four_modes,zz2, &
        			Wavefield%Px_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                ! Put back in dimensional space...
    			Wavefield%Px_I_bp(Gidx3)= ufactor*Wavefield%Px_I_bp(Gidx3)
    			Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
            ENDIF
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(k,i,j) = Gidx3
            !Wavefield%GidxTableBP(k,i-GhostGridX,j) = Gidx3
            ! Increment global index        
            Gidx3 = Gidx3+1     
        END DO
    END DO
    !
    IF (nfy>1) THEN
        ! South boundary
        j = 1 
        DO i = 1+GhostGridX, nfx-GhostGridX
            ! index of the first element (=2*Nz*(Ny-2*GhostY)+i)
            tmpIdx = Gidx3
            ! Wave elevation + derivative (Dimensional)
            Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(i,j+GhostGridY)
            Wavefield%Ex_I_bp(Gidx3) = FAC*Wavefield%Ex_I(i,j+GhostGridY)
            DO k = 1, nfz
                ! Wave elevation + derivative
                Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
                Wavefield%Ex_I_bp(Gidx3) = Wavefield%Ex_I_bp(tmpIdx)
                IF (South_refl==0) THEN
					Wavefield%Px_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
            	ELSE
                    ! transform on the real z-space... non-dimensional
                    ! FIXME: Here, use E_I or FAC*E_I ?? To check !
                    zbp=wavenum*((bottom(i,j+GhostGridY)+Wavefield%E(i,j+GhostGridY)+FAC*Wavefield%E_I(i,j+GhostGridY))*zp(k) &
                        -bottom(i,j+GhostGridY))
                    !zbp=bottom(i,j+GhostGridY)*zp(k)-bottom(i,j+GhostGridY)
                    !zbp=d*zp(k)-d
                    CALL phi_nonlinear(wavenum,xp(i,1),zbp,omega,tp_adim,n_four_modes,zz2, &
                        Wavefield%Px_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                    ! Put back in dimensional space...
    				Wavefield%Px_I_bp(Gidx3)= ufactor*Wavefield%Px_I_bp(Gidx3)
    				Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
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
            DO k = 1, nfz
                ! Wave elevation + derivative
                Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
                Wavefield%Ex_I_bp(Gidx3) = Wavefield%Ex_I_bp(tmpIdx)
                IF (North_refl==0) THEN
					Wavefield%Px_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
            	ELSE
                	! transform on the real z-space... non-dimensional
                    ! FIXME: Here, use E_I or FAC*E_I ?? To check !
                    zbp=wavenum*((bottom(i,j-GhostGridY)+Wavefield%E(i,j-GhostGridY)+FAC*Wavefield%E_I(i,j-GhostGridY))*zp(k) &
                        -bottom(i,j-GhostGridY))
                    !zbp=bottom(i,j-GhostGridY)*zp(k)-bottom(i,j-GhostGridY)
                    !zbp=d*zp(k)-d
                    CALL phi_nonlinear(wavenum,xp(i,1),zbp,omega,tp_adim,n_four_modes,zz2, &
                        Wavefield%Px_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                    ! Put back in dimensional space...
    				Wavefield%Px_I_bp(Gidx3)= ufactor*Wavefield%Px_I_bp(Gidx3)
    				Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
                ENDIF
                ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            	Wavefield%GidxTableBP(k,i,j) = Gidx3
            	!Wavefield%GidxTableBP(k,i,j-GhostGridY) = Gidx3
				! Increment global index 
                Gidx3 = Gidx3+1
            END DO
        END DO
    ENDIF
    ! FIXME: How to treat corners
    IF (West_refl==0) THEN
      Wavefield%SourceEx(1+GhostGridX,:) = zero
      Wavefield%SourcePx(1+GhostGridX,:) = zero
    ENDIF
    IF (East_refl==0) THEN
      Wavefield%SourceEx(nfx-GhostGridX,:) = zero
      Wavefield%SourcePx(nfx-GhostGridX,:) = zero
    ENDIF

!$$$$$$     ! wavefield along x, corners belong to x-boundaries
!$$$$$$     IF (nfy>1) THEN
!$$$$$$       IF (South_refl==0) THEN
!$$$$$$         Wavefield%SourceEx(:,1+GhostGridY) = zero
!$$$$$$         Wavefield%SourcePx(:,1+GhostGridY) = zero
!$$$$$$       ENDIF
!$$$$$$       IF (North_refl==0) THEN
!$$$$$$         Wavefield%SourceEx(:,nfy-GhostGridY) = zero
!$$$$$$         Wavefield%SourcePx(:,nfy-GhostGridY) = zero
!$$$$$$       ENDIF
!$$$$$$     ENDIF
ELSEIF((ABS(wf_type).EQ.2).AND.(curvilinearONOFF.EQ.0))THEN ! propagation along y
    !
    IF(nfy==1) THEN
		PRINT*, 'SWENSE Propagation along y with only one point 1 y-point'
        STOP
    ENDIF
	!
    ! Be careful call eta_nonlinear and phi_nonlinear with non-dimensional time !
    DO j=1,nfy
      	CALL eta_nonlinear(wavenum,yp(1,j),omega,tp_adim,n_four_modes,yy,&
        	Wavefield%E_I(1,j),Wavefield%Ey_I(1,j),Wavefield%Eyy_I(1,j),Wavefield%Et_I(1,j))
        ! Put back in dimensional space...
        Wavefield%E_I(1,j)   = kinv*Wavefield%E_I(1,j)
        Wavefield%Eyy_I(1,j) = wavenum*Wavefield%Eyy_I(1,j)
        Wavefield%Et_I(1,j)  = ufactor*Wavefield%Et_I(1,j)
    ENDDO
    ! Expand along x-direction if needed
	IF(nfx>1) THEN
        DO j=1,nfy
           Wavefield%E_I(:,j)   = Wavefield%E_I(1,j)
           Wavefield%Ey_I(:,j)  = Wavefield%Ey_I(1,j)
           Wavefield%Eyy_I(:,j) = Wavefield%Eyy_I(1,j)
           Wavefield%Et_I(:,j)  = Wavefield%Et_I(1,j)
           !
           Wavefield%Ex_I(:,j)   = zero
           Wavefield%Exx_I(:,j)  = zero
           Wavefield%Px_I_s(:,j) = zero
        ENDDO
    ENDIF
    !
    DO i=1,nfx
      DO j=1,nfy 
        ! Careful, eta^S can be 3D...
        ! Here phi, is computed on the free-surface z=eta^S+eta^I
        ! eta^I and eta^S are dimensional
        !
        zbp = wavenum*(Wavefield%E(i,j)+Wavefield%E_I(i,j))
        CALL phi_nonlinear(wavenum,yp(1,j),zbp,omega,tp_adim,n_four_modes,zz, &
            Wavefield%Py_I_s(i,j),Wavefield%Pz_I_s(i,j),Wavefield%P_I_s(i,j),Wavefield%Pt_I_s(i,j))
        !
        Wavefield%Py_I_s(i,j)= ufactor*Wavefield%Py_I_s(i,j)
        Wavefield%Pz_I_s(i,j)= ufactor*Wavefield%Pz_I_s(i,j)
        Wavefield%P_I_s(i,j) = pfactor*Wavefield%P_I_s(i,j)
        Wavefield%Pt_I_s(i,j)= ufactor**2*Wavefield%Pt_I_s(i,j)
      END DO
    END DO
    !
    ! Determine the incident wavefield on the boundary points...
    !
    IF(nfx>1) THEN
    	Wavefield%Ex_I_bp = zero
    	Wavefield%Px_I_bp = zero
        !
        Wavefield%SourceEx = zero
        Wavefield%SourcePx = zero
    ENDIF
    Wavefield%SourceEy = -FAC*Wavefield%Ey_I
    Wavefield%SourcePy = -FAC*Wavefield%Py_I_s
    ! Take into account the eventual ramp...
    zz2(1:2*n_four_modes+10) = zz(1:2*n_four_modes+10)
    DO j=1,n_four_modes
    	zz2(n_four_modes+j+10) = FAC*zz(n_four_modes+j+10)
    ENDDO
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
        	Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(i,j)
            ! transform on the real z-space... non-dimensional
            ! FIXME: Here, use E_I or FAC*E_I ?? To check !
            zbp=wavenum*((bottom(i,j)+Wavefield%E(i,j)+FAC*Wavefield%E_I(i,j))*zp(1+GhostGridZ)-bottom(i,j))
            !zbp=d*zp(2)-d
            CALL phi_nonlinear(wavenum,yp(1,j),zbp,omega,tp_adim,n_four_modes,zz2, &
        			Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
            ! Put back in dimensional space...
    		Wavefield%Py_I_bp(Gidx3)= ufactor*Wavefield%Py_I_bp(Gidx3)
    		Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(1,i,j) = Gidx3
            !Wavefield%GidxTableBP(1-GhostGridZ,i,j) = Gidx3
			! Increment global index 
            Gidx3 = Gidx3+1
		END DO
	END DO
    !
    IF (nfx>1) THEN
       ! West  boundary
       !
       i = 1
       !
       DO j = 1+GhostGridY, nfy-GhostGridY
           ! index of the first element
           tmpIdx = Gidx3
       	   ! Wave elevation + derivative (Dimensional)
       	   Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(i+GhostGridX,j)
       	   Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(i+GhostGridX,j)
           DO k = 1, nfz
               ! Wave elevation + derivative
               Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
               Wavefield%Ey_I_bp(Gidx3) = Wavefield%Ey_I_bp(tmpIdx)
               IF (West_refl==0) THEN
                   Wavefield%Py_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
               ELSE
                   ! transform on the real z-space... non-dimensional
                   ! FIXME: Here, use E_I or FAC*E_I ?? To check !
                   zbp=wavenum*((bottom(i+GhostGridX,j)+Wavefield%E(i+GhostGridX,j)+FAC*Wavefield%E_I(i+GhostGridX,j))*zp(k) &
                       -bottom(i+GhostGridX,j))
                   !zbp=bottom(i+GhostGridX,j)*zp(k)-bottom(i+GhostGridX,j)
                   !zbp=d*zp(k)-d
                   CALL phi_nonlinear(wavenum,yp(1,j),zbp,omega,tp_adim,n_four_modes,zz2, &
                       Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                   ! Put back in dimensional space...
                   Wavefield%Py_I_bp(Gidx3)= ufactor*Wavefield%Py_I_bp(Gidx3)
                   Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
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
       !
       DO j = 1+GhostGridY, nfy-GhostGridY
           ! index of the first element (=Nz*(Ny-2*GhostY))
           tmpIdx = Gidx3
           ! Wave elevation + derivative (Dimensional)
           Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(i-GhostGridX,j)
           Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(i-GhostGridX,j)
           DO k = 1, nfz
               ! Wave elevation + derivative
               Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
               Wavefield%Ey_I_bp(Gidx3) = Wavefield%Ey_I_bp(tmpIdx)
               IF (East_refl==0) THEN
                   Wavefield%Py_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
               ELSE
                   ! transform on the real z-space... non-dimensional
                   ! FIXME: Here, use E_I or FAC*E_I ?? To check !
                   zbp=wavenum*((bottom(i-GhostGridX,j)+Wavefield%E(i-GhostGridX,j)+FAC*Wavefield%E_I(i-GhostGridX,j))*zp(k) &
                       -bottom(i-GhostGridX,j))
                   !zbp=bottom(i-GhostGridX,j)*zp(k)-bottom(i-GhostGridX,j)
                   !zbp=d*zp(k)-d
                   CALL phi_nonlinear(wavenum,yp(1,j),zbp,omega,tp_adim,n_four_modes,zz2, &
                       Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                   ! Put back in dimensional space...
                   Wavefield%Py_I_bp(Gidx3)= ufactor*Wavefield%Py_I_bp(Gidx3)
                   Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
               ENDIF
               ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            	Wavefield%GidxTableBP(k,i,j) = Gidx3
            	!Wavefield%GidxTableBP(k,i-GhostGridX,j) = Gidx3
               ! Increment global index        
               Gidx3 = Gidx3+1     
           END DO
       END DO
    ENDIF  
    !
    ! South boundary
    j = 1 
    ! index of the first element (=2*Nz*(Ny-2*GhostY)+i)
    tmpIdx = Gidx3
    ! Wave elevation + derivative (Dimensional)
    Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(1,j+GhostGridY)
    Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(1,j+GhostGridY)
    !  
    DO i = 1+GhostGridX, nfx-GhostGridX
        DO k = 1, nfz
            ! Wave elevation + derivative
            Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
            Wavefield%Ey_I_bp(Gidx3) = Wavefield%Ey_I_bp(tmpIdx)
            IF (South_refl==0) THEN
                Wavefield%Py_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
            ELSE
                ! transform on the real z-space... non-dimensional
                ! FIXME: Here, use E_I or FAC*E_I ?? To check !
                zbp=wavenum*((bottom(i,j+GhostGridY)+Wavefield%E(i,j+GhostGridY)+FAC*Wavefield%E_I(i,j+GhostGridY))*zp(k) &
                    -bottom(i,j+GhostGridY))
                !zbp=bottom(i,j+GhostGridY)*zp(k)-bottom(i,j+GhostGridY)
                !zbp=d*zp(k)-d
                CALL phi_nonlinear(wavenum,yp(1,j+GhostGridY),zbp,omega,tp_adim,n_four_modes,zz2, &
                    Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                ! Put back in dimensional space...
                Wavefield%Py_I_bp(Gidx3)= ufactor*Wavefield%Py_I_bp(Gidx3)
                Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
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
    ! index of the first element (=2*Nz*(Ny-2*GhostY)+...+i)
    tmpIdx = Gidx3
    ! Wave elevation + derivative (Dimensional)
    Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(1,j-GhostGridY)
    Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(1,j-GhostGridY)
    !  
    DO i = 1+GhostGridX, nfx-GhostGridX
        DO k = 1, nfz
            ! Wave elevation + derivative
            Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
            Wavefield%Ey_I_bp(Gidx3) = Wavefield%Ey_I_bp(tmpIdx)
            IF (North_refl==0) THEN
                Wavefield%Py_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
            ELSE
                ! transform on the real z-space... non-dimensional
                ! FIXME: Here, use E_I or FAC*E_I ?? To check !
                zbp=wavenum*((bottom(i,j-GhostGridY)+Wavefield%E(i,j-GhostGridY)+FAC*Wavefield%E_I(i,j-GhostGridY))*zp(k) &
                    -bottom(i,j-GhostGridY))
                !zbp=bottom(i,j-GhostGridY)*zp(k)-bottom(i,j-GhostGridY)
                !zbp=d*zp(k)-d
                CALL phi_nonlinear(wavenum,yp(1,j-GhostGridY),zbp,omega,tp_adim,n_four_modes,zz2, &
                    Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                ! Put back in dimensional space...
                Wavefield%Py_I_bp(Gidx3)= ufactor*Wavefield%Py_I_bp(Gidx3)
                Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
            ENDIF
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(k,i,j) = Gidx3
            !Wavefield%GidxTableBP(k,i,j-GhostGridY) = Gidx3
            ! Increment global index 
            Gidx3 = Gidx3+1
        END DO
    END DO
    ! FIXME: How to treat boundaries...
    IF (South_refl==0) THEN
      Wavefield%SourceEy(:,1+GhostGridY) = zero
      Wavefield%SourcePy(:,1+GhostGridY) = zero
    ENDIF
    IF (North_refl==0) THEN
      Wavefield%SourceEy(:,nfy-GhostGridY) = zero
      Wavefield%SourcePy(:,nfy-GhostGridY) = zero
    ENDIF
!$$$$$$     ! wavefield along y, corners belong to y-boundaries
!$$$$$$     IF (nfy>1) THEN  
!$$$$$$       IF (West_refl==0) THEN
!$$$$$$         Wavefield%SourceEy(1+GhostGridX,:) = zero
!$$$$$$         Wavefield%SourcePy(1+GhostGridX,:) = zero
!$$$$$$       ENDIF
!$$$$$$       IF (East_refl==0) THEN
!$$$$$$         Wavefield%SourceEy(nfx-GhostGridX,:) = zero
!$$$$$$         Wavefield%SourcePy(nfx-GhostGridX,:) = zero
!$$$$$$       ENDIF
!$$$$$$     ENDIF    
ELSEIF((ABS(wf_type).EQ.3))THEN ! for boundary fitted ?
      	!FIXME : to do... (GD) calculate on each point of the FS...
        !voir comment sont formés les maillages mais a priori defini par
        !une matric pour chaque ligne ksi(x,y);eta(x,y) (WWWFB09) donc on calcule
        !pour chaque élement dnas espace physique de départ...
        write(*,*) 'case wf_type=3 not done yet...'
        stop
ELSEIF((ABS(wf_type).GT.3).OR.(curvilinearONOFF.EQ.1))THEN ! wavefield with an angle (straight boundaries) or curvilinear
    !FIXME: take into account fitted boundaries...
    !write(*,*) 'case wf_type>3 not done yet...'
    !stop
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
      ! Make it work with choice equal to 1 or 2 in curvilinear mode...
      IF(wf_type.EQ.1) THEN
        wf_type=360
      ELSEIF(wf_type.EQ.-1) THEN
        wf_type=180
      ELSEIF(wf_type.EQ.2) THEN
        wf_type=90
      ELSEIF(wf_type.EQ.-2) THEN
        wf_type=-90
      ENDIF
    ENDIF
    !
    angle = REAL(wf_type,long)*ACOS(-one)/180.d0
	!
    ! Be careful call eta_nonlinear and phi_nonlinear with non-dimensional time !
    DO i=1,nfx
       DO j=1,nfy
          CALL eta_nonlinear_3D(wavenum,angle,xp(i,j),yp(i,j),omega,tp_adim,n_four_modes,yy,&
              Wavefield%E_I(i,j),Wavefield%Ex_I(i,j),Wavefield%Ey_I(i,j),Wavefield%Exx_I(i,j),&
              Wavefield%Eyy_I(i,j),Wavefield%Et_I(i,j))
          ! Put back in dimensional space...
          Wavefield%E_I(i,j)   = kinv*Wavefield%E_I(i,j)
          Wavefield%Exx_I(i,j) = wavenum*Wavefield%Exx_I(i,j)
          Wavefield%Eyy_I(i,j) = wavenum*Wavefield%Eyy_I(i,j)
          Wavefield%Et_I(i,j)  = ufactor*Wavefield%Et_I(i,j)
          ! Careful, eta^S can be 3D...
          ! Here phi, is computed on the free-surface z=eta^S+eta^I
          ! eta^I and eta^S are dimensional
          !
          zbp = wavenum*(Wavefield%E(i,j)+Wavefield%E_I(i,j))
          !
          CALL phi_nonlinear_3D(wavenum,angle,xp(i,j),yp(i,j),zbp,omega,tp_adim,n_four_modes,zz, &
              Wavefield%Px_I_s(i,j),Wavefield%Py_I_s(i,j),Wavefield%Pz_I_s(i,j),Wavefield%P_I_s(i,j),Wavefield%Pt_I_s(i,j))
          !
          Wavefield%Px_I_s(i,j)= ufactor*Wavefield%Px_I_s(i,j)
          Wavefield%Py_I_s(i,j)= ufactor*Wavefield%Py_I_s(i,j)
          Wavefield%Pz_I_s(i,j)= ufactor*Wavefield%Pz_I_s(i,j)
          Wavefield%P_I_s(i,j) = pfactor*Wavefield%P_I_s(i,j)
          Wavefield%Pt_I_s(i,j)= ufactor**2*Wavefield%Pt_I_s(i,j)
      END DO
    END DO
    !
    ! Determine the incident wavefield on the boundary points...
    !
    Wavefield%SourceEx = -FAC*Wavefield%Ex_I
    Wavefield%SourcePx = -FAC*Wavefield%Px_I_s
    Wavefield%SourceEy = -FAC*Wavefield%Ey_I
    Wavefield%SourcePy = -FAC*Wavefield%Py_I_s
    ! Take into account the eventual ramp...
    zz2(1:2*n_four_modes+10) = zz(1:2*n_four_modes+10)
    DO j=1,n_four_modes
    	zz2(n_four_modes+j+10) = FAC*zz(n_four_modes+j+10)
    ENDDO
    !
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
            ! transform on the real z-space... non-dimensional
            ! FIXME: Here, use E_I or FAC*E_I ?? To check !
            zbp=wavenum*((bottom(i,j)+Wavefield%E(i,j)+FAC*Wavefield%E_I(i,j))*zp(1+GhostGridZ)-bottom(i,j))
            !zbp=d*zp(2)-d
            CALL phi_nonlinear_3D(wavenum,angle,xp(i,j),yp(i,j),zbp,omega,tp_adim,n_four_modes,zz2, &
            		Wavefield%Px_I_bp(Gidx3),Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
            ! Put back in dimensional space...
            Wavefield%Px_I_bp(Gidx3)= ufactor*Wavefield%Px_I_bp(Gidx3)
    		Wavefield%Py_I_bp(Gidx3)= ufactor*Wavefield%Py_I_bp(Gidx3)
    		Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
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
        ! index of the first element
        tmpIdx = Gidx3
        ! Wave elevation + derivative (Dimensional)
        Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(i+GhostGridX,j)
        Wavefield%Ex_I_bp(Gidx3) = FAC*Wavefield%Ex_I(i+GhostGridX,j)
        Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(i+GhostGridX,j)
        DO k = 1, nfz
            ! Wave elevation + derivative
            Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
            Wavefield%Ex_I_bp(Gidx3) = Wavefield%Ex_I_bp(tmpIdx)
            Wavefield%Ey_I_bp(Gidx3) = Wavefield%Ey_I_bp(tmpIdx)
            IF (West_refl==0) THEN
                Wavefield%Px_I_bp(Gidx3)=zero; Wavefield%Py_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
            ELSE
                ! transform on the real z-space... non-dimensional
                ! FIXME: Here, use E_I or FAC*E_I ?? To check !
                zbp=wavenum*((bottom(i+GhostGridX,j)+Wavefield%E(i+GhostGridX,j)+FAC*Wavefield%E_I(i+GhostGridX,j))*zp(k) &
                    -bottom(i+GhostGridX,j))
                !zbp=bottom(i+GhostGridX,j)*zp(k)-bottom(i+GhostGridX,j)
                !zbp=d*zp(k)-d
                CALL phi_nonlinear_3D(wavenum,angle,xp(i+GhostGridX,j),yp(i+GhostGridX,j),zbp,omega,tp_adim,n_four_modes,zz2, &
                 Wavefield%Px_I_bp(Gidx3),Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                 
                ! Put back in dimensional space...
                Wavefield%Px_I_bp(Gidx3)= ufactor*Wavefield%Px_I_bp(Gidx3)
                Wavefield%Py_I_bp(Gidx3)= ufactor*Wavefield%Py_I_bp(Gidx3)
                Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
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
    !
    DO j = 1+GhostGridY, nfy-GhostGridY
        ! index of the first element (=Nz*(Ny-2*GhostY))
        tmpIdx = Gidx3
        ! Wave elevation + derivative (Dimensional)
        Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(i-GhostGridX,j)
        Wavefield%Ex_I_bp(Gidx3) = FAC*Wavefield%Ex_I(i-GhostGridX,j)
        Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(i-GhostGridX,j)
        DO k = 1, nfz
            ! Wave elevation + derivative
            Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
            Wavefield%Ex_I_bp(Gidx3) = Wavefield%Ex_I_bp(tmpIdx)
            Wavefield%Ey_I_bp(Gidx3) = Wavefield%Ey_I_bp(tmpIdx)
            IF (East_refl==0) THEN
                Wavefield%Px_I_bp(Gidx3)=zero; Wavefield%Py_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
            ELSE
                ! transform on the real z-space... non-dimensional
                ! FIXME: Here, use E_I or FAC*E_I ?? To check !
                zbp=wavenum*((bottom(i-GhostGridX,j)+Wavefield%E(i-GhostGridX,j)+FAC*Wavefield%E_I(i-GhostGridX,j))*zp(k) &
                    -bottom(i-GhostGridX,j))
                !zbp=bottom(i-GhostGridX,j)*zp(k)-bottom(i-GhostGridX,j)
                !zbp=d*zp(k)-d
                CALL phi_nonlinear_3D(wavenum,angle,xp(i-GhostGridX,j),yp(i-GhostGridX,j),zbp,omega,tp_adim,n_four_modes,zz2, &
                 Wavefield%Px_I_bp(Gidx3),Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                ! Put back in dimensional space...
                Wavefield%Px_I_bp(Gidx3)= ufactor*Wavefield%Px_I_bp(Gidx3)
                Wavefield%Py_I_bp(Gidx3)= ufactor*Wavefield%Py_I_bp(Gidx3)
                Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
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
    j = 1 
    !  
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
                ! transform on the real z-space... non-dimensional
                ! FIXME: Here, use E_I or FAC*E_I ?? To check !
                zbp=wavenum*((bottom(i,j+GhostGridY)+Wavefield%E(i,j+GhostGridY)+FAC*Wavefield%E_I(i,j+GhostGridY))*zp(k) &
                    -bottom(i,j+GhostGridY))
                !zbp=bottom(i,j+GhostGridY)*zp(k)-bottom(i,j+GhostGridY)
                !zbp=d*zp(k)-d
                CALL phi_nonlinear_3D(wavenum,angle,xp(i,j+GhostGridY),yp(i,j+GhostGridY),zbp,omega,tp_adim,n_four_modes,zz2, &
            		Wavefield%Px_I_bp(Gidx3),Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                ! Put back in dimensional space...
                Wavefield%Px_I_bp(Gidx3)= ufactor*Wavefield%Px_I_bp(Gidx3)
                Wavefield%Py_I_bp(Gidx3)= ufactor*Wavefield%Py_I_bp(Gidx3)
                Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
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
    !  
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
                ! transform on the real z-space... non-dimensional
                ! FIXME: Here, use E_I or FAC*E_I ?? To check !
                zbp=wavenum*((bottom(i,j-GhostGridY)+Wavefield%E(i,j-GhostGridY)+FAC*Wavefield%E_I(i,j-GhostGridY))*zp(k) &
                    -bottom(i,j-GhostGridY))
                !zbp=bottom(i,j-GhostGridY)*zp(k)-bottom(i,j-GhostGridY)
                !zbp=d*zp(k)-d
                CALL phi_nonlinear_3D(wavenum,angle,xp(i,j-GhostGridY),yp(i,j-GhostGridY),zbp,omega,tp_adim,n_four_modes,zz2, &
            		Wavefield%Px_I_bp(Gidx3),Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                ! Put back in dimensional space...
                Wavefield%Px_I_bp(Gidx3)= ufactor*Wavefield%Px_I_bp(Gidx3)
                Wavefield%Py_I_bp(Gidx3)= ufactor*Wavefield%Py_I_bp(Gidx3)
                Wavefield%Pz_I_bp(Gidx3)= ufactor*Wavefield%Pz_I_bp(Gidx3)
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
ELSE
        write(*,*) 'Unknown wave propagation type'
        stop      
ENDIF
!
END SUBROUTINE incident_wf_finite
