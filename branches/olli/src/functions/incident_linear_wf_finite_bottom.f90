!>
!! Evaluate the linear incident wave field  from
!! linear theory with variable depth (eta, phi and derivatives)
!!
!! FIXME : not working for now ; velocity potential not OK !
!!
!!    Evaluation on the free surface
!!    Evaluation on boundaries
!<
!
SUBROUTINE incident_linear_wf_finite_bottom(wf_type, Wavefield, & !ramp_type,&
	 nfx,xp,nfy,yp,nfz,zp,bottom,bottomx,bottomy,bottomxx,bottomyy,tp,wavenum,grav,H, d)
!
USE PRECISION
USE GlobalVariables, ONLY: GhostGridX,GhostGridY,GhostGridZ, swenseTransientTime, &
     West_refl, East_refl, North_refl, South_refl, curvilinearONOFF
USE DataTypes
USE Constants
!
IMPLICIT NONE
INTERFACE
	SUBROUTINE phi_linear_bottom(H,k,x,z,d,hx,hxx,int_kx,omega,time,grav,phix,phiz,phi,phit)
    USE PRECISION
   	REAL(KIND=long), INTENT(IN)  :: H,k,x,z,d,omega,time,grav,int_kx,hx,hxx
   	REAL(KIND=long), INTENT(OUT) :: phix,phiz
   	REAL(KIND=long), INTENT(OUT), OPTIONAL :: phi,phit
   	END SUBROUTINE phi_linear_bottom
	!
	SUBROUTINE phi_linear_3D_bottom(H,k,angle,x,y,z,d,hx,hy,int_kx,omega,time,grav,phix,phiy,phiz,phi,phit)
    USE PRECISION
   	REAL(KIND=long), INTENT(IN)  :: H,k,angle,x,y,z,d,omega,time,grav,int_kx,hx,hy
   	REAL(KIND=long), INTENT(OUT) :: phix,phiy,phiz
   	REAL(KIND=long), INTENT(OUT), OPTIONAL :: phi,phit
   	END SUBROUTINE phi_linear_3D_bottom
END INTERFACE
!
INTEGER nfx, nfy, nfz, nbp, wf_type
REAL(KIND=long) tp, wavenum, grav, H, d
!REAL(KIND=long), DIMENSION(nfx) :: xp
!REAL(KIND=long), DIMENSION(nfy) :: yp
! For curvilinear
REAL(KIND=long), DIMENSION(nfx,nfy) :: xp,yp
REAL(KIND=long), DIMENSION(nfz) :: zp
REAL(KIND=long), DIMENSION(nfx,nfy) :: bottom,bottomx,bottomy,bottomxx,bottomyy
! For variable depth
REAL(KIND=long), DIMENSION(nfx,nfy) :: k_bottom, H_bottom, int_kx
REAL(KIND=long) :: disper

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
!$$$$$$ omega_t=omega*tp
!      wf_type=1  => Wavefield propagating along x>0,rectangular domain
!      wf_type=-1 => Wavefield propagating along x<0,rectangular domain
!      wf_type=2  => Wavefield propagating along y>0,rectangular domain
!      wf_type=-2 => Wavefield propagating along y<0,rectangular domain
!	   wf_type=3  => boundary-fitted ?
! FIXME: Is it necessary? (GD)
!
! FIXME: Think about the boundary fitted coordinate transformation...
! How to make the incident wave field in a clever way... Compute on all
! free surface points ? Going through an interpolation between 2 grids ?
!
! FIXME: linear case with non-flat bottom ?
omega=wf_type/abs(wf_type)*omega
!$$$$$$ omega_t=wf_type/abs(wf_type)*omega_t
!
!
! Determine the wavenumbers with respect to depth
! Careful: SFsol should give the correct wave frequency!
!
DO i=1,nfx
  DO j=1,nfy
    k_bottom(i,j) = disper(grav,two*pi/abs(omega),bottom(i,j),50)
    k_bottom(i,j) = k_bottom(i,j) / bottom(i,j) !disper gives dimensionless kh
    H_bottom(i,j) = H*SQRT((k_bottom(i,j)*(one+two*wavenum*d/(SINH(two*wavenum*d))))/ &
    	(wavenum*(one+two*k_bottom(i,j)*bottom(i,j)/(SINH(two*k_bottom(i,j)*bottom(i,j))))))
  ENDDO
ENDDO
!
! GD: FIXME: only 2D case here, see for the 3D case
! see the treatment of the integral...
! plus take into account bottomy
!
IF(nfy/=1) THEN
  print*,'3D case not available here...'
ENDIF
!
int_kx(1+GhostGridX,:) = zero
int_kx(1,:) = ((xp(1,1)-xp(2,1))*(k_bottom(1,:)+k_bottom(2,:))/two)
DO i=2+GhostGridX,nfx
    ! compute int(k(x) dx)/x
    int_kx(i,:) = (int_kx(i-1,:)+(xp(i,1)-xp(i-1,1))*(k_bottom(i,:)+k_bottom(i-1,:))/two)
ENDDO
!$$$$$$ DO i=1,nfx
!$$$$$$     print*, i, k_bottom(i,1), int_kx(i,1)
!$$$$$$ ENDDO
!$$$$$$ pause

IF((ABS(wf_type).EQ.1).AND.(curvilinearONOFF.EQ.0))THEN ! propagation along x
  	  IF(nfx==1) THEN
		PRINT*, 'SWENSE Propagation along x with only one point 1 x-point'
        STOP
      ENDIF
	  !     
      DO j=1,nfx
        CALL eta_linear_bottom(H_bottom(j,1),k_bottom(j,1),xp(j,1),int_kx(j,1),omega,tp, &
            Wavefield%E_I(j,1),Wavefield%Ex_I(j,1),Wavefield%Exx_I(j,1),Wavefield%Et_I(j,1))
        !
        ! Here phi, is computed on the free-surface z=0
        !    
        CALL phi_linear_bottom(H_bottom(j,1),k_bottom(j,1),xp(j,1),zero,bottom(j,1),bottomx(j,1),bottomxx(j,1),int_kx(j,1),&
        	omega,tp,grav,Wavefield%Px_I_s(j,1),Wavefield%Pz_I_s(j,1),Wavefield%P_I_s(j,1),Wavefield%Pt_I_s(j,1))
      END DO
      !
      IF(nfy>1) THEN
        DO j=1,nfx
           Wavefield%E_I(j,:)   = Wavefield%E_I(j,1)
           Wavefield%Ex_I(j,:)  = Wavefield%Ex_I(j,1)
           Wavefield%Exx_I(j,:) = Wavefield%Exx_I(j,1)
           Wavefield%Et_I(j,:)  = Wavefield%Et_I(j,1)
           !
           Wavefield%Px_I_s(j,:)= Wavefield%Px_I_s(j,1)
           Wavefield%Pz_I_s(j,:)= Wavefield%Pz_I_s(j,1)
           Wavefield%P_I_s(j,:) = Wavefield%P_I_s(j,1)
           Wavefield%Pt_I_s(j,:)= Wavefield%Pt_I_s(j,1)
           !
           Wavefield%Ey_I(j,:)   = zero
           Wavefield%Eyy_I(j,:)  = zero
           Wavefield%Py_I_s(j,:) = zero
        ENDDO
      ENDIF
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
            ! transform on the real z-space... Dimensional
            zbp=bottom(i,j)*zp(1+GhostGridZ)-bottom(i,j)
            !zbp=d*zp(2)-d
            CALL phi_linear_bottom(FAC*H_bottom(i,1),k_bottom(i,1),xp(i,1),zbp,bottom(i,1),bottomx(i,1),bottomxx(i,1),int_kx(i,1),&
            	omega,tp,grav,Wavefield%Px_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(1,i,j) = Gidx3
            !Wavefield%GidxTableBP(1+GhostGridZ,i,j) = Gidx3
			! Increment global index 
            Gidx3 = Gidx3+1
            print*,'bottom',i,bottomx(i,1)*Wavefield%Px_I_bp(Gidx3)+Wavefield%Pz_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3)
		END DO
	END DO
read*
!    pause
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
                ! transform on the real z-space... Dimensional
                zbp=bottom(i+GhostGridX,j)*zp(k)-bottom(i+GhostGridX,j)
                !zbp=d*zp(k)-d
                CALL phi_linear_bottom(FAC*H_bottom(i+GhostGridX,1),k_bottom(i+GhostGridX,1),bottomx(i+GhostGridX,1),&
                	bottomxx(i+GhostGridX,1),xp(i+GhostGridX,1),zbp,bottom(i+GhostGridX,1),&
                	int_kx(i+GhostGridX,1),omega,tp,grav,Wavefield%Px_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
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
                ! transform on the real z-space... Dimensional
                zbp=bottom(i-GhostGridX,j)*zp(k)-bottom(i-GhostGridX,j)
                !zbp=d*zp(k)-d
                CALL phi_linear_bottom(FAC*H_bottom(i-GhostGridX,1),k_bottom(i-GhostGridX,1),bottomx(i-GhostGridX,1), &
                	bottomxx(i-GhostGridX,1),xp(i-GhostGridX,1),zbp,bottom(i-GhostGridX,1), &
                	int_kx(i-GhostGridX,1),omega,tp,grav,Wavefield%Px_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
            ENDIF
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(k,i,j) = Gidx3
            !Wavefield%GidxTableBP(k,i-GhostGridX,j) = Gidx3
            ! Increment global index        
            Gidx3 = Gidx3+1     
        END DO
    END DO
    !
    ! FIXME: In the case of straight-sided geometries the following is not used...
    ! Keep this or not ?
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
                    ! transform on the real z-space... Dimensional
                    zbp=bottom(i,j+GhostGridY)*zp(k)-bottom(i,j+GhostGridY)
                    !zbp=d*zp(k)-d
                    CALL phi_linear_bottom(FAC*H_bottom(i,1),k_bottom(i,1),xp(i,1),zbp,bottom(i,1),bottomx(i,1),bottomxx(i,1),&
                    	int_kx(i,1),omega,tp,grav,Wavefield%Px_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
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
                  ! transform on the real z-space... Dimensional
                  zbp=bottom(i,j-GhostGridY)*zp(k)-bottom(i,j-GhostGridY)
                  !zbp=d*zp(k)-d
                  CALL phi_linear_bottom(FAC*H_bottom(i,1),k_bottom(i,1),xp(i,1),zbp,bottom(i,1),bottomx(i,1),bottomxx(i,1),&
                  	int_kx(i,1),omega,tp,grav,Wavefield%Px_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
                ENDIF
                ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            	Wavefield%GidxTableBP(k,i,j) = Gidx3
            	!Wavefield%GidxTableBP(k,i,j-GhostGridY) = Gidx3
				! Increment global index 
                Gidx3 = Gidx3+1
            END DO
        END DO
    ENDIF
    ! FIXME: check how to treat corners...
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
!$$$$$$         Wavefield%SourceEx(2+GhostGridX:nfx-1-GhostGridX,1+GhostGridY) = zero
!$$$$$$         Wavefield%SourcePx(2+GhostGridX:nfx-1-GhostGridX,1+GhostGridY) = zero
!$$$$$$       ENDIF
!$$$$$$       IF (North_refl==0) THEN
!$$$$$$         Wavefield%SourceEx(2+GhostGridX:nfx-1-GhostGridX,nfy-GhostGridY) = zero
!$$$$$$         Wavefield%SourcePx(2+GhostGridX:nfx-1-GhostGridX,nfy-GhostGridY) = zero
!$$$$$$       ENDIF
!$$$$$$     ENDIF
ELSEIF((ABS(wf_type).EQ.2).AND.(curvilinearONOFF.EQ.0))THEN ! propagation along y
    !
    IF(nfy==1) THEN
      PRINT*, 'SWENSE Propagation along x with only one point 1 x-point'
      STOP
    ENDIF
    !
    DO j=1,nfy
      CALL eta_linear_bottom(H_bottom(1,j),k_bottom(1,j),yp(1,j),int_kx(1,j),omega,tp, &
          Wavefield%E_I(1,j),Wavefield%Ey_I(1,j),Wavefield%Eyy_I(1,j),Wavefield%Et_I(1,j))
      !
      ! Here phi, is computed on the free-surface z=0
      !    
      CALL phi_linear_bottom(H_bottom(1,j),k_bottom(1,j),yp(1,j),zero,bottom(1,j),bottomy(1,j),bottomyy(1,j),int_kx(1,j), &
      	omega,tp,grav,Wavefield%Py_I_s(1,j),Wavefield%Pz_I_s(1,j),Wavefield%P_I_s(1,j),Wavefield%Pt_I_s(1,j))
    END DO
    !
    IF(nfx>1) THEN
      DO j=1,nfy
         Wavefield%E_I(:,j)   = Wavefield%E_I(1,j)
         Wavefield%Ey_I(:,j)  = Wavefield%Ey_I(1,j)
         Wavefield%Eyy_I(:,j) = Wavefield%Eyy_I(1,j)
         Wavefield%Et_I(:,j)  = Wavefield%Et_I(1,j)
         !
         Wavefield%Py_I_s(:,j)= Wavefield%Py_I_s(1,j)
         Wavefield%Pz_I_s(:,j)= Wavefield%Pz_I_s(1,j)
         Wavefield%P_I_s(:,j) = Wavefield%P_I_s(1,j)
         Wavefield%Pt_I_s(:,j)= Wavefield%Pt_I_s(1,j)
         !
         Wavefield%Ex_I(:,j)   = zero
         Wavefield%Exx_I(:,j)  = zero
         Wavefield%Px_I_s(:,j) = zero
      ENDDO
    ENDIF
	!
    ! Determine the incident wavefield on the boundary points...
    !
    IF (nfx>1) THEN
       Wavefield%Ex_I_bp = zero
       Wavefield%Px_I_bp = zero
       !
       Wavefield%SourceEx = zero
       Wavefield%SourcePx = zero
    ENDIF
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
        	Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(i,j)
            ! transform on the real z-space... Dimensional
            zbp=bottom(i,j)*zp(1+GhostGridZ)-bottom(i,j)
            !zbp=d*zp(k)-d
            CALL phi_linear_bottom(FAC*H_bottom(1,j),k_bottom(1,j),yp(1,j),zbp,bottom(1,j),bottomy(1,j),bottomyy(1,j),int_kx(1,j),&
            	omega,tp,grav,Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(1,i,j) = Gidx3
            !Wavefield%GidxTableBP(1+GhostGridZ,i,j) = Gidx3
            ! Increment global index 
            Gidx3 = Gidx3+1
		END DO
	END DO
    ! FIXME: In the case of straight-sided geometries the two first are not usefull...
    IF (nfx>1) THEN
       ! West  boundary
       !
       i = 1
       !
       DO j = 1+GhostGridY, nfy-GhostGridY
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
                 ! transform on the real z-space... Dimensional
                 zbp=bottom(i+GhostGridX,j)*zp(k)-bottom(i+GhostGridX,j)
                 !zbp=d*zp(k)-d
                 CALL phi_linear_bottom(FAC*H_bottom(1,j),k_bottom(1,j),yp(1,j),zbp,bottom(1,j),bottomy(1,j),bottomyy(1,j), &
                 	int_kx(1,j),omega,tp,grav, Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
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
                 ! transform on the real z-space... Dimensional
                 zbp=bottom(i-GhostGridX,j)*zp(k)-bottom(i-GhostGridX,j)
                 !zbp=d*zp(k)-d
                 CALL phi_linear_bottom(FAC*H_bottom(1,j),k_bottom(1,j),yp(1,j),zbp,int_kx(1,j),omega,bottom(1,j),bottomy(1,j), &
                 	bottomyy(1,j),tp,grav,Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
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
    ! index of the first element (=2*Nz*(Ny-2*GhostY)+...)
    tmpIdx = Gidx3
    ! Wave elevation + derivative (Dimensional)
    Wavefield%E_I_bp(Gidx3)  = FAC*Wavefield%E_I(1,j+GhostGridY)
    Wavefield%Ey_I_bp(Gidx3) = FAC*Wavefield%Ey_I(1,j+GhostGridY)
    DO i = 1+GhostGridX, nfx-GhostGridX
	  	DO k = 1, nfz
			! Wave elevation + derivative
            Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
            Wavefield%Ey_I_bp(Gidx3) = Wavefield%Ey_I_bp(tmpIdx)
            IF (South_refl==0) THEN
				Wavefield%Py_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
            ELSE
               ! transform on the real z-space... Dimensional
               zbp=bottom(i,j+GhostGridY)*zp(k)-bottom(i,j+GhostGridY)
               !zbp=d*zp(k)-d
               CALL phi_linear_bottom(FAC*H_bottom(1,j+GhostGridY),k_bottom(1,j+GhostGridY),yp(1,j+GhostGridY),zbp, &
               	bottom(1,j+GhostGridY), bottomy(1,j+GhostGridY),bottomyy(1,j+GhostGridY),int_kx(1,j+GhostGridY),&
                omega,tp,grav,Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
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
    DO i = 1+GhostGridX, nfx-GhostGridX
	  	DO k = 1, nfz
			! Wave elevation + derivative
            Wavefield%E_I_bp(Gidx3)  = Wavefield%E_I_bp(tmpIdx)
            Wavefield%Ey_I_bp(Gidx3) = Wavefield%Ey_I_bp(tmpIdx)
            IF (North_refl==0) THEN
				Wavefield%Py_I_bp(Gidx3)=zero; Wavefield%Pz_I_bp(Gidx3)=zero
            ELSE
               ! transform on the real z-space... Dimensional
               zbp=bottom(i,j-GhostGridY)*zp(k)-bottom(i,j-GhostGridY)
               !zbp=d*zp(k)-d
               CALL phi_linear_bottom(FAC*H_bottom(1,j-GhostGridY),H_bottom(1,j-GhostGridY),yp(1,j-GhostGridY),zbp,&
               	bottom(1,j-GhostGridY), bottomy(1,j-GhostGridY),bottomyy(1,j-GhostGridY),int_kx(1,j-GhostGridY),&
                omega,tp,grav,Wavefield%Py_I_bp(Gidx3),Wavefield%Pz_I_bp(Gidx3))
            ENDIF
            ! Global index for Boundary Points, which option ? for now cf NormalX,Y,Z
            Wavefield%GidxTableBP(k,i,j) = Gidx3
            !Wavefield%GidxTableBP(k,i,j-GhostGridY) = Gidx3
            ! Increment global index 
            Gidx3 = Gidx3+1
	  	END DO
    END DO
    ! FIXME: check how to treat corners...
    IF (South_refl==0) THEN
      Wavefield%SourceEy(:,1+GhostGridY) = zero
      Wavefield%SourcePy(:,1+GhostGridY) = zero
    ENDIF
    IF (North_refl==0) THEN
      Wavefield%SourceEy(:,nfy-GhostGridY) = zero
      Wavefield%SourcePy(:,nfy-GhostGridY) = zero
    ENDIF
!$$$$$$     IF (nfx>1) THEN
!$$$$$$       IF (West_refl==0) THEN
!$$$$$$         Wavefield%SourceEy(1+GhostGridX,2+GhostGridY:nfy-1-GhostGridY) = zero
!$$$$$$         Wavefield%SourcePy(1+GhostGridX,2+GhostGridY:nfy-1-GhostGridY) = zero
!$$$$$$       ENDIF
!$$$$$$       IF (East_refl==0) THEN
!$$$$$$         Wavefield%SourceEy(nfx-GhostGridX,2+GhostGridY:nfy-1-GhostGridY) = zero
!$$$$$$         Wavefield%SourcePy(nfx-GhostGridX,2+GhostGridY:nfy-1-GhostGridY) = zero
!$$$$$$       ENDIF
!$$$$$$     ENDIF
ELSEIF(ABS(wf_type).EQ.3)THEN ! for boundary fitted ?
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
      DO j=1,nfx
        DO k=1,nfy
	        CALL eta_linear_3D_bottom(H_bottom(j,k),k_bottom(j,k),angle,xp(j,k),yp(j,k),int_kx(j,k),omega,tp,Wavefield%E_I(j,k),&
            	Wavefield%Ex_I(j,k),Wavefield%Ey_I(j,k),Wavefield%Exx_I(j,k),Wavefield%Eyy_I(j,k),Wavefield%Et_I(j,k))
        	!
        	! Here phi, is computed on the free-surface z=0
        	!    
        	CALL phi_linear_3D_bottom(H_bottom(j,k),k_bottom(j,k),angle,xp(j,k),yp(j,k),zero,bottom(j,k),bottomx(j,k),&
            	bottomy(j,k),int_kx(j,k),omega,tp,grav, &
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
            CALL phi_linear_3D_bottom(FAC*H_bottom(i,j),k_bottom(i,j),angle,xp(i,j),yp(i,j),zbp,bottom(i,j),bottomx(i,j),&
            	bottomy(i,j),int_kx(i,j),omega,tp,grav, &
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
                CALL phi_linear_3D_bottom(FAC*H_bottom(i+GhostGridX,j),k_bottom(i+GhostGridX,j),angle,xp(i+GhostGridX,j),&
                	yp(i+GhostGridX,j),zbp,&
                	bottom(i+GhostGridX,j),bottomx(i+GhostGridX,j),bottomy(i+GhostGridX,j),int_kx(i+GhostGridX,j),omega,tp,grav, &
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
                CALL phi_linear_3D_bottom(FAC*H_bottom(i-GhostGridX,j),k_bottom(i-GhostGridX,j),angle,xp(i-GhostGridX,j),&
                	yp(i-GhostGridX,j),zbp,&
                	bottom(i-GhostGridX,j),bottomx(i-GhostGridX,j),bottomy(i-GhostGridX,j),int_kx(i-GhostGridX,j),omega,tp,grav, &
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
                CALL phi_linear_3D_bottom(FAC*H_bottom(i,j+GhostGridY),k_bottom(i,j+GhostGridY),angle,xp(i,j+GhostGridY),&
                	yp(i,j+GhostGridY),zbp, &
                	bottom(i,j+GhostGridY),bottomx(i,j+GhostGridY),bottomy(i,j+GhostGridY),int_kx(i,j+GhostGridY),omega,tp,grav, &
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
              CALL phi_linear_3D_bottom(FAC*H_bottom(i,j-GhostGridY),k_bottom(i,j-GhostGridY),angle,xp(i,j-GhostGridY), &
              	yp(i,j-GhostGridY),zbp, &
              	  bottom(i,j-GhostGridY),bottomx(i,j-GhostGridY),bottomy(i,j-GhostGridY),int_kx(i,j-GhostGridY),omega,tp,grav, &
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
ELSE
        write(*,*) 'Unknown wave propagation type'
        stop
ENDIF

RETURN
END SUBROUTINE incident_linear_wf_finite_bottom
