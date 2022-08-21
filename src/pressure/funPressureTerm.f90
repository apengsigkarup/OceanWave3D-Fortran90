!>
!! Evaluate a pressure forcing terms in the unsteady
!! kinematic and dynamic free surface conditions.
!!
!!    dE/dt = rhsE + Wavefield%NuD
!!    dP/dt = rhsP + Wavefield%P0 + Wavefield%Pd
!!
!! Wavefield%P0 represents an applied pressure patch on the free surface to model 
!! for example a simple ship-like geometry.  It is in general of the form
!!    
!!    Wavefield%P0 = -p/rho
!!
!! with rho the density of the fluid and p the applied free surface pressure.
!!
!! Pd and NuD are pressures corresponding to the application of friction damping 
!! to either the potential and elevation; or the velocity and elevation.  These 
!! pressures are only non-zero in the wave absorption regions defined by the 
!! Data Structure "PDampZones" which is set up in the subroutine 
!! PreprocessPDampingZones.f90. 
!! 
!! This function needs to be updated and recompiled for 
!! changes to take effect.
!!
!! NuD and Pd are friction damping terms applied to absorb waves.  
!<
SUBROUTINE funPressureTerm(t,g,Nx,Ny,FineGrid,Wavefield)
  USE Precision
  USE Constants
  USE DataTypes
  USE GlobalVariables, ONLY: PressureTermONOFF, Lz, NDampZones, PDampingOnOff, PDampZones, alpha, &
       CurrentFlux
  IMPLICIT NONE
  INTEGER :: Nx,Ny, i, j, k, nd, rank, id0, job, id, idebug, i0,i1,j0,j1
  REAL(KIND=long) :: t, g, x0, y0, sigma, sigmay, Lx, Ly, dx, Fr, V, rhs(Nx*Ny),            &
       tmp(Nx*Ny), magk, L, W, draft, a, b, U0, T0, xs, xp, yp, ffunc, gfunc
  !REAL(KIND=long), DIMENSION(Nx,Ny), INTENT(OUT) :: term ! includes ghost points
  TYPE (Level_def), INTENT(IN) :: FineGrid
  TYPE (Wavefield_FS), INTENT(OUT) :: Wavefield
  !
  ! The implemented types of pressure patches are:  
  ! 
  !  PressureTermOnOff=0  -> No pressure (default) 
  !  PressureTermOnOff=1  -> Stationary 2D Gaussian hump
  !  PressureTermOnOff=2  -> Moving 3D Gaussian hump
  !  PressureTermOnOff=3  -> Moving 3D tanh hump
  !  PressureTermOnOff=4  -> Moving 3D cos^2 ship-like form
  !
  !
  If (PressureTermOnOff==1) Then
     !
     ! Apply a 2D stationary Gaussian hump 
     !
     Lx=FineGrid%x(Nx-1,1)-FineGrid%x(2,1); 
     if(Lx .lt. 1.e-14) then 
        Lx=one
     endif
     !
     ! The center and variance of the Gaussian.
     !
     x0=Lx/2;  sigma=Lx/8;
     !
     DO j = 1, Ny ! ghost points implicitly assumed to be included
        DO i = 1 , Nx ! ghost points implicitly assumed to be included
           ! This following line defines the pressure forcing for each point
           Wavefield%P0(i,j) = -.01*g*Exp(-(FineGrid%x(i,j)-x0)**2/(2*sigma**2))
        END DO
     END DO
  ELSEIf (PressureTermOnOff==2) Then
     !
     ! Apply a moving Gaussian hump in 3D 
     !
     Lx=FineGrid%x(Nx-1,1)-FineGrid%x(2,1); dx=FineGrid%x(2,1)-FineGrid%x(1,1);
     Ly=FineGrid%y(1,Ny-1)-FineGrid%y(1,2);
     !
     if(Lx .lt. 1.e-14) then 
        Lx=one
     endif
     if(Ly .lt. 1.e-14) then 
        Ly=one
     endif
     !
     ! The center and variance of the Gaussian.
     !
     sigma=sqrt(3*dx); x0=3*sigma**2+2*t; 
     y0=0; sigmay=sigma/2;
     !
     DO j = 1, Ny ! ghost points implicitly assumed to be included
        DO i = 1 , Nx-1 ! ghost points implicitly assumed to be included
           !
           ! The pressure forcing for each point
           !
           Wavefield%P0(i,j) =-0.01*g*Exp(-(FineGrid%x(i,j)-x0)**2/(2*sigma**2)  &
                -(FineGrid%y(i,j)-y0)**2/(2*sigmay**2))
        END DO
     END DO
  ELSEIf (PressureTermOnOff==3) Then
     !
     ! Apply a moving tanh hump in 3D 
     !
     Lx=FineGrid%x(Nx-1,1)-FineGrid%x(2,1); dx=FineGrid%x(2,1)-FineGrid%x(1,1);
     Ly=FineGrid%y(1,Ny-1)-FineGrid%y(1,2);
     !
     if(Lx .lt. 1.e-14) then 
        Lx=one
     endif
     if(Ly .lt. 1.e-14) then 
        Ly=one
     endif
     !
     ! The center and variance of the function.
     !
     Fr=0.9;
     V=Fr*sqrt(g*Lz);
     sigma=sqrt(3*dx); x0=3*sigma**2+V*t; 
     y0=0; sigmay=sigma/2;
     !
     DO j = 1, Ny ! ghost points implicitly assumed to be included
        DO i = 1 , Nx-1 ! ghost points implicitly assumed to be included
           !
           ! The pressure forcing for each point
           !
           ! Wavefield%P0(i,j) =-0.07*g*Exp(-(FineGrid%x(i,j)-x0)**2/(2*sigma**2)  &
           !     -(FineGrid%y(i,j)-y0)**2/(2*sigmay**2))
           !
           !! Li and Sclavounos
           !If ((FineGrid%x(i,j)-x0-1.3*t) .lt. (2*10.25)) then 
           !            Wavefield%P0(i,j) =-g*0.735*((COS(Pi*(FineGrid%x(i,j)-x0-1.3*t)/(2*10.25)))**2) &
           !*((COS(Pi*(FineGrid%y(i,j)-y0)/(2*1.825)))**2)
           !else 
           !Wavefield%P0(i,j) =0
           !endif
           !
           !! Sung & Grilli (2005)
           !
           Wavefield%P0(i,j) =-g*0.01*(TANH(FineGrid%x(i,j)-x0+2)-TANH(FineGrid%x(i,j)-x0-2)) &
                *(TANH(FineGrid%y(i,j)-y0+2)-TANH(FineGrid%y(i,j)-y0-2))
        END DO
     END DO
  ELSEIf (PressureTermOnOff==4) Then
     !
     ! Apply a moving cos^2 ship-like form in 3D 
     !
     Lx=FineGrid%x(Nx-1,1)-FineGrid%x(2,1); dx=FineGrid%x(2,1)-FineGrid%x(1,1);
     Ly=FineGrid%y(1,Ny-1)-FineGrid%y(1,2);
     !
     if(Lx .lt. 1.e-14) then 
        Lx=one
     endif
     if(Ly .lt. 1.e-14) then 
        Ly=one
     endif
     !
     ! The center and distribution of the moving pressure function.
     !
!     U0=13.5; L=266.; W=42.; draft=6.8; a=0.5; b=0.25; T0=213.; x0=500.; y0=zero
     U0=13.5; L=266.; W=42.; draft=6.8; a=0.5; b=0.25; T0=50.; x0=500.; y0=zero

     if (t<pi/2*T0) then
        V=U0*sin(t/T0)
        xs=x0 + U0*T0*(1-cos(t/T0))
     else
        V=U0
        xs=x0+U0*T0*sin(pi/4)**2+V*(t-pi/2*T0);
     endif
     
     DO j = 1, Ny 
        yp=abs(FineGrid%y(1,j)-y0)
        if (half*b*W <= yp .and. yp <= half*W) then
           gfunc=cos(pi*(yp-half*b*W)/((1-b)*W))**2
        elseif (yp <= half*b*W) then
           gfunc=one
        else
           gfunc=zero
        endif
        DO i = 1 , Nx 
           xp=abs(FineGrid%x(i,j)-xs)
           
           if (half*a*L <= xp .and. xp <= L/2) then
              ffunc=cos(pi*(xp-half*a*L)/((1-a)*L))**2
           elseif (xp <= half*a*L)then
              ffunc=one
           else
              ffunc=zero
           endif
           !
           ! Xinshu Zhang's ship function. The hydrostatic term:
           !
           Wavefield%P0(i,j)=-g*draft*ffunc*gfunc
           !write(it,*)FineGrid%x(i,j),FineGrid%y(i,j),Wavefield%P0(i,j)
        END DO
     END DO
     !
     ! Add in the dynamic correction for t>0.
     !
     i0=nint(xs-half*L); i1=nint(xs+half*L);
     j0=1; j1=nint(half*W);
     if (t>zero)then
        Do j=j0,j1
           do i=i0,i1
              Wavefield%P0(i,j) = (half*(Wavefield%Px(i,j)**2+Wavefield%Py(i,j)**2-Wavefield%W(i,j)**2*(one + Wavefield%Ex(i,j)**2 + Wavefield%Ey(i,j)**2))+g*Wavefield%E(i,j))
              !Wavefield%P0(i,j) = (g*Wavefield%E(i,j))
           end do
        end do
        !
        ! Smooth the pressure patch
        !
        do j=2,Ny-1
              Wavefield%P0(i,j)=0.25_long*Wavefield%P0(i,j-1)+half*Wavefield%P0(i,j)+0.25_long*Wavefield%P0(i,j+1)
           do i=2,Nx-1
              Wavefield%P0(i,j)=0.25_long*Wavefield%P0(i-1,j)+half*Wavefield%P0(i,j)+0.25_long*Wavefield%P0(i+1,j)
           end do
        end do
     end if
  endif
  !
  ! Pressure damping terms.  
  !
  If (PDampingOnOff /= 0 .and. t>0) Then
     id=1
     IF(PDampZones(id)%type==0) THEN
        !
        ! Apply a friction damping pressure to the velocity and/or elevation.
        !
        !
        ! Solve the equation L^2 p_d = Div(gamma Grad Phi) over the 
        ! pressure damping region to get p_d there.
        !  
        nd=PDampZones(id)%nx; id0=PDampZones(id)%idx(1); rank=2*alpha+1
        !
        ! A test case for debugging...  
        idebug=0
        IF (idebug==1)THEN
           magk=10.*two*pi/(FineGrid%x(Nx,1)-FineGrid%x(1,1))
           do i=1,Nx
              Wavefield%P(i,1)=cos(magk*(FineGrid%x(i,1)-FineGrid%x(PDampZones(id)%idx(2),1)))
           end do
        end IF
        ! 
        ! Extract phi over the damping zone region
        !
        tmp(1:nd)=Wavefield%P(id0:id0+nd-1,1)  
        ! 
        ! Take Grad phi 
        !
        CALL DfDx_1D_Uneven(tmp,nd,PDampZones(id)%Grad,rank,rhs)
        !
        ! Multiply by gamma and take the divergence. Note here that we damp on the
        ! difference between the flow velocity and the possible uniform current. 
        !
        Do i=1,nd
           tmp(i)=PDampZones(id)%gamPhi(id0+i-1)*(rhs(i)-CurrentFlux%Ur)     
        END Do
        CALL DfDx_1D_Uneven(tmp,nd,PDampZones(id)%Grad,rank,rhs)
        !
        ! Add in the Dirichlet condition at the start of the zone and the Neumann 
        ! condition at the end.
        !
        rhs(1) = Wavefield%P(id0,1); rhs(nd)=zero 
        !
        ! Solve for p_d using the LU-factored Laplacian matrix.
        !
        JOB = 3;   ! Solution
        CALL MA41AD(JOB, nd, PDampZones(id)%Lop%nnz, PDampZones(id)%Lop%irn, PDampZones(id)%Lop%icn,    &
             PDampZones(id)%Lop%val, rhs, PDampZones(id)%Lop%COLSCA, PDampZones(id)%Lop%ROWSCA,         &
             PDampZones(id)%Lop%KEEP, PDampZones(id)%Lop%IS_HSL, PDampZones(id)%Lop%MAXIS,              &
             PDampZones(id)%Lop%SS, PDampZones(id)%Lop%MAXS, PDampZones(id)%Lop%CNTL,                   &
             PDampZones(id)%Lop%ICNTL, PDampZones(id)%Lop%INFOHSL, PDampZones(id)%Lop%RINFO)

        if (PDampZones(id)%Lop%INFOHSL(1) .LT. 0) then ! Check for problems
           print *, 'Problems with MA41, JOB = 3 (Solution), in funPressureTerm.f90'
        end if
        !
        ! Apply friction damping to the elevation (kinematic FSBC). 
        !
        j=1
        Do i=1,Nx
           k=(j-1)*Nx+i
           Wavefield%NuD(i,j)=-PDampZones(id)%gamEta(k)*Wavefield%E(i,j)
        END Do
        !
        ! Apply p_d to the dynamic FSBC to apply friction damping to the velocity.  
        !
        !
        ! Set the pressure at the first point of the zone to zero to remove the arbitrary 
        ! constant of integration.  Then apply the pressure over the damping zone.  
        !
        Wavefield%Pd(1:Nx,1:Ny)=zero;
        Do i=1,nd
           Wavefield%Pd(id0+i-1,j)=-(rhs(i)-rhs(1))
        END Do
        !
        !
        IF(idebug==1)THEN
           Do i=1,nx
              write(202,*)FineGrid%x(i,1),WaveField%P(i,1),Wavefield%Pd(i,1)
           end Do
           stop
        END IF
     ELSE
        !
        ! Apply a friction damping pressure to the potential and/or elevation.
        !
        Do j=1,Ny
           Do i=1,Nx
              k=(j-1)*Nx+i
              Wavefield%NuD(i,j)=-PDampZones(id)%gamEta(k)*Wavefield%E(i,j)
              Wavefield%Pd(i,j)=-PDampZones(id)%gamPhi(k)*Wavefield%P(i,j)
           END Do
        END Do
     END IF
  END If
END SUBROUTINE funPressureTerm

!>
!! Evaluate the initial free surface elevation and the applied pressure.  
!!
!!    dP/dt = rhsP + Wavefield%P0
!!
!! Wavefield%P0 is in general of the form
!!    
!!    Wavefield%P0 = -p/rho
!!
!! with rho the density of the fluid and p the free surface pressure.
!!
!! This function needs to be updated and recompiled for 
!! changes to take effect. 
!<
SUBROUTINE funInitialFreeSurfaceElevation(g,Nx,Ny,FineGrid,WaveField)
  USE Precision
  USE Constants
  USE DataTYpes
  USE GlobalVariables, ONLY: PressureTermONOFF
  IMPLICIT NONE
  INTEGER :: Nx,Ny, i, j
  REAL(KIND=long) :: g, x0, y0, sigma, sigmay, Lx, Ly, dx
  TYPE (Level_def), INTENT(IN) :: FineGrid
  TYPE (Wavefield_FS), INTENT(OUT) :: Wavefield
  !
  ! 
  !  This routine should only be called when PressureTerOnOff/=0.
  !
  ! Get the initial pressure and compute the elevation. 
  !
  !
  call funPressureTerm(zero,g,Nx,Ny,FineGrid,Wavefield)
  
  DO j = 1, Ny 
     DO i = 1 , Nx 
        Wavefield%E(i,j) = Wavefield%P0(i,j)/g
        Wavefield%P(i,j)=zero; Wavefield%W(i,j)=zero
     END DO
  END DO

END SUBROUTINE funInitialFreeSurfaceElevation


