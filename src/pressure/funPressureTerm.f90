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
  USE GlobalVariables, ONLY: PressureTermONOFF, Lz, NDampZones, PDampingOnOff, PDampZones, alpha
  IMPLICIT NONE
  INTEGER :: Nx,Ny, i, j, k, nd, rank, id0, job, id, idebug
  REAL(KIND=long) :: t, g, x0, y0, sigma, sigmay, Lx, Ly, dx, Fr, V, rhs(Nx*Ny),            &
       tmp(Nx*Ny), magk
  !REAL(KIND=long), DIMENSION(Nx,Ny), INTENT(OUT) :: term ! includes ghost points
  TYPE (Level_def), INTENT(IN) :: FineGrid
  TYPE (Wavefield_FS), INTENT(OUT) :: Wavefield
  !
  ! The implemented types of pressure patches are:  
  ! 
  !  PressureTermOnOff=0  -> No pressure (default) 
  !  PressureTermOnOff=1  -> Stationary 2D Gaussian hump
  !  PressureTermOnOff=2  -> Moving 3D Gaussian hump
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
     ! The center and variance of the Gaussian.
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
  endif
  !
  ! Pressure damping terms.  
  !
  If (PDampingOnOff /= 0) Then
     id=1
     IF(PDampZones(id)%type==0) THEN
        !
        ! Apply a friction damping pressure to the velocity and/or elevation.
        !
        !
        ! Solve the equation L^2 p_d = Div(gamma Grad Phi) over the 
        ! pressure damping region to get p_d there.
        !  
        nd=PDampZones(id)%nx; 
        id0=PDampZones(id)%idx(1); 
        rank=2*alpha+1

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
         
        ! Take Grad phi 
        !
        CALL DfDx_1D_Uneven(tmp,nd,PDampZones(id)%Grad,rank,rhs)
        
        ! Multiply by gamma and take the divergence
        !
        Do i=1,nd
           tmp(i)=PDampZones(id)%gamPhi(id0+i-1)*rhs(i)
        END Do
        CALL DfDx_1D_Uneven(tmp,nd,PDampZones(id)%Grad,rank,rhs)
        
        ! Add in the Dirichlet condition at the start of the zone and the Neumann 
        ! condition at the end.
        !
        rhs(1) = Wavefield%P(id0,1); 
        rhs(nd)=zero 
        !
        ! Solve for p_d using the LU-factored Laplacian matrix.
        !

        CALL LUSOL(nd,rhs,rhs,PDampZones(id)%Lop%alu,PDampZones(id)%Lop%jlu,PDampZones(id)%Lop%ju)
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
!! changes to take effect and it should be synched with the applied 
!! pressure subroutine "funPressureTerm".  
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
  !  PressureTermOnOff=0  -> No pressure (default) 
  !  PressureTermOnOff=1  -> Stationary 2D Gaussian hump
  !
  !
  If (PressureTermOnOff==1) Then
     !
     ! Apply a stationary Gaussian hump in 2D
     !
     Lx=FineGrid%x(Nx-1,1)-FineGrid%x(2,1); 
     if(Lx .lt. 1.e-14) then 
        Lx=one
     endif
     x0=Lx/2;  sigma=Lx/8;
     !
     DO j = 1, Ny ! ghost points implicitly assumed to be included
        DO i = 1 , Nx ! ghost points implicitly assumed to be included
           ! This following line defines the pressure forcing for each point
           Wavefield%P0(i,j) = -.01*g*Exp(-(FineGrid%x(i,j)-x0)**2/(2*sigma**2))
           Wavefield%E(i,j) = Wavefield%P0(i,j)/g
           Wavefield%P(i,j)=zero; Wavefield%W(i,j)=zero
        END DO
     END DO
  ELSEIf (PressureTermOnOff==2) Then
     !
     ! Apply a stationary Gaussian hump in 3D
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
     ! The variance and center of the Gaussian
     !
     sigma=sqrt(3*dx); x0=3*sigma**2; 
     y0=0; sigmay=sigma/2;
     !
     DO j = 1, Ny ! ghost points implicitly assumed to be included
        DO i = 1 , Nx-1 ! ghost points implicitly assumed to be included
           ! This line defines the pressure forcing for each point
           Wavefield%P0(i,j) =-0.01*g*Exp(-(FineGrid%x(i,j)-x0)**2/(2*sigma**2) &
                -(FineGrid%y(i,j)-y0)**2/(2*sigmay**2))
           Wavefield%E(i,j) = Wavefield%P0(i,j)/g
           Wavefield%P(i,j)=zero; Wavefield%W(i,j)=zero
        END DO
     END DO
   ELSEIf (PressureTermOnOff==3) Then
     !
     ! Apply a stationary tanh hump in 3D
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
     ! The variance and center of the Gaussian
     !
     sigma=sqrt(3*dx); x0=3*sigma**2; 
     y0=0; sigmay=sigma/2;
     !
     DO j = 1, Ny ! ghost points implicitly assumed to be included
        DO i = 1 , Nx-1 ! ghost points implicitly assumed to be included
           ! This line defines the pressure forcing for each point
           ! Wavefield%P0(i,j) =-0.05*g*Exp(-(FineGrid%x(i,j)-x0)**2/(2*sigma**2) &
           !     -(FineGrid%y(i,j)-y0)**2/(2*sigmay**2))
 !! Li and Sclavounos
!If ((FineGrid%x(i,j)-x0) .lt. (2*10.25)) then 
 !           Wavefield%P0(i,j) =-g*0.735*((COS(Pi*(FineGrid%x(i,j)-x0)/(2*10.25)))**2) &
!*((COS(Pi*(FineGrid%y(i,j)-y0)/(2*1.825)))**2)
!else 
!Wavefield%P0(i,j) =0
!endif
!
           ! Sung & Grilli (2005)
           !
           Wavefield%P0(i,j) =-g*0.01*(TANH(FineGrid%x(i,j)-x0+2)-TANH(FineGrid%x(i,j)-x0-2)) &
                *(TANH(FineGrid%y(i,j)-y0+2)-TANH(FineGrid%y(i,j)-y0-2))
           !
           Wavefield%E(i,j) = Wavefield%P0(i,j)/g
           Wavefield%P(i,j)=zero; Wavefield%W(i,j)=zero
        END DO
     END DO
 endif
END SUBROUTINE funInitialFreeSurfaceElevation


