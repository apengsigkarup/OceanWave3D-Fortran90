!>
!! Evaluate a pressure forcing term in the unsteady
!! dynamic free surface condition.
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
SUBROUTINE funPressureTerm(t,g,Nx,Ny,FineGrid,Wavefield)
  USE Precision
  USE Constants
  USE DataTypes
  USE GlobalVariables, ONLY: PressureTermONOFF, Lz
  IMPLICIT NONE
  INTEGER :: Nx,Ny, i, j
  REAL(KIND=long) :: t, g, x0, y0, sigma, sigmay, Lx, Ly, dx, Fr, V
  !REAL(KIND=long), DIMENSION(Nx,Ny), INTENT(OUT) :: term ! includes ghost points
  TYPE (Level_def), INTENT(IN) :: FineGrid
  TYPE (Wavefield_FS), INTENT(OUT) :: Wavefield
  !
  ! 
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


