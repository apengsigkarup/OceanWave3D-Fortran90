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
SUBROUTINE funPressureTerm(t,g,Nx,Ny,Nz,FineGrid,Wavefield)
  USE Precision
  USE Constants
  USE DataTypes
  USE GlobalVariables, ONLY: PressureTermONOFF,Lz,phi,alpha,GhostGridX,&
  GhostGridY,GhostGridZ,Lship,Uship,Dship,x0ship

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! Input:
  !-----------------------------------------------------------------------------
  INTEGER                       :: Nx,Ny,Nz
  REAL(KIND=long)               :: t, g
  TYPE (Level_def), INTENT(IN)  :: FineGrid

  !-----------------------------------------------------------------------------
  ! Output:
  !-----------------------------------------------------------------------------
  TYPE (Wavefield_FS), INTENT(OUT) :: Wavefield

  !-----------------------------------------------------------------------------
  ! Local parameters:
  !-----------------------------------------------------------------------------
  INTEGER           :: i,j
  REAL(KIND=long)   :: x0, y0,z0, sigma, sigmay, Lx, Ly, dx
  REAL(KIND=long)   :: dphidx(Nz,Nx,Ny)
   

  !-----------------------------------------------------------------------------
  ! Set the pressure distribution on the free surface:
  !-----------------------------------------------------------------------------
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

    !---------------------------------------------------------------------------
    ! Cosine function in x-direction: 
    !---------------------------------------------------------------------------
    DO j = 1, Ny ! ghost points implicitly assumed to be included
      DO i = 1 , Nx-1 ! ghost points implicitly assumed to be included
        if (abs(FineGrid%x(i,j)-x0ship) .lt. Lship/2.0)  then
          Wavefield%P0(i,j) = &
          -Dship*g*(cos(pi*(FineGrid%x(i,j)-x0ship)/Lship)**2.0)
        endif
      END DO
    END DO

  ELSEIf (PressureTermOnOff==3) Then
    !
    ! Apply a moving Gaussian hump in 3D 
    !
    Lx=FineGrid%x(Nx-1,1)-FineGrid%x(2,1); 
    dx=FineGrid%x(2,1)-FineGrid%x(1,1);
    Ly=FineGrid%y(1,Ny-1)-FineGrid%y(1,2);
    !
    if(Lx .lt. 1.e-14) then 
      Lx=one
    endif
    if(Ly .lt. 1.e-14) then 
      Ly=one
    endif

    !--------------------------------------------------------------------------
    ! Calculation of pressure term dphi/dt = - Uship*dphi/dx:
    !--------------------------------------------------------------------------
    ! Compute dphi/dx:
    CALL DiffXEven(phi,dphidx,1, FineGrid%Nx+2*GhostGridX,    &
    FineGrid%Ny+2*GhostGridY,         &
    FineGrid%Nz+GhostGridZ,           &
    FineGrid%DiffStencils,alpha)
    !
    sigma  =   sqrt(Lship) 
    sigmay =   sqrt(Lship)
    x0     =   40 + 0*2*t; 
    y0     =   0;
    z0     =   sin(2.0*3.14*t/1.0)
    !
    DO j = 1, Ny ! ghost points implicitly assumed to be included
      DO i = 1 , Nx-1 ! ghost points implicitly assumed to be included
        !--------------------------------------------------------------------
        ! The pressure forcing for each point
        !--------------------------------------------------------------------
        if ((abs(FineGrid%x(i,j)-x0) .lt. Lship/2.0) .AND. &
            (abs(FineGrid%y(i,j)) .lt. 4.0/2.0 )) then
            Wavefield%P0(i,j) = &
        -0.1*g*(cos(pi*(FineGrid%x(i,j)-x0)/Lship)**2.0) &
        *(cos(pi*FineGrid%y(i,j)/4.0)**2.0)
      endif
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
USE GlobalVariables, ONLY: PressureTermONOFF,Lship,Bship,Dship,x0ship,Fnship,&
Uship,shipfile,wgt_ship,GhostGridX,GhostGridY,alpha,beta,nx_ship,dzdxi1,dzdxi2

IMPLICIT NONE 

!-------------------------------------------------------------------------------
! Input / output:
!-------------------------------------------------------------------------------
INTEGER                             :: Nx,Ny
REAL(KIND=long)                     :: g
TYPE (Level_def), INTENT(IN)        :: FineGrid
TYPE (Wavefield_FS), INTENT(OUT)    :: Wavefield

!-------------------------------------------------------------------------------
! Local parameters:
!-------------------------------------------------------------------------------
INTEGER             :: i,j,i0ship,Nxship,Nyship,Nst,i1,j1,iter,flag_filter
REAL(KIND=long)     :: x0,y0,sigma,Lx,zship,scale,tol,dx1,dx2,dy1,dy2,w1,w2,dist
REAL(KIND=long)     :: dist_smooth,tmp1,tmp2,xmin,xmax,ymax,dxi1,dxi2,dxdxi1,dydxi1,dxdxi2,dydxi2,lng
CHARACTER(LEN=256)  :: dummy
INTEGER,         dimension(:),   ALLOCATABLE :: Nstpnt
REAL(kind=long) :: n1,n2,n3
REAL(KIND=long), dimension(:,:), ALLOCATABLE :: xst,yst,zst
REAL(KIND=long), dimension(:), ALLOCATABLE   :: xwl,ywl
REAL(KIND=long)     ::  nx_edge1,ny_edge1,nx_edge2,ny_edge2,dp1,dp2,y1,y2,z1,z2



!-------------------------------------------------------------------------------
! Set pressure on the free surface:
!-------------------------------------------------------------------------------
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

    !--------------------------------------------------------------------------
    ! Cosine function in x-direction: 
    !--------------------------------------------------------------------------
    Fnship  =   0.19_long
    Lship   =   0.25_long
    Dship   =   0.01_long
    x0ship  =   0.75_long
    Uship   =   Fnship*((g*Lship)**(0.5_long)) 
    
    DO j = 1,Ny ! ghost points implicitly assumed to be included
      DO i = 1,Nx ! ghost points implicitly assumed to be included
        ! This line defines the pressure forcing for each point
        if (abs(FineGrid%x(i,j)-x0ship) .lt. Lship/2.0) then
          Wavefield%P0(i,j) = &
          -Dship*g*(cos(pi*(FineGrid%x(i,j)-x0ship)/Lship)**2.0) 
        endif
        Wavefield%E(i,j) = Wavefield%P0(i,j)/g
        print*,i,j,Wavefield%E(i,j) 
        Wavefield%P(i,j) = zero 
        Wavefield%W(i,j) = zero
      END DO
    END DO

ELSEIf (PressureTermOnOff==3) Then

    !--------------------------------------------------------------------------
    ! Cosine function in both x and y directions:
    !--------------------------------------------------------------------------
    Fnship  =   0.30_long
    Lship   =   0.25_long
    Bship   =   0.05_long
   Dship   =   0.01_long
    x0ship  =   0.75_long
    Uship   =   Fnship*((g*Lship)**(0.5_long)) 
  DO j = 1, Ny ! ghost points implicitly assumed to be included
    DO i = 1 , Nx-1 ! ghost points implicitly assumed to be included
      ! This line defines the pressure forcing for each point
      if ((abs(FineGrid%x(i,j)-x0ship) .lt. Lship/2.0) .AND. &
          (abs(FineGrid%y(i,j))        .lt. Bship/2.0 )) then

      Wavefield%P0(i,j) = &
      -Dship*g*(cos(pi*(FineGrid%x(i,j)-x0)/Lship)**2.0) &
      *(cos(pi*FineGrid%y(i,j)/Bship)**2.0)

    endif
    Wavefield%E(i,j) = Wavefield%P0(i,j)/g
    Wavefield%P(i,j) = zero 
    Wavefield%W(i,j) = zero
  END DO
END DO

ELSEIf (PressureTermOnOff==4) Then
  
    !---------------------------------------------------------------------------
    ! Open ship file:
    !---------------------------------------------------------------------------
    OPEN(UNIT=9,FILE=shipfile)

    !---------------------------------------------------------------------------
    ! Read ship grid parameters:
    !---------------------------------------------------------------------------
    READ(UNIT=9,FMT='(I10)')    i0ship 
    READ(UNIT=9,FMT='(I10)')    Nxship 
    READ(UNIT=9,FMT='(I10)')    Nyship 
    READ(UNIT=9,FMT=101)        Lship 
    READ(UNIT=9,FMT=101)        Fnship 

    !---------------------------------------------------------------------------
    ! Calculate ship parameters:
    !---------------------------------------------------------------------------
    Uship   =   Fnship*((g*Lship)**(0.5_long)) 

    !---------------------------------------------------------------------------
    ! Write the ship grid parameters to screen:
    !---------------------------------------------------------------------------
    WRITE(*,'(A,I10)')  'i0ship: ',i0ship
    WRITE(*,'(A,I10)')  'Nxship: ',Nxship
    WRITE(*,'(A,I10)')  'Nyship: ',Nyship
    WRITE(*,'(A,F8.4)') 'Lship:  ',Lship 
    WRITE(*,'(A,F8.4)') 'Fnship: ',Fnship 
    WRITE(*,'(A,F8.4)') 'Uship:  ',Uship 

    !---------------------------------------------------------------------------
    ! Read ship elevations:
    !---------------------------------------------------------------------------
    DO i = 1,Nxship ! ghost points implicitly assumed to be included
      DO j = 1,Nyship ! ghost points implicitly assumed to be included
        READ(UNIT=9,FMT=101) zship 
        Wavefield%P0(i+i0ship,j)   = g*zship
        Wavefield%E( i+i0ship,j)   = Wavefield%P0(i+i0ship,j)/g
        Wavefield%P( i+i0ship,j)   = zero 
        Wavefield%W( i+i0ship,j)   = zero
      END DO
    END DO

    101 FORMAT (F8.4)
ELSEIf (PressureTermOnOff==5) Then

!--------------------------------------------------------------------------
Fnship  =   0.35_long
Lship   =   0.10_long
Bship   =   0.20_long
Dship   =   0.01_long
x0ship  =   0.5_long
Uship   =   Fnship*((g*Lship)**(0.5_long)) 
DO j = 1, Ny ! ghost points implicitly assumed to be included
  DO i = 1 , Nx-1 ! ghost points implicitly assumed to be included
    ! This line defines the pressure forcing for each point
!    if ((abs(FineGrid%x(i,j)-x0ship) .lt. Lship/2.0) .AND. &
!    (    abs(FineGrid%y(i,j)-x0ship) .lt. Bship/2.0)) then
    if ((abs(FineGrid%x(i,j)-x0ship) .lt. Lship/2.0+(1e-6)) .and. &
        (abs(FineGrid%y(i,j)-x0ship) .lt. Lship/2.0-(1e-6))) then
    Wavefield%P0(i,j) = -Dship*g
    wgt_ship(i,j) = 1
  else
    wgt_ship(i,j) = 0

  endif
  Wavefield%E(i,j) = Wavefield%P0(i,j)/g
  Wavefield%P(i,j) = zero 
  Wavefield%W(i,j) = zero
END DO
END DO

elseif (pressureTermOnOff==6) then

  
write(*,'(A)') 'Ship type: Series 69, CB 0.6'

    !---------------------------------------------------------------------------
    ! Open ship file:
    !---------------------------------------------------------------------------
    OPEN(UNIT=9,FILE='shipfile')

    !---------------------------------------------------------------------------
    ! Read ship grid parameters:
    !---------------------------------------------------------------------------
    READ(UNIT=9,FMT=*)        Fnship 
    READ(UNIT=9,FMT=*)        x0
    READ(UNIT=9,FMT=*)        y0 

    !---------------------------------------------------------------------------
    ! Write the ship grid parameters to screen:
    !---------------------------------------------------------------------------
    WRITE(*,'(A,F8.4)') 'Fnship: ',Fnship 
    WRITE(*,'(A,F8.4)') 'x0:     ',x0
    WRITE(*,'(A,F8.4)') 'y0:     ',y0 

    tol = 1e-6
  !-----------------------------------------------------------------------------
  ! Series 60 ship parameters:
  !-----------------------------------------------------------------------------
  Lship   =   1.0_long
  Uship   =   Fnship*((g*Lship)**(0.5_long)) 
  scale   =   1.0

  !-----------------------------------------------------------------------------
  ! Open series 60 file:
  !-----------------------------------------------------------------------------
  OPEN(UNIT=9,FILE='series60_06.pln')

  !---------------------------------------------------------------------------
  ! Read file:
  !---------------------------------------------------------------------------
  READ(UNIT=9,FMT='(A)')    dummy
  READ(UNIT=9,FMT=*) 
  READ(UNIT=9,FMT=*) Lship,x0ship
  READ(UNIT=9,FMT=*) Dship
  READ(UNIT=9,FMT=*) Nst

  !-----------------------------------------------------------------------------
  ! Write parameters to screen:
  !-----------------------------------------------------------------------------
  write(*,'(A,F8.4)')   'Lship              = ',Lship
  write(*,'(A,F8.4)')   'x0ship             = ',x0ship
  write(*,'(A,F8.4)')   'Dship              = ',Dship
  write(*,'(A,I4)')     'Number of stations = ',Nst

  !-----------------------------------------------------------------------------
  ! Allocate and read stations:
  !-----------------------------------------------------------------------------
  ALLOCATE(Nstpnt(Nst))
  ALLOCATE(xst(Nst,20))
  ALLOCATE(yst(Nst,20))
  ALLOCATE(zst(Nst,20))
  DO i=1,Nst
    READ(UNIT=9,FMT=*) Nstpnt(i)
    DO j=1,Nstpnt(i)
      READ(UNIT=9,FMT=*) xst(i,j),yst(i,j),zst(i,j)
!      write(*,FMT=*) xst(i,j),yst(i,j),zst(i,j)
    enddo
  ENDDO

  !-----------------------------------------------------------------------------
  ! Allocate and real waterline:
  !-----------------------------------------------------------------------------
  ALLOCATE(xwl(26))
  ALLOCATE(ywl(26))
  OPEN(UNIT=9,FILE='series60BC06waterline')
  DO i=1,26
      READ(UNIT=9,FMT=*) xwl(i),ywl(i)
      write(*,FMT=*) xwl(i),ywl(i)
  ENDDO

  !-----------------------------------------------------------------------------
  ! Make ship weight function based on distance to waterline:
  !-----------------------------------------------------------------------------
!  ALLOCATE(wgt_ship(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
!  dist_smooth = 0.0
  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
      wgt_ship(i,j)     = 0.0
      DO i1=1,25

        if (xwl(i1)-tol    .lt. FineGrid%x(i,j)-x0 .and. &
          FineGrid%x(i,j)-x0 .lt. xwl(i1+1)+tol) then 
          nx_edge1 = - (ywl(i1+1) - ywl(i1)) 
          ny_edge1 =   (xwl(i1+1) - xwl(i1)) 

          dx1 = 0.5*(xwl(i1)+xwl(i1+1))-(FineGrid%x(i,j)-x0)  
          dy1 = 0.5*(ywl(i1)+ywl(i1+1))-(FineGrid%y(i,j)-y0)

          dp1 = nx_edge1*dx1 + ny_edge1*dy1
          if (dp1 .gt. 0.0) then
            dist = dp1/sqrt(nx_edge1**2.0 + ny_edge1**2.0)
            if (dist .lt. 2.0*dist_smooth) then
!              wgt_ship(i,j)=0.5 + (dist-dist_smooth)/(2*dist_smooth) &
!              + 1/(2*pi)*sin(pi*(dist-dist_smooth)/dist_smooth) 
            else
              wgt_ship(i,j) = 1.0
            endif
          endif

      endif
      enddo
    enddo
  enddo

  !-----------------------------------------------------------------------------
  ! Loop all grid points:
  !-----------------------------------------------------------------------------
  DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
    DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX

      Wavefield%P0(i,j) = 0.0
!      wgt_ship(i,j)     = 0

      if ( x0-0.5 .lt. FineGrid%x(i,j) .and. FineGrid%x(i,j) .lt. x0+1.5 .and. &
           y0-tol .lt. FineGrid%y(i,j) .and. FineGrid%y(i,j) .lt. y0+tol) then

        !-----------------------------------------------------------------------
        ! The center line of the ship:
        !-----------------------------------------------------------------------

        !--------------------------------------------------------------------
        ! Loop stations (first and last station are bow and stern):
        !--------------------------------------------------------------------
        do i1 = 2,Nst-1
          if ( xst(i1+1,1)+x0ship-tol .lt. FineGrid%x(i,j)-x0  .and. &
               xst(i1,1)  +x0ship+tol .gt. FineGrid%x(i,j)-x0) then

            !-------------------------------------------------------------------
            ! Linear interpolation:
            !-------------------------------------------------------------------
            dx1 = (FineGrid%x(i,j)-x0) - (xst(i1+1,1)+x0ship)
            dx2 = (xst(i1,1)  +x0ship) - (FineGrid%x(i,j)-x0)
            w1 = dx2/(dx1+dx2)
            w2 = dx1/(dx1+dx2)
            Wavefield%P0(i,j) = (w1*zst(i1+1,1) + w2*zst(i1,1))*g
!            wgt_ship(i,j)     = 1

          endif
        enddo

        !-----------------------------------------------------------------------
        ! Stern station:
        !-----------------------------------------------------------------------
        i1 = Nst
        do j1 = 1,Nstpnt(i1)-1
          if ( xst(i1,j1+1)+x0ship-tol .lt. FineGrid%x(i,j)-x0  .and. &
               xst(i1,j1)  +x0ship+tol .gt. FineGrid%x(i,j)-x0) then

            !-------------------------------------------------------------------
            ! Linear interpolation:
            !-------------------------------------------------------------------
            dx1 = (FineGrid%x(i,j)-x0) - (xst(i1,j1+1)+x0ship)
            dx2 = (xst(i1,j1)  +x0ship) - (FineGrid%x(i,j)-x0)
            w1 = dx2/(dx1+dx2)
            w2 = dx1/(dx1+dx2)
            Wavefield%P0(i,j) = (w1*zst(i1,j1+1) + w2*zst(i1,j1))*g
!            wgt_ship(i,j)     = 1
          endif
        enddo

      elseif(x0-0.5.lt.FineGrid%x(i,j) .and. FineGrid%x(i,j) .lt. x0+1.5 .and. &
             y0-0.5 .lt. FineGrid%y(i,j) .and. FineGrid%y(i,j) .lt. y0+0.5) then

        !-----------------------------------------------------------------------
        ! Inside a box around ship with size 2*Lship x 1*Lship:
        !-----------------------------------------------------------------------

        !--------------------------------------------------------------------
        ! Loop stations (first and last station are bow and stern):
        !--------------------------------------------------------------------
        do i1 = 2,Nst-4
          if ( xst(i1+1,1)+x0ship-tol .lt. FineGrid%x(i,j)-x0  .and. &
               xst(i1,1)  +x0ship+tol .gt. FineGrid%x(i,j)-x0) then

          do j1 = 1,Nstpnt(i1)-1
            
            !-------------------------------------------------------------------
            ! Check if grid point is inside quadrilateral:
            !-------------------------------------------------------------------
            nx_edge1 =  yst(i1,j1)   - yst(i1+1,j1) 
            ny_edge1 =  xst(i1+1,j1) - xst(i1,j1) 

            nx_edge2 =  yst(i1+1,j1+1) - yst(i1,j1+1) 
            ny_edge2 =  xst(i1,j1+1)   - xst(i1+1,j1+1) 

            dx1 = xst(i1+1,j1)+x0ship - (FineGrid%x(i,j)-x0)  
            dy1 = yst(i1+1,j1)        - abs(FineGrid%y(i,j)-y0)

            dx2 = xst(i1+1,j1+1)+x0ship - (FineGrid%x(i,j)-x0)  
            dy2 = yst(i1+1,j1+1)        - abs(FineGrid%y(i,j)-y0)

            dp1 = nx_edge1*dx1 + ny_edge1*dy1
            dp2 = nx_edge2*dx2 + ny_edge2*dy2


            if (dp1 .gt. 0.0 .and. dp2 .gt. 0.0) then

              dx1 = (FineGrid%x(i,j)-x0)  - (xst(i1+1,j1)+x0ship)
              dx2 = (xst(i1,j1)  +x0ship) - (FineGrid%x(i,j)-x0)
              w1 = dx2/(dx1+dx2)
              w2 = dx1/(dx1+dx2)
            
              y1 = (w1*yst(i1+1,j1)   + w2*yst(i1,j1)  )
              y2 = (w1*yst(i1+1,j1+1) + w2*yst(i1,j1+1))

              z1 = (w1*zst(i1+1,j1)   + w2*zst(i1,j1))
              z2 = (w1*zst(i1+1,j1+1) + w2*zst(i1,j1+1))

              dy1 = abs(FineGrid%y(i,j)-y0)    - y1 
              dy2 = y2                      - abs(FineGrid%y(i,j)-y0)
              w1 = dy2/(dy1+dy2)
              w2 = dy1/(dy1+dy2)

              Wavefield%P0(i,j) = (w1*z1 + w2*z2)*g
!            wgt_ship(i,j)     = 1
            endif
          enddo
        endif
      enddo
    else
        Wavefield%P0(i,j) = 0.0
!        wgt_ship(i,j)     = 0
      endif
        !-----------------------------------------------------------------------
        ! If pressure is negative set it to zero (the case for the free board of
        ! the ship which is above the water line):
        !-----------------------------------------------------------------------
        if (Wavefield%P0(i,j) .gt. 0.0) then
            Wavefield%P0(i,j) = 0.0 
!            wgt_ship(i,j)     = 0
        endif
!Wavefield%P0(i,j) = wgt_ship(i,j)*Wavefield%P0(i,j)
! Wavefield%E(i,j) = Wavefield%P0(i,j)/g
!  Wavefield%P(i,j) = zero 
!  Wavefield%W(i,j) = zero
    END DO
  END DO

  !-----------------------------------------------------------------------------
  ! Do a light smoothing of the ship geometry. "Just to break the corners".
  ! Maybe use a even lighter filter with only three points
  !-----------------------------------------------------------------------------
  flag_filter = 1
  do iter = 1,4
    DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
      DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX

        if(x0-0.5.lt.FineGrid%x(i,j) .and. FineGrid%x(i,j) .lt. x0+1.5 .and. &
        y0-0.5 .lt. FineGrid%y(i,j) .and. FineGrid%y(i,j) .lt. y0+0.5) then

        !-----------------------------------------------------------------------
        ! Inside a box around ship with size 2*Lship x 1*Lship:
        !-----------------------------------------------------------------------

        !--------------------------------------------------------------------
        ! Filter:
        !--------------------------------------------------------------------
        if (flag_filter .eq. 1) then
        tmp1 = Wavefield%P0(i-1,j)*0.010867541574776 &
        + Wavefield%P0(i,j)*0.978264916850449   &
        + Wavefield%P0(i+1,j)*0.010867541574776 

        if (j==2) then
          tmp2 =  Wavefield%P0(i,j+1)*0.010867541574776 &
          + Wavefield%P0(i,j)*0.978264916850449   &
          + Wavefield%P0(i,j+1)*0.010867541574776 

        else
          tmp2 =  Wavefield%P0(i,j-1)*0.010867541574776 &
          + Wavefield%P0(i,j)*0.978264916850449   &
          + Wavefield%P0(i,j+1)*0.010867541574776 

        endif

        elseif (flag_filter .eq. 2) then

        tmp1 = Wavefield%P0(i-2,j)*0.000591270361240 &
        + Wavefield%P0(i-1,j)*0.125419103234394 &
        + Wavefield%P0(i,j)*0.747979252808734   &
        + Wavefield%P0(i+1,j)*0.125419103234393 &
        + Wavefield%P0(i+2,j)*0.000591270361240

        if (j==2) then
          tmp2 = Wavefield%P0(i,j+2)*0.000591270361240 &
          + Wavefield%P0(i,j+1)*0.125419103234394 &
          + Wavefield%P0(i,j)*0.747979252808734   &
          + Wavefield%P0(i,j+1)*0.125419103234393 &
          + Wavefield%P0(i,j+2)*0.000591270361240

          elseif (j==3) then
          tmp2 = Wavefield%P0(i,j)*0.000591270361240 &
          + Wavefield%P0(i,j-1)*0.125419103234394 &
          + Wavefield%P0(i,j)*0.747979252808734   &
          + Wavefield%P0(i,j+1)*0.125419103234393 &
          + Wavefield%P0(i,j+2)*0.000591270361240

        else
          tmp2 = Wavefield%P0(i,j-2)*0.000591270361240 &
          + Wavefield%P0(i,j-1)*0.125419103234394 &
          + Wavefield%P0(i,j)*0.747979252808734   &
          + Wavefield%P0(i,j+1)*0.125419103234393 &
          + Wavefield%P0(i,j+2)*0.000591270361240

        endif
      endif
        Wavefield%P0(i,j) = 0.5*(tmp1 + tmp2)
      endif
    END DO
  END DO
enddo

DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
  DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
  Wavefield%E(i,j) = Wavefield%P0(i,j)/g
  Wavefield%P(i,j) = zero 
  Wavefield%W(i,j) = zero
    END DO
  END DO

  
DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
  DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
    wgt_ship(i,j) = 0.0
  Wavefield%E(i,j) =  zero!Wavefield%P0(i,j)/g
  Wavefield%P0(i,j) = zero!Wavefield%P0(i,j)/g
  Wavefield%P(i,j) = zero 
  Wavefield%W(i,j) = zero
    END DO
  END DO

  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
  DO j = 1+GhostGridY,10
  Wavefield%E(i,j) =  0.05
    END DO
  END DO
   
  elseif (pressureTermOnOff==7) then
  
write(*,'(A)') 'Ship type: Kriso Container Ship (KCS)'

    !---------------------------------------------------------------------------
    ! Open ship file:
    !---------------------------------------------------------------------------
!    OPEN(UNIT=9,FILE=shipfile)

    !---------------------------------------------------------------------------
    ! Read ship grid parameters:
    !---------------------------------------------------------------------------
!    READ(UNIT=9,FMT=101)        Fnship 
!   READ(UNIT=9,FMT=101)        x0
!    READ(UNIT=9,FMT=101)        y0 

Fnship = 0.26
x0 = 1.5
y0 = 0.0 

    !---------------------------------------------------------------------------
    ! Write the ship grid parameters to screen:
    !---------------------------------------------------------------------------
    WRITE(*,'(A,F8.4)') 'Fnship: ',Fnship 
    WRITE(*,'(A,F8.4)') 'x0:     ',x0
    WRITE(*,'(A,F8.4)') 'y0:     ',y0 

  tol = 1e-6
  !-----------------------------------------------------------------------------
  ! Series 60 ship parameters:
  !-----------------------------------------------------------------------------
  Lship   =   1.0_long
  Uship   =   Fnship*((g*Lship)**(0.5_long)) 
  scale   =   1.0

  !-----------------------------------------------------------------------------
  ! Open KCS_hull_SVA.ow3 file:
  !-----------------------------------------------------------------------------
  OPEN(UNIT=9,FILE='KCS_hull_SVA.ow3')

  !---------------------------------------------------------------------------
  ! Read file:
  !---------------------------------------------------------------------------
  ALLOCATE(zst(81,1001))
  DO i=1,1001
    DO j=1,81
      READ(UNIT=9,FMT=*) zst(j,i)
    enddo
  ENDDO

  !-----------------------------------------------------------------------------
  ! Loop all grid points:
  !-----------------------------------------------------------------------------
        xmin = -0.043478260869565
        xmax =  1.043478260869565
        ymax =  0.086956521739130

  DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
    DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX

      Wavefield%P0(i,j) = 0.0
      wgt_ship(i,j)     = 0.0

      if (xmin.lt.FineGrid%x(i,j)-x0 .and. FineGrid%x(i,j)-x0.lt.xmax .and. &
          abs(FineGrid%y(i,j)-y0) .lt. ymax) then

        !-----------------------------------------------------------------------
        ! Inside a box around ship with size 2*Lship x 1*Lship:
        !-----------------------------------------------------------------------

        i1 = floor((FineGrid%x(i,j) - x0 + 0.043478260869565)/0.001086956521739)
        j1 = floor((FineGrid%y(i,j) - y0)/0.001086956521739)

dx1 =(FineGrid%x(i,j)-x0+0.043478260869565)/0.001086956521739 - real(i1)
dx2 = 1.0 - dx1

dy1 = (FineGrid%y(i,j) - y0)/0.001086956521739 - real(j1)
dy2 = 1.0 - dy1

              w1 = dx2/(dx1+dx2)
              w2 = dx1/(dx1+dx2)

              z1 = (w1*zst(j1+1,i1+1) + w2*zst(j1+1,i1+2))
              z2 = (w1*zst(j1+2,i1+1) + w2*zst(j1+2,i1+2))

              w1 = dy2/(dy1+dy2)
              w2 = dy1/(dy1+dy2)

Wavefield%P0(i,j) = (w1*z1 + w2*z2)*g
if (abs(Wavefield%P0(i,j)) .gt. 1e-10) then
  wgt_ship(i,j)     = 1.0
endif

      endif
    END DO
  END DO

  !-----------------------------------------------------------------------------
  ! Do a light smoothing of the ship geometry. "Just to break the corners".
  ! Maybe use a even lighter filter with only three points
  !-----------------------------------------------------------------------------
  flag_filter = 1
  do iter = 1,4
    DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
      DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX

        if(x0-0.5.lt.FineGrid%x(i,j) .and. FineGrid%x(i,j) .lt. x0+1.5 .and. &
        y0-0.5 .lt. FineGrid%y(i,j) .and. FineGrid%y(i,j) .lt. y0+0.5) then

        !-----------------------------------------------------------------------
        ! Inside a box around ship with size 2*Lship x 1*Lship:
        !-----------------------------------------------------------------------

        !--------------------------------------------------------------------
        ! Filter:
        !--------------------------------------------------------------------
        if (flag_filter .eq. 1) then
        tmp1 = Wavefield%P0(i-1,j)*0.010867541574776 &
        + Wavefield%P0(i,j)*0.978264916850449   &
        + Wavefield%P0(i+1,j)*0.010867541574776 

        if (j==2) then
          tmp2 =  Wavefield%P0(i,j+1)*0.010867541574776 &
          + Wavefield%P0(i,j)*0.978264916850449   &
          + Wavefield%P0(i,j+1)*0.010867541574776 

        else
          tmp2 =  Wavefield%P0(i,j-1)*0.010867541574776 &
          + Wavefield%P0(i,j)*0.978264916850449   &
          + Wavefield%P0(i,j+1)*0.010867541574776 

        endif

        elseif (flag_filter .eq. 2) then

        tmp1 = Wavefield%P0(i-2,j)*0.000591270361240 &
        + Wavefield%P0(i-1,j)*0.125419103234394 &
        + Wavefield%P0(i,j)*0.747979252808734   &
        + Wavefield%P0(i+1,j)*0.125419103234393 &
        + Wavefield%P0(i+2,j)*0.000591270361240

        if (j==2) then
          tmp2 = Wavefield%P0(i,j+2)*0.000591270361240 &
          + Wavefield%P0(i,j+1)*0.125419103234394 &
          + Wavefield%P0(i,j)*0.747979252808734   &
          + Wavefield%P0(i,j+1)*0.125419103234393 &
          + Wavefield%P0(i,j+2)*0.000591270361240

          elseif (j==3) then
          tmp2 = Wavefield%P0(i,j)*0.000591270361240 &
          + Wavefield%P0(i,j-1)*0.125419103234394 &
          + Wavefield%P0(i,j)*0.747979252808734   &
          + Wavefield%P0(i,j+1)*0.125419103234393 &
          + Wavefield%P0(i,j+2)*0.000591270361240

        else
          tmp2 = Wavefield%P0(i,j-2)*0.000591270361240 &
          + Wavefield%P0(i,j-1)*0.125419103234394 &
          + Wavefield%P0(i,j)*0.747979252808734   &
          + Wavefield%P0(i,j+1)*0.125419103234393 &
          + Wavefield%P0(i,j+2)*0.000591270361240

        endif
      endif
        Wavefield%P0(i,j) = 0.5*(tmp1 + tmp2)
      endif
    END DO
  END DO
enddo

  DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
    DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
  Wavefield%E(i,j) = Wavefield%P0(i,j)/g
  Wavefield%P(i,j) = zero 
  Wavefield%W(i,j) = zero
    END DO
  END DO

  endif
  
  ALLOCATE(dzdxi1(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
  ALLOCATE(dzdxi2(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))
  ALLOCATE(nx_ship(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY))

  dxi1 = 1.0_long/(FineGrid%Nx-1)
  dxi2 = 1.0_long/(FineGrid%Ny-1)

  CALL DetermineGenericStencils(FineGrid%CurvilinearStuff%DiffStencils,alpha)
  CALL DiffXuniform2D(WaveField%E(1+GhostGridX:FineGrid%Nx+GhostGridX,&
                                  1+GhostGridY:FineGrid%Ny+GhostGridY),&
                      dzdxi1(   1+GhostGridX:FineGrid%Nx+GhostGridX,&
                                1+GhostGridY:FineGrid%Ny+GhostGridY),1,&
  FineGrid%CurvilinearStuff%DiffStencils%StencilG,FineGrid%Nx,&
  FineGrid%Ny,alpha)
  dzdxi1 = 1.0_long/dxi1*dzdxi1

  CALL DetermineGenericStencils(FineGrid%CurvilinearStuff%DiffStencils,beta)
  CALL DiffYuniform2D(WaveField%E(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),&
                           dzdxi2(1+GhostGridX:FineGrid%Nx+GhostGridX,1+GhostGridY:FineGrid%Ny+GhostGridY),1,&
  FineGrid%CurvilinearStuff%DiffStencils%StencilG,FineGrid%Nx,&
  FineGrid%Ny,beta)
  dzdxi2 = 1.0_long/dxi2*dzdxi2

  !-----------------------------------------------------------------------------
  ! Normal vector to ship surface:
  !-----------------------------------------------------------------------------
    OPEN(UNIT=9,FILE='normalx')

  DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
    DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
      

      dxdxi1 = 1.0_long/dxi1*FineGrid%CurvilinearStuff%xe(i,j)
      dydxi1 = 1.0_long/dxi1*FineGrid%CurvilinearStuff%ye(i,j)
      dxdxi2 = 1.0_long/dxi2*FineGrid%CurvilinearStuff%xn(i,j)
      dydxi2 = 1.0_long/dxi2*FineGrid%CurvilinearStuff%yn(i,j)
      
      n1 = -dydxi2*dzdxi1(i,j)
      n2 = -dxdxi1*dzdxi2(i,j)
      n3 =  dxdxi1*dydxi2

      lng           = sqrt(n1**2.0 + n2**2.0 + n3**2.0)
      nx_ship(i,j)  = -n1/lng
!WRITE(unit=9,fmt='(F24.16,F24.16,F24.16)') FineGrid%X(i,j),FineGrid%y(i,j),nx_ship(i,j)

!      if (nx_ship(i,j)>0.00001) then
!      print *,nx_ship(i,j)
!    endif
    END DO
!WRITE(unit=9,fmt='(A)')'' 
  END DO
!  DEALLOCATE(dzdxi1)
!  DEALLOCATE(dzdxi2)


  END SUBROUTINE funInitialFreeSurfaceElevation
