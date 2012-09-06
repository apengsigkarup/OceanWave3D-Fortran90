!-------------------------------------------------------------------------------
! Finite difference WENO with Lax-Friedric flux:
!-------------------------------------------------------------------------------
SUBROUTINE FD2DWENOLF(u,dxi,A,Fx)
  
  !-----------------------------------------------------------------------------
  ! Global parameters:
  !-----------------------------------------------------------------------------
  USE GlobalVariables, ONLY: long,FineGrid,ghostgridx,ghostgridy,ghostgridz,&
  aw,bw,cw,uL,uR

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! External subroutines:
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! Input parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long)       :: dxi,A
  REAL(KIND=long),DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: u

  !-----------------------------------------------------------------------------
  ! Output parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long),DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: Fx

  !-----------------------------------------------------------------------------
  ! Local parameters:
  !-----------------------------------------------------------------------------
  INTEGER i,j
  REAL(KIND=long) :: epsilon = 1e-14_long
  REAL(KIND=long) :: aL,bL,cL,aC,bC,cC,aR,bR,cR,OIL,OIC,OIR,wL,wC,wR 
  REAL(KIND=long) :: sumw 
  REAL(KIND=long) :: C 
  REAL(KIND=long) :: fm,fp

  !-----------------------------------------------------------------------------
  ! WENO flux in x-direction:
  !-----------------------------------------------------------------------------
  DO i = 2+GhostGridX,FineGrid%Nx-1+GhostGridX
    DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY

      aL = 0.5/dxi**2*(u(i-2,j) - 2*u(i-1,j) +   u(i,j));
      bL = 0.5/dxi*   (u(i-2,j) - 4*u(i-1,j) + 3*u(i,j));
      cL =                                       u(i,j);

      aC = 0.5/dxi**2*( u(i-1,j) - 2*u(i,j) + u(i+1,j));
      bC = 0.5/dxi*   (-u(i-1,j)            + u(i+1,j));
      cC =                                    u(i  ,j);

      aR =  0.5/dxi**2*(  u(i,j) - 2*u(i+1,j) + u(i+2,j));
      bR = -0.5/dxi*   (3*u(i,j) - 4*u(i+1,j) + u(i+2,j));
      cR =                                      u(i  ,j);

      OIL = 13/3*aL**2*dxi**4 + bL**2*dxi**2;
      OIC = 13/3*aC**2*dxi**4 + bC**2*dxi**2;
      OIR = 13/3*aR**2*dxi**4 + bR**2*dxi**2;

      wL =   1/(OIL + epsilon)**8;
      wC = 1e6/(OIC + epsilon)**8;
      wR =   1/(OIR + epsilon)**8;

      sumw = wL + wC + wR;
      wL = wL/sumw;
      wC = wC/sumw;
      wR = wR/sumw;

      aw(i,j) =   wL*aL + wC*aC + wR*aR;
      bw(i,j) =   wL*bL + wC*bC + wR*bR;
      cw(i,j) =   wL*cL + wC*cC + wR*cR;

    END DO
  END DO

  !-----------------------------------------------------------------------------
  ! Boundary conditions:
  !-----------------------------------------------------------------------------
  DO i = GhostGridX+1,GhostGridX+2
    DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
      aw(i,j) =  0.0_long 
      bw(i,j) =  0.0_long 
      cw(i,j) =  u(i,j) 
    END DO
  END DO

  DO i = FineGrid%Nx-1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
      aw(i,j) =  0.0_long 
      bw(i,j) =  0.0_long 
      cw(i,j) =  u(i,j) 
    END DO
  END DO

  !-----------------------------------------------------------------------------
  ! WENO reconstruction:
  !-----------------------------------------------------------------------------
  DO i = 1,FineGrid%Nx+2*GhostGridX
    DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY

      uL(i,j)= aw(i,j)*((-0.5*dxi)**2-1/12*dxi**2) + bw(i,j)*(-0.5*dxi)+cw(i,j)
      uR(i,j)= aw(i,j)*(( 0.5*dxi)**2-1/12*dxi**2) + bw(i,j)*( 0.5*dxi)+cw(i,j)

    END DO
  END DO

  !-----------------------------------------------------------------------------
  ! Max eigenvalue of flux jacobian:
  !-----------------------------------------------------------------------------
  C = ABS(A)

  !-----------------------------------------------------------------------------
  ! Lax-Friedric flux and flux difference:
  !-----------------------------------------------------------------------------
  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY

      fm        = 0.5*A*(uR(i-1,j)+uL(i,j)  ) - 0.5*C*(uL(i  ,j) - uR(i-1,j))
      fp        = 0.5*A*(uR(i,j)  +uL(i+1,j)) - 0.5*C*(uL(i+1,j) - uR(i  ,j))
      Fx(i,j)   = 1.0/dxi*(fp - fm);

    END DO
  END DO

END SUBROUTINE FD2DWENOLF 


!-------------------------------------------------------------------------------
! Finite difference WENO with Lax-Friedric flux:
!-------------------------------------------------------------------------------
SUBROUTINE FD2DWENOReconXi2(u,dxi,uL,uR)
  
  !-----------------------------------------------------------------------------
  ! Global parameters:
  !-----------------------------------------------------------------------------
  USE GlobalVariables, ONLY: long,FineGrid,ghostgridx,ghostgridy,ghostgridz,&
  aw,bw,cw

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! External subroutines:
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! Input parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long)       :: dxi
  REAL(KIND=long),DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: u,uL,uR

  !-----------------------------------------------------------------------------
  ! Local parameters:
  !-----------------------------------------------------------------------------
  INTEGER i,j
  REAL(KIND=long) :: epsilon = 1e-14_long
  REAL(KIND=long) :: aL,bL,cL,aC,bC,cC,aR,bR,cR,OIL,OIC,OIR,wL,wC,wR 
  REAL(KIND=long) :: sumw 

  !-----------------------------------------------------------------------------
  ! WENO flux in y-direction:
  !-----------------------------------------------------------------------------
  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = 2+GhostGridY,FineGrid%Ny-1+GhostGridY

      aL = 0.5/dxi**2*(u(i,j-2) - 2*u(i,j-1) +   u(i,j));
      bL = 0.5/dxi*   (u(i,j-2) - 4*u(i,j-1) + 3*u(i,j));
      cL =                                       u(i,j);

      aC = 0.5/dxi**2*( u(i,j-1) - 2*u(i,j) + u(i,j+1));
      bC = 0.5/dxi*   (-u(i,j-1)            + u(i,j+1));
      cC =                                    u(i  ,j);

      aR =  0.5/dxi**2*(  u(i,j) - 2*u(i,j+1) + u(i,j+2));
      bR = -0.5/dxi*   (3*u(i,j) - 4*u(i,j+1) + u(i,j+2));
      cR =                                      u(i  ,j);

      OIL = 13/3*aL**2*dxi**4 + bL**2*dxi**2;
      OIC = 13/3*aC**2*dxi**4 + bC**2*dxi**2;
      OIR = 13/3*aR**2*dxi**4 + bR**2*dxi**2;

      wL =   1/(OIL + epsilon)**8;
      wC = 1e6/(OIC + epsilon)**8;
      wR =   1/(OIR + epsilon)**8;

      sumw = wL + wC + wR;
      wL = wL/sumw;
      wC = wC/sumw;
      wR = wR/sumw;

      bw(i,j) =   wL*bL + wC*bC + wR*bR;

    END DO
  END DO

  !-----------------------------------------------------------------------------
  ! Boundary conditions:
  !-----------------------------------------------------------------------------
  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = GhostGridY+1,GhostGridY+2
      aw(i,j) =  0.0_long 
      bw(i,j) =  0.0_long 
      cw(i,j) =  u(i,j) 
    END DO
  END DO

  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = FineGrid%Ny-1+GhostGridY,FineGrid%Ny+GhostGridY
      aw(i,j) =  0.0_long 
      bw(i,j) =  0.0_long 
      cw(i,j) =  u(i,j) 
    END DO
  END DO
  
  
  !-----------------------------------------------------------------------------
  ! WENO reconstruction:
  !-----------------------------------------------------------------------------
  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = FineGrid%Ny-1+GhostGridY,FineGrid%Ny+GhostGridY
 
      uL(i,j)= aw(i,j)*((-0.5*dxi)**2-1/12*dxi**2) + bw(i,j)*(-0.5*dxi)+cw(i,j)
      uR(i,j)= aw(i,j)*(( 0.5*dxi)**2-1/12*dxi**2) + bw(i,j)*( 0.5*dxi)+cw(i,j)

    END DO
  END DO

END SUBROUTINE FD2DWENOReconXi2 


!-------------------------------------------------------------------------------
! Finite difference WENO with Lax-Friedric flux:
!-------------------------------------------------------------------------------
SUBROUTINE FD2DWENOLFy(u,v,dxi,Fy)
  
  !-----------------------------------------------------------------------------
  ! Global parameters:
  !-----------------------------------------------------------------------------
  USE GlobalVariables, ONLY: long,FineGrid,ghostgridx,ghostgridy,ghostgridz,&
  aw,bw,cw,uL,uR,vL,vR

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! External subroutines:
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! Input parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long)       :: dxi
  REAL(KIND=long),DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: u
  REAL(KIND=long),DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: v

  !-----------------------------------------------------------------------------
  ! Output parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long),DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: Fy

  !-----------------------------------------------------------------------------
  ! Local parameters:
  !-----------------------------------------------------------------------------
  INTEGER i,j
  REAL(KIND=long) :: C 
  REAL(KIND=long) :: fm,fp

  !-----------------------------------------------------------------------------
  ! Weno reconstruction:
  !-----------------------------------------------------------------------------
  CALL FD2DWENOReconXi2(u,dxi,uL,uR)
  CALL FD2DWENOReconXi2(v,dxi,vL,vR)

  !-----------------------------------------------------------------------------
  ! Max eigenvalue of flux jacobian:
  !-----------------------------------------------------------------------------
  C = 0.0
  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
      C = MAX(ABS(v(i,j)),C)
    END DO
  END DO


  !-----------------------------------------------------------------------------
  ! Lax-Friedric flux and flux difference:
  !-----------------------------------------------------------------------------
  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY

      fm = 0.5*(vR(i,j-1)*uR(i,j-1) + vL(i,j)*uL(i,j)) &
         - 0.5*C*(uL(i,j)   - uR(i,j-1))

      fp = 0.5*(vR(i,j)*uR(i,j)     + vL(i,j+1)*uL(i,j+1)) &
         - 0.5*C*(uL(i,j+1) - uR(i,j))

      Fy(i,j)   = 1.0/dxi*(fp - fm);

    END DO
  END DO

END SUBROUTINE FD2DWENOLFy 

!-------------------------------------------------------------------------------
! Finite difference WENO Derivative:
!-------------------------------------------------------------------------------
SUBROUTINE FD2DWENODxi2(u,dxi,Dxi2)
  
  !-----------------------------------------------------------------------------
  ! Global parameters:
  !-----------------------------------------------------------------------------
  USE GlobalVariables, ONLY: long,FineGrid,ghostgridx,ghostgridy,ghostgridz,&
  aw,bw,cw,uL,uR

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! External subroutines:
  !-----------------------------------------------------------------------------
  
  !-----------------------------------------------------------------------------
  ! Input parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long)       :: dxi
  REAL(KIND=long),DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: u

  !-----------------------------------------------------------------------------
  ! Output parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long),DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: Dxi2

  !-----------------------------------------------------------------------------
  ! Local parameters:
  !-----------------------------------------------------------------------------
  INTEGER i,j
  REAL(KIND=long) :: epsilon = 1e-14_long
  REAL(KIND=long) :: aL,bL,cL,aC,bC,cC,aR,bR,cR,OIL,OIC,OIR,wL,wC,wR 
  REAL(KIND=long) :: sumw 

  !-----------------------------------------------------------------------------
  ! WENO flux in y-direction:
  !-----------------------------------------------------------------------------
  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = 2+GhostGridY,FineGrid%Ny-1+GhostGridY

      aL = 0.5/dxi**2*(u(i,j-2) - 2*u(i,j-1) +   u(i,j));
      bL = 0.5/dxi*   (u(i,j-2) - 4*u(i,j-1) + 3*u(i,j));
      cL =                                       u(i,j);

      aC = 0.5/dxi**2*( u(i,j-1) - 2*u(i,j) + u(i,j+1));
      bC = 0.5/dxi*   (-u(i,j-1)            + u(i,j+1));
      cC =                                    u(i  ,j);

      aR =  0.5/dxi**2*(  u(i,j) - 2*u(i,j+1) + u(i,j+2));
      bR = -0.5/dxi*   (3*u(i,j) - 4*u(i,j+1) + u(i,j+2));
      cR =                                      u(i  ,j);

      OIL = 13/3*aL**2*dxi**4 + bL**2*dxi**2;
      OIC = 13/3*aC**2*dxi**4 + bC**2*dxi**2;
      OIR = 13/3*aR**2*dxi**4 + bR**2*dxi**2;

      wL =   1/(OIL + epsilon)**8;
      wC = 1e6/(OIC + epsilon)**8;
      wR =   1/(OIR + epsilon)**8;

      sumw = wL + wC + wR;
      wL = wL/sumw;
      wC = wC/sumw;
      wR = wR/sumw;

      bw(i,j) =   wL*bL + wC*bC + wR*bR;

    END DO
  END DO

  !-----------------------------------------------------------------------------
  ! Boundary conditions:
  !-----------------------------------------------------------------------------
  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = GhostGridY+1,GhostGridY+2
      bw(i,j) =  0.0_long 
    END DO
  END DO

  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = FineGrid%Ny-1+GhostGridY,FineGrid%Ny+GhostGridY
      bw(i,j) =  0.0_long 
    END DO
  END DO

  !-----------------------------------------------------------------------------
  ! WENO derivative:
  !-----------------------------------------------------------------------------
  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = 1,FineGrid%Ny+2*GhostGridY
      Dxi2(i,j)=bw(i,j)
    END DO
  END DO

END SUBROUTINE FD2DWENODxi2

!-------------------------------------------------------------------------------
! ODE right hand side for Linear ship modelling:
!-------------------------------------------------------------------------------
SUBROUTINE RHSLinearShip2D(time,WaveField,RHSE,RHSP)
  
  !-----------------------------------------------------------------------------
  ! Global parameters:
  !-----------------------------------------------------------------------------
  USE GlobalVariables, ONLY: dt,long,RHS,FineGrid,Uship,alpha,beta,relaxonoff,&
  PHI,ghostgridx,ghostgridy,ghostgridz,alpha,gamma,g,WaveField_FS,Lship,dx,dy,&
  FEx,FPx,FEy,FPy,wgt_ship

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  ! External subroutines:
  !-----------------------------------------------------------------------------
  EXTERNAL BuildLinearSystem, BuildLinearSystemTransformedCurvilinear

  !-----------------------------------------------------------------------------
  ! Input parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long)       :: time, dxi1,dxi2,jac,jac2,Evisc,Pvisc
  TYPE(WaveField_FS)    :: WaveField

  !-----------------------------------------------------------------------------
  ! Output parameters:
  !-----------------------------------------------------------------------------
  REAL(KIND=long),DIMENSION(FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY)&
  :: RHSE,RHSP

  !-----------------------------------------------------------------------------
  ! Local parameters:
  !-----------------------------------------------------------------------------
  INTEGER i,j
  REAL(KIND=long) :: relax!,invdx2,invdy2

!  INTEGER :: rank,order, Diff
!  REAL(KIND=long), DIMENSION(2*alpha+1) :: Stencil
!  INTEGER, DIMENSION(2*alpha+1) :: idx

  relax = min(time/2.5,1.0_long)
relax = 0.0
!  invdx2=1.0/(dx*dx)
!  invdy2=1.0/(dy*dy)

  !-----------------------------------------------------------------------------
  ! Relax the free surface elevation and the downstream boundary:
  !-----------------------------------------------------------------------------
!  IF (relaxONOFF==1) THEN
!    CALL RelaxationModule_new(WaveField%E,WaveField%P,time)
!  ENDIF

  !-----------------------------------------------------------------------------
  ! Right hand side of Laplace equation (the potential on the free surface):
  !-----------------------------------------------------------------------------
  RHS(FineGrid%Nz+GhostGridZ,1+GhostGridX:FineGrid%Nx+GhostGridX,&
  1+GhostGridY:FineGrid%Ny+GhostGridY) = &
  Wavefield%P(1+GhostGridX:FineGrid%Nx+GhostGridX,&
  1+GhostGridY:FineGrid%Ny+GhostGridY)

  !-----------------------------------------------------------------------------
  ! Solve Laplace equation for the velocity potential:
  !-----------------------------------------------------------------------------
!  CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX) &
!  *(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
!  BuildLinearSystem,PHI,FineGrid)
  CALL iterative_solution(RHS,(FineGrid%Nx+2*GhostGridX) &
  *(FineGrid%Ny+2*GhostGridY)*(FineGrid%Nz+GhostGridZ),&
  BuildLinearSystemTransformedCurvilinear,PHI,FineGrid)

  !-----------------------------------------------------------------------------
  ! Vertical velocity at the free surface:
  !-----------------------------------------------------------------------------
  CALL VerticalFreeSurfaceVelocity(Wavefield%W,FineGrid%Nx+2*GhostGridX,&
  FineGrid%Ny+2*GhostGridY,FineGrid%Nz+GhostGridZ,PHI,FineGrid%DiffStencils,&
  FineGrid%dsigmanew(:,:,:,5), gamma)

  !-----------------------------------------------------------------------------
  ! Boundary conditions:
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! Derivatives at the free surface:
  !-----------------------------------------------------------------------------
  CALL DifferentiationsFreeSurfacePlane(Wavefield,GhostGridX,GhostGridY,&
  FineGrid,alpha,beta)

!  i = FineGrid%Nx+2*GhostGridX
!  DO j = 1,FineGrid%Ny+2*GhostGridY
!    WaveField%Pxx(i,j)      = 0.0
!    WaveField%Pxx(i-1,j)    = 0.0
!    WaveField%Pxx(i-2,j)    = 0.0
!    WaveField%Pxx(i-3,j)    = 0.0
!    WaveField%Pxx(i-4,j)    = 0.0
!    WaveField%Pxx(i-5,j)    = 0.0
!    WaveField%Pxx(i-6,j)    = 0.0
!  END DO


  !-----------------------------------------------------------------------------
  ! Onesided derivative of velocity potential at downstream boundary:
  !-----------------------------------------------------------------------------
!  rank = 2*alpha+1
!  DO i=-alpha,alpha,1
!    idx(i+alpha+1) = i
!  END DO
!  order=1
!  DO i = GhostGridX+1,alpha+GhostGridX
!    Diff = alpha+GhostGridX+1-i
!    CALL TaylorFDStencils1DArbitrary(alpha-Diff,alpha+Diff,order,Stencil,&
!    FineGrid%x((GhostGridX+1):(GhostGridX+rank),1))
!    DO j = 1,FineGrid%Ny+2*GhostGridY
!      WaveField%Px(i,j) = DOT_PRODUCT(Stencil,&
!      WaveField%P((GhostGridX+1):(GhostGridX+rank),j))
!      WaveField%Ex(i,j) = DOT_PRODUCT(Stencil,&
!      WaveField%E((GhostGridX+1):(GhostGridX+rank),j))
!     END DO
!  END DO
!
!  DO i = ceiling(GhostGridX+FineGrid%Nx*0.2),floor(FineGrid%Nx*0.9+GhostGridX)
!    DO j = 1,floor(FineGrid%Ny*0.8+GhostGridY)
!      WaveField%Exx(i,j) = 0.0_long
!      WaveField%Eyy(i,j) = 0.0_long
!      WaveField%Pxx(i,j) = 0.0_long
!      WaveField%Pyy(i,j) = 0.0_long
!     END DO
!  END DO


  !-----------------------------------------------------------------------------
  ! WENO flux:
  !-----------------------------------------------------------------------------
  dxi1 = 1.0_long/(FineGrid%Nx-1)
  CALL FD2DWENOLF(WaveField%E,dxi1,-relax*Uship,FEx)
  CALL FD2DWENOLF(WaveField%P,dxi1,-relax*Uship,FPx) 

  !-----------------------------------------------------------------------------
  ! WENO transverse velocity:
  !-----------------------------------------------------------------------------
  dxi2 = 1.0_long/(FineGrid%Ny-1)
  CALL FD2DWENODxi2(WaveField%P,dxi2,FPy)
!  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
!    DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
!    jac = FineGrid%CurvilinearStuff%yn(i,j)/dxi2
!     FPy(i,j) = 1.0_long/(jac)*FPy(i,j)
!    END DO
!  END DO
!  CALL FD2DWENOLFy(FPy,WaveField%E,dxi2,FEy) 

  !-----------------------------------------------------------------------------
  ! Upstream boundary condition:
  !-----------------------------------------------------------------------------
  i = FineGrid%Nx+GhostGridX
  DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY
    FEx(i,j) = 0.0
    FPx(i,j) = 0.0
  end do

  !-----------------------------------------------------------------------------
  ! Summation of right hand sides:
  !-----------------------------------------------------------------------------
  DO i = 1+GhostGridX,FineGrid%Nx+GhostGridX
    DO j = 1+GhostGridY,FineGrid%Ny+GhostGridY


      if( 1.0 .lt. FineGrid%x(i,j) .and. FineGrid%x(i,j) .lt. 3.0 .and. &
         -0.5 .lt. FineGrid%y(i,j) .and. FineGrid%y(i,j) .lt. 0.5) then
         Evisc = 0.2*min(dx,dy)**2.0/dt
         Pvisc = 0.2*min(dx,dy)**2.0/dt
       else 
         Evisc = 0.2*min(dx,dy)**2.0/dt
         Pvisc = 0.2*min(dx,dy)**2.0/dt
       endif




    jac = FineGrid%CurvilinearStuff%xe(i,j)/dxi1
!    jac2 = FineGrid%CurvilinearStuff%yn(i,j)/dxi2

!      RHSE(i,j) = - 1.0_long/jac*FEx(i,j) &
!                  - 1.0_long/jac2*FEy(i,j) + Wavefield%W(i,j)

      RHSE(i,j) = - 1.0_long/jac*FEx(i,j) + Wavefield%W(i,j) &
    + Evisc*(WaveField%Exx(i,j) + WaveField%Eyy(i,j)) 
          

      RHSP(i,j) = - 1.0_long/jac*FPx(i,j) - g*Wavefield%E(i,j) &
                  + Wavefield%P0(i,j) + wgt_ship(i,j)/jac*FPx(i,j) &
    + Pvisc*(WaveField%Pxx(i,j) + WaveField%Pyy(i,j)) 

!      if (j==10) then
!        write(*,'(I8,F8.4,F8.4,F8.4,F8.4)')i,FineGrid%x(i,j),WaveField%E(i,j),RHSE(i,j)
!      endif
    END DO
  END DO

  !-----------------------------------------------------------------------------
  ! Set right hand sides to zero in the upstream corners:
  !-----------------------------------------------------------------------------
  RHSE(1,:)                         = 0.0_long
  RHSP(1,:)                         = 0.0_long
  RHSE(FineGrid%Nx+2*GhostGridX,:)  = 0.0_long
  RHSP(FineGrid%Nx+2*GhostGridX,:)  = 0.0_long
  RHSE(:,1)                         = 0.0_long 
  RHSP(:,1)                         = 0.0_long
  RHSE(:,FineGrid%Ny+2*GhostGridY)  = 0.0_long
  RHSP(:,FineGrid%Ny+2*GhostGridY)  = 0.0_long

  END SUBROUTINE RHSLinearShip2D
