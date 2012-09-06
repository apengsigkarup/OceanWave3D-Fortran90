!-------------------------------------------------------------------------------
! Finite difference WENO with Lax-Friedric flux:
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
  ! WENO flux in x-direction:
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
