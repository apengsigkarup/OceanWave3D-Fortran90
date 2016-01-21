SUBROUTINE LocalSmoothing3D(nx,ny,nz,itime,fop)
  ! 
  ! This subroutine attempts to apply automatic localized smoothing to avoid any 
  ! breakdowns due to discontinuities in the free-surface.  We check the Lagrangian 
  ! downward vertical acceleration everywhere and when it exceeds a given factor times g 
  ! we smooth a 10-point region centered at the point with the classical 3-point 
  ! filter [.25 .5 .25].  The feature is turned off by switching the threshold 
  ! factor to something larger than 100.  
  !
  ! **  This is the initial 3D version with smoothing along x-lines only **
  !
  USE GlobalVariables
  IMPLICIT NONE
  INTEGER :: nx, ny, nz, itime, fop

  ! Local variables and automatic workspace
  INTEGER, SAVE :: i_first=0
  INTEGER ::  i, j, i_flag, N_smoothed, i_smoothed(nx*ny), j_smoothed(nx*ny), k
  REAL(KIND=long)  accel_tol, dt_inv, Max_accel
  REAL(KIND=long), POINTER :: et(:,:), ph(:,:), x(:,:), y(:,:), WHist(:,:,:), etax(:,:),      &
       etay(:,:), hx(:,:), hy(:,:), z(:), h(:,:), EtatHist(:,:,:)
  REAL(KIND=long) :: etatem(nx,ny), Phitem(nx,ny), etatt(nx,ny), alpha_s(nx,ny), Wt(nx,ny),   &
       dWdt(nx,ny)
  REAL(KIND=long) :: U(Nz,Nx,Ny), V(Nz,Nx,Ny), W(Nz,Nx,Ny), Wz(Nz,Nx,Ny), Wx(Nz,Nx,Ny),       &
       Wy(Nz,Nx,Ny), d(Nx,Ny)
  !
  ! Initial set up of the output file
  !
  IF(i_first==0)THEN
     i_first=1
     OPEN(fop,file='local.smoothing',status='unknown')
     WRITE(fop,10)accel_tol_fact
10   FORMAT('# Points where Lagrangian -dw/dt exceeded',e10.2, &
          ' times g and local smoothing was applied:  t, x, y, dw/dt')
  END IF
  !
  ! Do nothing for a factor larger than 100
  !
  IF(accel_tol_fact .le. 100. .and. itime .ge. 3)THEN
     !
     ! Search for high accelerations and smooth an area around points where they occur.  
     !
     !
     ! Initialize
     !
     x => FineGrid%x; y => FineGrid%y; et => WaveField%E; ph => WaveField%P; WHist => WaveField%WHist;
     etax => WaveField%Ex; etay => WaveField%Ey; hx => FineGrid%hx; hy => FineGrid%hy; 
     z => FineGrid%z; h => FineGrid%h; EtatHist => WaveField%EtatHist;
     etatem=et
     Phitem=ph
     alpha_s=zero
     Do j=1,Ny
        Do i=1,Nx
           d(i,j)=h(i,j)+etatem(i,j);
        end Do
     End Do
     !
     ! Take the computational-space derivatives of phi
     !
     CALL DiffXEven(phi,U,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
          FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,alpha)
     CALL DiffXEven(phi,V,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
          FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,beta)
     CALL DiffZArbitrary(phi,W,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
          FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
     !
     ! Form the Eulerian fluid velocities
     !
     Do j=1,Ny
        Do i=1,Nx
           Do k=1,Nz
              W(k,i,j) = W(k,i,j)/d(i,j) 
              U(k,i,j) = U(k,i,j) + ((1-z(k))*hx(i,j)-z(k)*etax(i,j))*W(k,i,j)   
              V(k,i,j) = V(k,i,j) + ((1-z(k))*hy(i,j)-z(k)*etay(i,j))*W(k,i,j)
           END Do
        END Do
     END Do
     !
     ! Take computational-space derivatives of w and convert to physical-space derivatives.
     !     
     CALL DiffZArbitrary(W,Wz,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
          FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
     CALL DiffXEven(W,Wx,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
          FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,alpha)
     CALL DiffYEven(W,Wy,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
          FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,beta)
     Do j=1,Ny
        Do i=1,Nx
           Do k=1,Nz
              Wz(k,i,j) = Wz(k,i,j)/d(i,j) 
              Wx(k,i,j) = Wx(k,i,j) + ((1-z(k))*hx(i,j)-z(k)*etax(i,j))*Wz(k,i,j)   
              Wy(k,i,j) = Wy(k,i,j) + ((1-z(k))*hy(i,j)-z(k)*etay(i,j))*Wz(k,i,j)   
           END Do
        END Do
     END DO
     !
     ! Take a 3-point one-sided difference to get the computational-space dw/dt at the 
     ! current time step and convert it to a physical-space dw/dt.   This is only done 
     ! on the free-surface, i.e. at vertical grid point k=nz.  
     !
     dt_inv=one/dt
     Do j=1,Ny
        DO i=1,nx
           Wt(i,j)=half*dt_inv*(three*WHist(i,j,1)-four*WHist(i,j,2)+WHist(i,j,3)) &
                -Wz(Nz,i,j)*EtatHist(i,j,1)
        END Do
     END DO
     !
     ! Form the vertical particle acceleration on the free-surface, (the Lagrangian dw/dt). 
     ! 
     Do j=1,Ny
        DO i=1,nx
           Wt(i,j)=Wt(i,j)+U(Nz,i,j)*Wx(Nz,i,j) + V(Nz,i,j)*Wy(Nz,i,j) + W(Nz,i,j)*Wz(Nz,i,j)
        END DO
     END Do
     !
     ! First build a coefficient vector which is .25 at those points which surround a point 
     ! of high accelerations.  
     !
     accel_tol=accel_tol_fact*g
     i_flag=0
     N_smoothed=0; Max_accel=zero;
     Do j=1,ny
        DO i=6,nx-5
           !
           ! If the downward vertical acceleration exceeds the tolerance, we smooth.
           !
           IF(-Wt(i,j) >= accel_tol)THEN
              i_flag=1
              N_smoothed=N_smoothed+1
              i_smoothed(N_smoothed)=i
              j_smoothed(N_smoothed)=j
              alpha_s(i-5:i+5,j)=0.25_long
           END IF
        END DO
     END Do
     IF(i_flag==1)THEN
        !
        ! Smooth the elevation and surface velocity.  
        !
        Do j=1,ny
           DO i=2,nx-1
              et(i,j)=alpha_s(i,j)*etatem(i-1,j)+(one-two*alpha_s(i,j))*etatem(i,j) &
                   +alpha_s(i,j)*etatem(i+1,j)
              ph(i,j)=alpha_s(i,j)*Phitem(i-1,j)+(one-two*alpha_s(i,j))*Phitem(i,j) &
                   +alpha_s(i,j)*Phitem(i+1,j)
           END DO
        END Do
        !
        ! Write the positions where smoothing was applied and the accelerations found.
        !
        WRITE(fop,11)
        DO i=1,N_smoothed
           WRITE(fop,12)(itime-1)*dt,x(i_smoothed(i),j_smoothed(i)),                &
                y(i_smoothed(i),j_smoothed(i)),Wt(i_smoothed(i),j_smoothed(i))
        END DO
11      FORMAT()
12      FORMAT(4e16.4)
     END IF

  END IF
  RETURN
END SUBROUTINE LocalSmoothing3D
