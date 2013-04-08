!
! ****************************************************************
!
SUBROUTINE StoreKinematicData(Nx,Ny,Nz,io,it)
!
! ****************************************************************
!
!>
!! Write 3D kinematics data to an unformatted binary file
!!
!!
!! By Allan P. Engsig-Karup.
!<
USE GlobalVariables
IMPLICIT NONE
! Input parameters
INTEGER :: Nx, Ny, Nz, io, it
! Local variables
INTEGER ::  i, j, k, i0, i1, is, j0, j1, js
REAL(KIND=long), DIMENSION(:,:), POINTER :: x, y, h, hx, hy, eta, etax, etay
REAL(KIND=long), DIMENSION(:), POINTER   :: z
! Automatic work space
REAL(KIND=long) :: U(Nz,Nx,Ny), V(Nz,Nx,Ny), W(Nz,Nx,Ny), d(Nx,Ny)
!
! Assign the local pointers
!
x => FineGrid%x; y => FineGrid%y; z => FineGrid%z; h => FineGrid%h; hx => FineGrid%hx
hy => FineGrid%hy; eta => WaveField%E; etax => WaveField%Ex; etay => WaveField%Ey
!
! Shift the horizontal grid point position to account for the ghost points.  
!
i0=Output(io)%xbeg+GhostGridX; i1=Output(io)%xend+GhostGridX; is=Output(io)%xstride; 
j0=Output(io)%ybeg+GhostGridY; j1=Output(io)%yend+GhostGridY; js=Output(io)%ystride; 
IF(it==0)THEN
   !
   ! Save the grid data on the first call
   !
   write (fileop(io+1)) Output(io)%xbeg,Output(io)%xend,Output(io)%xstride, &
        Output(io)%ybeg, Output(io)%yend, Output(io)%ystride,               &
        Output(io)%tbeg,Output(io)%tend,Output(io)%tstride, dt,             &
        FineGrid%Nz+GhostGridZ
   !
   WRITE (fileop(io+1)) ( ( x(i,j), y(i,j), h(i,j), hx(i,j), hy(i,j), &
        i=i0,i1,is ), j=j0,j1,js ) 
   WRITE (fileop(io+1)) (z(i),i=1,Nz)
   IF(curvilinearOnOff/=0)THEN
      Print *, 'StoreKinematicData:  Saving horizontal fluid velocities is not yet implemented for curvilinear grids.'
   END IF
ELSE
   !
   IF(curvilinearOnOff == 0)THEN
      !
      ! Dump this solution slice to the output file
      !
      ! First the free surface elevation and gradient at all points in this slice
      ! 
      WRITE (fileop(io+1)) ( ( eta(i,j), i=i0,i1,is ), j=j0,j1,js) 
      IF(Nx > 1) THEN
         WRITE (fileop(io+1)) ( ( etax(i,j), i=i0,i1,is ), j=j0,j1,js) 
      ELSE
         WRITE (fileop(io+1)) ( ( zero, i=i0,i1,is ), j=j0,j1,js) 
      END IF
      IF(Ny > 1) THEN
         WRITE (fileop(io+1)) ( ( etay(i,j), i=i0,i1,is ), j=j0,j1,js) 
      ELSE
         WRITE (fileop(io+1)) ( ( zero, i=i0,i1,is ), j=j0,j1,js) 
      END IF
      !
      ! The fluid thickness d=h+eta
      !
      Do j=1,Ny
         Do i=1,Nx
            d(i,j)=h(i,j)+eta(i,j);
         end Do
      End Do
      !
      ! Then the velocity potential at all points in this horizontal slice and at all sigma 
      ! (vertical) locations.  
      !
      WRITE (fileop(io+1)) ( ( ( phi(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js) 
      !
      ! Then the velocities at all points in this horizontal slice and at all sigma 
      ! (vertical) locations.  
      !
      !
      ! Compute dphi/dsigma
      !
      CALL DiffZArbitrary(phi,W,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
           FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
      IF (FineGrid%Nx>1) THEN
         !
         ! Compute dphi/dx 
         !	
         CALL DiffXEven(phi,U,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
              FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,alpha)
         IF ( LinearOnOff /= 0) THEN
            !
            ! Add in the chain rule contribution to get the velocity
            !
            Do j=1,Ny
               Do i=1,Nx
                  Do k=1,Nz
                     U(k,i,j) = U(k,i,j) + ((1-z(k))/d(i,j)*hx(i,j)-z(k)/d(i,j)*etax(i,j))*W(k,i,j)
                  END Do
               END Do
            END Do
         END IF
      ELSE
         U=zero
      END IF
      !
      WRITE (fileop(io+1)) ( ( ( U(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)
      ! 
      IF (FineGrid%Ny>1) THEN
         ! dphi/dy
         CALL DiffYEven(phi,V,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
              FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,beta)
         IF ( LinearOnOff /= 0) THEN
            Do j=1,Ny
               Do i=1,Nx
                  Do k=1,Nz
                     V(k,i,j) = V(k,i,j)+((1-z(k))/d(i,j)*hy(i,j)-z(k)/d(i,j)*etay(i,j))*W(k,i,j)
                  END Do
               END Do
            END Do
         END IF
      ELSE
         V=zero
      END IF
      WRITE (fileop(io+1)) ( ( ( V(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)
      !
      ! Write the vertical velocity
      !
      WRITE (fileop(io+1)) ( ( ( W(k,i,j)/d(i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js) 
      !
      ! Compute du/dsigma and write du/dz to disc
      !
      CALL DiffZArbitrary(U,W,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
           FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
      WRITE (fileop(io+1)) ( ( ( W(k,i,j)/d(i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js) 
      !
      CALL DiffZArbitrary(V,W,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
           FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
      WRITE (fileop(io+1)) ( ( ( W(k,i,j)/d(i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js) 
   END IF
End IF

END SUBROUTINE StoreKinematicData
