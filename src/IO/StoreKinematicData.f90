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
INTEGER :: FOUT
REAL(KIND=long), DIMENSION(:,:), POINTER :: x, y, h, hx, hy, eta, etax, etay
REAL(KIND=long), DIMENSION(:), POINTER   :: z
! Automatic work space
REAL(KIND=long) :: U(Nz,Nx,Ny), V(Nz,Nx,Ny), W(Nz,Nx,Ny), d(Nx,Ny)
CHARACTER(len=30) :: form
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
! Determine fileoutput
print*,'formattype=',formattype
IF(formattype==21)THEN
   WRITE(unit=filename, FMT="(A,I2.2,A,I5.5,A)") "Kinematics_",io,"_",it,".bin"
   form="unformatted" ! binary format chosen
   FOUT = 22 ! file handle
   OPEN (unit=FOUT, file=filename,form=form)
   WRITE(*,FMT='(A,A)') '  File output = ',filename
ELSE
   FOUT = FILEOP(io+1)
   WRITE(*,FMT='(A,I2)') '  File output = ',FOUT
END IF
IF(it==0)THEN
   !
   ! Save the grid data on the first call
   !
   write (FOUT) Output(io)%xbeg,Output(io)%xend,Output(io)%xstride, &
        Output(io)%ybeg, Output(io)%yend, Output(io)%ystride,               &
        Output(io)%tbeg,Output(io)%tend,Output(io)%tstride, dt,             &
        FineGrid%Nz+GhostGridZ
   !
   WRITE (FOUT) ( ( x(i,j), y(i,j), h(i,j), hx(i,j), hy(i,j), &
        i=i0,i1,is ), j=j0,j1,js ) 
   WRITE (FOUT) (z(i),i=1,Nz)
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
      WRITE (FOUT) ( ( eta(i,j), i=i0,i1,is ), j=j0,j1,js) 
      IF(Nx > 1) THEN
         WRITE (FOUT) ( ( etax(i,j), i=i0,i1,is ), j=j0,j1,js) 
      ELSE
         WRITE (FOUT) ( ( zero, i=i0,i1,is ), j=j0,j1,js) 
      END IF
      IF(Ny > 1) THEN
         WRITE (FOUT) ( ( etay(i,j), i=i0,i1,is ), j=j0,j1,js) 
      ELSE
         WRITE (FOUT) ( ( zero, i=i0,i1,is ), j=j0,j1,js) 
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
      WRITE (FOUT) ( ( ( phi(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js) 
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
      WRITE (FOUT) ( ( ( U(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)
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
      WRITE (FOUT) ( ( ( V(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)
      !
      ! Write the vertical velocity
      !
      WRITE (FOUT) ( ( ( W(k,i,j)/d(i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js) 
      !
      ! Compute du/dsigma and write du/dz to disc
      !
      CALL DiffZArbitrary(U,W,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
           FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
      WRITE (FOUT) ( ( ( W(k,i,j)/d(i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js) 
      !
      CALL DiffZArbitrary(V,W,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
           FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
      WRITE (FOUT) ( ( ( W(k,i,j)/d(i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js) 
   END IF
END IF

IF(formattype==21)THEN
   ! close file access
   CLOSE(FOUT)
END IF
END SUBROUTINE StoreKinematicData
