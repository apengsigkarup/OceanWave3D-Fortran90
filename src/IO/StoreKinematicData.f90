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
USE hl_hdf5
use kinematicsArray
IMPLICIT NONE
! Input parameters
INTEGER :: Nx, Ny, Nz, io, it
! Local variables
INTEGER ::  i, j, k, i0, i1, is, j0, j1, js
INTEGER :: FOUT, i3
REAL(KIND=long), DIMENSION(:,:), POINTER :: x, y, h, hx, hy, eta, etax, etay
REAL(KIND=long), DIMENSION(:), POINTER   :: z, tmpxval
! Automatic work space
REAL(KIND=long) :: U(Nz,Nx,Ny), V(Nz,Nx,Ny), W(Nz,Nx,Ny), d(Nx,Ny)
REAL(KIND=long) :: tmpx(Nx), tmpy(Ny)
CHARACTER(len=30) :: form
REAL(KIND=long) :: hint, etaint, dint
REAL(KIND=long) :: Uint(Nz), Vint(Nz)
REAL(KIND=long) :: Ux(Nz,Nx,Ny),Uy(Nz,Nx,Ny),Uz(Nz,Nx,Ny)
REAL(KIND=long) :: Vx(Nz,Nx,Ny),Vy(Nz,Nx,Ny),Vz(Nz,Nx,Ny)
REAL(KIND=long) ::                           Wz(Nz,Nx,Ny)
CHARACTER(LEN=30) :: h5file
INTEGER(HID_T) :: extended_dimension_id
INTEGER(HSIZE_T), ALLOCATABLE :: dims_ext(:)
INTEGER(HSIZE_T), SAVE :: maxdims1(1), maxdims2(3), maxdims3(4), &
                    chunkdims1(1), chunkdims2(3), chunkdims3(4), &
                    extdims1(1), extdims2(3), extdims3(4)
INTEGER(HSIZE_T):: nx_save, ny_save, nz_save, onei = 1, zeroi = 0
REAL(KIND=long) :: x3d(Nz, Ny, Nx), y3d(Nz, Ny, Nx), z3d(Nz, Ny, Nx)
! Assign the local pointers
!
x => FineGrid%x; y => FineGrid%y; z => FineGrid%z; h => FineGrid%h; hx => FineGrid%hx
hy => FineGrid%hy; eta => WaveField%E; etax => WaveField%Ex; etay => WaveField%Ey
!
! Shift the horizontal grid point position to account for the ghost points.  
!
IF(FORMATTYPE/=22)THEN
  i0=Output(io)%xbeg+GhostGridX; i1=Output(io)%xend+GhostGridX; is=Output(io)%xstride; 
  j0=Output(io)%ybeg+GhostGridY; j1=Output(io)%yend+GhostGridY; js=Output(io)%ystride; 


   ! Some parameters that are necessary for h5 saving
  nx_save = size((/(i, i=i0,i1,is)/))
  ny_save = size((/(j, j=j0,j1,js)/))
  nz_save = Nz-GhostGridZ

  ! Set max dimensions for h5 writing
  ! FP 20190705: According to my tests, this syntax works on HDF5/1.8.17 and HDF5/1.8.5
  maxdims1 = (/-1/)
  maxdims2 = (/-1, -1, -1/)
  maxdims3 = (/-1, -1, -1, -1/)
  ! FP 20190705: According to my tests, this syntax works only on HDF5/1.8.5
  ! maxdims1 = (/H5S_UNLIMITED_F/)
  ! maxdims2 = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
  ! maxdims3 = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
  
  ! Chunk dims
  ! FP 20190705: The chunk dims is quite big in the time-axis.
  ! This is done to facilitate the access of the data, which are usually read
  ! as whole time series: usually, all the time steps of one dataset are read at one 
  ! particular x,y,z location
  
  chunkdims1 = (/30000*onei/)
  chunkdims2 = (/ny_save, nx_save, 30000*onei/)
  chunkdims3 = (/nz_save, ny_save, nx_save, 30000*onei/)
  
  ! Dimensions of the extended dataset, when appending to an existing dataset.
  extdims1 = (/onei/)
  extdims2 = (/ny_save, nx_save, onei/)
  extdims3 = (/nz_save, ny_save, nx_save, onei/)

ELSE
 ! PRINT*,'Storing kinematics data...'
  ! determine indices for stencils
  tmpx = x(1:Nx,1)-Output(io)%x
  ! search
  DO i=1,Nx
    IF(tmpx(i)>0)THEN
       Output(io)%idx(1) = i-alpha
       Output(io)%idx(2) = i+alpha
       EXIT ! Out of DO loop
    END IF
  END DO
  IF(Ny>1)THEN
    tmpy = y(1,1:Ny)-Output(io)%y
    DO j=1,Ny
      IF(tmpy(j)>0)THEN
         Output(io)%idx(3) = j-beta
         Output(io)%idx(4) = j+beta
         EXIT ! Out of DO loop
      END IF
    END DO
  ELSE
    Output(io)%idx(3) = 1
    Output(io)%idx(4) = 1
  ENDIF
  ! determine stencil weights for the interpolation
  ALLOCATE( Output(io)%stencilx(2*alpha+1) )
  ALLOCATE( tmpxval(2*alpha+1) )
  tmpxval = x(Output(io)%idx(1) : Output(io)%idx(2),1)
  tmpxval(alpha+1) = Output(io)%x
  CALL TaylorFDStencils1DArbitrary(alpha,alpha,0,Output(io)%stencilx,tmpxval)
!  Output(io)%stencilx = zero
!  print*,'stencilx = ',Output(io)%stencilx
  IF(Ny>1)THEN
     ALLOCATE( Output(io)%stencily(2*beta+1) )
     ! weights in stream_func_wave_finite.f
     CALL TaylorFDStencils1DArbitrary(beta,beta,0,Output(io)%stencily,Output(io)%y)
  ENDIF
  ! formattype=22
  i0=Output(io)%idx(1) ! xmin
  i1=Output(io)%idx(2) ! xmax
  j0=Output(io)%idx(3) ! ymin
  j1=Output(io)%idx(4) ! ymax
!  print*,'i0=',i0
!  print*,'i1=',i1
!  print*,'j0=',j0
!  print*,'j1=',j1
END IF


! Determine fileoutput
IF(formattype==21)THEN
   WRITE(unit=filename, FMT="(A,I2.2,A,I5.5,A)") "Kinematics_",io,"_",it,".bin"
   form="unformatted" ! binary format chosen
   FOUT = 22 ! file handle
   OPEN (unit=FOUT, file=filename,form=form)
!   WRITE(*,FMT='(A,A)') '  File output = ',filename
ELSE IF(formattype==22)THEN
   WRITE(unit=filename, FMT="(A,I2.2,A,I5.5,A)") "Kinematics_",io,"_",it,".bin"
   form="unformatted" ! binary format chosen
   FOUT = 22 ! file handle
   OPEN (unit=FOUT, file=filename,form=form)
   WRITE(*,FMT='(A,A)') '  File output = ',filename
ELSE IF(formattype==30)THEN   
   WRITE(h5file, "(A18,I3.3,A3)") "WaveKinematicsZone",fntH5(io),".h5"
   WRITE(*,FMT='(A,A)') '  File output of h5 file number = ',h5file
ELSE
   FOUT = FILEOP(io+1)
!   WRITE(*,FMT='(A,I2)') '  File output unit number = ',FOUT
END IF

IF(it==0)THEN
   !
   ! Save the grid data on the first call
   !
   IF(formattype==22)THEN
!     WRITE (FOUT) Nx,Ny,Nz
     WRITE (FOUT) Output(io)%x, Output(io)%y, Output(io)%tbeg,Output(io)%tend, &
          Output(io)%tstride, dt, FineGrid%Nz+GhostGridZ
     IF(FineGrid%Nx>1 .AND. FineGrid%Ny>1) THEN
        PRINT*,'3D setup for formattype=22 in StoreKinematicData.f90 not setup yet.'
        STOP
!        hint   = DOT_PRODUCT( Output(io)%stencilx,h(i0:i1,j0:j1)   )
!        etaint = DOT_PRODUCT( Output(io)%stencilx,eta(i0:i1,j0:j1) )
!        dint   = hint + etaint
        WRITE (FOUT) hint, etaint, dint
     ELSE IF(FineGrid%Nx>1) THEN
        hint   = DOT_PRODUCT( Output(io)%stencilx,h(i0:i1,1)   )
        print*,'hint = ',hint
        etaint = DOT_PRODUCT( Output(io)%stencilx,eta(i0:i1,1) )
        print*,'eta = ',etaint
        dint   = hint + etaint
        print*,'dint = ',dint
        WRITE (FOUT) hint, etaint, dint
     ELSE IF(FineGrid%Ny>1) THEN
        PRINT*,'2D setup in y-direction for formattype=22 in StoreKinematicData.f90 not setup yet.'
        STOP
     END IF
!     WRITE (FOUT) ( Output(io)%x, Output(io%y, hint )
!     WRITE (FOUT) ( ( x(i,j), y(i,j), h(i,j), hx(i,j), hy(i,j), &
!          i=i0,i1,is ), j=j0,j1,js )
     ! are these z values the sigma values??
     WRITE (FOUT) (z(i),i=1,Nz)
     print*,'z=',z
     IF(curvilinearOnOff/=0)THEN
        Print *, 'StoreKinematicData:  Saving horizontal fluid velocities is not yet implemented for curvilinear grids.'
     END IF
   ELSEIF (formattype == 30) THEN !hdf5 file
         ! FabioPierella 20190319
         ! Output in HDF5 files. We need some preparation on the dataset before
         ! we can output it to file.
         ! 1. Fortran is column-major. We want to therefore transpose the arrays before output to have a row-major HDF5 file structure.
         ! 2. This HDF5 file output wants  to be consistent with OW3D-GPU version, therefore some variables need to be made 3D before being output (e.g. the position arrays)

         ! The kinematicsArray will hold 5 kinematics timesteps.
         ! They will be used to compute the kinematic acceleration at the location
         ! so that we can save it in the .h5 file.
         call allocatePointers(Zones(io), nx_save, ny_save, nz_save)
         call increment_timestep_counter(Zones(io))

         ! Initialize all datasets
         WRITE(h5file, "(A18,I3.3,A3)") "WaveKinematicsZone",fntH5(io),".h5"
         ! Create
         call h5_dataset_create_chunked(h5file, 'time', INT(1, HID_T), &
                  & extdims1, maxdims1, chunkdims1) 
         ! Surface elevation variables
         call h5_dataset_create_chunked(h5file, 'surface_elevation', INT(3, HID_T), &
                  & extdims2, maxdims2, chunkdims2)
         call h5_dataset_create_chunked(h5file, 'surface_elevation_derivative_etax', INT(3, HID_T), &
                  & extdims2, maxdims2, chunkdims2)
         call h5_dataset_create_chunked(h5file, 'surface_elevation_derivative_etay', INT(3, HID_T), &
                  & extdims2, maxdims2, chunkdims2)
         ! Position variables 
         call h5_dataset_create_chunked(h5file, 'position_x', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)
         call h5_dataset_create_chunked(h5file, 'position_y', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)
         call h5_dataset_create_chunked(h5file, 'position_z', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)
         ! Velocity variables 
         call h5_dataset_create_chunked(h5file, 'velocity_u', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)
         call h5_dataset_create_chunked(h5file, 'velocity_v', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)
         call h5_dataset_create_chunked(h5file, 'velocity_w', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)
         ! Velocity x-gradient variables 
         call h5_dataset_create_chunked(h5file, 'velocity_derivative_ux', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)
         call h5_dataset_create_chunked(h5file, 'velocity_derivative_vx', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)
         ! Velocity y-gradient variables 
         call h5_dataset_create_chunked(h5file, 'velocity_derivative_uy', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)
         call h5_dataset_create_chunked(h5file, 'velocity_derivative_vy', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)      
         ! Velocity z-gradient variables 
         call h5_dataset_create_chunked(h5file, 'velocity_derivative_uz', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)
         call h5_dataset_create_chunked(h5file, 'velocity_derivative_vz', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)      
         call h5_dataset_create_chunked(h5file, 'velocity_derivative_wz', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)      

         ! Velocity kinematic accelerations
         call h5_dataset_create_chunked(h5file, 'velocity_derivative_ut', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)
         call h5_dataset_create_chunked(h5file, 'velocity_derivative_vt', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)      
         call h5_dataset_create_chunked(h5file, 'velocity_derivative_wt', INT(4, HID_T), &
                  & extdims3, maxdims3, chunkdims3)      

         ! Write the first timestep
         ! call h5_write(h5file, 'time', (/it*dt/))
         ! allocate the position arrays

         ! Fluid thickness
         DO j=1,Ny
            DO i=1,Nx
                d(i,j)=h(i,j)+eta(i,j);
            END DO
         END DO

         do k=1,Nz
            do j=1,Ny
               do i=1,Nx
                  x3d(k,j,i) = x(i,j)
                  y3d(k,j,i) = y(i,j)
                  z3d(k,j,i)  = z(k)*d(i,j)-h(i,j)
               end do 
            end do 
         end do

         extended_dimension_id = 4
         call h5_write(h5file, 'position_x', x3d(1+GhostGridZ:, j0:j1:js, i0:i1:is))                  
         call h5_write(h5file, 'position_y', y3d(1+GhostGridZ:, j0:j1:js, i0:i1:is))                          
         call h5_write(h5file, 'position_z', z3d(1+GhostGridZ:, j0:j1:js, i0:i1:is))

      IF(curvilinearOnOff/=0)THEN
      Print *, 'StoreKinematicData:  Saving horizontal fluid velocities is not yet implemented for curvilinear grids.'
      END IF
   ELSE
     ! formattype /= 22
     write (FOUT) Output(io)%xbeg,Output(io)%xend,Output(io)%xstride, &
          Output(io)%ybeg, Output(io)%yend, Output(io)%ystride,               &
          Output(io)%tbeg,Output(io)%tend,Output(io)%tstride, dt, Nz
     !
     WRITE (FOUT) ( ( x(i,j), y(i,j), h(i,j), hx(i,j), hy(i,j), &
          i=i0,i1,is ), j=j0,j1,js ) 
     WRITE (FOUT) (z(i),i=1,Nz)
     IF(curvilinearOnOff/=0)THEN
        Print *, 'StoreKinematicData:  Saving horizontal fluid velocities is not yet implemented for curvilinear grids.'
     END IF
   END IF

ELSE !IF(it==0)THEN

   !
   IF(curvilinearOnOff == 0)THEN
      IF(formattype==22)THEN
         !
         ! Dump free surface elevation, still water depth and kinematics to the output file.
         ! Minimal no of variables and storage is output.
         !
         ! WRITE (FOUT) ( ( eta(i,j), i=i0,i1,is ), j=j0,j1,js)
         ! WRITE (FOUT) ( ( h(i,j), i=i0,i1,is ), j=j0,j1,js) 
         !                                                                  
         ! The fluid thickness d=h+eta
         !            
         DO j=1,Ny
            DO i=1,Nx
                d(i,j)=h(i,j)+eta(i,j);
            END DO
         END DO
        !  CALL DiffZArbitrary(phi,W,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
         !  FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
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
         ! Interpolate to point in question
         DO k = 1, FineGrid%Nz+GhostGridZ
            DO j=j0,j1
               Uint(k) = DOT_PRODUCT( Output(io)%stencilx,U(k,i0:i1,j) )
            END DO
         END DO
         WRITE (FOUT) Uint
         ! print*,'Uint:stencilx=',Output(io)%stencilx
         ! print*,'Uint = ',Uint
         ! WRITE (FOUT) ( ( ( U(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)
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
            ! WRITE (FOUT) ( ( ( V(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)
            ! Interpolate to point in question
            DO k = 1, FineGrid%Nz+GhostGridZ
               DO j=j0,j1
                  Vint(k) = DOT_PRODUCT( Output(io)%stencilx,V(k,i0:i1,j) )
               END DO
            END DO
            WRITE (FOUT) Vint
         ELSE
            V=zero
         END IF
      !
      ! Write the vertical velocity
      !
      ! WRITE (FOUT) ( ( ( W(k,i,j)/d(i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js) 
      ELSE 
         !
         ! Dump this solution slice to the output file
         !
         ! First the free surface elevation and gradient at all points in this slice
         ! 

         !
         ! The fluid thickness d=h+eta
         !
         DO j=1,Ny
            DO i=1,Nx
               d(i,j)=h(i,j)+eta(i,j);
            END DO
         END DO

         do k=1,Nz
            do j=1,Ny
               do i=1,Nx
                  x3d(k,j,i) = x(i,j)
                  y3d(k,j,i) = y(i,j)
                  z3d(k,j,i)  = z(k)*d(i,j)-h(i,j)
               end do 
            end do 
         end do
         
         !
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

         !
         !--------------------------------------- d /dsigma
         !
         CALL DiffZArbitrary(U,Uz,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
            FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
         CALL DiffZArbitrary(V,Vz,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
            FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
         CALL DiffZArbitrary(W,Wz,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
         FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,gamma)
         ! ---- See below: we still need to divide all three by fluid thickness
         !


         !--------------------------------------- d /dy
         IF (FineGrid%Ny>1) THEN
         CALL DiffYEven(U,Uy,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
               FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,beta)
         CALL DiffYEven(V,Vy,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
               FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,beta)

         ! ---- See below: we still need to divide dw/dy by fluid thickness   
            
         ELSE
            Uy=zero
            Vy=zero
         END IF

         !--------------------------------------- d /dx
         IF (FineGrid%Nx>1) THEN
         CALL DiffXEven(U,Ux,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
               FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,alpha)
         CALL DiffXEven(V,Vx,1,FineGrid%Nx+2*GhostGridX,FineGrid%Ny+2*GhostGridY,  &
               FineGrid%Nz+GhostGridZ,FineGrid%DiffStencils,alpha)               
         ! ---- See below: we still need to divide dw/dx by fluid thickness                        
         ELSE
            Ux=zero
            Vx=zero
         END IF

      ! sigs ends lines

         ! After computation, store according to user choice
         IF (formattype==30)THEN

            ! Save the velocity information in the array of Kinematics, for later
            ! computation of the kinematic acceleration

            call cycle(Zones(io)%Kinematics)
            Zones(io)%Kinematics(5)%U = U(1+GhostGridZ:, i0:i1:is, j0:j1:js);
            Zones(io)%Kinematics(5)%V = V(1+GhostGridZ:, i0:i1:is, j0:j1:js);
            DO K=1+GhostGridZ,NZ
               Zones(io)%Kinematics(5)%W(K-GhostGridZ,:,:) = W(K, i0:i1:is, j0:j1:js)/d(i0:i1:is, j0:j1:js);
               Zones(io)%Kinematics(5)%Uz(K-GhostGridZ,:,:) = Uz(K, i0:i1:is, j0:j1:js)/d(i0:i1:is, j0:j1:js);
               Zones(io)%Kinematics(5)%Vz(K-GhostGridZ,:,:) = Vz(K, i0:i1:is, j0:j1:js)/d(i0:i1:is, j0:j1:js);               
               !FP20190430: We need to divide Wz by twice the fluid thickness. Before today it was mistakenly divided only once.
               Zones(io)%Kinematics(5)%Wz(K-GhostGridZ,:,:) = Wz(K, i0:i1:is, j0:j1:js)/d(i0:i1:is, j0:j1:js)**2; 

            END DO
            Zones(io)%Kinematics(5)%Eta = Eta(i0:i1:is, j0:j1:js);

            call increment_timestep_counter(Zones(io)) ! signal that we have added another timestep

            if (Zones(io)%number_of_saved_timesteps == 5) then 
               ! if we have saved more than 5 timesteps  then we calculate the centered 5-points derivative
               call calculateKinAcceleration(Zones(io), dt*Output(io)%tstride, z(1+GhostGridZ:))
            end if 

            WRITE(h5file, "(A18,I3.3,A3)") "WaveKinematicsZone",fntH5(io),".h5"

            if (it == Output(io)%tbeg) then

               ! FP20190420 Unfortunately this is a little cumbersome.
               ! At step 1 we need to overwrite what was written at step 0, where the dataset
               ! was created. This is because we cannot create a dataset without writing in it, so
               ! when we get to step 1 ther is already what was written in step 0.
               extended_dimension_id = 1
               ! call h5_extend(h5file, 'time', extended_dimension_id, extdims1, (/it*dt/))
               call h5_write_at_step(h5file, 'time', extended_dimension_id, 0, extdims1, (/it*dt/))

               extended_dimension_id = 4
               call h5_write_at_step(h5file, 'position_x', extended_dimension_id, 0, extdims3, &
                  & x3d(1+GhostGridZ:,j0:j1:js, i0:i1:is))         
               call h5_write_at_step(h5file, 'position_y', extended_dimension_id, 0, extdims3, &
                  & y3d(1+GhostGridZ:,j0:j1:js, i0:i1:is)) 

               call h5_write_at_step(h5file, 'position_z', extended_dimension_id, 0, extdims3, &
                  & z3d(1+GhostGridZ:, j0:j1:js, i0:i1:is))

               extended_dimension_id = 3
               call h5_write_at_step(h5file, 'surface_elevation', extended_dimension_id,  0,extdims2, &
                  & transpose(eta(i0:i1:is, j0:j1:js)))
               
               ! Write only if there is more than one point in the x direction
               if (FineGrid%Nx>1) then
                  call h5_write_at_step(h5file, 'surface_elevation_derivative_etax', extended_dimension_id, 0, extdims2, &
                  & transpose(etax(i0:i1:is, j0:j1:js)))            
               else ! just write 0s
                  call h5_write_at_step(h5file, 'surface_elevation_derivative_etax', extended_dimension_id, 0, extdims2, &
                  & reshape((/(zero, i=1,nx_save*ny_save)/), shape=(/nx_save, ny_save/)))     
               end if

               ! Write only if there is more than one point in the y direction
               if (FineGrid%Ny>1) then
                  call h5_write_at_step(h5file, 'surface_elevation_derivative_etay', extended_dimension_id, 0, extdims2, &
                  & transpose(etay(i0:i1:is, j0:j1:js)))            
               else ! just write 0s
                  call h5_write_at_step(h5file, 'surface_elevation_derivative_etay', extended_dimension_id, 0, extdims2, &
                  & reshape((/(zero, i=1,nx_save*ny_save)/), shape=(/nx_save, ny_save/))) 
               end if
            
               extended_dimension_id = 4
               ! velocities
               call h5_write_at_step(h5file, 'velocity_u', extended_dimension_id, 0, extdims3, &
            & reshape(U(1+GhostGridZ:, i0:i1:is, j0:j1:js), shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))
               call h5_write_at_step(h5file, 'velocity_v', extended_dimension_id, 0, extdims3, &
               & reshape(V(1+GhostGridZ:, i0:i1:is, j0:j1:js), shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))
               call h5_write_at_step(h5file, 'velocity_w', extended_dimension_id, 0, extdims3, &
               & reshape(Zones(io)%Kinematics(5)%W, shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))

               ! Kinematics accelerations
               ! Writes a zero with the shape of (Nz_save, Ny_save, Nx_save), since we don't have enough
               ! information to compute a correct version of the acceleration.
               call h5_write_at_step(h5file, 'velocity_derivative_ut', extended_dimension_id, 0, extdims3, &
            & reshape(Zones(io)%Kinematics(3)%Ut * 0., shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_write_at_step(h5file, 'velocity_derivative_vt', extended_dimension_id, 0, extdims3, &
               & reshape(Zones(io)%Kinematics(3)%Vt * 0., shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_write_at_step(h5file, 'velocity_derivative_wt', extended_dimension_id, 0, extdims3, &
               & reshape(Zones(io)%Kinematics(3)%Wt * 0., shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))

               ! velocity z gradients
               call h5_write_at_step(h5file, 'velocity_derivative_uz', extended_dimension_id, 0, extdims3, &
            & reshape(Zones(io)%Kinematics(5)%Uz, shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_write_at_step(h5file, 'velocity_derivative_vz', extended_dimension_id, 0, extdims3, &
               & reshape(Zones(io)%Kinematics(5)%Vz, shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_write_at_step(h5file, 'velocity_derivative_wz', extended_dimension_id, 0, extdims3, &
               & reshape(Zones(io)%Kinematics(5)%Wz, shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))                            


               ! velocity x gradients
               call h5_write_at_step(h5file, 'velocity_derivative_ux', extended_dimension_id, 0, extdims3, &
            & reshape(Ux(1+GhostGridZ:, i0:i1:is, j0:j1:js), shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_write_at_step(h5file, 'velocity_derivative_vx', extended_dimension_id, 0, extdims3, &
               & reshape(Vx(1+GhostGridZ:, i0:i1:is, j0:j1:js), shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    

               ! velocity y gradients
               call h5_write_at_step(h5file, 'velocity_derivative_uy', extended_dimension_id, 0, extdims3, &
            & reshape(Uy(1+GhostGridZ:, i0:i1:is, j0:j1:js), shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_write_at_step(h5file, 'velocity_derivative_vy', extended_dimension_id, 0, extdims3, &
               & reshape(Vy(1+GhostGridZ:, i0:i1:is, j0:j1:js), shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/))) 
            
            else ! if (it == 1) then

               extended_dimension_id = 1
               ! call h5_extend(h5file, 'time', extended_dimension_id, extdims1, (/it*dt/))
               call h5_extend(h5file, 'time', extended_dimension_id, extdims1, (/it*dt/))

               extended_dimension_id = 4
               call h5_extend(h5file, 'position_x', extended_dimension_id, extdims3, &
                  & x3d(1+GhostGridZ:,j0:j1:js, i0:i1:is))         
               call h5_extend(h5file, 'position_y', extended_dimension_id, extdims3, &
                  & y3d(1+GhostGridZ:,j0:j1:js, i0:i1:is)) 

               call h5_extend(h5file, 'position_z', extended_dimension_id, extdims3, &
                  & z3d(1+GhostGridZ:, j0:j1:js, i0:i1:is))

               extended_dimension_id = 3
               call h5_extend(h5file, 'surface_elevation', extended_dimension_id, extdims2, &
                  & transpose(eta(i0:i1:is, j0:j1:js)))
               
               ! Write only if there is more than one point in the x direction
               if (FineGrid%Nx>1) then
                  call h5_extend(h5file, 'surface_elevation_derivative_etax', extended_dimension_id, extdims2, &
                  & transpose(etax(i0:i1:is, j0:j1:js)))            
               else ! just write 0s
                  call h5_extend(h5file, 'surface_elevation_derivative_etax', extended_dimension_id, extdims2, &
                  & reshape((/(zero, i=1,nx_save*ny_save)/), shape=(/nx_save, ny_save/)))     
               end if

               ! Write only if there is more than one point in the y direction
               if (FineGrid%Ny>1) then
                  call h5_extend(h5file, 'surface_elevation_derivative_etay', extended_dimension_id, extdims2, &
                  & transpose(etay(i0:i1:is, j0:j1:js)))            
               else ! just write 0s
                  call h5_extend(h5file, 'surface_elevation_derivative_etay', extended_dimension_id, extdims2, &
                  & reshape((/(zero, i=1,nx_save*ny_save)/), shape=(/nx_save, ny_save/))) 
               end if
            
               extended_dimension_id = 4
               ! velocities
               call h5_extend(h5file, 'velocity_u', extended_dimension_id, extdims3, &
            & reshape(U(1+GhostGridZ:, i0:i1:is, j0:j1:js), shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))
               call h5_extend(h5file, 'velocity_v', extended_dimension_id, extdims3, &
               & reshape(V(1+GhostGridZ:, i0:i1:is, j0:j1:js), shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))
               call h5_extend(h5file, 'velocity_w', extended_dimension_id, extdims3, &
               & reshape(Zones(io)%Kinematics(5)%W, shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))

               ! Kinematics accelerations
               ! They are a little different from the others. Since we compute them with a centered scheme with 5 points
               ! We need always to first extend the file, writing a "zero" value, so that it has the same nr. of points
               ! as the other quantities.
               
               call h5_extend(h5file, 'velocity_derivative_ut', extended_dimension_id, extdims3, &
            & reshape(Zones(io)%Kinematics(5)%Ut *0., shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_extend(h5file, 'velocity_derivative_vt', extended_dimension_id, extdims3, &
               & reshape(Zones(io)%Kinematics(5)%Vt *0., shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_extend(h5file, 'velocity_derivative_wt', extended_dimension_id, extdims3, &
               & reshape(Zones(io)%Kinematics(5)%Wt *0., shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))

               ! But then, if we are beyond it = 3, we start writing the centered scheme to (it-3)
               ! reason behind MAX((/it-3, 0/): it-3 < 0 until we are at step 3, where we can start writing at the 
               ! right place (0 to N)
               
               i3 = (((it-2*Output(io)%tstride)-Output(io)%tbeg)/Output(io)%tstride)
               call h5_write_at_step(h5file, 'velocity_derivative_ut', extended_dimension_id, MAXVAL((/i3, 0/)), extdims3, &
            & reshape(Zones(io)%Kinematics(3)%Ut, shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_write_at_step(h5file, 'velocity_derivative_vt', extended_dimension_id, MAXVAL((/i3, 0/)), extdims3, &
               & reshape(Zones(io)%Kinematics(3)%Vt, shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_write_at_step(h5file, 'velocity_derivative_wt', extended_dimension_id, MAXVAL((/i3, 0/)), extdims3, &
               & reshape(Zones(io)%Kinematics(3)%Wt, shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))

               ! velocity z gradients
               call h5_extend(h5file, 'velocity_derivative_uz', extended_dimension_id, extdims3, &
            & reshape(Zones(io)%Kinematics(5)%Uz, shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_extend(h5file, 'velocity_derivative_vz', extended_dimension_id, extdims3, &
               & reshape(Zones(io)%Kinematics(5)%Vz, shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_extend(h5file, 'velocity_derivative_wz', extended_dimension_id, extdims3, &
               & reshape(Zones(io)%Kinematics(5)%Wz, shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))                            


               ! velocity x gradients
               call h5_extend(h5file, 'velocity_derivative_ux', extended_dimension_id, extdims3, &
            & reshape(Ux(1+GhostGridZ:, i0:i1:is, j0:j1:js), shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_extend(h5file, 'velocity_derivative_vx', extended_dimension_id, extdims3, &
               & reshape(Vx(1+GhostGridZ:, i0:i1:is, j0:j1:js), shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    

               ! velocity y gradients
               call h5_extend(h5file, 'velocity_derivative_uy', extended_dimension_id, extdims3, &
            & reshape(Uy(1+GhostGridZ:, i0:i1:is, j0:j1:js), shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))    
               call h5_extend(h5file, 'velocity_derivative_vy', extended_dimension_id, extdims3, &
               & reshape(Vy(1+GhostGridZ:, i0:i1:is, j0:j1:js), shape=(/nz_save, ny_save, nx_save/), order=(/1,3,2/)))                
            end if ! if (it == 1) then
         else

            ! Dump the surface elevation and its gradients
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

            ! Then the velocity potential at all points in this horizontal slice and at all sigma 
            ! (vertical) locations.  
            !
            WRITE (FOUT) ( ( ( phi(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js) 

      !  
            ! Then the U, V, W velocity   
            WRITE (FOUT) ( ( ( U(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)
            WRITE (FOUT) ( ( ( V(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)
            WRITE (FOUT) ( ( ( W(k,i,j)/d(i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js) 

            ! Then the U, V, W velocity   z  gradients
            WRITE (FOUT) ( ( ( Uz(k,i,j)/d(i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js) 
            WRITE (FOUT) ( ( ( Vz(k,i,j)/d(i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)             
            WRITE (FOUT) ( ( ( Wz(k,i,j)/d(i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)             

            ! Then the U, V velocity   y gradients
            WRITE (FOUT) ( ( ( Uy(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)                  
            WRITE (FOUT) ( ( ( Vy(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)            

            ! Then the U, V velocity   x gradients
            WRITE (FOUT) ( ( ( Ux(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)
            WRITE (FOUT) ( ( ( Vx(k,i,j), k=1,FineGrid%Nz+GhostGridZ), i=i0,i1,is), j=j0,j1,js)

         end if    


      END IF

   END IF
END IF

if (formattype == 22 .or. formattype == 21) then 
   CLOSE(FOUT)
end if

contains

subroutine chainRuleContribution(InOutArray, h_gradient, eta_gradient)

   real(kind=long),intent(INOUT) :: InOutArray(:,:,:)
   real(kind=long),intent(IN) :: eta_gradient(:,:), h_gradient(:,:)

   Do j=1,Ny
      Do i=1,Nx
         Do k=1,Nz
            InOutArray(k,i,j) = InOutArray(k,i,j) + ((1-z(k))/d(i,j)*h_gradient(i,j)-z(k)/d(i,j)*eta_gradient(i,j))*W(k,i,j)
         END Do
      END Do
   END Do

end subroutine chainRuleContribution


END SUBROUTINE StoreKinematicData

