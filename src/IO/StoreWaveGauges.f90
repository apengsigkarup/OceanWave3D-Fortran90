!
! ****************************************************************
!
SUBROUTINE StoreWaveGauges(Nx,Ny,Nz,io,it)
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
REAL(KIND=long), DIMENSION(:,:), POINTER :: x, y, h, eta
REAL(KIND=long), DIMENSION(:), POINTER   :: z
! Automatic work space
REAL(KIND=long) :: d(Nx,Ny)
!
! Assign the local pointers
!
x => FineGrid%x; y => FineGrid%y; z => FineGrid%z; h => FineGrid%h; 
eta => WaveField%E
!
! Shift the horizontal grid point position to account for the ghost points.  
!


!; i1=Output(io)%xend+GhostGridX; is=Output(io)%xstride; 
!j0=Output(io)%ybeg+GhostGridY; j1=Output(io)%yend+GhostGridY; js=Output(io)%ystride; 

IF(it==0)THEN
   !
   ! Save the grid data on the first call
   WRITE(fileop(io),'(A)',ADVANCE='no') 'x: '
   DO i = 1,nOutFiles
           WRITE(fileop(io),'(F10.2)',advance='no') (x(Output(i)%xbeg+GhostGridX,Output(i)%ybeg+GhostGridY))
   END DO
   WRITE(fileop(io),*) ''

   WRITE(fileop(io),'(A)',ADVANCE='no') 'y: '
   DO i = 1,nOutFiles
           WRITE(fileop(io),'(F10.2)',advance='no') (y(Output(i)%xbeg+GhostGridX,Output(i)%ybeg+GhostGridY))
   END DO

   WRITE(fileop(io),*) ''
   WRITE(fileop(io),'(A)') 'time, wg1, wg2, ..., wgN '
   
   IF(curvilinearOnOff/=0)THEN
      Print *, 'StoreKinematicData:  Saving horizontal fluid velocities is not yet implemented for curvilinear grids.'
   END IF
ELSE
   !
   IF(curvilinearOnOff == 0)THEN
           write(fileop(io),'(F10.3 )',advance='no') time
           DO i = 1,nOutFiles
                   WRITE(fileop(io),'(F10.3)',advance='no') (eta(Output(i)%xbeg+GhostGridX,Output(i)%ybeg+GhostGridY))
           END DO
   WRITE(fileop(io),*) ''
    ENDIF
End IF
FLUSH(fileop(io))
END SUBROUTINE StoreWaveGauges
