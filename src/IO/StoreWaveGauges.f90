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

IF(it==0)THEN
   !
   ! Save the grid data on the first call
   WRITE(fileop(io),'(A)',ADVANCE='no') 'Time'
   DO i = 1,nOutFiles
           WRITE(fileop(io),'(A)',advance='no') ' WG'
   END DO
   WRITE(fileop(io),*) ''

   WRITE(fileop(io),'(A)',ADVANCE='no') '-1 '
   DO i = 1,nOutFiles
           WRITE(fileop(io),'(F10.6)',advance='no') (x(Output(i)%xbeg+GhostGridX,Output(i)%ybeg+GhostGridY))
   END DO
   WRITE(fileop(io),*) ''


   WRITE(fileop(io),'(A)',ADVANCE='no') '-2 '
   DO i = 1,nOutFiles
           WRITE(fileop(io),'(F10.6)',advance='no') (y(Output(i)%xbeg+GhostGridX,Output(i)%ybeg+GhostGridY))
   END DO
   WRITE(fileop(io),*) ''


   WRITE(fileop(io),'(A)',ADVANCE='no') '-3 '
   DO i = 1,nOutFiles
           WRITE(fileop(io),'(F10.6)',advance='no') (h(Output(i)%xbeg+GhostGridX,Output(i)%ybeg+GhostGridY))
   END DO
   WRITE(fileop(io),*) ''


   
   IF(curvilinearOnOff/=0)THEN
      Print *, 'StoreKinematicData:  Saving horizontal fluid velocities is not yet implemented for curvilinear grids.'
   END IF
ELSE
   !
   IF(curvilinearOnOff == 0)THEN
           write(fileop(io),'(F14.6 )',advance='no') time
           DO i = 1,nOutFiles
                   WRITE(fileop(io),'(F10.6)',advance='no') (eta(Output(i)%xbeg+GhostGridX,Output(i)%ybeg+GhostGridY))
           END DO
   WRITE(fileop(io),*) ''
    ENDIF
End IF
FLUSH(fileop(io))
END SUBROUTINE StoreWaveGauges
