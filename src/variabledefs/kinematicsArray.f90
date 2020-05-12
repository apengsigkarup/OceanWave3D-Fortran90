module kinematicsArray

! (c) Fabio Pirella 20190402
! In this module we create an array of pointers that will contain
! the velocities for the Wave Kinematics zone data storage (i == 20 or i == 30),
! see ReadInputFileParameters.f90
! We will save 5 time instants and then take their derivative in order to calculate 
! the kinematic acceleration and store it into the WaveKinematicsZone**.h5 file

use hdf5
implicit none

type kinArray
    real(kind=8), allocatable :: U(:,:,:)
    real(kind=8), allocatable :: V(:,:,:)
    real(kind=8), allocatable :: W(:,:,:)
    real(kind=8), allocatable :: Uz(:,:,:)
    real(kind=8), allocatable :: Vz(:,:,:)
    real(kind=8), allocatable :: Wx(:,:,:)
    real(kind=8), allocatable :: Wy(:,:,:)
    real(kind=8), allocatable :: Wz(:,:,:)
    real(kind=8), allocatable :: Ut(:,:,:)
    real(kind=8), allocatable :: Vt(:,:,:)
    real(kind=8), allocatable :: Wt(:,:,:)
    real(kind=8), allocatable :: Eta(:,:)
end type kinArray

type zoneKin
    type(kinArray), dimension(5) :: Kinematics !5: newest, 1:oldest
    integer :: number_of_saved_timesteps = 0
    integer :: id
end type zoneKin

contains

subroutine increment_timestep_counter(inZone)
    ! This subroutine increases a counter
    ! that says how many samples have been already saved
    ! The counter cannot be over 5 (max 5 timesteps are saved)

    implicit none
    type(zoneKin) :: inZone

    if (inZone%number_of_saved_timesteps < 5) then
        inZone%number_of_saved_timesteps = inZone%number_of_saved_timesteps +1
    end if

end subroutine increment_timestep_counter

subroutine allocatePointers(inZone, nx_save, ny_save, nz_save)

    implicit none
    type(zoneKin)  :: inZone
    integer(kind=8):: nx_save, ny_save, nz_save ! total of points to be saved
    integer        :: i

    DO I=1,5
        ALLOCATE(inZone%Kinematics(i)%U(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%V(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%W(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Ut(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Vt(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Wt(nz_save, nx_save, ny_save))            
        ALLOCATE(inZone%Kinematics(i)%Uz(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Vz(nz_save, nx_save, ny_save))
        ALLOCATE(inZone%Kinematics(i)%Wz(nz_save, nx_save, ny_save))            
        ALLOCATE(inZone%Kinematics(i)%Eta(nx_save, ny_save))

        inZone%Kinematics(i)%U = 0.;
        inZone%Kinematics(i)%V = 0.;
        inZone%Kinematics(i)%W = 0.;
        inZone%Kinematics(i)%Ut = 0.;
        inZone%Kinematics(i)%Vt = 0.;
        inZone%Kinematics(i)%Wt = 0.;
        inZone%Kinematics(i)%Uz = 0.;
        inZone%Kinematics(i)%Vz = 0.;
        inZone%Kinematics(i)%Wz = 0.;
        inZone%Kinematics(i)%Eta = 0.

    END DO

end subroutine allocatePointers

subroutine cycle(inKinArray)
    ! This subroutine cycles the pointers so that 
    ! we can replace the oldest timestep with the newest
    implicit none
    type(kinArray) :: inKinArray(:)

    inKinArray = cshift(inKinArray, 1)
end subroutine cycle

subroutine calculateKinAcceleration(inZone, dt, sigma)
    ! This subroutine calculate the derivative
    ! of eta in time for a certain zone

    implicit none
    type(zoneKin)  :: inZone
    real(kind=8)   :: dt
    real(kind=8)   :: sigma(:)
    real(kind=8)   :: time(5), Eta_t, &
                      Ut_nocorr, Vt_nocorr, Wt_nocorr        ! Accelerations without the correction
    integer        :: i, j, k, ii, it, &
                      time_steps, &     ! size of the stencil we use for time differenciation
                      alpha, &          ! half-width of the stencil
                      nz_save, nx_save, ny_save ! size of the zone
    real(kind=8), allocatable, save   :: FDstencil(:,:)


    time_steps = inZone%number_of_saved_timesteps
    alpha = 2
    
    ! Construct the time array
    Do i=1,5
        time(i)=dt*(i-1)
    END DO

    if (.NOT.(allocated(FDStencil))) allocate(FDStencil(time_steps, time_steps))

    CALL BuildStencilsGridX(alpha,1,time,FDStencil,5,1)

    nz_save = size(inZone%Kinematics(1)%U, 1) ! How many points in the z-direction we are saving
    nx_save = size(inZone%Kinematics(1)%U, 2) ! How many points in the x-direction we are saving
    ny_save = size(inZone%Kinematics(1)%U, 3) ! How many points in the y-direction we are saving
    

    DO j=1,ny_save
        DO i=1,nx_save
            Eta_t = &
            Dot_Product(FDStencil(3,:), (/(inZone%Kinematics(it)%Eta(i,j),it=1,5)/))
            ! Central, 5-points derivative of the surface elevation
            DO k=1,nz_save
                Ut_nocorr = &
                Dot_Product(FDStencil(3,:), (/(inZone%Kinematics(it)%U(k,i,j),it=1,5)/))
                Vt_nocorr = &
                Dot_Product(FDStencil(3,:), (/(inZone%Kinematics(it)%V(k,i,j),it=1,5)/))                        
                Wt_nocorr = &
                Dot_Product(FDStencil(3,:), (/(inZone%Kinematics(it)%W(k,i,j),it=1,5)/))

                inZone%Kinematics(3)%Ut(k,i,j) = Ut_nocorr - &
                    sigma(k)*inZone%Kinematics(3)%Uz(k,i,j)*Eta_t
                inZone%Kinematics(3)%Vt(k,i,j) = Vt_nocorr - &
                    sigma(k)*inZone%Kinematics(3)%Vz(k,i,j)*Eta_t
                inZone%Kinematics(3)%Wt(k,i,j) = Wt_nocorr - &
                    sigma(k)*inZone%Kinematics(3)%Wz(k,i,j)*Eta_t                
            END DO
        END DO
    END Do

end subroutine calculateKinAcceleration

end module kinematicsArray