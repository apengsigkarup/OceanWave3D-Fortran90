module HL_HDF5 !High level HDF5 interface


    USE HDF5 ! This module contains all necessary modules 
    USE ISO_C_BINDING
    IMPLICIT NONE
    INTEGER     ::   hdferr ! Error flag

    interface h5_write
        module procedure h5_write_1d
        module procedure H5_write_2d
        module procedure H5_write_3d
    end interface 

    interface h5_write_at_step
        module procedure h5_write_1d_at_step
        module procedure H5_write_2d_at_step
        module procedure H5_write_3d_at_step
    end interface 

    interface h5_extend
        module procedure h5_extend_1d
        module procedure h5_extend_2d
        module procedure h5_extend_3d
    end interface     

    contains

    subroutine h5_check_version()
        integer :: major = 0;
        integer :: minor = 0;
        integer :: patch = 0;
        integer :: error = 0;

        call h5check_version_f(major, minor, patch, error)

        print*, "This version of h5 is: ", major, ".", minor, ".", patch

    end subroutine h5_check_version

    subroutine h5_file_create(file_name, file_id)
       
        CHARACTER(*), intent(IN) :: file_name
        INTEGER(HID_T), intent(INOUT) :: file_id       ! File identifier
        logical :: dummy

        CALL h5open_f(hdferr)
        dummy = check_return_value(hdferr, "h5_file_create", "h5fopen_f")
        CALL h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, hdferr) !H5F_ACC_TRUNC_F overwrite existing file
        dummy = check_return_value(hdferr, "h5_file_create", "h5fcreate_f")
        CALL h5fclose_f(file_id, hdferr)
        dummy = check_return_value(hdferr, "h5_file_create", "h5fclose_f")
    end subroutine h5_file_create

    function h5_dataset_exists(file_name, dataset_name)

        character(*) :: file_name, dataset_name
        integer(HID_T) :: file
        logical :: isExisting, h5_dataset_exists, dummy

        call h5fopen_f(file_name, H5F_ACC_RDWR_F, file, hdferr);
        isExisting = .FALSE.

        call H5LExists_f(file, dataset_name, isExisting, hdferr)
     
        call H5Fclose_f(file, hdferr)

        dummy = check_return_value(hdferr, "h5_dataset_exists", "h5fclose_f")

        h5_dataset_exists = isExisting

    end function h5_dataset_exists

    subroutine h5_dataset_dimension(file_name, dataset_name, dimId, dataset_dimension)
        ! Returns the dataset dimension in a certain direction
        character(*) :: file_name, dataset_name
        integer(HID_T) :: file_id, dataset_id, dataspace_id, rank
        logical :: dummy
        integer(HSIZE_T),allocatable :: dims(:), maxdims(:)
        integer(HSIZE_T) :: dataset_dimension, dimId

        call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr);
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5fopen")

        call h5dopen_f(file_id, dataset_name, dataset_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5dopen")

        call h5dget_space_f(dataset_id, dataspace_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5dget_space")

        call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sget_simple_extent_ndims")

        allocate(dims(rank))
        call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sget_simple_extent_dims")

        call h5dclose_f(dataset_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5dclose")

        call h5sclose_f(dataspace_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sclose")

        call h5fclose_f(file_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5fclose")

        dataset_dimension = dims(dimId)

    end subroutine h5_dataset_dimension

    subroutine h5_dataset_create_chunked(file_name, dataset_name, rank, dims, maxdims, chunk_dims)

        character(*) :: file_name, dataset_name
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id, prop_id
        logical :: dummy
        integer(HSIZE_T) :: dims(:), maxdims(:), chunk_dims(:)

        integer :: rank
        ! real(kind=8) :: data(:,:,:,:)

        call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr);
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5fopen")


        call h5screate_simple_f(rank, dims, dataspace_id, hdferr, maxdims)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5screate_simple_f")

        call h5pcreate_f(H5P_DATASET_CREATE_F, prop_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5p_create_f")        
        call h5pset_chunk_f(prop_id, size(chunk_dims), chunk_dims, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5pset_chunk")   

        call h5dcreate_f(file_id, dataset_name, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr, &
                    &prop_id, H5P_DEFAULT_F, H5P_DEFAULT_F)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5dcreate_f")

        call h5dclose_f(dataset_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5dclose")

        call h5pclose_f(prop_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_create_chunked", "h5pclose")

        call h5sclose_f(dataspace_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sclose")

        call h5fclose_f(file_id, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5fclose")

        
    end subroutine h5_dataset_create_chunked


    ! 3D routines
    subroutine h5_write_3d(file_name, dataset_name, data)

        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:,:,:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id
        logical :: dummy

        include "h5_write.f90"

    end subroutine h5_write_3d
       
    subroutine h5_write_3d_at_step(file_name, dataset_name, extended_dimension_id, &
            step, dims, data)

        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:,:,:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id, extended_dimension_id, memspace_id
        integer(HSIZE_T) :: dims(:)
        integer(HSIZE_T),allocatable :: maxdims(:), &
                dims_old(:), dims_new(:), offset(:)        
        logical :: dummy
        integer :: rank, i, step

        include "h5_write_at_step.f90"

    end subroutine h5_write_3d_at_step

    subroutine h5_extend_3d(file_name, dataset_name, extended_dimension_id, &
            dims_ext, data)

        ! Data must have the same bounds of dims_ext
        ! Routine to insert a 3D array at a location

        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:,:,:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id, extended_dimension_id, memspace_id
        logical :: dummy
        integer(HSIZE_T) :: dims_ext(:)
        integer(HSIZE_T),allocatable :: maxdims(:), &
                dims_old(:), dims_new(:), offset(:)
        integer :: rank

        include "h5_extend.f90"

    end subroutine h5_extend_3d

    ! 2D routines
    subroutine h5_write_2d(file_name, dataset_name, data)

        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:,:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id
        logical :: dummy

        include "h5_write.f90"

    end subroutine h5_write_2d

    subroutine h5_write_2d_at_step(file_name, dataset_name, extended_dimension_id, &
        step, dims, data)

    character(*) :: file_name, dataset_name
    real(kind=8) :: data(:,:)
    integer(HID_T) :: file_id, dataset_id, &
        dataspace_id, extended_dimension_id, memspace_id
    integer(HSIZE_T) :: dims(:)
    integer(HSIZE_T),allocatable :: maxdims(:), &
            dims_old(:), dims_new(:), offset(:)        
    logical :: dummy
    integer :: rank, i, step

    include "h5_write_at_step.f90"

end subroutine h5_write_2d_at_step

    subroutine h5_extend_2d(file_name, dataset_name, extended_dimension_id, &
            dims_ext, data)

        ! Data must have the same bounds of dims_ext
        ! Routine to insert a 2D array at a location

        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:,:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id, extended_dimension_id, memspace_id
        logical :: dummy
        integer(HSIZE_T) :: dims_ext(:)
        integer(HSIZE_T),allocatable :: maxdims(:), &
                dims_old(:), dims_new(:), offset(:)
        integer :: rank

        include "h5_extend.f90"

    end subroutine h5_extend_2d

    ! 1D routines
    subroutine h5_write_1d(file_name, dataset_name, data)
        ! We can also write a scalar as 1D array
        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id
        logical :: dummy

        include "h5_write.f90"

    end subroutine h5_write_1d
       

    subroutine h5_write_1d_at_step(file_name, dataset_name, extended_dimension_id, &
        step, dims, data)

    character(*) :: file_name, dataset_name
    real(kind=8) :: data(:)
    integer(HID_T) :: file_id, dataset_id, &
        dataspace_id, extended_dimension_id, memspace_id
    integer(HSIZE_T) :: dims(:)
    integer(HSIZE_T),allocatable :: maxdims(:), &
            dims_old(:), dims_new(:), offset(:)        
    logical :: dummy
    integer :: rank, i, step

    include "h5_write_at_step.f90"

end subroutine h5_write_1d_at_step

    subroutine h5_extend_1d(file_name, dataset_name, extended_dimension_id, &
            dims_ext, data)

        ! Data must have the same bounds of dims_ext

        character(*) :: file_name, dataset_name
        real(kind=8) :: data(:)
        integer(HID_T) :: file_id, dataset_id, &
            dataspace_id, extended_dimension_id, memspace_id
        logical :: dummy
        integer(HSIZE_T) :: dims_ext(:)
        integer(HSIZE_T),allocatable :: maxdims(:), &
                dims_old(:), dims_new(:), offset(:)
        integer :: rank

        include "h5_extend.f90"

    end subroutine h5_extend_1d    

    function check_return_value(status, calling_function, returning_function)
        character(*) :: calling_function, returning_function
        integer(HID_T) :: status
        logical :: check_return_value

        if (status < 0) then
            print *, "Error using ", calling_function
            print *, "Failed during call to ", returning_function
            print *, "Status is ", status
            check_return_value = .FALSE.
        end if

        check_return_value = .TRUE.

    end function check_return_value

end module HL_HDF5