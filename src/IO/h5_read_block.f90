! call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr);
! dummy = check_return_value(hdferr, "h5_read_block", "h5fopen")

! call h5dopen_f(file_id, dataset_name, dataset_id, hdferr)
! dummy = check_return_value(hdferr, "h5_read_block", "h5dopen")

! call h5dget_space_f(dataset_id, dataspace_id, hdferr)
! dummy = check_return_value(hdferr, "h5_read_block", "h5dget_space")

! call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
! dummy = check_return_value(hdferr, "h5_read_block", "h5sget_simple_extent_ndims")

! ! Allocate the array to hold the dims
! allocate(dimsr(rank), maxdims(rank))

! call h5sget_simple_extent_dims_f(dataspace_id, dimsr, maxdims, hdferr)
! dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sget_simple_extent_dims")

! call H5Dget_create_plist_f(dataset_id, creation_id, hdferr)
! dummy = check_return_value(hdferr, "h5_read_block", "H5Dget_create_plist")

! ! call h5pget_layout_f (creation_id, layout, hdferr)
! ! dummy = check_return_value(hdferr, "h5_read_block", "h5pget_layout")

! ! if (H5D_CHUNKED == layout) then
! ! ... (see OW3D-CUDA code)
! ! end if

! allocate(h5count(rank), h5stride(rank), h5block(rank), h5start(rank))

! do i=1,rank
!     h5stride(i) = 1
!     h5count(i) = 1
!     h5block(i) = dimsr(i)
!     h5start(i) = block_offset(i)
! end do

! h5block(rank) = 1 ! We want to read one timestep at a time.
! allocate(data(dimsr(1), dimsr(2), dimsr(3), dimsr(4)))

! call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, int(shape(data),HSIZE_T), hdferr)
! dummy = check_return_value(hdferr, "h5_read_block", "h5read")

! call h5dclose_f(dataset_id, hdferr)
! dummy = check_return_value(hdferr, "h5_read_block", "h5dclose")

! call h5fclose_f(file_id, hdferr)
! dummy = check_return_value(hdferr, "h5_read_block", "h5fclose")

        ! Check the dimension
        ! call h5_dataset_dimension_all(file_name, dataset_name, dataset_dimension)
        ! allocate(data(dataset_dimension))

        call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr);
        dummy = check_return_value(hdferr, "h5_read_block", "h5fopen")

        call h5dopen_f(file_id, dataset_name, dataset_id, hdferr)
        dummy = check_return_value(hdferr, "h5_read_block", "h5dopen")

        call h5dget_space_f(dataset_id, dataspace_id, hdferr)
        dummy = check_return_value(hdferr, "h5_read_block", "h5dget_space")

        call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
        dummy = check_return_value(hdferr, "h5_read_block", "h5sget_simple_extent_ndims")

        ! Allocate the array to hold the dims
        allocate(dimsr(rank), maxdims(rank))

        call h5sget_simple_extent_dims_f(dataspace_id, dimsr, maxdims, hdferr)
        dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sget_simple_extent_dims")

        call H5Dget_create_plist_f(dataset_id, creation_id, hdferr)
        dummy = check_return_value(hdferr, "h5_read_block", "H5Dget_create_plist")

        allocate(h5count(rank), h5stride(rank), h5block(rank), h5start(rank))

        do i=1,rank
            h5stride(i) = 1
            h5count(i) = 1
            h5block(i) = dimsr(i)
            h5start(i) = block_offset(i)
        end do

        h5block(rank) = 1 ! We want to read one timestep at a time.
        ! allocate buffer for reading the file
        ! allocate(data(h5block(1), h5block(2), h5block(3), h5block(4)))

        ! Create a memory space to read to
        call h5screate_simple_f(rank, h5block, memspace, hdferr)
        dummy = check_return_value(hdferr, "h5_read_block", "h5screate_simple_f")

        ! Select an hyoerslab from the file with dimension h5block and 
        ! stride, count, and start
        call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, h5start, h5count, hdferr, h5stride, h5block)
        dummy = check_return_value(hdferr, "h5_read_block", "h5sselect_hyperslab_f")

        ! Read to a memory space
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, int(shape(data),HSIZE_T), hdferr, memspace, dataspace_id)
        dummy = check_return_value(hdferr, "h5_read_block", "h5read")

        call h5dclose_f(dataset_id, hdferr)
        dummy = check_return_value(hdferr, "h5_read_block", "h5dclose")

        call h5fclose_f(file_id, hdferr)
        dummy = check_return_value(hdferr, "h5_read_block", "h5fclose")

        
        deallocate(dimsr, maxdims, h5count, h5stride, h5block, h5start)
