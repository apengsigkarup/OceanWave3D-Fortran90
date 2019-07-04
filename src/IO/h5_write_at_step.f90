if (allocated(offset)) deallocate(offset)
if (allocated(dims_new)) deallocate(dims_new)
if (allocated(dims_old)) deallocate(dims_old)
if (allocated(maxdims)) deallocate(maxdims)

call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr);
dummy = check_return_value(hdferr, "h5_extend", "h5fopen")

call h5dopen_f(file_id, dataset_name, dataset_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5dopen")

call h5dget_space_f(dataset_id, dataspace_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5dget_space")   

! find the old dimensions
! allocate(dims_old, mold=dims)
call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5sget_simple_extent_ndims")

allocate(dims_old(rank))
allocate(maxdims(rank))
call h5sget_simple_extent_dims_f(dataspace_id, dims_old, maxdims, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5sget_simple_extent_dims")

allocate(offset(rank))
offset= 0 ; offset(extended_dimension_id) = step

call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, dims, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5sselect_hyperslab_f")

! print*, 'offset ',offset
! print*, 'dims_ext', dims_ext

! Create a simple memspace where to store the new data
call h5screate_simple_f(rank, dims, memspace_id, hdferr, maxdims)
dummy = check_return_value(hdferr, "h5_extend", "h5screate_simple_f")

! Write the new data to the file
call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr, memspace_id, dataspace_id)
dummy = check_return_value(hdferr, "h5_extend", "h5dwrite")        

! close resources
call h5dclose_f(dataset_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5dclose")
call h5sclose_f(dataspace_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5sclose_dataspace")
call h5sclose_f(memspace_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5sclose_memspace")
call h5fclose_f(file_id, hdferr)
dummy = check_return_value(hdferr, "h5_extend", "h5fclose")    

