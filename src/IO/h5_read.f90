call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr);
dummy = check_return_value(hdferr, "h5_write", "h5fopen")

call h5dopen_f(file_id, dataset_name, dataset_id, hdferr)
dummy = check_return_value(hdferr, "h5_write", "h5dopen")

call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, data, int(shape(data),HSIZE_T), hdferr)
dummy = check_return_value(hdferr, "h5_read", "h5read")

call h5dclose_f(dataset_id, hdferr)
dummy = check_return_value(hdferr, "h5_write", "h5dclose")

!call h5sclose_f(dataspace_id, hdferr)
!dummy = check_return_value(hdferr, "h5_dataset_dimensions", "h5sclose")

call h5fclose_f(file_id, hdferr)
dummy = check_return_value(hdferr, "h5_write", "h5fclose")