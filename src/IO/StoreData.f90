SUBROUTINE StoreData(nx,ny,E,P,FineGrid,nr,formattype)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: nx, ny, i, j, nr, GIDX, formattype
REAL(KIND=long), DIMENSION(nx*ny) :: E,P
TYPE (Level_def) :: FineGrid
CHARACTER(len=30) :: filename,form
WRITE(unit=filename, FMT="(A,I5.5,A)") "EP_",nr,".bin"
print *, '!!!!!!!!! DB3'
print *, 'Format type:',formattype
SELECT CASE (formattype)
	CASE (1)
		! Relatively efficient for large data storage
		! information about records stored
        form="unformatted"
	CASE DEFAULT
		! Most efficient for large data storage
		! no information about records stored
        form="binary"
END SELECT
OPEN (unit=22, file=filename,form=form)
! Size of free surface
WRITE(22) nx,ny
! Grid information for free surface
WRITE(22) FineGrid%x, FineGrid%y
! Computed solution
WRITE(22) E, P
CLOSE(22)
!
! READ IN MATLAB USING:
!
!    switch formattype
!        case 1 % unformatted file
!            fid=fopen('EP00000.bin','r');
!            [nrec,count] = fread(fid,1,'int32');
!            Nx=fread(fid,1,'int32');
!            Ny=fread(fid,1,'int32');
!            [nrec,count] = fread(fid,1,'int32');
!
!            [nrec,count] = fread(fid,1,'int32');
!            X=fread(fid,[Nx Ny],'float64');
!            Y=fread(fid,[Nx Ny],'float64');
!            [nrec,count] = fread(fid,1,'int32');
!
!            [nrec,count] = fread(fid,1,'int32');
!            E=fread(fid,[Nx Ny],'float64');
!            P=fread(fid,[Nx,Ny],'float64');
!            fclose(fid);
!
!        otherwise  % binary file
!            fid=fopen('EP00000.bin','r');
!            Nx=fread(fid,1,'int32');
!            Ny=fread(fid,1,'int32');
!            X=fread(fid,[Nx Ny],'float64');
!            Y=fread(fid,[Nx Ny],'float64');
!            E=fread(fid,[Nx Ny],'float64');
!            P=fread(fid,[Nx,Ny],'float64');
!            fclose(fid);
!    end
!
!
END SUBROUTINE StoreData
