SUBROUTINE GatherData(nx,ny,nz,err,gridtypeZ,r,kh,formattype)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: nx, ny, nz, gridtypeZ, r, formattype
REAL(KIND=long) :: kh, err
CHARACTER(len=20) :: filename, form
WRITE(unit=filename, FMT="(A)") "convergresults.bin"
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
OPEN (unit=22, file=filename,form=form, position="append")
! Size of free surface
WRITE(22) nx,ny,nz,gridtypeZ,r,kh,err
CLOSE(22)
END SUBROUTINE GatherData
