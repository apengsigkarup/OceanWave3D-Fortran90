SUBROUTINE StoreRealArray3D(Array,r1,r2,r3,filename,formattype)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: r1, r2, r3, formattype
REAL(KIND=long), DIMENSION(r1,r2,r3) :: Array
CHARACTER(len=40) :: filename,access,form
SELECT CASE (formattype)
	CASE (1)
		! Relatively efficient for large data storage
		! information about records stored
        form="unformatted"
		OPEN (unit=22, file=filename,form=form)
	CASE DEFAULT
		! Most efficient for large data storage
		! no information about records stored
        form="binary"
		OPEN (unit=22, file=filename,form=form)
END SELECT
WRITE(22) r1, r2, r3, Array
CLOSE(22)
END SUBROUTINE StoreRealArray3D
