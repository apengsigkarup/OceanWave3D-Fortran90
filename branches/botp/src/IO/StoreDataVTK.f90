SUBROUTINE StoreDataVTK(nx,ny,E,P,FineGrid,nr,formattype)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
USE LIB_VTK_IO
IMPLICIT NONE
INTEGER(4) :: nx, ny, i, j, nr, GIDX, formattype
REAL(KIND=long), DIMENSION(nx*ny) :: E,P
REAL(8) :: X(nx*ny), Y(nx*ny), Z(nx*ny)
TYPE (Level_def) :: FineGrid
CHARACTER(len=30) :: filename,form
CHARACTER(len=1) :: varname
INTEGER(4) :: E_IO, nz , ntotal, NC_NN
WRITE(unit=filename, FMT="(A,I5.5,A)") "E_",nr,".vtk"
! E_IO = VTK_INI('Binary',filename,'Free surface Field','STRUCTURED_GRID')
X = zero !FineGrid%x
Y = zero !FineGrid%y
Z = zero
nz=1
ntotal = nx*ny
!E_IO = VTK_GEO(nx,ny,nz,ntotal,X,Y,Z) 
!E_IO = VTK_VAR(NC_NN = ntotal, varname='E',var=E)
!E_IO = VTK_END()
END SUBROUTINE StoreDataVTK
