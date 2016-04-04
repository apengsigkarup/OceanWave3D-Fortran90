Subroutine PreprocessPDampingZones
!
! Set up for the pressure damping zones to absorb waves at the boundaries.  
! 
! Only implemented for a +x-directed zone at this point.  
!
! By Allan P. Engsig-Karup.
USE Precision
USE GlobalVariables
USE pdamp_CSR
Implicit none
Integer:: i, j, i1, i2, i3, i4, nx, ny, nxd, nyd, nd, HSL_job
REAL(kind=long) :: x0, x1, y0, y1, magk, idebug
REAL(KIND=long), DIMENSION(FineGrid%Nx+2*GhostGridX) :: tmpx
REAL(KIND=long), DIMENSION(FineGrid%Ny+2*GhostGridY) :: tmpy
Real(kind=long) :: gtmp((FineGrid%Nx+2*GhostGridX)*(FineGrid%Ny+2*GhostGridY))
REAL(kind=long), allocatable :: tmp(:), test_rhs(:), test_lhs(:), x(:)
!
IF (curvilinearONOFF==1) THEN 
   print *, 'PreprocessPDampingZones:  Not implemented for curvilinear yet.'
   stop
END IF

tmpx(:) = FineGrid%x(:,1)
tmpy(:) = FineGrid%y(1,:)
nx=FineGrid%Nx+2*GhostGridX; ny=FineGrid%Ny+2*GhostGridY
Do i=1,NDampZones
   IF (FineGrid%Nx>1) THEN
      i1 = TracePosition(tmpx(:),PDampZones(i)%BBox(1),nx)
      i2 = TracePosition(tmpx(:),PDampZones(i)%BBox(2),nx)
   ELSE
      i1=1; i2=1;
   ENDIF
   PDampZones(i)%idx(1) = i1!-GhostGridX
   PDampZones(i)%idx(2) = i2!+GhostGridX
   nxd=i2-i1+1; PDampZones(i)%nx=nxd
   IF (FineGrid%Ny>1) THEN
      i3 = TracePosition(tmpy(:),PDampZones(i)%BBox(3),ny)
      i4 = TracePosition(tmpy(:),PDampZones(i)%BBox(4),ny)
   ELSE
      i3=1; i4=1;
   END IF
   PDampZones(i)%idx(3) = i3!-GhostGridY
   PDampZones(i)%idx(4) = i4!+GhostGridY
   nyd=i3-i4+1; PDampZones(i)%ny=nyd; nd=nxd*nyd
   print *, ' '
   print *, 'Pressure damping zone ',i,' will be applied from x=',tmpx(i1),' to x=',tmpx(i2)
   print *, '                                        and from y=',tmpy(i3),' to y=',tmpy(i4)
   print *, ' '
   !
   ! Compute the gamma function
   !
   Allocate(PDampZones(i)%gamPhi(nx*ny),PDampZones(i)%gamEta(nx*ny))
   PDampZones(i)%gamPhi(1:nx*ny)=zero; PDampZones(i)%gamEta(1:nx*ny)=zero; 
   !
   CAll GammaFunctions(tmpx(i1:i2),nd,1,11,gtmp,3.5d00)
   PDampZones(i)%gamPhi(i1:i2)=two*pi*PDampZones(i)%g0Phi*gtmp(1:nxd)
   PDampZones(i)%gamEta(i1:i2)=two*pi*PDampZones(i)%g0Eta*gtmp(1:nxd)
   ! 
   ! Copy the function along all y-grid lines.
   !
   Do j=2,nyd
      PDampZones(i)%gamPhi((j-1)*nxd+1:j*nxd)=PDampZones(i)%gamPhi(1:nxd)  
      PDampZones(i)%gamEta((j-1)*nxd+1:j*nxd)=PDampZones(i)%gamEta(1:nxd)  
   end Do
   !
   ! If we're damping the velocity, build and LU-factor the 2D Laplacian over the 
   ! damping zone.  Also save the Gradient operator.  
   !
   If (PDampZones(i)%type==0) THEN
      CALL BuildPDampMatrices(nxd,nyd,PDampZones(i),FineGrid,alpha,beta)

      PRINT *, ' '
      PRINT *, 'PreprocessPDampingZones:  Factoring the Laplacian matrix.'
      
      !      
      ! 
      ! Create the LU matrix using SPARSKIT
      !
      ifil = 2*alpha ! declare
      IWK = PDampZones(i)%Lop%nnz + ifil + 1
      droptol = 0.000
      workspaceSize =(nd+3)*(1+2) +1 
      ALLOCATE(PDampZones(i)%Lop%alu(IWK),PDampZones(i)%Lop%jlu(IWK),PDampZones(i)%Lop%ju(nd),PDampZones(i)%Lop%w(workspaceSize),PDampZones(i)%Lop%jw(3*nd))

      CALL ILUT(nd,PDampZones(i)%Lop%val,PDampZones(i)%Lop%icn,PDampZones(i)%Lop%irn,ifil,droptol,PDampZones(i)%Lop%alu,PDampZones(i)%Lop%jlu,PDampZones(i)%Lop%ju,IWK,PDampZones(i)%Lop%w,PDampZones(i)%Lop%jw,ierrSPK)
   
   
   
   END IF
END DO
!
! Write out the damping functions for inspection.  
!
Open(fileop(15),file='Pdamp.chk',status='unknown')
write(fileop(15),790)'% Pressure damping coefficients:  x,y,gammaPhi,gammaEta'
Do i=1,NDampZones
   Do j=1,nx
      write(fileop(15),792)tmpx(j),tmpy(1),PDampZones(i)%gamPhi(j),PDampZones(i)%gamEta(j)
   END Do
END Do
790 FORMAT(A)
792 FORMAT(4E16.6)
close(fileop(15))


CONTAINS

INTEGER FUNCTION TracePosition(x,x0,n) RESULT(ipos)
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: n, i
REAL(KIND=long) :: x(n), x0
DO i = 1 , n
    ipos = i
	IF (x0<=x(i)) THEN
		RETURN
	ENDIF
END DO
END FUNCTION TracePosition
END Subroutine PreprocessPDampingZones
