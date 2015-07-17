Subroutine PreprocessPDampingZones_spk
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
Integer:: i, j, i1, i2, i3, i4, nx, ny, nxd, nyd, nd, HSL_job, k
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
      i1 = TracePosition(tmpx(:),PDampZones_csr(i)%BBox(1),nx)
      i2 = TracePosition(tmpx(:),PDampZones_csr(i)%BBox(2),nx)
   ELSE
      i1=1; i2=1;
   ENDIF
   PDampZones_csr(i)%idx(1) = i1!-GhostGridX
   PDampZones_csr(i)%idx(2) = i2!+GhostGridX
   nxd=i2-i1+1; PDampZones_csr(i)%nx=nxd
   IF (FineGrid%Ny>1) THEN
      i3 = TracePosition(tmpy(:),PDampZones_csr(i)%BBox(3),ny)
      i4 = TracePosition(tmpy(:),PDampZones_csr(i)%BBox(4),ny)
   ELSE
      i3=1; i4=1;
   END IF
   PDampZones_csr(i)%idx(3) = i3!-GhostGridY
   PDampZones_csr(i)%idx(4) = i4!+GhostGridY
   nyd=i3-i4+1; PDampZones_csr(i)%ny=nyd; nd=nxd*nyd
   print *, ' '
   print *, 'Pressure damping zone ',i,' will be applied from x=',tmpx(i1),' to x=',tmpx(i2)
   print *, '                                        and from y=',tmpy(i3),' to y=',tmpy(i4)
   print *, ' '
   !
   ! Compute the gamma function
   !
   Allocate(PDampZones_csr(i)%gamPhi(nx*ny),PDampZones_csr(i)%gamEta(nx*ny))
   PDampZones_csr(i)%gamPhi(1:nx*ny)=zero; PDampZones_csr(i)%gamEta(1:nx*ny)=zero; 
   !
   CAll GammaFunctions(tmpx(i1:i2),nd,1,11,gtmp,3.5d00)
   PDampZones_csr(i)%gamPhi(i1:i2)=two*pi*PDampZones_csr(i)%g0Phi*gtmp(1:nxd)
   PDampZones_csr(i)%gamEta(i1:i2)=two*pi*PDampZones_csr(i)%g0Eta*gtmp(1:nxd)
   ! 
   ! Copy the function along all y-grid lines.
   !
   Do j=2,nyd
      PDampZones_csr(i)%gamPhi((j-1)*nxd+1:j*nxd)=PDampZones_csr(i)%gamPhi(1:nxd)  
      PDampZones_csr(i)%gamEta((j-1)*nxd+1:j*nxd)=PDampZones_csr(i)%gamEta(1:nxd)  
   end Do
   !
   ! If we're damping the velocity, build and LU-factor the 2D Laplacian over the 
   ! damping zone.  Also save the Gradient operator.  
   !
   If (PDampZones_csr(i)%type==0) THEN
      CALL BuildPDampMatrices_spk(nxd,nyd,PDampZones_csr(i),FineGrid,alpha,beta)
      !
      ! Temporary work space
      allocate(tmp(nd), test_rhs(nd), test_lhs(nd), x(nd) )
      idebug=0  ! Flag for testing the operators
      If (idebug==1) THEN
         !Write out the matrix
         write(203,*),'row pointer = ',PDampZones_csr(i)%Lop%irn
         write(203,*),'col = ',PDampZones_csr(i)%Lop%icn
         write(203,*),'val  = ',PDampZones_csr(i)%Lop%val
         !do j=1,PDampZones_csr(i)%Lop%nnz
         !   write(203,*)PDampZones_csr(i)%Lop%irn(j),PDampZones_csr(i)%Lop%icn(j),PDampZones_csr(i)%Lop%val(j)
         !end do
      end If
      
      ! Set up the Harwell workspace and factor the matrix
      !
      !PDampZones_csr(i)%Lop%MAXIS=5*PDampZones_csr(i)%Lop%nnz+14*nd+1
      !Allocate( PDampZones_csr(i)%Lop%COLSCA(nd), PDampZones_csr(i)%Lop%ROWSCA(nd), PDampZones_csr(i)%Lop%RINFO(20),   &
      !     PDampZones_csr(i)%Lop%INFOHSL(20), PDampZones_csr(i)%Lop%IS_HSL(PDampZones_csr(i)%Lop%MAXIS) )
      !PDampZones_csr(i)%Lop%CNTL(1:10) = (/ 0.01_long, 1e-10_long, zero, zero, zero, zero,               &
      !     zero, zero, zero, zero  /)
      !PDampZones_csr(i)%Lop%ICNTL(1:20) = (/ 6, 0, -1, 2, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,               &
      !     0, 0, 0, 0, 0  /)
      ! Analyse and factor the Laplacian matrix
      PRINT *, ' '
      PRINT *, 'PreprocessPDampingZones:  Factoring the Laplacian matrix.'

      ifil = 8 ! declare
      IWK = PDampZones_csr(i)%Lop%nnz + ifil + 1
      droptol = 0.001
      workspaceSize =(nd+3)*(1+2) +1 
      ALLOCATE(PDampZones_csr(i)%Lop%alu(IWK),PDampZones_csr(i)%Lop%jlu(IWK),PDampZones_csr(i)%Lop%ju(nd),PDampZones_csr(i)%Lop%w(workspaceSize),PDampZones_csr(i)%Lop%jw(3*nd))

      CALL ILUT(nd,PDampZones_csr(i)%Lop%val,PDampZones_csr(i)%Lop%icn,PDampZones_csr(i)%Lop%irn,ifil,droptol,PDampZones_csr(i)%Lop%alu,PDampZones_csr(i)%Lop%jlu,PDampZones_csr(i)%Lop%ju,IWK,PDampZones_csr(i)%Lop%w,PDampZones_csr(i)%Lop%jw,ierrSPK)
      
      If (idebug==1) THEN
         !
         ! Test the integration with a sine wave.  
         !
         x=FineGrid%x(PDampZones_csr(i)%idx(1):PDampZones_csr(i)%idx(1)+nd-1,1); 
         magk=two*pi/(x(nd)-x(1)); 
         do j=1,nd
            test_rhs(j) = -magk**2*cos(magk*x(j))
            test_lhs(j) = cos(magk*x(j))
         end do
         ! Test the matrix by applying it directly to the lhs to take the derivative:  
         !tmp=zero
         !tmp(1)=test_rhs(1)  ! Dirichlet point
         !do j=2,nd-1
         !   do k = 1,PDampZones_csr(i)%Lop%irn(3)-PDampZones_csr(i)%Lop%irn(2)
         !           tmp(j) = tmp(j)+PDampZones_csr(i)%Lop%alu(j)*test_lhs(PDampZones_csr(i)%Lop%icn(j))
         !   end do
         !end do
         tmp=test_rhs; 
         tmp(1)=test_lhs(1); ! Dirichlet condition to the left
         tmp(nd)=zero        ! Neumann condition to the right
         CALL LUSOL(nd,tmp,tmp,PDampZones_csr(i)%Lop%alu,PDampZones_csr(i)%Lop%jlu,PDampZones_csr(i)%Lop%ju)
         do j=1,nd
            write(205,*) test_rhs(j), tmp(j), (tmp(j)-test_rhs(j))/magk**2
         end do


         print *, 'stop'
         stop
      !   ! Now a solve with the Laplacian
      !   tmp=test_rhs; 
      !   tmp(1)=test_lhs(1); ! Dirichlet condition to the left
      !   tmp(nd)=zero        ! Neumann condition to the right
      !   JOB = 3  ! Solution
      !   CALL MA41AD(JOB, nd, PDampZones(i)%Lop%nnz, PDampZones(i)%Lop%irn, PDampZones(i)%Lop%icn,            &
      !        PDampZones(i)%Lop%val, tmp, PDampZones(i)%Lop%COLSCA, PDampZones(i)%Lop%ROWSCA,                 &
      !        PDampZones(i)%Lop%KEEP, PDampZones(i)%Lop%IS_HSL, PDampZones(i)%Lop%MAXIS,                      &
      !        PDampZones(i)%Lop%SS, PDampZones(i)%Lop%MAXS, PDampZones(i)%Lop%CNTL, PDampZones(i)%Lop%ICNTL,  &
      !        PDampZones(i)%Lop%INFOHSL, PDampZones(i)%Lop%RINFO)
      !
      !   if (PDampZones(i)%Lop%INFOHSL(1) .LT. 0) then ! Check for problems
      !      print *, 'Problems with MA41, JOB = 3 (Solution), in PreprocessPDampingZones.f90'
      !   end if
      !   do j=1,nd
      !      write(204,*) test_lhs(j), test_rhs(j), tmp(j), (tmp(j)-test_lhs(j))
      !   end do
      END If
      deallocate(tmp, test_rhs, test_lhs)
!      
   END If
END Do
!
! Write out the damping functions for inspection.  
!
Open(fileop(15),file='Pdamp.chk',status='unknown')
write(fileop(15),790)'% Pressure damping coefficients:  x,y,gammaPhi,gammaEta'
Do i=1,NDampZones
   Do j=1,nx
      write(fileop(15),792)tmpx(j),tmpy(1),PDampZones_csr(i)%gamPhi(j),PDampZones_csr(i)%gamEta(j)
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
END Subroutine PreprocessPDampingZones_spk
