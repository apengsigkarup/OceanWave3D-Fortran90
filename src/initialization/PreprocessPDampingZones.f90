Subroutine PreprocessPDampingZones
!
! Set up for the pressure damping zones to absorb waves at the boundaries.  
! 
! Only implemented for a +x-directed zone at this point.  
!
! By Allan P. Engsig-Karup.
USE Precision
USE GlobalVariables
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
   !hbb   CAll GammaFunctions(tmpx(i1:i2),nd,1,11,gtmp,3.5d00)
   ! I replaced the exponential here with a Gaussian. The choice should be
   ! moved to the input file. 
   CAll GammaFunctions(tmpx(i1:i2),nd,1,12,gtmp,0.25d00)
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
      !
      ! Temporary work space
      allocate(tmp(nd), test_rhs(nd), test_lhs(nd), x(nd) )
      idebug=0  ! Flag for testing the operators
      If (idebug==1) THEN
         !Write out the matrix
         do j=1,PDampZones(i)%Lop%nnz
            write(203,*)PDampZones(i)%Lop%irn(j),PDampZones(i)%Lop%icn(j),PDampZones(i)%Lop%val(j)
         end do
      end If
      !
      ! Set up the Harwell workspace and factor the matrix
      !
      PDampZones(i)%Lop%MAXIS=5*PDampZones(i)%Lop%nnz+14*nd+1
      Allocate( PDampZones(i)%Lop%COLSCA(nd), PDampZones(i)%Lop%ROWSCA(nd), PDampZones(i)%Lop%RINFO(20),   &
           PDampZones(i)%Lop%INFOHSL(20), PDampZones(i)%Lop%IS_HSL(PDampZones(i)%Lop%MAXIS) )
      PDampZones(i)%Lop%CNTL(1:10) = (/ 0.01_long, 1e-10_long, zero, zero, zero, zero,               &
           zero, zero, zero, zero  /)
      PDampZones(i)%Lop%ICNTL(1:20) = (/ 6, 0, -1, 2, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,               &
           0, 0, 0, 0, 0  /)
      ! Analyse and factor the Laplacian matrix
      PRINT *, ' '
      PRINT *, 'PreprocessPDampingZones:  Factoring the Laplacian matrix.'
      HSL_JOB = 1 ! Analysis
      CALL MA41ID(PDampZones(i)%Lop%CNTL,PDampZones(i)%Lop%ICNTL,PDampZones(i)%Lop%KEEP) ! Sets default values
      tmp=one
      CALL MA41AD(HSL_JOB, nd, PDampZones(i)%Lop%nnz, PDampZones(i)%Lop%irn, PDampZones(i)%Lop%icn,            &
           PDampZones(i)%Lop%val, tmp, PDampZones(i)%Lop%COLSCA, PDampZones(i)%Lop%ROWSCA,                 &
           PDampZones(i)%Lop%KEEP, PDampZones(i)%Lop%IS_HSL, PDampZones(i)%Lop%MAXIS,                      &
           PDampZones(i)%Lop%SS, PDampZones(i)%Lop%MAXS, PDampZones(i)%Lop%CNTL, PDampZones(i)%Lop%ICNTL,  &
           PDampZones(i)%Lop%INFOHSL, PDampZones(i)%Lop%RINFO )
      if (PDampZones(i)%Lop%INFOHSL(1) .LT. 0) then ! Error
         print *, 'Problems with MA41, HSL_JOB = 1 (Analysis)'
         print *, 'INFOHSL(1) = ', PDampZones(i)%Lop%INFOHSL(1)
         stop
      end if
      if (PDampZones(i)%Lop%INFOHSL(1) .GT. 0) THEN ! Warning
         print *, 'Warning from MA41, HSL_JOB = 1 (Analysis)'
         print *, 'INFOHSL(1) = ', PDampZones(i)%Lop%INFOHSL(1)
      end if
      PDampZones(i)%Lop%MAXS = PDampZones(i)%Lop%INFOHSL(8) ! Minimum size allowable from analysis
 !     print *, 'PreprocessPDampingZones:  MAXS_Laplacian = ', PDampZones(i)%Lop%MAXS
      ALLOCATE ( PDampZones(i)%Lop%SS(PDampZones(i)%Lop%MAXS) )
      !
      HSL_JOB = 2 ! Factorization
      CALL MA41AD(HSL_JOB, nd, PDampZones(i)%Lop%nnz, PDampZones(i)%Lop%irn, PDampZones(i)%Lop%icn,        &
           PDampZones(i)%Lop%val, tmp, PDampZones(i)%Lop%COLSCA, PDampZones(i)%Lop%ROWSCA,                 &
           PDampZones(i)%Lop%KEEP, PDampZones(i)%Lop%IS_HSL, PDampZones(i)%Lop%MAXIS,                      &
           PDampZones(i)%Lop%SS, PDampZones(i)%Lop%MAXS, PDampZones(i)%Lop%CNTL, PDampZones(i)%Lop%ICNTL,  &
           PDampZones(i)%Lop%INFOHSL, PDampZones(i)%Lop%RINFO )
      if (PDampZones(i)%Lop%INFOHSL(1) .LT. 0) then
         print *, 'Problems with MA41, HSL_JOB = 2 (Numerical factorization)'
         stop
      end if
      If (idebug==1) THEN
         !
         ! Test the integration with a sine wave.  
         !
         x=FineGrid%x(PDampZones(i)%idx(1):PDampZones(i)%idx(1)+nd-1,1); magk=two*pi/(x(nd)-x(1)); 
         do j=1,nd
            test_rhs(j) = -magk**2*cos(magk*x(j))
            test_lhs(j) = cos(magk*x(j))
         end do
         ! Test the matrix by applying it directly to the lhs to take the derivative:  
         tmp=zero
         tmp(1)=test_rhs(1)  ! Dirichlet point
         do j=2,PDampZones(i)%Lop%nnz-1
            tmp(PDampZones(i)%Lop%irn(j))=tmp(PDampZones(i)%Lop%irn(j))  &
                 +PDampZones(i)%Lop%val(j)*test_lhs(PDampZones(i)%Lop%icn(j))
         end do
         tmp(nd)=test_rhs(nd) ! Neumann point
         do j=1,nd
            write(205,*) test_rhs(j), tmp(j), (tmp(j)-test_rhs(j))/magk**2
         end do
         ! Now a solve with the Laplacian
         tmp=test_rhs; 
         tmp(1)=test_lhs(1); ! Dirichlet condition to the left
         tmp(nd)=zero        ! Neumann condition to the right
         JOB = 3  ! Solution
         CALL MA41AD(JOB, nd, PDampZones(i)%Lop%nnz, PDampZones(i)%Lop%irn, PDampZones(i)%Lop%icn,            &
              PDampZones(i)%Lop%val, tmp, PDampZones(i)%Lop%COLSCA, PDampZones(i)%Lop%ROWSCA,                 &
              PDampZones(i)%Lop%KEEP, PDampZones(i)%Lop%IS_HSL, PDampZones(i)%Lop%MAXIS,                      &
              PDampZones(i)%Lop%SS, PDampZones(i)%Lop%MAXS, PDampZones(i)%Lop%CNTL, PDampZones(i)%Lop%ICNTL,  &
              PDampZones(i)%Lop%INFOHSL, PDampZones(i)%Lop%RINFO)
      
         if (PDampZones(i)%Lop%INFOHSL(1) .LT. 0) then ! Check for problems
            print *, 'Problems with MA41, JOB = 3 (Solution), in PreprocessPDampingZones.f90'
         end if
         do j=1,nd
            write(204,*) test_lhs(j), test_rhs(j), tmp(j), (tmp(j)-test_lhs(j))
         end do
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
