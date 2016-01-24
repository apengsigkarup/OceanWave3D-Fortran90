SUBROUTINE PreprocessRelaxationZones
!
! Prepares relaxation zones
!
! FIXME: Take into account ghost layers
! GD: FIXME for curvilinear...
!
! By Allan P. Engsig-Karup.
USE Precision
USE GlobalVariables
IMPLICIT NONE
INTEGER, DIMENSION(FineGrid%Nx,FineGrid%Ny) :: L
INTEGER :: i, j, k, i1, i2, i3, i4
REAL(KIND=long), DIMENSION(FineGrid%Nx+2*GhostGridX) :: tmpx
REAL(KIND=long), DIMENSION(FineGrid%Ny+2*GhostGridY) :: tmpy
! GD: Initialization (necessary in 2D cases) + ghosts...
i1=1; i2=1; i3=1; i4=1
! Global numbering
DO i = 1, FineGrid%Nx
	DO j = 1, FineGrid%Ny
		L(i,j) = (j-1)*FineGrid%Ny+FineGrid%Nx
	END DO
END DO
WRITE(*,*) 'Relaxation zone setup:'
WRITE(fileop(1),*) 'Relaxation zone setup:'
DO i = 1, relaxNo
	RelaxZones(i)%dir = SIGN(1,RelaxZones(i)%ftype)
	! Define indexmaps
	RelaxZones(i)%idx = 1
    IF (curvilinearONOFF==1) THEN ! curvilinear coordinates... FIXME for make it work whatever transformation
      write(*,*) 'This has to be done...'
      !stop
      ! In the curvilinear case, BBox have to be the indexes of the corresponding relaxation zones
      ! FIXME: use Lx and Ly? Depends how x and y are defined...
      DO j=1,FineGrid%Nx+2*GhostGridX
      	tmpx(j) = REAL(j-2,long)
      ENDDO
      DO j=1,FineGrid%Ny+2*GhostGridY
      	tmpy(j) = REAL(j-2,long)
      ENDDO
    ELSE
      tmpx(:) = FineGrid%x(:,1)
      tmpy(:) = FineGrid%y(1,:)
      !IF (FineGrid%Nx>1) THEN
      !  i1 = TracePosition(FineGrid%x(:,1),RelaxZones(i)%BBox(1),FineGrid%Nx+2*GhostGridX)
      !  i2 = TracePosition(FineGrid%x(:,1),RelaxZones(i)%BBox(2),FineGrid%Nx+2*GhostGridX)
      !  RelaxZones(i)%idx(1) = i1!-GhostGridX
      !  RelaxZones(i)%idx(2) = i2!+GhostGridX
      !ENDIF
      !IF (FineGrid%Ny>1) THEN
      !  i3 = TracePosition(FineGrid%y(1,:),RelaxZones(i)%BBox(3),FineGrid%Ny+2*GhostGridY)
      !  i4 = TracePosition(FineGrid%y(1,:),RelaxZones(i)%BBox(4),FineGrid%Ny+2*GhostGridY)
      !  RelaxZones(i)%idx(3) = i3!-GhostGridY
      !  RelaxZones(i)%idx(4) = i4!+GhostGridY
      !ENDIF
    ENDIF
    IF (FineGrid%Nx>1) THEN
      i1 = TracePosition(tmpx(:),RelaxZones(i)%BBox(1),FineGrid%Nx+2*GhostGridX)
      i2 = TracePosition(tmpx(:),RelaxZones(i)%BBox(2),FineGrid%Nx+2*GhostGridX)
      RelaxZones(i)%idx(1) = i1!-GhostGridX
      RelaxZones(i)%idx(2) = i2!+GhostGridX
    ENDIF
    IF (FineGrid%Ny>1) THEN
      i3 = TracePosition(tmpy(:),RelaxZones(i)%BBox(3),FineGrid%Ny+2*GhostGridY)
      i4 = TracePosition(tmpy(:),RelaxZones(i)%BBox(4),FineGrid%Ny+2*GhostGridY)
      RelaxZones(i)%idx(3) = i3!-GhostGridY
      RelaxZones(i)%idx(4) = i4!+GhostGridY
    ENDIF
	! Determine relaxation function (1D)
	IF (RelaxZones(i)%XorY=='X') THEN
	    ALLOCATE( RelaxZones(i)%gam(i2-i1+1) )
        CALL GammaFunctions(tmpx(i1:i2),i2-i1+1,RelaxZones(i)%dir,RelaxZones(i)%ftype,RelaxZones(i)%gam,&
        	 RelaxZones(i)%param)
        !CALL GammaFunctions(FineGrid%x(i1:i2,1),i2-i1+1,RelaxZones(i)%dir,RelaxZones(i)%ftype,RelaxZones(i)%gam,&
        !	 RelaxZones(i)%param)
	ELSE IF (RelaxZones(i)%XorY=='Y') THEN
	    ALLOCATE( RelaxZones(i)%gam(i4-i3+1) )
        CALL GammaFunctions(tmpy(i3:i4),i4-i3+1,RelaxZones(i)%dir,RelaxZones(i)%ftype,RelaxZones(i)%gam,&
        	 RelaxZones(i)%param)
        !CALL GammaFunctions(FineGrid%y(1,i3:i4),i4-i3+1,RelaxZones(i)%dir,RelaxZones(i)%ftype,RelaxZones(i)%gam,&
        !	 RelaxZones(i)%param)
	ELSE
	   WRITE (*,*) 'Error: Define coordinate direction for ',i,' th relaxation zone. Should be "X" or "Y".'
	   STOP
	ENDIF
	IF (RelaxZones(i)%XorYgen=='X') THEN
    	ALLOCATE( RelaxZones(i)%Ea(i2-i1+1), RelaxZones(i)%Pa(i2-i1+1) )
        !ALLOCATE( RelaxZones(i)%Ea(RelaxZones(i)%idx(2)-RelaxZones(i)%idx(1)+1), RelaxZones(i)%Pa(RelaxZones(i)%idx(2)-&
		!RelaxZones(i)%idx(1)+1) )
	ELSE IF (RelaxZones(i)%XorYgen=='Y') THEN
    	ALLOCATE( RelaxZones(i)%Ea(i4-i3+1), RelaxZones(i)%Pa(i4-i3+1) )
        !ALLOCATE( RelaxZones(i)%Ea(RelaxZones(i)%idx(4)-RelaxZones(i)%idx(3)+1), RelaxZones(i)%Pa(RelaxZones(i)%idx(4)-&
		!RelaxZones(i)%idx(3)+1) )
	ELSE
	   WRITE (*,*) 'Error: Define coordinate direction for ',i,' th relaxation zone. Should be "X" or "Y".'
	ENDIF
    ! OUTPUT SETUP
	WRITE(*,789) '  Zone ',i,': idx=(',RelaxZones(i)%idx(1),',',RelaxZones(i)%idx(2),',',RelaxZones(i)%idx(3),',',&
	    RelaxZones(i)%idx(4),')'
	WRITE(fileop(1),789) '  Zone ',i,': idx=(',RelaxZones(i)%idx(1),',',RelaxZones(i)%idx(2),',',RelaxZones(i)%idx(3),',',&
	    RelaxZones(i)%idx(4),')'
   	789 FORMAT(A,I2,A,I6,A,I6,A,I6,A,I6,A)
    !stop
END DO
!
! Write out the relaxation functions for inspection.  
!
Open(fileop(13),file='relax.chk',status='unknown')
Do i=1,relaxNo
IF (RelaxZones(i)%XorY=='X') THEN
   k=0
   Do j=RelaxZones(i)%idx(1),RelaxZones(i)%idx(2)
      k=k+1
      write(fileop(13),792)j,RelaxZones(i)%idx(3),RelaxZones(i)%gam(k)
   END Do
ELSE
   k=0
   Do j=RelaxZones(i)%idx(3),RelaxZones(i)%idx(4)
      k=k+1
      write(fileop(13),792)RelaxZones(i)%idx(1),j,RelaxZones(i)%gam(k)
   END Do
END IF
END Do
792 FORMAT(2i8,1E16.6)
Close(fileop(13))

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
END SUBROUTINE PreprocessRelaxationZones
