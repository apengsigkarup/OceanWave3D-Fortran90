!> 
!! Carry out initial preprocessing for 3D Multigrid solver
!!
!! Allan P. Engsig-Karup, 23 Aug 2007.
!<
SUBROUTINE MGPreProcess (FineGrid,GhostGridX,GhostGridY,GhostGridZ,MGCoarseningStrategy,alpha,beta,gamma,Precond,&
	MGmaxgrids,CurvilinearONOFF)
USE Precision
USE Constants
USE MGLevels
USE GlobalVariables, ONLY: filename, formattype
IMPLICIT NONE
INTEGER :: k, GhostGridX, GhostGridY, GhostGridZ, alpha, beta, gamma, Precond, MGmaxgrids
INTEGER :: xcoarsen, ycoarsen, zcoarsen, MGCoarseningStrategy, CurvilinearONOFF
INTEGER :: NFx, NFy, NFz, XFAC, YFAC, ZFAC
TYPE (Level_def) :: FineGrid
REAL(KIND=long), DIMENSION(:,:), POINTER :: x_fine, y_fine
REAL(KIND=long), DIMENSION(:),   POINTER :: z_fine, dx, dy
INTEGER :: Nxf, Nyf, Nzf, Nxg, Nyg, Nzg, STAT, ii, jj, kk, Count, Gidx

! number of points in each direction on fine grid (without ghost points included)
Nxf = FineGrid%Nx !SIZE(FineGrid%x,2)
Nyf = FineGrid%Ny !SIZE(FineGrid%y,3)
Nzf = FineGrid%Nz !SIZE(FineGrid%z,1)-GhostGridZ ! Subtract GhostGridZ from length

Nxg = Nxf + 2*GhostGridX
Nyg = Nyf + 2*GhostGridY
Nzg = Nzf + GhostGridZ

! Determine number of grid levels
MG_N_levels = NoGrids(Nxf,Nyf,Nzf)
IF (MG_N_levels > MGmaxgrids .AND. MGmaxgrids>0) THEN
	MG_N_levels = MGmaxgrids
	PRINT*,'  WARNING: The number of grids is user-defined to be: ', MG_N_levels,' grids'
	PRINT*,'           although algorithm suggest that more levels can be used.'
ENDIF

! Allocate array in memory
ALLOCATE ( arrLevels( MG_N_levels),STAT=STAT )
CALL CheckError(STAT,3)

! FIXME: ghost points not included here
ALLOCATE(x_fine(Nxg,Nyg), y_fine(Nxg,Nyg), z_fine(Nzg),STAT=STAT)
CALL CheckError(STAT,2)
x_fine = FineGrid%x
y_fine = FineGrid%y
z_fine = FineGrid%z

!print*,'x_fine=',x_fine
!read*

! Determine grids using combined standard coarsening and semi-coarsening
! FIXME: should be based on spatial resolution parameters, e.g. dx, dy, dz 
! rather than Nx, Ny, Nz as employed here
!
! COARSENING STRATEGY: COARSEN IN DIRECTION WITH MOST POINTS
arrLevels(MG_N_levels)%Nx = Nxf;
arrLevels(MG_N_levels)%Ny = Nyf;
arrLevels(MG_N_levels)%Nz = Nzf;
ALLOCATE ( arrLevels(MG_N_levels)%x(Nxg,Nyg),STAT=STAT )
CALL CheckError(STAT,2)
arrLevels(MG_N_levels)%x = x_fine
! FIXME: too much memory allocated below (APEK)
IF (Nxf>1) THEN
	ALLOCATE ( arrLevels(MG_N_levels)%hx(Nxg,Nyg ) ,STAT=STAT)
	CALL CheckError(STAT,2)
	ALLOCATE ( arrLevels(MG_N_levels)%hxx(Nxg,Nyg) ,STAT=STAT)
	CALL CheckError(STAT,2)
	arrLevels(MG_N_levels)%hx  = FineGrid%hx
	arrLevels(MG_N_levels)%hxx = FineGrid%hxx
ENDIF
ALLOCATE ( arrLevels(MG_N_levels)%y(Nxg,Nyg) ,STAT=STAT)
CALL CheckError(STAT,2)
arrLevels(MG_N_levels)%y = y_fine
IF (Nyf>1) THEN
	ALLOCATE ( arrLevels(MG_N_levels)%hy(Nxg,Nyg ) ,STAT=STAT)
	CALL CheckError(STAT,2)
	ALLOCATE ( arrLevels(MG_N_levels)%hyy(Nxg,Nyg) ,STAT=STAT)
	CALL CheckError(STAT,2)
	arrLevels(MG_N_levels)%hy  = FineGrid%hy
	arrLevels(MG_N_levels)%hyy = FineGrid%hyy
ENDIF
ALLOCATE ( arrLevels(MG_N_levels)%h(Nxg,Nyg  ) )
arrLevels(MG_N_levels)%h   = FineGrid%h
ALLOCATE ( arrLevels(MG_N_levels)%z(Nzg),STAT=STAT )
CALL CheckError(STAT,2)
arrLevels(MG_N_levels)%z = z_fine

WRITE (*,*)   '  Multigrid preproccesing...'
WRITE (*,100)   '   A total of ',MG_N_levels,' grids will be used:'
WRITE (*,100) '     (Nx,Ny,Nz) = (',arrLevels(MG_N_levels)%Nx,',',arrLevels(MG_N_levels)%Ny,',',arrLevels(MG_N_levels)%Nz,')'
100 FORMAT (A,I5,A,I5,A,I5,A)

! SETUP COARSENING STATEGY & SOME GRID INFORMATION
DO k = MG_N_levels-1, 1, -1
	! points at fine grid level
	NFx = arrLevels(k+1)%Nx
	NFy = arrLevels(k+1)%Ny
	NFz = arrLevels(k+1)%Nz
	xcoarsen = 0
	ycoarsen = 0
	zcoarsen = 0
	IF (MGCoarseningStrategy==1) THEN ! Ole's coarsening strategy
		! Coarsen by comparing Nx with both Ny,Nz and
		!                      Ny with both Nx,Nz
		IF (NFx == NFy .AND. NFx == NFz) THEN
			! 3D standard coarsening
			xcoarsen = 1; ycoarsen = 1;	zcoarsen = 1
		ELSE IF (NFx > NFy .AND. NFx > NFz) THEN
			! 3D x-semi-coarsening
			xcoarsen = 1; 
		ELSE IF (NFy > NFx .AND. NFy > NFz) THEN
			! 3D y-semi-coarsening
			ycoarsen = 1; 
		ELSE IF (NFz > NFx .AND. NFz > NFy) THEN
			! 3D z-semi-coarsening
			zcoarsen = 1; 
		ELSE IF (NFx == NFy .AND. NFx > NFz) THEN
			! 3D xy-semi-coarsening
			xcoarsen = 1; ycoarsen = 1
		ELSE IF (NFx == NFz .AND. NFx > NFy) THEN
			! 3D xz-semi-coarsening
			xcoarsen = 1; zcoarsen = 1
		ELSE IF (NFy == NFz .AND. NFy > NFx) THEN
			! 3D yz-semi-coarsening
			ycoarsen = 1; zcoarsen = 1
		ENDIF
	ELSE IF (MGCoarseningStrategy==2) THEN ! Allan's second coarsening strategy
		! Coarsen by comparing Nx with both Ny,Nz and
		!                      Ny with both Nx,Nz
		IF (NFx>3) THEN
		    xcoarsen = 1
		ENDIF 
		IF (NFy>3) THEN
		    ycoarsen = 1
		ENDIF 
		IF (NFz>3) THEN
		    zcoarsen = 1
		ENDIF 
	ELSE ! STANDARD STRATEGY: Allan's coarsening strategy
		! Coarsen by comparing Nx with Nz and
		!                      Ny with Nz
		! FIXME: This strategy should be based on spatial resolution rather than on total points

		! FIRST CHECK IF THE SETUP IS ONLY 2D
		IF (NFx == 1) THEN
			IF (NFy == NFz) THEN
				! 2D standard coarsening
				ycoarsen = 1; zcoarsen = 1;
			ELSE ! "1D" SEMI-COARSENING IN Y-DIRECTION
				ycoarsen = 1;
			END IF
		ELSE IF (NFy == 1) THEN
			IF (NFx == NFz) THEN
				! 2D standard coarsening
				xcoarsen = 1; zcoarsen = 1;
			ELSE ! "1D" SEMI-COARSENING IN Y-DIRECTION
				xcoarsen = 1;
			END IF
		ELSE
			IF (NFx == NFy .AND. NFx == NFz) THEN
				! 3D standard coarsening
				xcoarsen = 1; ycoarsen = 1;	zcoarsen = 1
			ELSE ! "2D" SEMI-COARSENING IN XZ- AND YZ-PLANES
				IF (NFx>NFz) THEN
					! x-semi-coarsening
					xcoarsen = 1
				END IF
				IF (NFy>NFz) THEN
					! y-semi-coarsening
					ycoarsen = 1
				END IF
				IF (NFx<NfZ .AND. NFy<NFz) THEN
					! z-semi-coarsening
					zcoarsen = 1
				END IF
			ENDIF
		ENDIF
	ENDIF
	! X-COARSENING
	XFAC = 1
	IF (xcoarsen==1) THEN
		arrLevels(k)%Nx = (arrLevels(k+1)%Nx-1)/2+1
		IF (MOD(arrLevels(k+1)%Nx,2)==0) THEN
		   PRINT *,'Error: cannot coarsen grid in x-direction (multigrid). Nx should be odd-numbered. Nx=',MOD(arrLevels(k)%Nx,2)
		   STOP
		ENDIF
		XFAC = 2
	ELSE
		arrLevels(k)%Nx = arrLevels(k+1)%Nx
	END IF
	! Y-COARSENING
	YFAC = 1
	IF (ycoarsen==1) THEN
		arrLevels(k)%Ny = (arrLevels(k+1)%Ny-1)/2+1
		IF (MOD(arrLevels(k+1)%Ny,2)==0) THEN
		   PRINT * , 'Error: cannot coarsen grid in y-direction (multigrid). Ny should be odd-numbered. Ny=',arrLevels(k)%Ny
		   STOP
		ENDIF
		YFAC = 2
	ELSE
		arrLevels(k)%Ny = arrLevels(k+1)%Ny
	END IF
	! Z-COARSENING
	ZFAC = 1
	IF (zcoarsen==1) THEN
		arrLevels(k)%Nz = (arrLevels(k+1)%Nz-1)/2+1
		IF (MOD(arrLevels(k)%Nz,2)==0) THEN
		   PRINT *,'Error: cannot coarsen grid in z-direction (multigrid). Nz should be odd-numbered. Nz=',MOD(arrLevels(k)%Nz,2)
		   STOP
		ENDIF
		ZFAC = 2
	ELSE
		arrLevels(k)%Nz = arrLevels(k+1)%Nz
	END IF
	! DETERMINE COARSE GRID FROM FINE GRID
	Nxg = arrLevels(k)%Nx + 2*GhostGridX
	Nyg = arrLevels(k)%Ny + 2*GhostGridY
	IF (Nxg>1) THEN
		ALLOCATE ( arrLevels(k)%x(Nxg,Nyg)   )
		ALLOCATE ( arrLevels(k)%hx(Nxg,Nyg)  )
		ALLOCATE ( arrLevels(k)%hxx(Nxg,Nyg) )
		arrLevels(k)%x(1+GhostGridX:Nxg-GhostGridX,1+GhostGridY:Nyg-GhostGridY)   = &
        	arrLevels(k+1)%x(1+GhostGridX:arrLevels(k+1)%Nx+GhostGridX:XFAC,1+GhostGridY:arrLevels(k+1)%Ny+GhostGridY:YFAC)
		arrLevels(k)%hx(1+GhostGridX:Nxg-GhostGridX,1+GhostGridY:Nyg-GhostGridY)  = &
        	arrLevels(k+1)%hx(1+GhostGridX:arrLevels(k+1)%Nx+GhostGridX:XFAC,1+GhostGridY:arrLevels(k+1)%Ny+GhostGridY:YFAC)
		arrLevels(k)%hxx(1+GhostGridX:Nxg-GhostGridX,1+GhostGridY:Nyg-GhostGridY) = &
        	arrLevels(k+1)%hxx(1+GhostGridX:arrLevels(k+1)%Nx+GhostGridX:XFAC,1+GhostGridY:arrLevels(k+1)%Ny+GhostGridY:YFAC)
        IF (GhostGridX==1) THEN
            arrLevels(k)%x(1,1+GhostGridY:Nyg-GhostGridY) = arrLevels(k)%x(2,1+GhostGridY:Nyg-GhostGridY) - &
                                    (arrLevels(k)%x(3,1+GhostGridY:Nyg-GhostGridY) - arrLevels(k)%x(2,1+GhostGridY:Nyg-GhostGridY))
            arrLevels(k)%x(Nxg,1+GhostGridY:Nyg-GhostGridY) = arrLevels(k)%x(Nxg-GhostGridX,1+GhostGridY:Nyg-GhostGridY) + &
                                    (arrLevels(k)%x(Nxg-GhostGridX,1+GhostGridY:Nyg-GhostGridY) - &
                                    arrLevels(k)%x(Nxg-2*GhostGridX,1+GhostGridY:Nyg-GhostGridY))
        END IF
!print*,'TESTER lige x-grid i MGPreProcess...'
!print*,'arrLevels(k)%x=',arrLevels(k)%x
!read*        
		IF (Nyg>1) THEN ! FIXME: should be GhostGridY==1?
			ALLOCATE(dx(arrLevels(k)%Ny))
			dx = arrLevels(k)%x(3,1+GhostGridY:arrLevels(k)%Ny+GhostGridY) - &
            	arrLevels(k)%x(2,1+GhostGridY:arrLevels(k)%Ny+GhostGridY)
			arrLevels(k)%x(1,1+GhostGridY:arrLevels(k)%Ny+GhostGridY) = &
            	arrLevels(k)%x(2,1+GhostGridY:arrLevels(k)%Ny+GhostGridY) - dx
			dx = arrLevels(k)%x(arrLevels(k)%Nx+2*GhostGridX,1+GhostGridY:arrLevels(k)%Ny+GhostGridY) - &
            	arrLevels(k)%x(arrLevels(k)%Nx+GhostGridX,1+GhostGridY:arrLevels(k)%Ny+GhostGridY)
			arrLevels(k)%x(arrLevels(k)%Nx+2*GhostGridX,1+GhostGridY:arrLevels(k)%Ny+GhostGridY) = &
            	arrLevels(k)%x(arrLevels(k)%Nx+GhostGridX,1+GhostGridY:arrLevels(k)%Ny+GhostGridY) + dx
			! set top and bottom
			arrLevels(k)%x(1:arrLevels(k)%Nx+2*GhostGridX,1) = arrLevels(k)%x(1:arrLevels(k)%Nx+2*GhostGridX,2) & 
				- (arrLevels(k)%x(1:arrLevels(k)%Nx+2*GhostGridX,3) - arrLevels(k)%x(1:arrLevels(k)%Nx+2*GhostGridX,2) )
			arrLevels(k)%x(1:arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY) = &
            	arrLevels(k)%x(1:arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY-1) & 
				+ (arrLevels(k)%x(1:arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY-1) - &
                arrLevels(k)%x(1:arrLevels(k)%Nx+2*GhostGridX,arrLevels(k)%Ny+2*GhostGridY-2) )
			DEALLOCATE(dx)
		END IF
	ENDIF
	IF (Nyg>1) THEN
		ALLOCATE ( arrLevels(k)%y(Nxg,Nyg)   )
		ALLOCATE ( arrLevels(k)%hy(Nxg,Nyg)  )
		ALLOCATE ( arrLevels(k)%hyy(Nxg,Nyg) )
		arrLevels(k)%y(1+GhostGridX:Nxg-GhostGridX,1+GhostGridY:Nyg-GhostGridY)   = &
        	arrLevels(k+1)%y(1+GhostGridX:arrLevels(k+1)%Nx+GhostGridX:XFAC,1+GhostGridY:arrLevels(k+1)%Ny+GhostGridY:YFAC)
		arrLevels(k)%hy(1+GhostGridX:Nxg-GhostGridX,1+GhostGridY:Nyg-GhostGridY)  = &
        	arrLevels(k+1)%hy(1+GhostGridX:arrLevels(k+1)%Nx+GhostGridX:XFAC,1+GhostGridY:arrLevels(k+1)%Ny+GhostGridY:YFAC)
		arrLevels(k)%hyy(1+GhostGridX:Nxg-GhostGridX,1+GhostGridY:Nyg-GhostGridY) = &
        	arrLevels(k+1)%hyy(1+GhostGridX:arrLevels(k+1)%Nx+GhostGridX:XFAC,1+GhostGridY:arrLevels(k+1)%Ny+GhostGridY:YFAC)
        IF (GhostGridY==1) THEN
            arrLevels(k)%y(1+GhostGridX:Nxg-GhostGridX,1) = arrLevels(k)%y(1+GhostGridX:Nxg-GhostGridX,2) - &
                                    (arrLevels(k)%y(1+GhostGridX:Nxg-GhostGridX,3) - arrLevels(k)%y(1+GhostGridX:Nxg-GhostGridx,2))
            arrLevels(k)%y(1+GhostGridX:Nxg-GhostGridX,Nyg) = arrLevels(k)%y(1+GhostGridX:Nxg-GhostGridX,Nyg-GhostGridY) + &
                                    (arrLevels(k)%y(1+GhostGridX:Nxg-GhostGridX,Nyg-GhostGridY) - &
                                    arrLevels(k)%y(1+GhostGridX:Nxg-GhostGridx,Nyg-2*GhostGridY))
        END IF
		IF (Nxg>1) THEN ! FIXME: should be GhostGridX==1?
			ALLOCATE(dy(arrLevels(k)%Nx))
			dy = arrLevels(k)%y(1+GhostGridX:arrLevels(k)%Nx+GhostGridX,3) - &
            	arrLevels(k)%y(1+GhostGridX:arrLevels(k)%Nx+GhostGridX,2)
			arrLevels(k)%y(1+GhostGridX:arrLevels(k)%Nx+GhostGridX,1) = &
            	arrLevels(k)%y(1+GhostGridX:arrLevels(k)%Nx+GhostGridX,2) - dy
			dy = arrLevels(k)%y(1+GhostGridX:arrLevels(k)%Nx+GhostGridX,arrLevels(k)%Ny+2*GhostGridY) - &
            	arrLevels(k)%y(1+GhostGridX:arrLevels(k)%Nx+GhostGridX,arrLevels(k)%Ny+GhostGridY)
			arrLevels(k)%y(1+GhostGridX:arrLevels(k)%Nx+GhostGridX,arrLevels(k)%Ny+2*GhostGridY) = &
            	arrLevels(k)%y(1+GhostGridX:arrLevels(k)%Nx+GhostGridX,arrLevels(k)%Ny+GhostGridY) + dy
			! set left and right boundaries
			arrLevels(k)%y(1,1:arrLevels(k)%Ny+2*GhostGridY) = arrLevels(k)%y(2,1:arrLevels(k)%Ny+2*GhostGridY) & 
				- (arrLevels(k)%y(3,1:arrLevels(k)%Ny+2*GhostGridY) - arrLevels(k)%y(2,1:arrLevels(k)%Ny+2*GhostGridY) )
			arrLevels(k)%y(arrLevels(k)%Nx+2*GhostGridX,1:arrLevels(k)%Ny+2*GhostGridY) = &
            	arrLevels(k)%x(arrLevels(k)%Nx+2*GhostGridX-1,1:arrLevels(k)%Ny+2*GhostGridY) & 
				+ (arrLevels(k)%y(arrLevels(k)%Nx+2*GhostGridX-1,1:arrLevels(k)%Ny+2*GhostGridY) - &
                arrLevels(k)%y(arrLevels(k)%Nx+2*GhostGridX-2,1:arrLevels(k)%Ny+2*GhostGridY) )			
			DEALLOCATE(dy)
		END IF
	ENDIF
	ALLOCATE ( arrLevels(k)%h(Nxg,Nyg),STAT=STAT )
	CALL CheckError(STAT,2)

	ALLOCATE ( arrLevels(k)%z(arrLevels(k)%Nz+GhostGridZ), STAT=STAT )
	CALL CheckError(STAT,2)
	IF (GhostGridZ==1) THEN
		arrLevels(k)%z(1)  = -arrLevels(k+1)%z(2+ZFAC)
		arrLevels(k)%z(2:) = arrLevels(k+1)%z(2:arrLevels(k+1)%Nz+GhostGridZ:ZFAC)
	ELSE
		arrLevels(k)%z = arrLevels(k+1)%z(1:arrLevels(k+1)%Nz:ZFAC)
	END IF
	! FIXME: h is not defined on ghost points since it is not needed there
	arrLevels(k)%h(1+GhostGridX:arrLevels(k)%Nx+GhostGridX,1+GhostGridY:arrLevels(k)%Ny+GhostGridY)   = &
    	arrLevels(k+1)%h(1+GhostGridX:arrLevels(k+1)%Nx+GhostGridX:XFAC,1+GhostGridY:arrLevels(k+1)%Ny+GhostGridY:YFAC)

	WRITE (*,100) '     (Nx,Ny,Nz) = (',arrLevels(k)%Nx,',',arrLevels(k)%Ny,',',arrLevels(k)%Nz,')'

END DO

! SETUP LINEAR LOW-ORDER SYSTEM OPERATORS FOR EACH GRID LEVEL
DO k = MG_N_levels, 1, -1

!print*,'GhostGridX=',GhostGridX
!print*,'GhostGridY=',GhostGridY
!print*,'GhostGridZ=',GhostGridZ
!print*,'alpha=',alpha
!print*,'beta=',beta
!print*,'gamma=',gamma
!print*,'CurvilinearONOFF=',CurvilinearONOFF
	CALL PreparePreconditioner(arrLevels(k)%PreconditioningMatrix,arrLevels(k),GhostGridX, GhostGridY, GhostGridZ, &
		alpha, beta, gamma, 1, CurvilinearONOFF)

IF (0==1) THEN
    WRITE(filename,'(a,i4.4,a4)') 'P',k,'.bin'        
    print*,'filename=',filename
    print*,'nnz=',arrLevels(k)%PreconditioningMatrix%nnz
    CALL StoreSparseMatrix(arrLevels(k)%PreconditioningMatrix,filename,formattype)
    print*,'Store preconditioning matrix for level k=',k
!READ*
END IF

	IF (k==1) THEN
		! Prepare for using direct LU-FACTORIZATION at coarsest gridlevel
		! we can only factor on one level with current setup due to global variables in LUFactor subroutine
		CALL FactorPreconditioner(arrLevels(k)%PreconditioningMatrix, &
    		 (arrLevels(k)%Nx+2*GhostGridX)*(arrLevels(k)%Ny+2*GhostGridY)*(arrLevels(k)%Nz+GhostGridZ))
    !         END IF
!             print*,'LU factorization of coarsest ceoffciient matrix in multigrid preconditinoing strategy done...'
!             print*,'k=',k
!print*,'arrLevels(k)%PreconditioningMatrix%val     = ',arrLevels(k)%PreconditioningMatrix%val
!print*,'arrLevels(k)%PreconditioningMatrix%col_ind = ',arrLevels(k)%PreconditioningMatrix%col_ind
!print*,'arrLevels(k)%PreconditioningMatrix%row_ptr = ',arrLevels(k)%PreconditioningMatrix%row_ptr
!             read*
	ELSE
		! Convert COO matrix to CSR matrix
		CALL ConvertCOOtoCSR(arrLevels(k)%PreconditioningMatrix,arrLevels(k)%IterationMatrix)
                arrLevels(k)%PreconditioningMatrix_CSR = arrLevels(k)%IterationMatrix

		! reorder CSR matrix to prepare for gauss-seidel operations
		CALL CSRdiaREORDER(arrLevels(k)%IterationMatrix%nrow,arrLevels(k)%IterationMatrix%nnz,&
			arrLevels(k)%IterationMatrix%val,arrLevels(k)%IterationMatrix%col_ind,arrLevels(k)%IterationMatrix%row_ptr)
	END IF
	! Determine stencils for use with the BuildLinearSystem subroutines
	IF (CurvilinearONOFF==1) THEN
	    CALL DetermineGenericStencils(arrLevels(k)%CurvilinearStuff%DiffStencils,gamma) ! FIXME: chosen gamma, assuming that alpha=beta=gamma
        CALL PreProcessDiffStencilsZ(arrLevels(k),arrLevels(k)%DiffStencils,GhostGridZ,gamma)
	ENDIF
END DO

!k=1
!             print*,'k=',k
!print*,'arrLevels(k)%PreconditioningMatrix%val     = ',arrLevels(k)%PreconditioningMatrix%val
!print*,'arrLevels(k)%PreconditioningMatrix%col_ind = ',arrLevels(k)%PreconditioningMatrix%col_ind
!print*,'arrLevels(k)%PreconditioningMatrix%row_ptr = ',arrLevels(k)%PreconditioningMatrix%row_ptr
!read*

! SETUP indexmaps for ghost points
DO k = MG_N_levels, 2, -1
	Nxg = arrLevels(k)%Nx+2*GhostGridX 
	Nyg = arrLevels(k)%Ny+2*GhostGridY
	Nzg = arrLevels(k)%Nz+GhostGridZ
!print*,'=================================='
!print*,'k=',k
!print*,'(Nxg,Nyg,Nzg)=(',Nxg,',',Nyg,',',Nzg,')'
    arrLevels(k)%mapNp = 0
    IF (Nyg>1 .AND. Nxg>1) THEN
        ! 3D
        arrLevels(k)%mapNp = Nxg*Nyg+2*Nxg*(Nzg-1)+2*(Nyg-2)*(Nzg-1)
    ELSE IF (Nxg>1) THEN   ! 2D
        arrLevels(k)%mapNp = Nxg+2*(Nzg-1)
    ELSE ! Nyg>1   ! 2D
        arrLevels(k)%mapNp = Nyg+2*(Nzg-1)
    END IF
!print*,'mapNp = ',arrLevels(k)%mapNp
	ALLOCATE(arrLevels(k)%iGhost(arrLevels(k)%mapNp),STAT=STAT)
    CALL CheckError(STAT,9)
	arrLevels(k)%iGhost(arrLevels(k)%mapNp) = 0
	Gidx  = 0
	Count = 0
	! column major order where first index varies faster than subsequent indexes
    IF (Nyg>1 .AND. Nxg>1) THEN
        ! 3D
	DO jj = 1, Nyg
		DO ii = 1, Nxg
			DO kk = 1, Nzg
				Gidx = Gidx + 1
				IF (kk==1 .OR. ii==1 .OR. ii==Nxg .OR. jj==1 .OR. jj==Nyg) THEN
				    Count = Count + 1				
					arrLevels(k)%iGhost(Count) = Gidx
				END IF
			END DO
		END DO
	END DO
    ELSE IF (Nxg>1) THEN   ! 2D
		DO ii = 1, Nxg
			DO kk = 1, Nzg
				Gidx = Gidx + 1
				IF (kk==1 .OR. ii==1 .OR. ii==Nxg) THEN
				    Count = Count + 1				
					arrLevels(k)%iGhost(Count) = Gidx
				END IF
			END DO
		END DO
    ELSE ! Nyg>1   ! 2D
  	    DO jj = 1, Nyg
			DO kk = 1, Nzg
				Gidx = Gidx + 1
				IF (kk==1 .OR. jj==1 .OR. jj==Nyg) THEN
				    Count = Count + 1				
					arrLevels(k)%iGhost(Count) = Gidx
				END IF
		    END DO
	    END DO    
    END IF
!	if (arrLevels(k)%mapNp/=Count) then
!	print*,'mapNp=',arrLevels(k)%mapNp
!	print*,'Count=',Count
!	print * ,'ighost not correct length for grid level.'
!	STOP
!	end if
END DO
print*,'Multigrid preprocessing step completed.'

CONTAINS

!> Determine number of grid levels based on 
!! number of points in each of the Cartesian directions
!!
!! FIXME: Include support for other grids than 2^p+1 grids
!<
FUNCTION NoGrids(nx,ny,nz) RESULT (K)
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: nx,ny,nz,K
REAL    :: Nmax, tmp
Nmax = MAX(nx,ny,nz)
tmp  = LOG10(Nmax-1)/LOG10(two)
K    = NINT(tmp)
END FUNCTION NoGrids

END SUBROUTINE MGPreProcess

