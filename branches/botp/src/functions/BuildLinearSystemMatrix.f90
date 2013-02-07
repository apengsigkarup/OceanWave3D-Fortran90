SUBROUTINE BuildLinearSystemMatrix(Nx,Ny,Nz,GlobalMatrix,GhostGridX,GhostGridY,GhostGridZ,dsigma,FineGrid,&
		   FullRankStencils,alpha,beta,gamma)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: FineGrid
TYPE (Diff_def) :: FullRankStencils
TYPE (SparseArray_COO) :: GlobalMatrix
INTEGER::Gidx,Gidx2, Nx, Ny, Nz, i, j, k, GhostGridX, GhostGridY, GhostGridZ, alpha, beta, gamma, nnz, stencilrank, dummy
REAL(KIND=long), DIMENSION(Nx,Ny) :: h, hx, hxx, hy, hyy
REAL(KIND=long), DIMENSION(Nx*Ny*Nz,5) :: dsigma
REAL(KIND=long), DIMENSION(:), ALLOCATABLE    :: A_val
INTEGER, DIMENSION(:), ALLOCATABLE :: A_colind, A_rowptr

h = FineGrid%h
IF (Nx>1) THEN
	hx  = FineGrid%hx
	hxx = FineGrid%hxx
ENDIF
IF (Ny>1) THEN
	hy  = FineGrid%hy
	hyy = FineGrid%hyy
ENDIF
stencilrank = (2*alpha+1)*(2*beta+1)*(2*gamma+1)

! ALLOCATE TEMPORARY VECTORS FOR SPARSE OPERATOR ASSEMBLY
! GD: modification to take into account the cross derivatives correctly...
!dummy = Nx*Ny*Nz*stencilrank
dummy = 3*Nx*Ny*Nz*stencilrank !3 times due to different treatment of XZ and YZ derivatives
ALLOCATE( A_val(dummy)    )
ALLOCATE( A_colind(dummy) )
ALLOCATE( A_rowptr(dummy) )
A_val =  zero; A_colind = 0; A_rowptr = 0;
Gidx = 0; nnz = 0
! FIXME: Not most efficient way of generating sparse matrix, however,
!        it only has to be done a single time as a part of pre-processing.
DO j = 1+GhostGridY, Ny-GhostGridY ! FIXME: Problem here if GhostGridY is zero, because same value will be overwritten below
	DO i = 1+GhostGridX, Nx-GhostGridX
		DO k = 1, Nz
			Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
!			Gidx = Gidx + 1
            IF (k==Nz) THEN
				! FREE SURFACE
				nnz = nnz + 1
				A_val(nnz) = one
				A_colind(nnz) = Gidx
  			    A_rowptr(nnz) = Gidx
			ELSE IF (k==1) THEN
				! BOTTOM, DIRECT KIN. BOTTOM CONDITION
				IF (Nx>1 .AND. Ny>1) THEN
					A_val(nnz+1:nnz+stencilrank) = (dsigma(Gidx+GhostGridZ,5) + &
					   hx(i,j)*dsigma(Gidx+GhostGridZ,2) + &
					   hy(i,j)*dsigma(Gidx+GhostGridZ,4) )*FullRankStencils%StencilZ(Gidx+GhostGridZ,:,1)+&
					   hx(i,j)*FullRankStencils%StencilX(Gidx+GhostGridZ,:,1)+&
					   hy(i,j)*FullRankStencils%StencilY(Gidx+GhostGridZ,:,1)
				ELSE IF (Nx>1) THEN
					A_val(nnz+1:nnz+stencilrank) = (dsigma(Gidx+GhostGridZ,5) + &
					   hx(i,j)*dsigma(Gidx+GhostGridZ,2) )*FullRankStencils%StencilZ(Gidx+GhostGridZ,:,1)+&
					   hx(i,j)*FullRankStencils%StencilX(Gidx+GhostGridZ,:,1)
				ELSE IF (Ny>1) THEN
					A_val(nnz+1:nnz+stencilrank) = (dsigma(Gidx+GhostGridZ,5) + &
					   hy(i,j)*dsigma(Gidx+GhostGridZ,4) )*FullRankStencils%StencilZ(Gidx+GhostGridZ,:,1)+&
					   hy(i,j)*FullRankStencils%StencilY(Gidx+GhostGridZ,:,1)
				ENDIF
				! impose boundary condition at what bottom boundary, so shift Gidx-index for stencil
				A_colind(nnz+1:nnz+stencilrank)  = FullRankStencils%Indexes(Gidx+GhostGridZ,:)
			    A_rowptr(nnz+1:nnz+stencilrank)  = Gidx
				nnz = nnz + stencilrank
			ELSE
				! INTERIOR POINTS
				!IF (Nx-2*GhostGridX>1 .AND. Ny-2*GhostGridY>1) THEN
				!	A_val(nnz+1:nnz+stencilrank) = FullRankStencils%StencilX(Gidx,:,2) + &
				!	   FullRankStencils%StencilY(Gidx,:,2) + dsigma(Gidx,3)*FullRankStencils%StencilZ(Gidx,:,1) +&
				!	   two*(dsigma(Gidx,2)*FullRankStencils%StencilXZorYZ(Gidx,:,1) + &
				!	   dsigma(Gidx,4)*FullRankStencils%StencilXZorYZ(Gidx,:,2) ) + &
				!	   (dsigma(Gidx,2)**2+dsigma(Gidx,4)**2+dsigma(Gidx,5)**2)*FullRankStencils%StencilZ(Gidx,:,2)
				!ELSE IF (Nx-2*GhostGridX>1) THEN
				!	A_val(nnz+1:nnz+stencilrank) = FullRankStencils%StencilX(Gidx,:,2) + &
				!	   dsigma(Gidx,3)*FullRankStencils%StencilZ(Gidx,:,1) + &
				!	   two*(dsigma(Gidx,2)*FullRankStencils%StencilXZorYZ(Gidx,:,1) ) + &
				!	   (dsigma(Gidx,2)**2+dsigma(Gidx,5)**2)*FullRankStencils%StencilZ(Gidx,:,2)
				!ELSE IF (Ny-2*GhostGridY>1) THEN
				!	A_val(nnz+1:nnz+stencilrank) = FullRankStencils%StencilY(Gidx,:,2) + &
				!	   dsigma(Gidx,3)*FullRankStencils%StencilZ(Gidx,:,1) + &
				!	   two*(dsigma(Gidx,4)*FullRankStencils%StencilXZorYZ(Gidx,:,2) ) + &
				!	   (dsigma(Gidx,4)**2+dsigma(Gidx,5)**2)*FullRankStencils%StencilZ(Gidx,:,2)
				!ENDIF
				!A_colind(nnz+1:nnz+stencilrank) = FullRankStencils%Indexes(Gidx,:)
			    !A_rowptr(nnz+1:nnz+stencilrank) = Gidx
				!nnz = nnz + stencilrank
                !
                ! GD: modification to take into account the cross derivatives correctly...
                ! INTERIOR POINTS first terms excluding cross derivatives...
                IF (Nx-2*GhostGridX>1 .AND. Ny-2*GhostGridY>1) THEN
                    A_val(nnz+1:nnz+stencilrank) = FullRankStencils%StencilX(Gidx,:,2) + &
                       FullRankStencils%StencilY(Gidx,:,2) + dsigma(Gidx,3)*FullRankStencils%StencilZ(Gidx,:,1) +&
                       (dsigma(Gidx,2)**2+dsigma(Gidx,4)**2+dsigma(Gidx,5)**2)*FullRankStencils%StencilZ(Gidx,:,2)
                ELSE IF (Nx-2*GhostGridX>1) THEN
                    A_val(nnz+1:nnz+stencilrank) = FullRankStencils%StencilX(Gidx,:,2) + &
                       dsigma(Gidx,3)*FullRankStencils%StencilZ(Gidx,:,1) + &
                       (dsigma(Gidx,2)**2+dsigma(Gidx,5)**2)*FullRankStencils%StencilZ(Gidx,:,2)
                ELSE IF (Ny-2*GhostGridY>1) THEN
                    A_val(nnz+1:nnz+stencilrank) = FullRankStencils%StencilY(Gidx,:,2) + &
                       dsigma(Gidx,3)*FullRankStencils%StencilZ(Gidx,:,1) + &
                       (dsigma(Gidx,4)**2+dsigma(Gidx,5)**2)*FullRankStencils%StencilZ(Gidx,:,2)
                ENDIF
!        print*,'FullRankStencils%StencilX(Gidx,:,2) =',FullRankStencils%StencilX(Gidx,:,2) 
!        print*,'FullRankStencils%StencilZ(Gidx,:,1) =',FullRankStencils%StencilZ(Gidx,:,1)
!        print*,'dsigma(Gidx,3)=',dsigma(Gidx,3)
!        print*,'FullRankStencils%StencilZ(Gidx,:,2) =',FullRankStencils%StencilZ(Gidx,:,2)
!        print*,'dsigma(Gidx,2) =',dsigma(Gidx,2)
!        print*,'dsigma(Gidx,5) =',dsigma(Gidx,5)
!        print*,'A_val=',A_val(nnz+1:nnz+stencilrank)
!        read*
                A_colind(nnz+1:nnz+stencilrank) = FullRankStencils%Indexes(Gidx,:)
                A_rowptr(nnz+1:nnz+stencilrank) = Gidx
                nnz = nnz + stencilrank
                ! INTERIOR POINTS XZ derivative
                IF (Nx-2*GhostGridX>1 .AND. Ny-2*GhostGridY>1) THEN
                    A_val(nnz+1:nnz+stencilrank) = two*(dsigma(Gidx,2)*FineGrid%DiffStencils%FullRankStencilXZ(Gidx,:)) !to test: remove the zero factor... FullRankStencilXZ(Gidx,:) !
                    !
                    A_colind(nnz+1:nnz+stencilrank) = FineGrid%DiffStencils%FullRankIndexXZ(Gidx,:)
                    A_rowptr(nnz+1:nnz+stencilrank) = Gidx
                    nnz = nnz + stencilrank
                ELSE IF (Nx-2*GhostGridX>1) THEN
                    A_val(nnz+1:nnz+stencilrank) = two*(dsigma(Gidx,2)*FineGrid%DiffStencils%FullRankStencilXZ(Gidx,1:stencilrank))!to test: remove the zero factor... FullRankStencilXZ(Gidx,:) !
                    !
                    A_colind(nnz+1:nnz+stencilrank) = FineGrid%DiffStencils%FullRankIndexXZ(Gidx,1:stencilrank)
                    A_rowptr(nnz+1:nnz+stencilrank) = Gidx
                    nnz = nnz + stencilrank
                ELSE IF (Ny-2*GhostGridY>1) THEN
                    !A_val(nnz+1:nnz+stencilrank) = zero
                    !XZ cross derivative not defined if Nx=1
                ENDIF
                ! INTERIOR POINTS YZ derivative
                IF (Nx-2*GhostGridX>1 .AND. Ny-2*GhostGridY>1) THEN
                    A_val(nnz+1:nnz+stencilrank) =  two*(dsigma(Gidx,4)*FineGrid%DiffStencils%FullRankStencilYZ(Gidx,:)) !to test: remove the zero factor...
                    !
                    A_colind(nnz+1:nnz+stencilrank) = FineGrid%DiffStencils%FullRankIndexYZ(Gidx,:)
                    A_rowptr(nnz+1:nnz+stencilrank) = Gidx
                    nnz = nnz + stencilrank
                ELSE IF (Nx-2*GhostGridX>1) THEN
                    !A_val(nnz+1:nnz+stencilrank) = zero
                    !YZ cross derivative not defined if Ny=1
                ELSE IF (Ny-2*GhostGridY>1) THEN
                    A_val(nnz+1:nnz+stencilrank) =  two*(dsigma(Gidx,4)*FineGrid%DiffStencils%FullRankStencilYZ(Gidx,1:stencilrank)) !to test: remove the zero factor... FullRankStencilYZ(Gidx,:) !
                    !
                    A_colind(nnz+1:nnz+stencilrank) = FineGrid%DiffStencils%FullRankIndexYZ(Gidx,1:stencilrank)
                    A_rowptr(nnz+1:nnz+stencilrank) = Gidx
                    nnz = nnz + stencilrank
                ENDIF
			END IF
		END DO
	END DO
END DO
! Impose Neumann type boundary conditions at vertical ghost layer points (in the plane!)
! We have deliberately expressed the top and bottom layer ghost points explicitly in terms of the
! interior values. Make sure that this is also done in the Direct matrix-vector routine
! for matching the equations in the gmres.
! FIXME: Assumption here is that the plane is aligned with the x- and y-axes. Extend to
!        general curvilinear coordinates.
IF (GhostGridX==1) THEN
	! West  boundary
	i = 1
	DO j = 1+GhostGridY, Ny-GhostGridY ! omit corners
		!DO k = 1, Nz ! top to bottom
        !GD: no Neumann condition on bottom points (in horizontal plane)
        DO k = 1+GhostGridZ, Nz
			Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz   ! ghost point
			Gidx2 = k + (i+1-1)*Nz + (j-1)*Nx*Nz ! connecting boundary point
			A_val(nnz+1:nnz+stencilrank)    = FullRankStencils%StencilX(Gidx2,:,1)
			A_colind(nnz+1:nnz+stencilrank) = FullRankStencils%Indexes(Gidx2,:)
		    A_rowptr(nnz+1:nnz+stencilrank) = Gidx
			nnz = nnz + stencilrank
		END DO
	END DO
	! East  boundary
	i = Nx
	DO j = 1+GhostGridY, Ny-GhostGridY ! omit corners
		!DO k = 1, Nz ! top to bottom
        !GD: no Neumann condition on bottom points (in horizontal plane)
        DO k = 1+GhostGridZ, Nz
			Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz    ! ghost point
			Gidx2 = k + (i-1-1)*Nz + (j-1)*Nx*Nz ! connecting boundary point
			A_val(nnz+1:nnz+stencilrank)    = FullRankStencils%StencilX(Gidx2,:,1)
			A_colind(nnz+1:nnz+stencilrank) = FullRankStencils%Indexes(Gidx2,:)
		    A_rowptr(nnz+1:nnz+stencilrank) = Gidx
			nnz = nnz + stencilrank
		END DO
	END DO
ENDIF
IF (GhostGridY==1) THEN
	! South boundary
	j = 1
	DO i = 1+GhostGridX, Nx-GhostGridX ! omit corners
		!DO k = 1, Nz ! top to bottom
        !GD: no Neumann condition on bottom points (in horizontal plane)
        DO k = 1+GhostGridZ, Nz
			Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz   ! ghost point
			Gidx2 = k + (i-1)*Nz + (j-1+1)*Nx*Nz ! connecting boundary point
			A_val(nnz+1:nnz+stencilrank)    = FullRankStencils%StencilY(Gidx2,:,1)
			A_colind(nnz+1:nnz+stencilrank) = FullRankStencils%Indexes(Gidx2,:)
		    A_rowptr(nnz+1:nnz+stencilrank) = Gidx
			nnz = nnz + stencilrank
		END DO
	END DO
	! North boundary
	j = Ny
	DO i = 1+GhostGridX, Nx-GhostGridX ! omit corners
		!DO k = 1, Nz ! top to bottom
        !GD: no Neumann condition on bottom points (in horizontal plane)
        DO k = 1+GhostGridZ, Nz
			Gidx  = k + (i-1)*Nz + (j-1)*Nx*Nz   ! ghost point
			Gidx2 = k + (i-1)*Nz + (j-1-1)*Nx*Nz ! connecting boundary point
			A_val(nnz+1:nnz+stencilrank)    = FullRankStencils%StencilY(Gidx2,:,1)
			A_colind(nnz+1:nnz+stencilrank) = FullRankStencils%Indexes(Gidx2,:)
		    A_rowptr(nnz+1:nnz+stencilrank) = Gidx
			nnz = nnz + stencilrank
		END DO
	END DO
ENDIF
! Four corner type points in domain (ghost points)
IF (GhostGridX+GhostGridY==2) THEN
	i = 1; j = 1;
	DO k = 1, Nz
		Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
		nnz           = nnz + 1
		A_val(nnz)    = one
		A_colind(nnz) = Gidx
		A_rowptr(nnz) = Gidx
	END DO
	i = Nx; j = 1;
	DO k = 1, Nz
		Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
		nnz           = nnz + 1
		A_val(nnz)    = one
		A_colind(nnz) = Gidx
		A_rowptr(nnz) = Gidx
	END DO
	i = 1; j = Ny;
	DO k = 1, Nz
		Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
		nnz           = nnz + 1
		A_val(nnz)    = one
		A_colind(nnz) = Gidx
		A_rowptr(nnz) = Gidx
	END DO
	i = Nx; j = Ny;
	DO k = 1, Nz
		Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
		nnz           = nnz + 1
		A_val(nnz)    = one
		A_colind(nnz) = Gidx
		A_rowptr(nnz) = Gidx
	END DO
    ! GD: FIXME the ghost corners on the bottom not treated ?
    ! FIXME: 2D case have to be treated also... (remove the two corner points at the bottom)
    k = 1; j = 1
    DO i = 1+GhostGridX, Nx-GhostGridX
      Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
      nnz           = nnz + 1
      A_val(nnz)    = one
      A_colind(nnz) = Gidx
      A_rowptr(nnz) = Gidx
    ENDDO
    !
    k = 1; j = Ny
    DO i = 1+GhostGridX, Nx-GhostGridX
      Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
      nnz           = nnz + 1
      A_val(nnz)    = one
      A_colind(nnz) = Gidx
      A_rowptr(nnz) = Gidx
    ENDDO
    k = 1; i = 1
    DO j = 1+GhostGridY, Ny-GhostGridY
      Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
      nnz           = nnz + 1
      A_val(nnz)    = one
      A_colind(nnz) = Gidx
      A_rowptr(nnz) = Gidx
    ENDDO
    k = 1; i = Nx
    DO j = 1+GhostGridY, Ny-GhostGridY
      Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
      nnz           = nnz + 1
      A_val(nnz)    = one
      A_colind(nnz) = Gidx
      A_rowptr(nnz) = Gidx
    ENDDO
!ENDIF
! GD: 2D cases bottom corners have to be treated
ELSE IF ((GhostGridX==1).AND.(Ny==1)) THEN
    k = 1; j = 1; i=1
    Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
    nnz           = nnz + 1
    A_val(nnz)    = one
    A_colind(nnz) = Gidx
    A_rowptr(nnz) = Gidx
    !
    k = 1; j = 1; i=Nx
    Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
    nnz           = nnz + 1
    A_val(nnz)    = one
    A_colind(nnz) = Gidx
    A_rowptr(nnz) = Gidx
!ENDIF
ELSE IF ((GhostGridY==1).AND.(Nx==1)) THEN
    k = 1; i = 1; j=1
    Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
    nnz           = nnz + 1
    A_val(nnz)    = one
    A_colind(nnz) = Gidx
    A_rowptr(nnz) = Gidx
    !
    k = 1; i = 1; j=Ny
    Gidx = k + (i-1)*Nz + (j-1)*Nx*Nz
    nnz           = nnz + 1
    A_val(nnz)    = one
    A_colind(nnz) = Gidx
    A_rowptr(nnz) = Gidx
ENDIF

ALLOCATE(GlobalMatrix%val(nnz),GlobalMatrix%row_ptr(nnz),GlobalMatrix%col_ind(nnz))

GlobalMatrix%val     = A_val(1:nnz)
GlobalMatrix%col_ind = A_colind(1:nnz)
GlobalMatrix%row_ptr = A_rowptr(1:nnz)
GlobalMatrix%nnz     = nnz
GlobalMatrix%nrow    = Nx*Ny*Nz

DEALLOCATE( A_val )
DEALLOCATE( A_colind )
DEALLOCATE( A_rowptr )

END SUBROUTINE BuildLinearSystemMatrix
