!>
!! User-defined data types.
!!
!! FIXME: For better compilations, all POINTER's should be changed to ALLOCATABLE.
!! By Allan P. Engsig-Karup.
!<
MODULE DataTypes
USE Precision
IMPLICIT NONE

! Used for interpolation between OceanWave3D and OpenFOAM, botp
TYPE Interpolation_def

        ! order = 2*alpha/beta + 1 derivatives the x/y-direction respectively.  
        !
        REAL(KIND=long), DIMENSION(:,:), POINTER :: dx
        REAL(KIND=long), DIMENSION(:,:), POINTER :: dy
        
        ! order = 2 * gamma + 1 derivatives in the z-direction.
        !
        REAL(KIND=long), DIMENSION(:,:,:), POINTER :: dz

        ! Local interpolation stencils 
        !
        REAL(KIND=long), DIMENSION(:), POINTER :: stencilX
        REAL(KIND=long), DIMENSION(:), POINTER :: stencilY
        REAL(KIND=long), DIMENSION(:), POINTER :: stencilZ

        ! Vector containing nearest neighbour 
        !
        INTEGER, DIMENSION(:), POINTER :: NN

        ! Integer indicating however a cell is above or below the free surface
        ! 0=below, 1=above
        INTEGER, POINTER :: inOrOut

END TYPE Interpolation_def

TYPE SparseArray_CSR
! Sparse storage using: Compressed Storage Row (CSR) format
	REAL(KIND=long), DIMENSION(:), POINTER :: val ! Array values
	INTEGER, DIMENSION(:), POINTER :: col_ind     ! Column indice vector
	INTEGER, DIMENSION(:), POINTER :: row_ptr     ! Row pointer
	INTEGER :: nnz  ! number of nonzero elements
	INTEGER :: nrow ! number of rows (matrix rank)
END TYPE SparseArray_CSR

TYPE SparseArray_COO
! Sparse storage using: Coordinate (COO) format
! Used by the harwell library routines
	REAL(KIND=long), DIMENSION(:), POINTER :: val ! Array values
	INTEGER, DIMENSION(:), POINTER :: col_ind  ! Column pointer
	INTEGER, DIMENSION(:), POINTER :: row_ptr  ! Row pointer
	INTEGER :: nnz ! number of nonzero elements
	INTEGER :: nrow ! rank of matrix
END TYPE SparseArray_COO

TYPE SparseArray_CSR_pdamp
	REAL(KIND=long),DIMENSION(:), ALLOCATABLE :: val(:) ! Array values
	REAL(KIND=long), DIMENSION(:), ALLOCATABLE :: alu(:),w(:)
    INTEGER, ALLOCATABLE, DIMENSION(:) :: jlu,ju,jw 
	INTEGER, allocatable :: irn(:), icn(:)! Row and Column pointers
	INTEGER :: nnz ! number of nonzero elements
	INTEGER :: nrow ! rank of matrix
END TYPE SparseArray_CSR_pdamp
TYPE SparseArray_COO_HBB
! Sparse storage using: Coordinate (COO) format
! Used by the harwell library routines.  This version without pointers 
! and with all the variables.  
	REAL(KIND=long), allocatable:: val(:) ! Array values
	INTEGER, allocatable :: irn(:), icn(:)! Row and Column pointers
	INTEGER :: nnz ! number of nonzero elements
	INTEGER :: nrow ! rank of matrix
        REAL(KIND=long) :: CNTL(10)
        INTEGER :: ICNTL(20), KEEP(50), MAXS, MAXIS
        INTEGER, ALLOCATABLE :: INFOHSL(:), IS_HSL(:)
        REAL(KIND=long), ALLOCATABLE :: COLSCA(:), ROWSCA(:), SS(:), RINFO(:)
END TYPE SparseArray_COO_HBB

! DERIVED DATA TYPE FOR DIFFERENTIAL OPERATIONS
!	DIMENSIONS: Global index, Stencil, order
TYPE Diff_def
	REAL(KIND=long), DIMENSION(:,:,:), POINTER :: StencilG  ! Generic stencil
	REAL(KIND=long), DIMENSION(:,:,:), POINTER :: StencilX  ! Either (Gidx,Coefficients,order) or (i,Coefficients,order)
	REAL(KIND=long), DIMENSION(:,:,:), POINTER :: StencilY  ! Either (Gidx,Coefficients,order) or (j,Coefficients,order)
	REAL(KIND=long), DIMENSION(:,:,:), POINTER :: StencilZ  ! Either (Gidx,Coefficients,order) or (k,Coefficients,order)
	REAL(KIND=long), DIMENSION(:,:,:), POINTER :: StencilXZorYZ 
    INTEGER, DIMENSION(:,:), POINTER :: Indexes
	! Below for ARRAY STORAGE 
    INTEGER, DIMENSION(:,:,:,:), POINTER :: Indexesnew   ! (k,i,j,indexes)
    INTEGER, DIMENSION(:,:,:,:), POINTER :: IndexesnewXZ ! (k,i,j,indexes)
    INTEGER, DIMENSION(:,:,:,:), POINTER :: IndexesnewYZ ! (k,i,j,indexes)
	REAL(KIND=long), DIMENSION(:,:,:,:,:), POINTER :: StencilXZorYZnew ! (k,i,j,indexes,Xz or YZ)
    INTEGER, DIMENSION(:,:,:,:), POINTER :: IndexesG ! (k,i,j,indexes) Generic stencil
    INTEGER, DIMENSION(:,:,:,:), POINTER :: IndexesX ! (k,i,j,indexes)
    INTEGER, DIMENSION(:,:,:,:), POINTER :: IndexesY
    INTEGER, DIMENSION(:,:,:,:), POINTER :: IndexesZ
    ! GD: STORAGE of Indexes and Stencils needed for the treatment of cross derivatives
    INTEGER, DIMENSION(:,:,:,:,:), POINTER :: IndexesX_XZorXY ! (k,i,j,rank,2) 1:XZ ; 2:XY
    INTEGER, DIMENSION(:,:,:,:,:), POINTER :: IndexesY_YZorXY ! (k,i,j,rank,2) 1:YZ ; 2:XY
    INTEGER, DIMENSION(:,:,:,:,:), POINTER :: IndexesZ_XZorYZ ! (k,i,j,rank,2) 1:XZ ; 2:YZ
    REAL(KIND=long), DIMENSION(:,:,:,:,:), POINTER :: StencilsX_XZorXY ! (k,i,j,rank,2) 1:XZ ; 2:XY
    REAL(KIND=long), DIMENSION(:,:,:,:,:), POINTER :: StencilsY_YZorXY ! (k,i,j,rank,2) 1:YZ ; 2:XY
    REAL(KIND=long), DIMENSION(:,:,:,:,:), POINTER :: StencilsZ_XZorYZ ! (k,i,j,rank,2) 1:XZ ; 2:YZ
    ! Full rank Stencils used for straight boundaries preconditionning matrices...
    INTEGER, DIMENSION(:,:), POINTER :: FullRankIndexXZ ! (Gidx,fullrank)
    INTEGER, DIMENSION(:,:), POINTER :: FullRankIndexYZ ! (Gidx,fullrank)
    REAL(KIND=long), DIMENSION(:,:), POINTER :: FullRankStencilXZ ! (Gidx,fullrank)
    REAL(KIND=long), DIMENSION(:,:), POINTER :: FullRankStencilYZ ! (Gidx,fullrank)
    !INTEGER, DIMENSION(:,:,:,:,:), POINTER :: IndexesXZ ! (k,i,j,rank,2) 1:IndexX ; 2:IndexZ : assumed that alpha=gamma
    !INTEGER, DIMENSION(:,:,:,:,:), POINTER :: IndexesYZ ! (k,i,j,rank,2) 1:IndexY ; 2:IndexZ : assumed that beta=gamma
    !INTEGER, DIMENSION(:,:,:,:,:), POINTER :: IndexesXY ! (k,i,j,rank,2) 1:IndexX ; 2:IndexY : assumed that alpha=beta
    !REAL(KIND=long), DIMENSION(:,:,:,:,:), POINTER :: StencilsXZ ! (k,i,j,rank,2) 1:StencilX ; 2:StencilZ : assumed that alpha=gamma
    !REAL(KIND=long), DIMENSION(:,:,:,:,:), POINTER :: StencilsYZ ! (k,i,j,rank,2) 1:StencilY ; 2:StencilZ : assumed that beta=gamma
    !REAL(KIND=long), DIMENSION(:,:,:,:,:), POINTER :: StencilsXY ! (k,i,j,rank,2) 1:StencilX ; 2:StencilY : assumed that alpha=beta
END TYPE Diff_def

! FOR CURVILINEAR PARTS
TYPE Level_def_curvilinear
	REAL(KIND=long) , DIMENSION(:,:), POINTER :: xe, xn, ye, yn, J, J3, xee, xen, xnn, yee, yen, ynn
	REAL(KIND=long) , DIMENSION(:,:), POINTER :: nx, ny, ex, ey, exy, nxy, exx, nxx, eyy, nyy
	REAL(KIND=long) , DIMENSION(:,:,:), POINTER :: NormalX, NormalY, NormalZ ! normal vector components ! FIXME: needs to be stored more compactly
	TYPE (Diff_def) :: DiffStencils         ! stencil and coefficient tables
	TYPE (Diff_def) :: DiffStencilsPrecond  ! to be used for preconditioner
END TYPE Level_def_curvilinear

! DEFINE DERIVED DATA TYPE FOR GRID
! FIXME: Reduce storage needs... only store information ones (in the right format)
TYPE Level_def  
	INTEGER :: Nx , Ny, Nz              ! Grid points in each Cartesian direction
	REAL(KIND=long) , DIMENSION(:,:), POINTER :: x   ! x coordinate array
	REAL(KIND=long) , DIMENSION(:,:), POINTER :: y   ! y coordinate array
	REAL(KIND=long) , DIMENSION(:),   POINTER :: z   ! z coordinate array
	REAL(KIND=long) , DIMENSION(:,:), POINTER :: h, hx, hxx, hy, hyy
	TYPE (Diff_def) :: DiffStencils  ! stencil and coefficient tables
	TYPE (LeveL_def_curvilinear) :: CurvilinearStuff
	TYPE (Diff_def) :: FullRankStencils ! stencil and coefficient tables
	TYPE (Diff_def) :: FullRankStencilsNEW ! stencil and coefficient tables
	REAL(KIND=long), DIMENSION(:,:,:,:), POINTER :: dsigmanew ! array with coordinate index, (k,i,j,:)
	TYPE (SparseArray_COO) :: PreconditioningMatrix
	TYPE (SparseArray_CSR) :: PreconditioningMatrix_CSR
	TYPE (SparseArray_CSR) :: IterationMatrix ! for use with jacobi or gauss-seidel
	INTEGER :: mapNp ! number of points in iGhost
	INTEGER, DIMENSION(:), ALLOCATABLE :: iGhost ! global indexes for ghost points
END TYPE Level_def

! DEFINE STRUCT FOR STREAM FUNCTION SOLUTION PARAMETERS

TYPE SFparam
	REAL(KIND=long) :: g, T, L, k, h, HH, c, e_or_s_vel 
	INTEGER :: i_wavel_or_per, i_euler_or_stokes, i_deep_or_finite, n_h_steps, n_four_modes, nwrk
	REAL(KIND=long), POINTER :: zz(:), yy(:)
END TYPE SFparam

! DEFINE a Structure for random and mono-chromatic wave generation parameters

TYPE RandomWaveParam
	REAL(KIND=long) :: Tp, Hs, h0, kh_max, dx, x0, y0, beta0, S0, gamma
	INTEGER :: ispec, seed, seed2, nf, nx, ny
        CHARACTER(len=30)  inc_wave_file
        REAL(KIND=long), allocatable :: eta(:,:), Phis(:,:), eta0(:), Phis0(:), beta(:)
END TYPE RandomWaveParam
! DEFINE a structure for 3D random wave generation by flux condition, botp
TYPE wave3DFluxStruct
	REAL(KIND=long) :: dt, rampTime
	INTEGER :: n2, order
    CHARACTER(len=30)  inc_wave_file
    REAL(KIND=long), allocatable :: flux(:,:), y(:), time(:), etax(:,:)
END TYPE wave3DFluxStruct
! DEFINE a Structure for the wave breaking model parameters

TYPE BreakingModelParam
   Integer :: i_breaking, n_rollers, i_break_time
   Integer, allocatable :: i_roller(:,:)
   REAL(KIND=long) :: T_half, tan_phi_b, tan_phi_0, del_fac, cel_fac, gamma_break
   REAL(KIND=long),allocatable :: x_roller(:), eta_roller(:), roller_thickness(:), &
        tan_phi(:), alpha_roller(:), beta_roller(:), c_roller(:)
END TYPE BreakingModelParam


TYPE OutputParam
	INTEGER :: xbeg, xend, xstride, ybeg, yend, ystride, tbeg, tend, tstride
END TYPE OutputParam


TYPE RelaxZone
	REAL(KIND=long)                :: BBox(4)        ! Bounding box [xmin xmax ymin ymax]
	REAL(KIND=long)                :: param          ! Parameter
	INTEGER                        :: dir            ! Direction of relaxation function -1 or 1)
	INTEGER                        :: ftype          ! relaxation function to be used
	INTEGER                        :: idx(4)         ! index list start stop [xmin xmax ymin ymax]
	REAL(KIND=long), DIMENSION(:), POINTER :: gam    ! relaxation function values
	CHARACTER(len=1)               :: XorY           ! coordinate direction for relaxation funtion
	INTEGER                        :: WavegenOnOff   ! turn on wavegeneration in the zone (0=off,1=on)
	INTEGER                        :: PhiOnOff       ! relax phi as well as eta (1) or just eta (0).  
	CHARACTER(len=1)               :: XorYgen        ! coordinate direction for generation
    REAL(KIND=long), DIMENSION(:), POINTER :: Ea, Pa ! Storage for analytical solution
    REAL(KIND=long)                :: degrees        ! coordinate rotation angle in degrees
END TYPE RelaxZone

! Pressure Damping Zones.  
! Friction damping on eta, i.e. -gamEta*eta is added to the kinematic free-surface boundary 
! condition.  
! Friction damping on either phi or grad phi:  -gamPhi*phi or L^-1(Div gamPhi*Grad phi) 
! are added to the dynamic free-surface boundary condition. 
!

! DELETE the following /botp
TYPE PDampZone
	REAL(KIND=long)                :: BBox(4)        ! Bounding box [xmin xmax ymin ymax]
	REAL(KIND=long)                :: g0Phi, g0Eta   
	INTEGER                        :: idx(4)         ! index list start stop [xmin xmax ymin ymax]
	INTEGER                        :: nx, ny         ! number of points in the zone
	INTEGER                        :: type           ! Damp GradPhi (0) or Phi (1)
	REAL(KIND=long), allocatable :: gamPhi(:), gamEta(:)  ! Damper function values
        TYPE(SparseArray_COO_HBB) :: Lop                 ! The 2D Laplacian operator in the damping zone
        REAL(kind=long), allocatable :: Grad(:,:,:)      ! The 2D gradient operator in the damping zone
END TYPE PDampZone

TYPE PDampZone_CSR
	REAL(KIND=long)                :: BBox(4)        ! Bounding box [xmin xmax ymin ymax]
	REAL(KIND=long)                :: g0Phi, g0Eta   
	INTEGER                        :: idx(4)         ! index list start stop [xmin xmax ymin ymax]
	INTEGER                        :: nx, ny         ! number of points in the zone
	INTEGER                        :: ierr
	INTEGER                        :: type           ! Damp GradPhi (0) or Phi (1)
	REAL(KIND=long), allocatable :: gamPhi(:), gamEta(:)  ! Damper function values
    TYPE(SparseArray_CSR_pdamp) :: Lop                 ! The 2D Laplacian operator in the damping zone
    REAL(kind=long), allocatable :: Grad(:,:,:)      ! The 2D gradient operator in the damping zone
END TYPE PDampZone_CSR

! New type for wavefield definition on Free Surface (scattered and incident)
TYPE Wavefield_FS
	REAL(KIND=long), DIMENSION(:,:), POINTER :: E, Ex, Exx, Ey, Eyy, P, Px, Py, W     ! Scattered wavefield
        REAL(KIND=long), DIMENSION(:,:), POINTER :: P0, NuD, Pd    ! Pressure terms on the FS
	REAL(KIND=long), DIMENSION(:,:), POINTER :: Qr_x, Mr_t     ! Breaking model terms
	REAL(KIND=long), DIMENSION(:,:), POINTER :: E_I, Ex_I, Exx_I, Ey_I, Eyy_I, Et_I   ! Incident wavefield (free surface)
	REAL(KIND=long), DIMENSION(:,:,:), POINTER :: EtatHist, WHist
	REAL(KIND=long), DIMENSION(:,:), POINTER :: P_I_s, Pz_I_s, Px_I_s, Py_I_s, Pt_I_s ! Incident wavefield (velocity potential)
	INTEGER :: nbp ! number boundary points
	REAL(KIND=long), DIMENSION(:), POINTER   :: E_I_bp, Ex_I_bp, Ey_I_bp                ! Incident wavefield on boundary, FIXME: can be removed
	REAL(KIND=long), DIMENSION(:), POINTER   :: Px_I_bp, Py_I_bp, Pz_I_bp               ! Incident wavefield on boundary
	REAL(KIND=long), DIMENSION(:,:), POINTER :: SourceEx, SourceEy, SourcePx, SourcePy  ! Incident wavefield for boundary treatment
	INTEGER, DIMENSION(:,:,:), POINTER       :: GidxTableBP ! Needed for the curvilinear part
END TYPE Wavefield_FS

END MODULE
