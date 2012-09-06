!>
!! This module contains definitions for all global variables.
!<
MODULE GlobalVariables
USE Precision
USE Constants
USE DataTypes
USE HSL_LU
IMPLICIT NONE 

! I/O File handle arrays
INTEGER           :: FILEIP(4), FILEOP(14)
CHARACTER(len=2) fnt(10)
CHARACTER(LEN=40) :: filenameINPUT, filename, fname_bottom
INTEGER           :: STAT

! DUMMY VARIABLES
!INTEGER :: i,j,k  !*hbb not a good idea to make these global!  
REAL(KIND=long) :: dum

CHARACTER(len=72) :: HEAD(4)

TYPE (Level_def) :: FineGrid

INTEGER :: alpha, beta, gamma, kappa
INTEGER :: alphaprecond, betaprecond, gammaprecond
INTEGER :: GridX, GridY, GridZ ! Type of grid in each Cartesian direction (0=even)
INTEGER :: GhostGridX, GhostGridY, GhostGridZ, DetermineBottomGradients
REAL(KIND=long) :: Lx, Ly, Lz  ! Length of computational domain in each direction
REAL(KIND=long) :: g, accel_tol_fact=1000.           ! gravitational acceleration constant

! SOLUTION VARIABLES
REAL(KIND=long), DIMENSION(:,:,:), POINTER :: PHI, RHS, PHI2 ! Solution variable in transformed domain
REAL(KIND=long), DIMENSION(:,:,:,:), POINTER :: LASTPHI ! Solution variable in transformed domain
TYPE (Wavefield_FS) :: Wavefield, Wavefield_tmp !tmp_wavefield
REAL(KIND=long), DIMENSION(:,:),   POINTER :: EA,WA
REAL(KIND=long), DIMENSION(:,:),   POINTER :: tmp2D, rhsE, rhsP

! Initial condition and wave generation flags
Integer :: IC, IncWaveType

! FIXME: TO BE DISCARDED
REAL(KIND=long), DIMENSION(:,:),   ALLOCATABLE :: A
REAL(KIND=long), DIMENSION(:),   ALLOCATABLE :: ee, tm

! TIME INTEGRATION PARAMETERS
INTEGER :: Nsteps, tstep, timemethod, RKSTAGES, extrapolationONOFF
REAL(KIND=long) :: dt, time, time0, CFL, c, dxmin, dymin, dx, dy, dsigmamin

INTEGER :: Precond, MGCoarseningStrategy, MGmaxgrids, GMRESmaxiterations

REAL(KIND=long), DIMENSION(:,:,:), ALLOCATABLE :: PHIfine, PHIcoarse

! ITERATIVE SCHEMES
REAL(KIND=long) :: reltol, abstol
INTEGER :: maxit
CHARACTER(len=1) :: cyclet
INTEGER :: MINITER, MAXITER, TOTALITER, TOTALITERFS, MINITERNOFS, MAXITERNOFS, TOTALITEROLD

! HARWELL ROUTINES
INTEGER :: ipar(16)
REAL(KIND=long) :: fpar(16)
REAL(KIND=long), DIMENSION(:), ALLOCATABLE :: workspace
INTEGER :: RESETSOLVER

! VARIABLES FOR MEASURING WALL CLOCK TIME
INTEGER :: count_0, count_1, count_rate, count_max, CPUinitial

CHARACTER(len=40) :: text

! VARIABLES FOR STREAM FUNCTION SOLUTION
TYPE (SFparam) :: SFsol

! VARIABLES FOR Random waves
INTEGER :: n_fft, n_wavem, j_wavem
TYPE (RandomWaveParam) :: RandomWave

! VARIABLES FOR the breaking model
TYPE (BreakingModelParam) :: BreakMod

! FILTERING
INTEGER :: filteringONOFF, filterNP, filterORDER, filterALPHA
REAL(KIND=long), DIMENSION(:), ALLOCATABLE :: filtercoefficients, tmpfilter
! GD addition (filtering with boundaries)
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: filtercoefficients2
REAL(KIND=long) :: sigma_filt(3)

! RELAXATION ZONES
INTEGER          :: relaxONOFF, relaxNo
CHARACTER(len=1) :: relaxXorY
REAL(KIND=long)  :: relaxTransientTime, relaxDegrees
TYPE (RelaxZone), DIMENSION(:), ALLOCATABLE :: RelaxZones ! Array with relaxation zone definitions

! Output Data Parameters
INTEGER :: StoreDataONOFF, formattype, iKinematics, nOutFiles
TYPE (OutputParam), DIMENSION(:), ALLOCATABLE :: Output

! Linear/Nonlinear
INTEGER :: LinearONOFF

! Pressure Term
INTEGER :: PressureTermONOFF

! SWENSE additions (incident wavefield is inside the wavefield_FS TYPE)
INTEGER    :: swenseONOFF, swenseDir, West_refl, East_refl, North_refl, South_refl
REAL(KIND=long) :: swenseTransientTime

! Curvilinear
INTEGER :: curvilinearONOFF

!-------------------------------------------------------------------------------
! Moving ship / pressure distribution parameters / surface integration:
!-------------------------------------------------------------------------------
REAL(KIND=long)     :: Uship,Lship,Bship,Dship,x0ship,Fnship
CHARACTER(len=72)   :: shipfile 
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: wgt_ship
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: nx_ship
REAL(KIND=long), dimension(:,:), ALLOCATABLE :: dzdxi1,dzdxi2

!-------------------------------------------------------------------------------
! WENO: 
!-------------------------------------------------------------------------------
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: FEx,FPx 
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: FEy,FPy
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: aw,bw,cw
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: uL,uR
REAL(KIND=long), DIMENSION(:,:), ALLOCATABLE :: vL,vR

END MODULE GlobalVariables
