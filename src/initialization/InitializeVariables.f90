SUBROUTINE InitializeVariables
! By Allan P. Engsig-Karup.
USE GlobalVariables
IMPLICIT NONE
INTEGER:: Nx, Ny, Nz
! TEMPORARY POINTERS/VARIABLES USED FOR ALLOCATION
Nx = FineGrid%Nx + 2*GhostGridX
Ny = FineGrid%Ny + 2*GhostGridY
Nz = FineGrid%Nz + GhostGridZ
! Remark the storage format of the variables;
!    x first, then y in free surface variables
! GD: use of a subroutine for initialization of wavefield_FS type... (usefull after)
! variables are set to zero inside
CALL ALLOCATE_Wavefield_Type(Wavefield, FineGrid%Nx, FineGrid%Ny, FineGrid%Nz, GhostGridX, GhostGridy, GhostGridZ, swenseONOFF)
!GD: To nullify pointers of Wavefield_tmp
CALL ALLOCATE_Wavefield_Type(Wavefield_tmp, 1, 1, 1, 0, 0, 0, 0)
CALL DEALLOCATE_Wavefield_Type(Wavefield_tmp, 1, 1, 1, 0)
! FIXME: is it still usefull this temporary array?... to check
ALLOCATE(tmp2D(Nx,Ny))
! FIXME: add this inside the wavefield_FS TYPE ?
ALLOCATE(rhsE(Nx,Ny), rhsP(Nx,Ny))
rhsE=zero; rhsP=zero
! Unknown (semi-potential) PHI's inside computational domain of Transformed Laplace problem
ALLOCATE(PHI(Nz,Nx,Ny), RHS(Nz,Nx,Ny), PHI2(Nz,Nx,Ny), LASTPHI(Nz,Nx,Ny,2))
PHI  = zero; RHS = zero; LASTPHI = zero
PHI2 = zero
time = zero
! HARWELL ROUTINES
RESETSOLVER = 0
! ITERATIONS COUNTERS
MINITER = 100; MAXITER=0; TOTALITER=0
MINITERNOFS = 100; MAXITERNOFS=0; TOTALITERFS=0
!
END SUBROUTINE InitializeVariables
