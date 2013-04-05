SUBROUTINE BuildPDampMatrices(Nx,Ny,PDamp,FineGrid,alpha,beta)
!
! Set up the FD schemes in the damping zone(s) and build the Laplacian matrix 
! in sparse format for the Harwell libs.  
!
! Only implemented in 1D (horizontal) and for one damping zone.  -HBB
!
  USE Precision
  USE Constants
  USE DataTypes
  IMPLICIT NONE
  TYPE (Level_def) :: FineGrid
  TYPE (PDampZone) :: PDamp
  INTEGER ::  Nx, Ny, i, ip, alpha, beta, nnz, rank, dummy, isten, index
  REAL(KIND=long), ALLOCATABLE    :: A_val(:)
  INTEGER, DIMENSION(:), ALLOCATABLE :: A_colind, A_rowptr
  !
  rank = (2*alpha+1)
  !
  ! Build the 1D, derivative operators
  !
  ALLOCATE( PDamp%Grad(rank,rank,Nx) ) 
  !CAll BuildStencil_1D_Uneven(Nx,rank,FineGrid%x(PDamp%idx(1),1),PDamp%Grad)
  ! GD: everything is done on Nx points afterward: don't see where the PDamp%idx(1) comes from? A vector has to be specified
  CAll BuildStencil_1D_Uneven(Nx,rank,FineGrid%x(:,1),PDamp%Grad)
  ! 
  ! ALLOCATE TEMPORARY VECTORS FOR SPARSE OPERATOR ASSEMBLY
  dummy = Nx*Ny*rank 
  ALLOCATE( A_val(dummy), A_colind(dummy), A_rowptr(dummy) ) 
    A_val =  zero; A_colind = 0; A_rowptr = 0; nnz = 0
  ! 
  ! 
  !
  If (Ny>1) Then
     print *, 'BuildPDampMatrices:  Only implemented in 1D at this point...'
     stop
  End If
  !
  ! Build the Laplacian matrix over the zone.  
  ! 
  ! Left boundary points
  i=1 ! Dirichlet point at the left boundary
  index=1
  A_val(index)=one; A_rowptr(index)=i; A_colind(index)=i;
  ! Off-centered Laplacian towards the left.
  Do i=2,alpha
     Do isten=1,rank
        index=index+1
        A_val(index)=PDamp%Grad(isten,3,i); A_rowptr(index)=i; 
        A_colind(index)=isten
     End Do
  END Do
  ! Centered Laplacian 
  Do i=alpha+1,Nx-alpha
     Do isten=1,rank
        index=index+1
        A_val(index)=PDamp%Grad(isten,3,i); A_rowptr(index)=i
        A_colind(index)=i-alpha-1+isten
     End Do
  END Do
  ! Right boundary points, Laplacian
  Do i=Nx-alpha+1,Nx-1
     Do isten=1,rank
        index=index+1
        A_val(index)=PDamp%Grad(isten,3,i); A_rowptr(index)=i 
        A_colind(index)=Nx-rank+isten
     End Do
  END Do
  ! Neumann condition for the last point
  i=Nx
  Do isten=1,rank
     index=index+1
     A_val(index)=PDamp%Grad(isten,2,i); A_rowptr(index)=i 
     A_colind(index)=Nx-rank+isten
  End Do

  nnz=index

  ALLOCATE(PDamp%Lop%val(nnz),PDamp%Lop%irn(nnz),PDamp%Lop%icn(nnz))

  PDamp%Lop%val     = A_val(1:nnz)
  PDamp%Lop%icn = A_colind(1:nnz)
  PDamp%Lop%irn = A_rowptr(1:nnz)
  PDamp%Lop%nnz     = nnz
  PDamp%Lop%nrow    = Nx*Ny

  DEALLOCATE( A_val )
  DEALLOCATE( A_colind )
  DEALLOCATE( A_rowptr )

END SUBROUTINE BuildPDampMatrices
