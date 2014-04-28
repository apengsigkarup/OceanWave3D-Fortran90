SUBROUTINE DfDx_1D_Uneven(f,nx,coeffs,rank,fx)
!
! A subroutine to take the first derivative of the function f on an 
! arbitrary grid with one-sided derivatives towards the two boundaries 
! using the coefficients coeffs(1:rank,2,1:nx).  The coefficients are 
! computed in BuildStencil_1D_Uneven.  
!
! By Allan P. Engsig-Karup.
  USE Precision
  IMPLICIT NONE
  INTEGER nx, i, rank, alpha
  REAL(kind=long) :: f(nx), coeffs(rank,rank,nx), fx(nx)

  alpha=(rank-1)/2


  ! One-sided derivatives towards the left end
  do i=1,alpha
     fx(i)=Dot_Product(f(1:rank),coeffs(1:rank,2,i))
  end do
  ! Centered schemes for the interior
  do i=alpha+1,Nx-alpha
     fx(i)=Dot_Product(f(i-alpha:i+alpha),coeffs(1:rank,2,i))
  end do
  ! One-sided derivatives towards the right end
  do i=Nx-alpha+1,Nx
     fx(i)=Dot_Product(f(Nx-rank+1:Nx),coeffs(1:rank,2,i))
  end do

  RETURN
END SUBROUTINE DfDx_1D_Uneven
