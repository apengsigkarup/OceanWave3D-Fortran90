SUBROUTINE DetermineCurvilinearTransform2D(CompGrid,alpha,beta,gamma,GhostGridX,GhostGridY,GhostGridZ)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
TYPE (Level_def) :: CompGrid
INTEGER :: alpha, beta, gamma, kappa, GhostGridX, GhostGridY, GhostGridZ, Nxg, Nyg, i1, i2
REAL(KIND=long) , DIMENSION(:,:), POINTER :: xe, xn, ye, yn, J, J3, xee, xen, xnn, yee, yen, ynn, &
												 nx, ny, ex, ey, exx, nxx, eyy, nyy, tmp2
Nxg = CompGrid%Nx + 2*GhostGridX
Nyg = CompGrid%Ny + 2*GhostGridY

! Allocate arrays
ALLOCATE(xe(Nxg,Nyg) , xn(Nxg,Nyg) , ye(Nxg,Nyg) , yn(Nxg,Nyg) )
ALLOCATE(J(Nxg,Nyg)  , J3(Nxg,Nyg) , xee(Nxg,Nyg), xen(Nxg,Nyg))
ALLOCATE(xnn(Nxg,Nyg), yee(Nxg,Nyg), yen(Nxg,Nyg), ynn(Nxg,Nyg))
ALLOCATE(exx(Nxg,Nyg), nxx(Nxg,Nyg), eyy(Nxg,Nyg), nyy(Nxg,Nyg))
ALLOCATE(ex(Nxg,Nyg) , ey(Nxg,Nyg) , nx(Nxg,Nyg) , ny(Nxg,Nyg) )
ALLOCATE(tmp2(Nxg,Nyg)                                         )

! Initialize
xe  = one;  xn  = zero; ye  = zero; yn  = one  ! note that two of these are initially set to one (this is convenient for 1D setup)
J   = zero; J3  = zero; xee = zero; xen = zero
xnn = zero; yee = zero; yen = zero; ynn = zero
exx = zero; nxx = zero; eyy = zero; nyy = zero
ex  = zero; ey  = zero; nx  = zero; ny  = zero
tmp2= zero

kappa = alpha
IF (alpha/=beta) THEN
	! FIXME: just picking the largest of alpha and beta here... perhaps check that they are equal in 3D
	kappa = MAX(alpha,beta)
END IF
!ALLOCATE(CompGrid%CurvilinearStuff%DiffStencils%StencilG(2*alpha+1,2*alpha+1,2))
!CALL DetermineGenericStencilsUniform(CompGrid%CurvilinearStuff%DiffStencils%StencilG,kappa)
CALL DetermineGenericStencils(CompGrid%CurvilinearStuff%DiffStencils,kappa)

! First derivatives
IF (CompGrid%Nx>1) THEN
  CALL DiffXuniform2D(CompGrid%x(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),&
     tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),1,&
     CompGrid%CurvilinearStuff%DiffStencils%StencilG,CompGrid%Nx,CompGrid%Ny,alpha)
!  CALL DiffXuniform2D(CompGrid%x,&
!     tmp2,1,CompGrid%CurvilinearStuff%DiffStencils%StencilG,&
!     CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,alpha)
!  xe = tmp2
! GD: to prevent from J=zero on boundary...
  xe(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)= &
  	tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)
  IF (CompGrid%Ny>1) THEN
	CALL DiffYuniform2D(CompGrid%x(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),&
     tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),1,&
     CompGrid%CurvilinearStuff%DiffStencils%StencilG,CompGrid%Nx,CompGrid%Ny,beta)
!	CALL DiffYuniform2D(CompGrid%x,&
!	  tmp2,1,CompGrid%CurvilinearStuff%DiffStencils%StencilG,&
!	  CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,beta)
!	xn = tmp2
! GD: to prevent from J=zero on boundary...
  xn(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)= &
  	tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)
  END IF
END IF
IF (CompGrid%Ny>1) THEN
	CALL DiffYuniform2D(CompGrid%y(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),&
     tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),1,&
     CompGrid%CurvilinearStuff%DiffStencils%StencilG,CompGrid%Nx,CompGrid%Ny,beta)
!   CALL DiffYuniform2D(CompGrid%y,&
 !    tmp2,1,CompGrid%CurvilinearStuff%DiffStencils%StencilG,&
  !   CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,beta)
   !yn = tmp2
   ! GD: to prevent from J=zero on boundary...
   yn(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)= &
  		tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)
   IF (CompGrid%Nx>1) THEN
     CALL DiffXuniform2D(CompGrid%y(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),&
       tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),1,&
       CompGrid%CurvilinearStuff%DiffStencils%StencilG,CompGrid%Nx,CompGrid%Ny,alpha)
!	  CALL DiffXuniform2D(CompGrid%y,&
!		tmp2,1,CompGrid%CurvilinearStuff%DiffStencils%StencilG,&
!		CompGrid%Nx+2*GhostGridX,CompGrid%Ny+2*GhostGridY,alpha)
	  !ye = tmp2
      ! GD: to prevent from J=zero on boundary...
  	  ye(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)= &
  		tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)
   END IF
END IF

! second derivatives
IF (CompGrid%Nx>1) THEN
  CALL DiffXuniform2D(CompGrid%x(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),&
     tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),2,&
     CompGrid%CurvilinearStuff%DiffStencils%StencilG,CompGrid%Nx,CompGrid%Ny,alpha)
  !xee = tmp2
  ! GD: to prevent from assigning non-correct value... tmp2 not evaluated on boundary (not necessary since ghost points not used?!)
  xee(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)= &
  	tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)
  IF (CompGrid%Ny>1) THEN
	CALL DiffYuniform2D(CompGrid%x(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),&
     tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),2,&
     CompGrid%CurvilinearStuff%DiffStencils%StencilG,CompGrid%Nx,CompGrid%Ny,beta)
    !xnn = tmp2
    ! GD: to prevent from assigning non-correct value... tmp2 not evaluated on boundary (not necessary since ghost points not used?!)
  	xnn(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)= &
  		tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)
  END IF
END IF
IF (CompGrid%Ny>1) THEN
  CALL DiffYuniform2D(CompGrid%y(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),&
     tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),2,&
     CompGrid%CurvilinearStuff%DiffStencils%StencilG,CompGrid%Nx,CompGrid%Ny,beta)
  !ynn = tmp2
  ! GD: to prevent from assigning non-correct value... tmp2 not evaluated on boundary (not necessary since ghost points not used?!)
  ynn(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)= &
  	tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)
  IF (CompGrid%Nx>1) THEN
	  CALL DiffXuniform2D(CompGrid%y(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),&
       tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),2,&
       CompGrid%CurvilinearStuff%DiffStencils%StencilG,CompGrid%Nx,CompGrid%Ny,alpha)
     !yee = tmp2
     ! GD: to prevent from assigning non-correct value... tmp2 not evaluated on boundary (not necessary since ghost points not used?!)
  	yee(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)= &
  		tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)
  END IF
END IF

! cross derivatives
!! FIXME: choice has been made with order of differentiation for cross-products here:
IF (CompGrid%Nx>1 .AND. CompGrid%Ny>1) THEN
!  CALL DiffYuniform2D(xe(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),&
!     tmp2,1,CompGrid%CurvilinearStuff%DiffStencils%StencilG,&
!     CompGrid%Nx,CompGrid%Ny,beta)
!  xen = tmp2
!  CALL DiffYuniform2D(ye(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),&
!     tmp2,1,CompGrid%CurvilinearStuff%DiffStencils%StencilG,&
!     CompGrid%Nx,CompGrid%Ny,beta)
!  yen = tmp2
! GD: error on indexes of tmp2
! GD: construct tmp2 on (Nx,Ny) points or (Nx+2*GhostGridX,Ny+2*GhostGridY)? 1st option after...
  CALL DiffYuniform2D(xe(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),&
     tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),1,&
     CompGrid%CurvilinearStuff%DiffStencils%StencilG,CompGrid%Nx,CompGrid%Ny,beta)
  !xen = tmp2
  ! GD: to prevent from assigning non-correct value... tmp2 not evaluated on boundary (not necessary since ghost points not used?!)
  xen(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)= &
  	tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)
  CALL DiffYuniform2D(ye(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),&
     tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY),1,&
     CompGrid%CurvilinearStuff%DiffStencils%StencilG,CompGrid%Nx,CompGrid%Ny,beta)
  !yen = tmp2
  ! GD: to prevent from assigning non-correct value... tmp2 not evaluated on boundary (not necessary since ghost points not used?!)
  yen(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)= &
  	tmp2(1+GhostGridX:CompGrid%Nx+GhostGridX,1+GhostGridY:CompGrid%Ny+GhostGridY)
END IF

J   = xe*yn-ye*xn ! 2D operator

nx  = -ye/J ! d\eta dx
ny  =  xe/J ! d\eta dy
ex  =  yn/J ! d\xi dx
ey  = -xn/J ! d\xi dy

J3  = J**3
exx = one/J3*(-xee*yn**3+two*xen*yn**2*ye-xnn*yn*ye**2+xn*ynn*ye**2-two*xn*yn*ye*yen+xn*yn**2*yee)
nxx = one/J3*(xee*yn**2*ye-two*xen*yn*ye**2-xe*ynn*ye**2+xnn*ye**3+two*xe*yn*ye*yen-xe*yn**2*yee)
eyy = one/J3*(-xnn*xe**2*yn+xn*(two*xe*xen*yn-xn*xee*yn+xe**2*ynn-two*xn*xe*yen+xn**2*yee))
nyy = one/J3*(-xe**3*ynn+xn**2*xee*ye+xe**2*(xnn*ye+two*xn*yen)-xn*xe*(two*xen*ye+xn*yee))
! not needed - but kept here for potential future need
!exy = one/J3*(xe*(-xen*yn**2+xnn*yn*ye-xn*ynn*ye+xn*yn*yen)+xn*(xee*yn**2-xen*yn*ye+xn*ye*yen-xn*yn*yee))
!nxy = one/J3*(xn*ye*(-xee*yn+xen*ye)+xe**2*(ynn*ye-yn*yen)+xe*(xen*yn*ye-xnn*ye**2-xn*ye*yen+xn*yn*yee))


CompGrid%CurvilinearStuff%xe  => xe
CompGrid%CurvilinearStuff%xn  => xn
CompGrid%CurvilinearStuff%xen => xen
CompGrid%CurvilinearStuff%xee => xee
CompGrid%CurvilinearStuff%xnn => xnn
CompGrid%CurvilinearStuff%ye  => ye
CompGrid%CurvilinearStuff%yn  => yn
CompGrid%CurvilinearStuff%yen => yen
CompGrid%CurvilinearStuff%yee => yee
CompGrid%CurvilinearStuff%ynn => ynn
CompGrid%CurvilinearStuff%J   => J
CompGrid%CurvilinearStuff%J3  => J3
CompGrid%CurvilinearStuff%ex  => ex
CompGrid%CurvilinearStuff%ey  => ey
CompGrid%CurvilinearStuff%nx  => nx
CompGrid%CurvilinearStuff%ny  => ny
CompGrid%CurvilinearStuff%exx => exx
CompGrid%CurvilinearStuff%nxx => nxx
CompGrid%CurvilinearStuff%eyy => eyy
CompGrid%CurvilinearStuff%nyy => nyy

END SUBROUTINE DetermineCurvilinearTransform2D
