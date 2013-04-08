SUBROUTINE MultiGridSolveLinearSystem(r,u,Nx,Ny,Nz,GhostGridX,GhostGridY,GhostGridZ,cyclet,maxit,reltol,alpha,beta,gamma)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE GlobalVariables, ONLY: FineGrid
USE MGLevels
IMPLICIT NONE
INTEGER          :: maxit, iter, gam, GhostGridX, GhostGridY, GhostGridZ, Nx, Ny, Nz, alpha, beta, gamma
INTEGER          :: Nxg, Nyg, Nzg
REAL(KIND=long)  :: res, resold, res_init, reltol
!REAL(KIND=long), DIMENSION(maxit+1) :: resvec
REAL(KIND=long), DIMENSION((Nx+2*GhostGridX)*(Ny+2*GhostGridY)*(Nz+GhostGridZ)), INTENT(INOUT) :: u
REAL(KIND=long), DIMENSION((Nx+2*GhostGridX)*(Ny+2*GhostGridY)*(Nz+GhostGridZ)) :: Au, tmp
REAL(KIND=long), DIMENSION((Nx+2*GhostGridX)*(Ny+2*GhostGridY)*(Nz+GhostGridZ)), INTENT(IN) :: r
CHARACTER(len=1) :: cyclet

Nxg = Nx+2*GhostGridX
Nyg = Ny+2*GhostGridY
Nzg = Nz+GhostGridZ

!print*,'r=',r
!print*,'u=',u
!read*

! Initial Guess on solution
u = zero

! Cycle Type
IF (cyclet=='W' .OR. cyclet=='F') THEN
	gam = 2
ELSE
	gam = 1
END IF

!IF (0) THEN
!  tmp = r !-Au
!  res_init = MAXNORM(Nx*Ny*(Nz+GhostGridZ),tmp) ! initial residual res = r-Mu = r
!  res       = one; !ten*reltol;
!  resvec    = zero;
!  resvec(1) = one;  ! normalized
!ENDIF

iter      = 0
DO
	iter = iter + 1
	IF (iter>maxit) EXIT
!
!	IF (res<reltol) EXIT
!
	IF (MG_N_levels==1 .AND. iter>1) EXIT

	IF (MG_N_levels==1) THEN
		CALL MultiGridSolver(u,r,MG_N_levels,gam,cyclet,arrLevels(MG_N_levels)%Nx,arrLevels(MG_N_levels)%Ny,&
		   arrLevels(MG_N_levels)%Nz,arrLevels(MG_N_levels)%Nx,arrLevels(MG_N_levels)%Ny,arrLevels(MG_N_levels)%Nz,&
		   GhostGridX, GhostGridY, GhostGridZ, MIN(1,alpha), MIN(1,beta), MIN(1,gamma))
	ELSE
		CALL MultiGridSolver(u,r,MG_N_levels,gam,cyclet,arrLevels(MG_N_levels)%Nx,arrLevels(MG_N_levels)%Ny,&
		   arrLevels(MG_N_levels)%Nz,arrLevels(MG_N_levels-1)%Nx,arrLevels(MG_N_levels-1)%Ny,&
		   arrLevels(MG_N_levels-1)%Nz,GhostGridX, GhostGridY, GhostGridZ, MIN(1,alpha), MIN(1,beta), MIN(1,gamma))
	END IF

!IF (0) THEN
!
!	! Determine relative change of residual
!	resold = res
!
	! DETERMINE RESIDUAL error vector (FIXME: DO NOT NEED IT FOR FIXED ITERATIONS)
!	tmp = u
!!	! Update Ghost layer
!	CALL BuildLinearSystemGhost(Nx,Ny,Nz+GhostGridZ,tmp,tmp,arrLevels(MG_N_levels)%dsigma,arrLevels(MG_N_levels)%h,arrLevels(MG_N_levels)%hx,arrLevels(MG_N_levels)%hxx,arrLevels(MG_N_levels)%hy,arrLevels(MG_N_levels)%hyy,arrLevels(MG_N_levels)%DiffStencils,
!!	! Do Matrix-Vector product A*u to determine new residual
!	CALL BuildLinearSystem(Nx,Ny,Nz+GhostGridZ,tmp,Au,arrLevels(MG_N_levels),alpha,beta,gamma)
!!	! Make sure values in ghost layer are zero corresponding to eliminated conditions
!    CALL DiscardGhostLayer(Nx,Ny,Nz+GhostGridZ,Au)
!!
!	tmp = r-Au
!	res = MAXNORM(Nx*Ny*(Nz+GhostGridZ),tmp)
!	res = res/res_init
!!
!	resvec(iter+1) = res
!ENDIF

!	print*,'u=',u
!	print*,'r=',r
!	print*,'Au=',Au
!	print*,'iter=',iter
!	print*,'Iter=',iter,', Residual=',res
!	print*,'maxit=',maxit
! print*,'Residual=',res
! read*
END DO
!print*,u
!read*
RETURN
END SUBROUTINE MultiGridSolveLinearSystem
