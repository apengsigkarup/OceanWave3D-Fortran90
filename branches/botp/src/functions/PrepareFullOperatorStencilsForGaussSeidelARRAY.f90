SUBROUTINE PrepareFullOperatorStencilsForGaussSeidelARRAY(DiffStencils,FullRankStencils,alpha,beta,gamma,Nx,Ny,Nz)
!
! Prepare all stencils for differential operators needed
! such that they can be added together as full-sized stencils.
!
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: alpha, beta, gamma, Nx, Ny, Nz, Gidx, i,j,k, Diffa, Diffb, Diffg, ii, jj, kk, Lidx, q, count
INTEGER :: ranka, rankb, rankg, fullrank
TYPE (Diff_def), INTENT(IN)  :: DiffStencils
TYPE (Diff_def), INTENT(OUT) :: FullRankStencils
INTEGER, DIMENSION(Nz,Nx,Ny) :: L
INTEGER, DIMENSION(2*gamma+1,2*alpha+1,2*beta+1) :: IndexStencil
INTEGER, DIMENSION((2*gamma+1),(2*alpha+1)) :: IndexStencilXZ
INTEGER, DIMENSION((2*gamma+1),(2*beta+1))  :: IndexStencilYZ
INTEGER, DIMENSION((2*gamma+1)*(2*alpha+1)*(2*beta+1)) :: StencilI, reordervector
INTEGER, DIMENSION((2*gamma+1)*(2*alpha+1)) :: reordervectorXZ
INTEGER, DIMENSION((2*gamma+1)*(2*beta+1))  :: reordervectorYZ
REAL(KIND=long), DIMENSION(2*gamma+1,2*alpha+1) :: Stencil ! FIXME: ASSUMED HERE THAT alpha=beta
!REAL(KIND=long), DIMENSION(2*gamma+1,2*alpha+1,2*beta+1) :: Stencil
INTEGER, DIMENSION(2*alpha+1) :: idxa, idxa2, reordervectorX
INTEGER, DIMENSION(2*beta+1)  :: idxb, idxb2, reordervectorY
INTEGER, DIMENSION(2*gamma+1) :: idxg, idxg2, reordervectorZ
ranka    = 2*alpha+1
rankb    = 2*beta+1
rankg    = 2*gamma+1
fullrank = ranka*rankb*rankg
! Generate array with global indices for every point in the domain
Gidx = 0
DO j = 1,Ny
	DO i = 1, Nx
		DO k = 1,Nz
			Gidx = Gidx + 1
			L(k,i,j) = Gidx
		END DO
	END DO
END DO
! Generate relative index list for a point in the domain for a full rank stencil
idxa=0; idxb=0; idxg=0
DO i = -alpha,alpha
	idxa(i+alpha+1) = i
END DO
DO j = -beta,beta
	idxb(j+beta+1)  = j
END DO
DO k = -gamma,gamma
	idxg(k+gamma+1) = k
END DO

IF (Nx>1) THEN
	ALLOCATE( FullRankStencils%StencilX(Nx,ranka,2) )
	ALLOCATE( FullRankStencils%IndexesX(Nz,Nx,Ny,ranka) )
	ALLOCATE( FullRankStencils%IndexesnewXZ(Nz,Nx,Ny,ranka*rankg) )
ENDIF
IF (Ny>1) THEN
	ALLOCATE( FullRankStencils%StencilY(Ny,rankb,2) )
	ALLOCATE( FullRankStencils%IndexesY(Nz,Nx,Ny,rankb) )
	ALLOCATE( FullRankStencils%IndexesnewYZ(Nz,Nx,Ny,rankb*rankg) )
ENDIF

ALLOCATE( FullRankStencils%StencilZ(Nz,rankg,2) )
ALLOCATE( FullRankStencils%IndexesZ(Nz,Nx,Ny,rankg) )

!ALLOCATE( FullRankStencils%StencilXZorYZnew(Nz,Nx,Ny,fullrank,2) ) ! DXZ, DYZ
ALLOCATE( FullRankStencils%StencilXZorYZnew(Nz,Nx,Ny,ranka*rankg,2) ) ! DXZ, DYZ ! FIXME: ASSUMED HERE THAT ranka = rankb
ALLOCATE( FullRankStencils%Indexesnew(Nz,Nx,Ny,fullrank) )

DO j = 1,Ny
	idxb2 = idxb + j
	Diffb=0
	IF (j<beta+1) THEN
		Diffb = beta+1-j
	ELSE IF (j>Ny-beta) THEN
		Diffb = Ny-beta-j
	ENDIF
	idxb2 = idxb2 + Diffb
	! Reorder stencil such that diagonal element is first and remaining elements
	! are ordered contiuously for current stencil
	IF (Ny>1) THEN
		Lidx = (beta+1)-Diffb ! index for diagonal element
		reordervectorY(1) = Lidx
		count = 1
		DO q = 1, rankb
			IF (q/=Lidx) THEN
				count = count + 1
				reordervectorY(count) = q
			ENDIF
		END DO
		! DY
		FullRankStencils%StencilY(j,:,1) = CSHIFT(DiffStencils%StencilY(j,:,1),SHIFT=Diffb)
		! DYY
		FullRankStencils%StencilY(j,:,2) = CSHIFT(DiffStencils%StencilY(j,:,2),SHIFT=Diffb)

		! REORDER
		FullRankStencils%StencilY(j,:,1) = FullRankStencils%StencilY(j,reordervectorY,1)
		FullRankStencils%StencilY(j,:,2) = FullRankStencils%StencilY(j,reordervectorY,2)
	ELSE
		reordervectorY = 1
	ENDIF

	DO i = 1, Nx
		idxa2 = idxa + i
		Diffa=0
		IF (i<alpha+1) THEN
			Diffa = alpha+1-i
		ELSE IF (i>Nx-alpha) THEN
			Diffa = Nx-alpha-i
		ENDIF
		idxa2 = idxa2 + Diffa

		! Reorder stencil such that diagonal element is first and remaining elements
		! are ordered contiuously for current stencil
		IF (Nx>1) THEN
			Lidx = (alpha+1)-Diffa ! index for diagonal element
			reordervectorX(1) = Lidx
			count = 1
			DO q = 1, ranka
				IF (q/=Lidx) THEN
					count = count + 1
					reordervectorX(count) = q
				ENDIF
			END DO
			! DX
			FullRankStencils%StencilX(i,:,1) = CSHIFT(DiffStencils%StencilX(i,:,1),SHIFT=Diffa)
			! DXX
			FullRankStencils%StencilX(i,:,2) = CSHIFT(DiffStencils%StencilX(i,:,2),SHIFT=Diffa)	

			! REORDER
			FullRankStencils%StencilX(i,:,1) = FullRankStencils%StencilX(i,reordervectorX,1)
			FullRankStencils%StencilX(i,:,2) = FullRankStencils%StencilX(i,reordervectorX,2)
		ELSE
			reordervectorX = 1
		ENDIF

		DO k = 1,Nz

			idxg2 = idxg + k
			Diffg=0
			IF (k<gamma+1) THEN
				Diffg = gamma+1-k
			ELSE IF (k>Nz-gamma) THEN
				Diffg = Nz-gamma-k
			ENDIF
			idxg2 = idxg2+Diffg

			! Determine list og global indices
			IndexStencil = L(idxg2,idxa2,idxb2)
			FullRankStencils%Indexesnew(k,i,j,1:fullrank) = PACK( IndexStencil ,.TRUE. )
			IF (Nx>1) THEN
				IndexStencilXZ = L(idxg2,idxa2,j)
				FullRankStencils%IndexesnewXZ(k,i,j,1:ranka*rankg) = PACK( IndexStencilXZ ,.TRUE. )
				FullRankStencils%IndexesX(k,i,j,1:ranka) = L(k,idxa2,j)

				! XZ
				Lidx = (((alpha+1)-Diffa)-1)*rankg + (gamma+1)-Diffg ! index for diagonal element
				reordervectorXZ(1) = Lidx
				count = 1
				DO q = 1, ranka*rankg
					IF (q/=Lidx) THEN
						count = count + 1
						reordervectorXZ(count) = q
					ENDIF
				END DO

				! DXZ
				Stencil=zero
		        DO ii = 1, ranka
				  DO kk = 1, rankg
					Stencil(kk,ii) = DiffStencils%StencilX(i,ii,1)*DiffStencils%StencilZ(k,kk,1)
				  END DO
				END DO
				FullRankStencils%StencilXZorYZnew(k,i,j,1:ranka*rankg,1) = RESHAPE2(ranka*rankg,Stencil)

				! REORDERING! (SUCH THAT DIAGONAL ELEMENT IS ALWAYS FIRST ELEMENT IN EACH STENCIL)
				FullRankStencils%IndexesX(k,i,j,:) = FullRankStencils%IndexesX(k,i,j,reordervectorX)
				FullRankStencils%StencilXZorYZnew(k,i,j,:,1) = FullRankStencils%StencilXZorYZnew(k,i,j,reordervectorXZ,1)
				FullRankStencils%IndexesnewXZ(k,i,j,:)       = FullRankStencils%IndexesnewXZ(k,i,j,reordervectorXZ)

			ENDIF
			IF (Ny>1) THEN
				IndexStencilYZ = L(idxg2,i,idxb2)
				FullRankStencils%IndexesnewYZ(k,i,j,1:rankb*rankg) = PACK( IndexStencilYZ ,.TRUE. )
				FullRankStencils%IndexesY(k,i,j,1:rankb) = L(k,i,idxb2)

				! YZ
				Lidx = (((beta+1)-Diffb)-1)*rankg + (gamma+1)-Diffg ! index for diagonal element
				reordervectorYZ(1) = Lidx
				count = 1
				DO q = 1, rankb*rankg
					IF (q/=Lidx) THEN
						count = count + 1
						reordervectorYZ(count) = q
					ENDIF
				END DO

				! DYZ
				Stencil=zero
				DO jj = 1,rankb
				  DO kk = 1, rankg
					Stencil(kk,jj) = DiffStencils%StencilY(j,jj,1)*DiffStencils%StencilZ(k,kk,1)
				  END DO
				END DO
				FullRankStencils%StencilXZorYZnew(k,i,j,1:rankb*rankg,2) = RESHAPE2(rankb*rankg,Stencil)

				! REORDERING! (SUCH THAT DIAGONAL ELEMENT IS ALWAYS FIRST ELEMENT IN EACH STENCIL)
				FullRankStencils%IndexesY(k,i,j,:) = FullRankStencils%IndexesY(k,i,j,reordervectorY)
				FullRankStencils%StencilXZorYZnew(k,i,j,:,2) = FullRankStencils%StencilXZorYZnew(k,i,j,reordervectorYZ,2)
				FullRankStencils%IndexesnewYZ(k,i,j,:)       = FullRankStencils%IndexesnewYZ(k,i,j,reordervectorYZ)
			ENDIF

			FullRankStencils%IndexesZ(k,i,j,1:rankg) = L(idxg2,i,j)

			! Reorder stencil such that diagonal element is first and remaining elements
			! are ordered contiuously for current stencil
			Lidx = (((alpha+1)-Diffa)-1)*rankg + (((beta+1)-Diffb)-1)*ranka*rankg + (gamma+1)-Diffg ! index for diagonal element
			reordervector(1) = Lidx
			count = 1
			DO q = 1, fullrank
				IF (q/=Lidx) THEN
					count = count + 1
					reordervector(count) = q
				ENDIF
			END DO

			Lidx = (gamma+1)-Diffg ! index for diagonal element
			reordervectorZ(1) = Lidx
			count = 1
			DO q = 1, rankg
				IF (q/=Lidx) THEN
					count = count + 1
					reordervectorZ(count) = q
				ENDIF
			END DO

			! DZ
			FullRankStencils%StencilZ(k,:,1) = DiffStencils%StencilZ(k,:,1)
			! DZZ
			FullRankStencils%StencilZ(k,:,2) = DiffStencils%StencilZ(k,:,2)

			! REORDERING! (SUCH THAT DIAGONAL ELEMENT IS ALWAYS FIRST ELEMENT IN EACH STENCIL)
			FullRankStencils%IndexesZ(k,i,j,:) = FullRankStencils%IndexesZ(k,i,j,reordervectorZ)
			FullRankStencils%StencilZ(k,:,1)   = FullRankStencils%StencilZ(k,reordervectorZ,1)
			FullRankStencils%StencilZ(k,:,2)   = FullRankStencils%StencilZ(k,reordervectorZ,2)
			FullRankStencils%Indexesnew(k,i,j,:) = FullRankStencils%Indexesnew(k,i,j,reordervector)

		END DO
	END DO
END DO

CONTAINS

FUNCTION RESHAPE2(nnz,A) RESULT(B)
! Convert array with nnz elements to one-dimensional array of nnz elements
USE Precision
IMPLICIT NONE
INTEGER:: nnz
REAL(KIND=long), DIMENSION(nnz) :: A, B
B=A
END FUNCTION RESHAPE2
END SUBROUTINE PrepareFullOperatorStencilsForGaussSeidelARRAY