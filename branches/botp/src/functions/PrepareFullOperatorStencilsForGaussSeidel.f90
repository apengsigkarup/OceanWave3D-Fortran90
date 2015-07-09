SUBROUTINE PrepareFullOperatorStencilsForGaussSeidel(DiffStencils,FullRankStencils,alpha,beta,gamma,Nx,Ny,Nz)
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
TYPE (Diff_def) :: DiffStencils, FullRankStencils
INTEGER, DIMENSION(Nz,Nx,Ny) :: L
INTEGER, DIMENSION(2*gamma+1,2*alpha+1,2*beta+1) :: IndexStencil
INTEGER, DIMENSION((2*gamma+1)*(2*alpha+1)*(2*beta+1)) :: StencilI, reordervector
REAL(KIND=long), DIMENSION(2*gamma+1,2*alpha+1,2*beta+1) :: Stencil
REAL(KIND=long), DIMENSION(2*alpha+1) :: StencilX
REAL(KIND=long), DIMENSION(2*beta+1) :: StencilY
INTEGER, DIMENSION(2*alpha+1) :: idxa, idxa2
INTEGER, DIMENSION(2*beta+1)  :: idxb, idxb2
INTEGER, DIMENSION(2*gamma+1) :: idxg, idxg2
INTEGER :: STAT
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

ALLOCATE( FullRankStencils%StencilX(Nx*Ny*Nz,fullrank,2),STAT=STAT)
CALL CheckError(STAT,2)
ALLOCATE( FullRankStencils%StencilY(Nx*Ny*Nz,fullrank,2),STAT=STAT )
CALL CheckError(STAT,2)
ALLOCATE( FullRankStencils%StencilZ(Nx*Ny*Nz,fullrank,2),STAT=STAT )
CALL CheckError(STAT,2)
ALLOCATE( FullRankStencils%StencilXZorYZ(Nx*Ny*Nz,fullrank,2),STAT=STAT ) ! DXZ, DYZ
CALL CheckError(STAT,2)
ALLOCATE( FullRankStencils%Indexes(Nx*Ny*Nz,fullrank),STAT=STAT )
CALL CheckError(STAT,2)

Gidx = 0
DO j = 1,Ny
	idxb2 = idxb + j
	Diffb=0
	IF (j<beta+1) THEN
		Diffb = beta+1-j
	ELSE IF (j>Ny-beta) THEN
		Diffb = Ny-beta-j
	ENDIF
	idxb2 = idxb2 + Diffb
	DO i = 1, Nx
		idxa2 = idxa + i
		Diffa=0
		IF (i<alpha+1) THEN
			Diffa = alpha+1-i
		ELSE IF (i>Nx-alpha) THEN
			Diffa = Nx-alpha-i
		ENDIF
		idxa2 = idxa2 + Diffa
		DO k = 1,Nz
			Gidx = Gidx + 1
			idxg2 = idxg + k
			Diffg=0
			IF (k<gamma+1) THEN
				Diffg = gamma+1-k
			ELSE IF (k>Nz-gamma) THEN
				Diffg = Nz-gamma-k
			ENDIF
			idxg2 = idxg2+Diffg
			! Determine list of global indices
			IndexStencil = L(idxg2,idxa2,idxb2)
 			FullRankStencils%Indexes(Gidx,:)    = PACK( IndexStencil ,.TRUE. )

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

			! DX
			Stencil=zero
			Stencil(gamma+1-Diffg,:,beta+1-Diffb) = CSHIFT(DiffStencils%StencilX(i,:,1),SHIFT=Diffa)
			FullRankStencils%StencilX(Gidx,:,1) = RESHAPE2(fullrank,Stencil)

			! DXX
			Stencil=zero
			Stencil(gamma+1-Diffg,:,beta+1-Diffb) = CSHIFT(DiffStencils%StencilX(i,:,2),SHIFT=Diffa)
			FullRankStencils%StencilX(Gidx,:,2) = RESHAPE2(fullrank,Stencil)

			! DY
			Stencil=zero
			Stencil(gamma+1-Diffg,alpha+1-Diffa,:) = CSHIFT(DiffStencils%StencilY(j,:,1),SHIFT=Diffb)
			FullRankStencils%StencilY(Gidx,:,1) = RESHAPE2(fullrank,Stencil)

			! DYY
			Stencil=zero
			Stencil(gamma+1-Diffg,alpha+1-Diffa,:) = CSHIFT(DiffStencils%StencilY(j,:,2),SHIFT=Diffb)
			FullRankStencils%StencilY(Gidx,:,2) = RESHAPE2(fullrank,Stencil)

			! DZ
			Stencil=zero
			Stencil(:,alpha+1-Diffa,beta+1-Diffb) = DiffStencils%StencilZ(k,:,1)
			FullRankStencils%StencilZ(Gidx,:,1) = RESHAPE2(fullrank,Stencil)

			! DZZ
			Stencil=zero
			Stencil(:,alpha+1-Diffa,beta+1-Diffb) = DiffStencils%StencilZ(k,:,2)
			FullRankStencils%StencilZ(Gidx,:,2) = RESHAPE2(fullrank,Stencil)

			! DXZ
			StencilX = CSHIFT(DiffStencils%StencilX(i,:,1),SHIFT=Diffa);
			Stencil=zero
			DO jj = beta+1-Diffb,beta+1-Diffb
	  	      DO ii = 1, ranka
				DO kk = 1, rankg
					Stencil(kk,ii,jj) = StencilX(ii)*DiffStencils%StencilZ(k,kk,1)
				END DO
			  END DO
			END DO
			FullRankStencils%StencilXZorYZ(Gidx,:,1) = RESHAPE2(fullrank,Stencil)

			! DYZ
			StencilY = CSHIFT(DiffStencils%StencilY(j,:,1),SHIFT=Diffb)
			Stencil=zero
			DO jj = 1,rankb
	  	      DO ii = alpha+1-Diffa,alpha+1-Diffa
				DO kk = 1, rankg
					Stencil(kk,ii,jj) = StencilY(jj)*DiffStencils%StencilZ(k,kk,1)
				END DO
			  END DO
			END DO
			FullRankStencils%StencilXZorYZ(Gidx,:,2) = RESHAPE2(fullrank,Stencil)

			! REORDERING! (SUCH THAT DIAGONAL ELEMENT IS ALWAYS FIRST ELEMENT IN EACH STENCIL)
			FullRankStencils%StencilX(Gidx,:,1) = FullRankStencils%StencilX(Gidx,reordervector,1)
			FullRankStencils%StencilX(Gidx,:,2) = FullRankStencils%StencilX(Gidx,reordervector,2)
			FullRankStencils%StencilY(Gidx,:,1) = FullRankStencils%StencilY(Gidx,reordervector,1)
			FullRankStencils%StencilY(Gidx,:,2) = FullRankStencils%StencilY(Gidx,reordervector,2)
			FullRankStencils%StencilZ(Gidx,:,1) = FullRankStencils%StencilZ(Gidx,reordervector,1)
			FullRankStencils%StencilZ(Gidx,:,2) = FullRankStencils%StencilZ(Gidx,reordervector,2)
			FullRankStencils%StencilXZorYZ(Gidx,:,1) = FullRankStencils%StencilXZorYZ(Gidx,reordervector,1)
			FullRankStencils%StencilXZorYZ(Gidx,:,2) = FullRankStencils%StencilXZorYZ(Gidx,reordervector,2)
			FullRankStencils%Indexes(Gidx,:)    = FullRankStencils%Indexes(Gidx,reordervector)

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
END SUBROUTINE PrepareFullOperatorStencilsForGaussSeidel