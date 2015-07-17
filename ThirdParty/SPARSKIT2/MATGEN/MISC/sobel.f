      subroutine sobel(n,nrowc,ncolc,c,jc,ic,a,ja,ia,b,jb,ib,nzmax,ierr)
      integer i, n, ia(*), ja(*), ib(*), jb(*)
      integer nrowa, ncola, nrowb, ncolb, nrowc, ncolc, ipos
      integer ic(*), jc(*), offset, ierr
      real*8  a(*), b(*), c(*)
c-----------------------------------------------------------------------
c     This subroutine generates a matrix used in the statistical problem
c     presented by Prof. Sobel. The matrix is formed by a series of
c     small submatrix on or adjancent to the diagonal. The submatrix on
c     the diagonal is square and the size goes like 1, 1, 2, 2, 3, 3,...
c     Each of the diagonal block is a triadiagonal matrix, each of the
c     off-diagonal block is a bidiagonal block. The values of elements
c     in the off-diagonal block are all -1. So are the values of the
c     elements on the sub- and super-diagonal of the blocks on the
c     diagonal. The first element(1,1) of the diagonal block is alternating
c     between 3 and 5, the rest of the diagonal elements (of the block
c     on the diagonal) are 6.
c-----------------------------------------------------------------------
c     This subroutine calls following subroutines to generate the three
c     thypes of submatrices:
c     diagblk -- generates diagonal block.
c     leftblk -- generates the block left of the diagonal one.
c     rightblk-- generates the block right of the diagonal one.
c-----------------------------------------------------------------------
      if (n.lt.2) return

      ipos = 1
      offset = 1
      call diagblk(1, nrowc, ncolc, c, jc, ic)
      do 10 i=2, n-2
         nrowa = nrowc
         ncola = ncolc
         call copmat (nrowc,c,jc,ic,a,ja,ia,1,1)
         call rightblk(i-1, nrowb, ncolb, b,jb,ib)
         call addblk(nrowa,ncola,a,ja,ia,ipos,ipos+offset,1,
     $        nrowb,ncolb,b,jb,ib,nrowc,ncolc,c,jc,ic,nzmax,ierr)
         call leftblk(i,nrowb,ncolb,b,jb,ib)
         call addblk(nrowc,ncolc,c,jc,ic,ipos+offset,ipos,1,
     $        nrowb,ncolb,b,jb,ib,nrowa,ncola,a,ja,ia,nzmax,ierr)
         ipos = ipos + offset
         call diagblk(i,nrowb,ncolb,b,jb,ib)
         call addblk(nrowa,ncola,a,ja,ia,ipos,ipos,1,
     $        nrowb,ncolb,b,jb,ib,nrowc,ncolc,c,jc,ic,nzmax,ierr)
         offset = 1 + (i-1)/2
 10   continue
      end
c-----------------------------------------------------------------------
      subroutine diagblk(n, nrow, ncol, a, ja, ia)
      implicit none
      integer n, nrow, ncol, ia(1:*), ja(1:*)
      real*8  a(1:*)
c-----------------------------------------------------------------------
c     generates the diagonal block for the given problem.
c-----------------------------------------------------------------------
      integer i, k
      nrow = 1 + (n-1)/2
      ncol = nrow
      k = 1
      ia(1) = 1
      ja(1) = 1
      if (mod(n, 2) .eq. 1) then
         a(1) = 3
      else
         a(1) = 5
      end if
      k = k + 1
      if (ncol.gt.1) then
         ja(2) = 2
         a(2) = -1.0
         k = k + 1
      end if
      
      do 10 i = 2, nrow
         ia(i) = k
         ja(k) = i-1
         a(k) = -1.0
         k = k + 1
         ja(k) = i
         a(k) = 6.0
         k = k + 1
         if (i.lt.nrow) then
            ja(k) = i + 1
            a(k) = -1.0
            k = k+1
         end if
 10   continue
      ia(nrow+1) = k
      return
c---------end-of-diagblk------------------------------------------------ 
      end
c-----------------------------------------------------------------------
      subroutine leftblk(n, nrow, ncol, a, ja, ia)
      implicit none
      integer n, nrow, ncol, ja(1:*), ia(1:*)
      real*8  a(1:*)
c-----------------------------------------------------------------------
c     Generate the subdiagonal block for the problem given.
c-----------------------------------------------------------------------
      integer i, k
      nrow = 1 + (n-1)/2
      ncol = n/2
      k = 1
      do 10 i = 1, nrow
         ia(i) = k
         if (nrow.ne.ncol) then
            if (i.gt.1) then
               ja(k) = i-1
               a(k) = -1.0
               k = k+1
            end if
         end if
         if (i.le.ncol) then
            ja(k) = i
            a(k) = -1.0
            k = k+1
         end if
         if (nrow.eq.ncol) then
            if (i.lt.ncol) then
               ja(k) = i+1
               a(k) = -1.0
               k = k+1
            end if
         end if
 10   continue
      ia(nrow+1) = k
      return
c---------end-of-leftblk------------------------------------------------ 
      end
c-----------------------------------------------------------------------
      subroutine rightblk(n, nrow, ncol, a, ja, ia)
      implicit none
      integer n, nrow, ncol, ja(1:*), ia(1:*)
      real*8  a(1:*)
      integer i, k
      nrow = 1 + (n-1)/2
      ncol = 1 + n/2
      k = 1
      do 10 i = 1, nrow
         ia(i) = k
         if (nrow.eq.ncol) then
            if (i.gt.1) then
               ja(k) = i-1
               a(k) = -1.0
               k = k+1
            end if
         end if
         ja(k) = i
         a(k) = -1.0
         k = k+1
         if (nrow.ne.ncol) then
            ja(k) = i+1
            a(k) = -1.0
            k = k+1
         end if
 10   continue
      ia(nrow+1) = k
c---------end-of-rightblk----------------------------------------------- 
      end
