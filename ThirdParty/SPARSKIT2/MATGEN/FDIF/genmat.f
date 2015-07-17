c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c    MATRIX GENERATION ROUTINES  -- FINITE DIFFERENCE MATRICES         c
c----------------------------------------------------------------------c
c contents:                                                            c
c----------                                                            c
c gen57pt  : generates 5-point and 7-point matrices.                   c
c gen57bl  : generates block 5-point and 7-point matrices.             c
c                                                                      c
c supporting routines:                                                 c
c---------                                                             c
c gensten  : generate the stencil (point version)                      c
c bsten    : generate the stencil (block version)                      c
c fdaddbc  : finite difference add boundary conditions                 c
c fdreduce : reduce the system to eliminate node with known values     c
c clrow    : clear a row of a CSR matrix                               c
c lctcsr   : locate the position of A(i,j) in CSR format               c
c----------------------------------------------------------------------c
      subroutine gen57pt(nx,ny,nz,al,mode,n,a,ja,ia,iau,rhs)
      integer ja(*),ia(*),iau(*), nx, ny, nz, mode, n
      real*8 a(*), rhs(*), al(6)
c-----------------------------------------------------------------------
c On entry:
c
c nx      = number of grid points in x direction
c ny	  = number of grid points in y direction
c nz	  = number of grid points in z direction
c al      = array of size 6, carries the coefficient alpha of the
c           boundary conditions
c mode    = what to generate:
c           < 0 : generate the graph only,
c           = 0 : generate the matrix,
c           > 0 : generate the matrix and the right-hand side.
c
c On exit:
c
c n       = number of nodes with unknown values, ie number of rows
c           in the matrix
c
c a,ja,ia = resulting matrix in row-sparse format
c
c iau     = integer*n, containing the poisition of the diagonal element
c           in the a, ja, ia structure
c
c rhs     = the right-hand side
c
c External functions needed (must be supplied by caller)
c     afun, bfun, cfun, dfun, efun, ffun, gfun, hfun
c     betfun, gamfun
c They have the following prototype:
c     real*8 function xfun(x, y, z)
c     real*8 x, y, z
c-----------------------------------------------------------------------
c This subroutine computes the sparse matrix in compressed sparse row
c format for the elliptic equation:
c       d    du    d    du    d    du      du     du     du
c L u = --(A --) + --(B --) + --(C --) + D -- + E -- + F -- + G u = H u
c       dx   dx    dy   dy    dz   dz      dx     dy     dz
c
c with general Mixed Boundary conditions, on a rectangular 1-D,
c 2-D or 3-D grid using 2nd order centered difference schemes.
c
c The functions a, b, ..., g, h are known through the
c as afun, bfun, ..., gfun, hfun in this subroutine.
c NOTE: To obtain the correct matrix, any function that is not
c needed should be set to zero.  For example for two-dimensional
c problems, nz should be set to 1 and the functions cfun and ffun
c should be zero functions.
c
c The Boundary condition is specified in the following form:
c           du
c     alpha -- + beta u = gamma
c           dn
c Where alpha is constant at each side of the boundary surfaces.  Alpha
c is represented by parameter al.  It is expected to an array that
c contains enough elements to specify the boundaries for the problem,
c 1-D case needs two elements, 2-D needs 4 and 3-D needs 6.  The order
c of the boundaries in the array is left(west), right(east),
c bottom(south), top(north), front, rear.  Beta and gamma are functions
c of type real with three arguments x, y, z.  These two functions are
c known subroutine 'addbc' as betfun and gamfun.  They should following
c the same notion as afun ... hfun.  For more restriction on afun ...
c hfun, please read the documentation follows the subroutine 'getsten',
c and, for more on betfun and gamfun, please refer to the documentation
c under subroutine 'fdaddbc'.
c
c The nodes are ordered using natural ordering, first x direction, then
c y, then z.  The mesh size h is uniform and determined by grid points
c in the x-direction.
c
c The domain specified for the problem is [0 .ge. x .ge. 1],
c [0 .ge. y .ge. (ny-1)*h] and [0 .ge. z .ge. (nz-1)*h], where h is
c 1 / (nx-1).  Thus if non-Dirichlet boundary condition is specified,
c the mesh will have nx points along the x direction, ny along y and
c nz along z.  For 1-D case, both y and z value are assumed to zero
c when calling relavent functions that have three parameters.
c Similarly, for 2-D case, z is assumed to be zero.
c
c About the expectation of nx, ny and nz:
c nx is required to be .gt. 1 always;
c if the second dimension is present in the problem, then ny should be
c .gt. 1, else 1;
c if the third dimension is present in the problem, nz .gt. 1, else 1.
c when ny is 1, nz must be 1.
c-----------------------------------------------------------------------
c
c     stencil [1:7] has the following meaning:
c
c     center point = stencil(1)
c     west point = stencil(2)
c     east point = stencil(3)
c     south point = stencil(4)
c     north point = stencil(5)
c     front point = stencil(6)
c     back point = stencil(7)
c
c     al[1:6] carry the coefficient alpha in the similar order
c
c     west  side = al(1)
c     east  side = al(2)
c     south side = al(3)
c     north side = al(4)
c     front side = al(5)
c     back  side = al(6)
c
c                           al(4)
c                           st(5)
c                            |
c                            |
c                            |           al(6)
c                            |          .st(7)
c                            |     .
c         al(1)              | .             al(2)
c         st(2) ----------- st(1) ---------- st(3)
c                       .    |
c                   .        |
c               .            |
c            st(6)           |
c            al(5)           |
c                            |
c                           st(4)
c                           al(3)
c
c-------------------------------------------------------------------
c     some constants
c
      real*8 one
      parameter (one=1.0D0)
c
c     local variables
c
      integer ix, iy, iz, kx, ky, kz, node, iedge
      real*8  r, h, stencil(7)
      logical value, genrhs
c
c     nx has to be larger than 1
c
      if (nx.le.1) return
      h = one / dble(nx-1)
c
c     the mode
c
      value = (mode.ge.0)
      genrhs = (mode.gt.0)
c
c     first generate the whole matrix as if the boundary condition does
c     not exist
c
      kx = 1
      ky = nx
      kz = nx*ny
      iedge = 1
      node = 1
      do 100 iz = 1,nz
         do 90 iy = 1,ny
            do 80 ix = 1,nx
               ia(node) = iedge
c
c     compute the stencil at the current node
c
               if (value) call
     &              getsten(nx,ny,nz,mode,ix-1,iy-1,iz-1,stencil,h,r)
c     west
               if (ix.gt.1) then
                  ja(iedge)=node-kx
		  if (value) a(iedge) = stencil(2)
                  iedge=iedge + 1
               end if
c     south
               if (iy.gt.1) then
                  ja(iedge)=node-ky
		  if (value) a(iedge) = stencil(4)
                  iedge=iedge + 1
               end if
c     front plane
               if (iz.gt.1) then
                  ja(iedge)=node-kz
		  if (value) a(iedge) = stencil(6)
                  iedge=iedge + 1
               endif
c     center node
               ja(iedge) = node
               iau(node) = iedge
               if (value) a(iedge) = stencil(1)
               iedge = iedge + 1
c     east
               if (ix.lt.nx) then
                  ja(iedge)=node+kx
		  if (value) a(iedge) = stencil(3)
                  iedge=iedge + 1
               end if
c     north
               if (iy.lt.ny) then
                  ja(iedge)=node+ky
		  if (value) a(iedge) = stencil(5)
                  iedge=iedge + 1
               end if
c     back plane
               if (iz.lt.nz) then
                  ja(iedge)=node+kz
                  if (value) a(iedge) = stencil(7)
                  iedge=iedge + 1
               end if
c     the right-hand side
               if (genrhs) rhs(node) = r
               node=node+1
 80         continue
 90      continue
 100  continue
      ia(node)=iedge
c
c     Add in the boundary conditions
c
      call fdaddbc(nx,ny,nz,a,ja,ia,iau,rhs,al,h)
c
c     eliminate the boudary nodes from the matrix
c
      call fdreduce(nx,ny,nz,al,n,a,ja,ia,iau,rhs,stencil)
c
c     done
c
      return
c-----end-of-gen57pt----------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine getsten (nx,ny,nz,mode,kx,ky,kz,stencil,h,rhs)
      integer nx,ny,nz,mode,kx,ky,kz
      real*8 stencil(*),h,rhs,afun,bfun,cfun,dfun,efun,ffun,gfun,hfun
      external afun,bfun,cfun,dfun,efun,ffun,gfun,hfun
c-----------------------------------------------------------------------
c     This subroutine calculates the correct stencil values for
c     centered difference discretization of the elliptic operator
c     and the right-hand side
c
c L u = delx( A delx u ) + dely ( B dely u) + delz ( C delz u ) +
c	delx ( D u ) + dely (E u) + delz( F u ) + G u = H
c
c   For 2-D problems the discretization formula that is used is:
c
c h**2 * Lu == A(i+1/2,j)*{u(i+1,j) - u(i,j)} +
c	       A(i-1/2,j)*{u(i-1,j) - u(i,j)} +
c              B(i,j+1/2)*{u(i,j+1) - u(i,j)} +
c              B(i,j-1/2)*{u(i,j-1) - u(i,j)} +
c              (h/2)*D(i,j)*{u(i+1,j) - u(i-1,j)} +
c              (h/2)*E(i,j)*{u(i,j+1) - u(i,j-1)} +
c              (h/2)*E(i,j)*{u(i,j+1) - u(i,j-1)} +
c              (h**2)*G(i,j)*u(i,j)
c-----------------------------------------------------------------------
c     some constants
c
      real*8 zero, half
      parameter (zero=0.0D0,half=0.5D0)
c
c     local variables
c
      integer k
      real*8 hhalf,cntr, x, y, z, coeff
c
c     if mode < 0, we shouldn't have come here
c
      if (mode .lt. 0) return
c
      do 200 k=1,7
         stencil(k) = zero
 200  continue
c
      hhalf = h*half
      x = h*dble(kx)
      y = h*dble(ky)
      z = h*dble(kz)
      cntr = zero
c     differentiation wrt x:
      coeff = afun(x+hhalf,y,z)
      stencil(3) = stencil(3) + coeff
      cntr = cntr + coeff
c
      coeff = afun(x-hhalf,y,z)
      stencil(2) = stencil(2) + coeff
      cntr = cntr + coeff
c
      coeff = dfun(x,y,z)*hhalf
      stencil(3) = stencil(3) + coeff
      stencil(2) = stencil(2) - coeff
      if (ny .le. 1) goto 99
c
c     differentiation wrt y:
c
      coeff = bfun(x,y+hhalf,z)
      stencil(5) = stencil(5) + coeff
      cntr = cntr + coeff
c
      coeff = bfun(x,y-hhalf,z)
      stencil(4) = stencil(4) + coeff
      cntr = cntr + coeff
c
      coeff = efun(x,y,z)*hhalf
      stencil(5) = stencil(5) + coeff
      stencil(4) = stencil(4) - coeff
      if (nz .le. 1) goto 99
c
c differentiation wrt z:
c
      coeff = cfun(x,y,z+hhalf)
      stencil(7) = stencil(7) + coeff
      cntr = cntr + coeff
c
      coeff = cfun(x,y,z-hhalf)
      stencil(6) = stencil(6) + coeff
      cntr = cntr + coeff
c
      coeff = ffun(x,y,z)*hhalf
      stencil(7) = stencil(7) + coeff
      stencil(6) = stencil(6) - coeff
c
c contribution from function G:
c
 99   coeff = gfun(x,y,z)
      stencil(1) = h*h*coeff - cntr
c
c     the right-hand side
c
      if (mode .gt. 0) rhs = h*h*hfun(x,y,z)
c
      return
c------end-of-getsten---------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine gen57bl (nx,ny,nz,nfree,na,n,a,ja,ia,iau,stencil)
c     implicit real*8 (a-h,o-z)
      integer ja(*),ia(*),iau(*),nx,ny,nz,nfree,na,n
      real*8 a(na,1), stencil(7,1)
c--------------------------------------------------------------------
c This subroutine computes the sparse matrix in compressed
c format for the elliptic operator
c
c L u = delx( a . delx u ) + dely ( b . dely u) + delz ( c . delz u ) +
c	delx ( d . u ) + dely (e . u) + delz( f . u ) + g . u
c
c Here u is a vector of nfree componebts and each of the functions
c a, b, c, d, e, f, g   is an (nfree x nfree) matrix depending of
c the coordinate (x,y,z).
c with Dirichlet Boundary conditions, on a rectangular 1-D,
c 2-D or 3-D grid using centered difference schemes.
c
c The functions a, b, ..., g are known through the
c subroutines  afunbl, bfunbl, ..., gfunbl. (user supplied) .
c
c uses natural ordering, first x direction, then y, then z
c mesh size h is uniform and determined by grid points
c in the x-direction.
c 
c The output matrix is in Block -- Sparse Row format. 
c
c--------------------------------------------------------------------
c parameters:
c-------------
c Input:
c ------
c nx      = number of points in x direction
c ny	  = number of points in y direction
c nz	  = number of points in z direction
c nfree   = number of degrees of freedom per point
c na	  = first dimension of array a as declared in calling
c           program. Must be .ge. nfree**2
c
c Output: 
c ------ 
c n	  = dimension of matrix (output)
c
c a, ja, ia = resulting matrix in  Block Sparse Row format
c           a(1:nfree**2, j ) contains a nonzero block and ja(j) 
c           contains the (block) column number of this block.
c           the block dimension of the matrix is n (output) and 
c           therefore the total number of (scalar) rows is n x nfree.
c     
c iau     = integer*n containing the position of the diagonal element
c           in the a, ja, ia structure
c
c Work space:
c------------ 
c stencil =  work array of size (7,nfree**2) [stores local stencils]
c
c--------------------------------------------------------------------
c
c     stencil (1:7,*) has the following meaning:
c
c     center point = stencil(1)
c     west point   = stencil(2)
c     east point   = stencil(3)
c     south point  = stencil(4)
c     north point  = stencil(5)
c     front point  = stencil(6)
c     back point   = stencil(7)
c
c
c                           st(5)
c                            |
c                            |
c                            |
c                            |          .st(7)
c                            |     .
c                            | .
c         st(2) ----------- st(1) ---------- st(3)
c                       .    |
c                   .        |
c               .            |
c            st(6)           |
c                            |
c                            |
c                           st(4)
c
c-------------------------------------------------------------------
c     some constants
c
      real*8 one
      parameter (one=1.0D0)
c
c     local variables
c
      integer iedge,ix,iy,iz,k,kx,ky,kz,nfree2,node
      real*8  h
c
      h = one/dble(nx+1)
      kx = 1
      ky = nx
      kz = nx*ny
      nfree2 = nfree*nfree
      iedge = 1
      node = 1
      do 100 iz = 1,nz
         do 90 iy = 1,ny
            do 80 ix = 1,nx
               ia(node) = iedge
               call bsten(nx,ny,nz,ix,iy,iz,nfree,stencil,h)
c     west
               if (ix.gt.1) then
                  ja(iedge)=node-kx
	          do 4 k=1,nfree2
		     a(k,iedge) = stencil(2,k)
 4		  continue
                  iedge=iedge + 1
               end if
c     south
               if (iy.gt.1) then
                  ja(iedge)=node-ky
	          do 5 k=1,nfree2
		     a(k,iedge) = stencil(4,k)
 5		  continue
                  iedge=iedge + 1
               end if
c     front plane
               if (iz.gt.1) then
                  ja(iedge)=node-kz
	          do 6 k=1,nfree2
		     a(k,iedge) = stencil(6,k)
 6		  continue
                  iedge=iedge + 1
               endif
c     center node
               ja(iedge) = node
               iau(node) = iedge
               do 7 k=1,nfree2
                  a(k,iedge) = stencil(1,k)
 7             continue
               iedge = iedge + 1
c     -- upper part
c     east
               if (ix.lt.nx) then
                  ja(iedge)=node+kx
	          do 8 k=1,nfree2
		     a(k,iedge) = stencil(3,k)
 8		  continue
                  iedge=iedge + 1
               end if
c     north
               if (iy.lt.ny) then
                  ja(iedge)=node+ky
	          do 9 k=1,nfree2
		     a(k,iedge) = stencil(5,k)
 9		  continue
                  iedge=iedge + 1
               end if
c     back plane
               if (iz.lt.nz) then
                  ja(iedge)=node+kz
	          do 10 k=1,nfree2
                     a(k,iedge) = stencil(7,k)
 10		  continue
                  iedge=iedge + 1
               end if
c------next node -------------------------
               node=node+1
 80         continue
 90      continue
 100  continue
c     
c     -- new version of BSR -- renumbering removed. 
c     change numbering of nodes so that each ja(k) will contain the
c     actual column number in the original matrix of entry (1,1) of each
c     block (k).
c      do 101 k=1,iedge-1
c         ja(k) = (ja(k)-1)*nfree+1
c 101  continue
c
c      n = (node-1)*nfree
      n = node-1 
      ia(node)=iedge
      return
c--------------end-of-gen57bl-------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine bsten (nx,ny,nz,kx,ky,kz,nfree,stencil,h)
c-----------------------------------------------------------------------
c     This subroutine calcultes the correct block-stencil values for
c     centered difference discretization of the elliptic operator
c     (block version of stencil)
c
c L u = delx( a delx u ) + dely ( b dely u) + delz ( c delz u ) +
c       d delx ( u ) + e dely (u) + f delz( u ) + g u
c
c   For 2-D problems the discretization formula that is used is:
c
c h**2 * Lu == a(i+1/2,j)*{u(i+1,j) - u(i,j)} +
c	       a(i-1/2,j)*{u(i-1,j) - u(i,j)} +
c              b(i,j+1/2)*{u(i,j+1) - u(i,j)} +
c              b(i,j-1/2)*{u(i,j-1) - u(i,j)} +
c              (h/2)*d(i,j)*{u(i+1,j) - u(i-1,j)} +
c              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} +
c              (h/2)*e(i,j)*{u(i,j+1) - u(i,j-1)} +
c              (h**2)*g(i,j)*u(i,j)
c-----------------------------------------------------------------------
c     some constants
c
      real*8  zero,half
      parameter(zero=0.0D0,half=0.5D0)
c
c     local variables
c
      integer i,k,kx,ky,kz,nfree,nfree2,nx,ny,nz
      real*8 stencil(7,*), cntr(225), coeff(225),h,h2,hhalf,x,y,z
c------------
      if (nfree .gt. 15) then
         print *, ' ERROR ** nfree too large '
         stop
      endif
c
      nfree2 = nfree*nfree
      do 200 k=1, nfree2
         cntr(k) = zero
         do 199 i=1,7
            stencil(i,k) = zero
 199     continue
 200  continue
c------------
      hhalf = h*half
      h2 = h*h
      x = h*dble(kx)
      y = h*dble(ky)
      z = h*dble(kz)
c differentiation wrt x:
      call afunbl(nfree,x+hhalf,y,z,coeff)
      do 1 k=1, nfree2
      stencil(3,k) = stencil(3,k) + coeff(k)
      cntr(k) = cntr(k) + coeff(k)
 1    continue
c
      call afunbl(nfree,x-hhalf,y,z,coeff)
      do 2 k=1, nfree2
         stencil(2,k) = stencil(2,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 2    continue
c
      call dfunbl(nfree,x,y,z,coeff)
      do 3 k=1, nfree2
         stencil(3,k) = stencil(3,k) + coeff(k)*hhalf
         stencil(2,k) = stencil(2,k) - coeff(k)*hhalf
 3    continue
      if (ny .le. 1) goto 99
c
c differentiation wrt y:
c
      call bfunbl(nfree,x,y+hhalf,z,coeff)
      do 4 k=1,nfree2
         stencil(5,k) = stencil(5,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 4    continue
c
      call bfunbl(nfree,x,y-hhalf,z,coeff)
      do 5 k=1, nfree2
         stencil(4,k) = stencil(4,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 5    continue
c
      call efunbl(nfree,x,y,z,coeff)
      do 6 k=1, nfree2
         stencil(5,k) = stencil(5,k) + coeff(k)*hhalf
         stencil(4,k) = stencil(4,k) - coeff(k)*hhalf
 6    continue
      if (nz .le. 1) goto 99
c
c differentiation wrt z:
c
      call cfunbl(nfree,x,y,z+hhalf,coeff)
      do 7 k=1, nfree2
         stencil(7,k) = stencil(7,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 7    continue
c
      call cfunbl(nfree,x,y,z-hhalf,coeff)
      do 8 k=1, nfree2
         stencil(6,k) = stencil(6,k) + coeff(k)
         cntr(k) = cntr(k) + coeff(k)
 8    continue
c
      call ffunbl(nfree,x,y,z,coeff)
      do 9 k=1, nfree2
         stencil(7,k) = stencil(7,k) + coeff(k)*hhalf
         stencil(6,k) = stencil(6,k) - coeff(k)*hhalf
 9    continue
c
c discretization of  product by g:
c
 99   call gfunbl(nfree,x,y,z,coeff)
      do 10 k=1, nfree2
         stencil(1,k) = h2*coeff(k) - cntr(k)
 10   continue
c
      return
c------------end of bsten-----------------------------------------------
c-----------------------------------------------------------------------
      end
      subroutine fdreduce(nx,ny,nz,alpha,n,a,ja,ia,iau,rhs,stencil)
      implicit none
      integer nx,ny, nz, n, ia(*), ja(*), iau(*)
      real*8  alpha(*), a(*), rhs(*), stencil(*)
c-----------------------------------------------------------------------
c This subroutine tries to reduce the size of the matrix by looking
c for Dirichlet boundary conditions at each surface and solve the boundary
c value and modify the right-hand side of related nodes, then clapse all
c the boundary nodes.
c-----------------------------------------------------------------------
c     parameters
c
      real*8   zero
      parameter(zero=0.0D0)
c
c     local variables
c
      integer  i,j,k,kx,ky,kz,lx,ux,ly,uy,lz,uz,node,nbnode,lk,ld,iedge
      real*8   val
      integer  lctcsr
      external lctcsr
c
c     The first half of this subroutine will try to change the right-hand
c     side of all the nodes that has a neighbor with Dirichlet boundary
c     condition, since in this case the value of the boundary point is
c     known.
c     Then in the second half, we will try to eliminate the boundary
c     points with known values (with Dirichlet boundary condition).
c
      kx = 1
      ky = nx
      kz = nx*ny
      lx = 1
      ux = nx
      ly = 1
      uy = ny
      lz = 1
      uz = nz
c
c     Here goes the first part. ----------------------------------------
c
c     the left (west) side
c
      if (alpha(1) .eq. zero) then
         lx = 2
         do 10 k = 1, nz
            do 11 j = 1, ny
               node = (k-1)*kz + (j-1)*ky + 1
               nbnode = node + kx
               lk = lctcsr(nbnode, node, ja, ia)
               ld = iau(node)
               val = rhs(node)/a(ld)
c     modify the rhs
               rhs(nbnode) = rhs(nbnode) - a(lk)*val
 11         continue
 10      continue
      endif
c
c     right (east) side
c
      if (alpha(2) .eq. zero) then
         ux = nx - 1
         do 20 k = 1, nz
            do 21 j = 1, ny
               node = (k-1)*kz + (j-1)*ky + nx
               nbnode = node - kx
               lk = lctcsr(nbnode, node, ja, ia)
               ld = iau(node)
               val = rhs(node)/a(ld)
c     modify the rhs
               rhs(nbnode) = rhs(nbnode) - a(lk)*val
 21         continue
 20      continue
      endif
c
c     if it's only 1-D, skip the following part
c
      if (ny .le. 1) goto 100
c
c     the bottom (south) side
c
      if (alpha(3) .eq. zero) then
         ly = 2
         do 30 k = 1, nz
            do 31 i = lx, ux
               node = (k-1)*kz + i
               nbnode = node + ky
               lk = lctcsr(nbnode, node, ja, ia)
               ld = iau(node)
               val = rhs(node)/a(ld)
c     modify the rhs
               rhs(nbnode) = rhs(nbnode) - a(lk)*val
 31         continue
 30      continue
      endif
c
c     top (north) side
c
      if (alpha(4) .eq. zero) then
         uy = ny - 1
         do 40 k = 1, nz
            do 41 i = lx, ux
               node = (k-1)*kz + i + (ny-1)*ky
               nbnode = node - ky
               lk = lctcsr(nbnode, node, ja, ia)
               ld = iau(node)
               val = rhs(node)/a(ld)
c     modify the rhs
               rhs(nbnode) = rhs(nbnode) - a(lk)*val
 41         continue
 40      continue
      endif
c
c     if only 2-D skip the following section on z
c
      if (nz .le. 1) goto 100
c
c     the front surface
c
      if (alpha(5) .eq. zero) then
         lz = 2
         do 50 j = ly, uy
            do 51 i = lx,  ux
               node = (j-1)*ky + i
               nbnode = node + kz
               lk = lctcsr(nbnode, node, ja, ia)
               ld = iau(node)
               val = rhs(node)/a(ld)
c     modify the rhs
               rhs(nbnode) = rhs(nbnode) - a(lk)*val
 51         continue
 50      continue
      endif
c
c     rear surface
c
      if (alpha(6) .eq. zero) then
         uz = nz - 1
         do 60 j = ly, uy
            do 61 i = lx, ux
               node = (nz-1)*kz + (j-1)*ky + i
               nbnode = node - kz
               lk = lctcsr(nbnode, node, ja, ia)
               ld = iau(node)
               val = rhs(node)/a(ld)
c     modify the rhs
               rhs(nbnode) = rhs(nbnode) - a(lk)*val
 61         continue
 60      continue
      endif
c
c     now the second part ----------------------------------------------
c
c     go through all the actual nodes with unknown values, collect all
c     of them to form a new matrix in compressed sparse row format.
c
 100  kx = 1
      ky = ux - lx + 1
      kz = (uy - ly + 1) * ky
      node = 1
      iedge = 1
      do 80 k = lz, uz
         do 81 j = ly, uy
            do 82 i = lx, ux
c
c     the corresponding old node number
               nbnode = ((k-1)*ny + j-1)*nx + i
c
c     copy the row into local stencil, copy is done is the exact
c     same order as the stencil is written into array a
               lk = ia(nbnode)
               if (i.gt.1) then
                  stencil(2) = a(lk)
                  lk = lk + 1
               end if
               if (j.gt.1) then
                  stencil(4) = a(lk)
                  lk = lk + 1
               end if
               if (k.gt.1) then
                  stencil(6) = a(lk)
                  lk = lk + 1
               end if
               stencil(1) = a(lk)
               lk = lk + 1
               if (i.lt.nx) then
                  stencil(3) = a(lk)
                  lk = lk + 1
               endif
               if (j.lt.ny) then
                  stencil(5) = a(lk)
                  lk = lk + 1
               end if
               if (k.lt.nz) stencil(7) = a(lk)
c
c     first the ia pointer -- points to the beginning of each row
               ia(node) = iedge
c
c     move the values from the local stencil to the new matrix
c
c     the neighbor on the left (west)
               if (i.gt.lx) then
                  ja(iedge)=node-kx
                  a(iedge) =stencil(2)
                  iedge=iedge + 1
               end if
c     the neighbor below (south)
               if (j.gt.ly) then
                  ja(iedge)=node-ky
                  a(iedge)=stencil(4)
                  iedge=iedge + 1
               end if
c     the neighbor in the front
               if (k.gt.lz) then
                  ja(iedge)=node-kz
                  a(iedge)=stencil(6)
                  iedge=iedge + 1
               endif
c     center node (itself)
               ja(iedge) = node
               iau(node) = iedge
               a(iedge) = stencil(1)
               iedge = iedge + 1
c     the neighbor to the right (east)
               if (i.lt.ux) then
                  ja(iedge)=node+kx
                  a(iedge)=stencil(3)
                  iedge=iedge + 1
               end if
c     the neighbor above (north)
               if (j.lt.uy) then
                  ja(iedge)=node+ky
                  a(iedge)=stencil(5)
                  iedge=iedge + 1
               end if
c     the neighbor at the back
               if (k.lt.uz) then
                  ja(iedge)=node+kz
                  a(iedge)=stencil(7)
                  iedge=iedge + 1
               end if
c     the right-hand side
               rhs(node) = rhs(nbnode)
c------next node -------------------------
               node=node+1
c
 82         continue
 81      continue
 80   continue
c
      ia(node) = iedge
c
c     the number of nodes in the final matrix is stored in n
c
      n = node - 1
      return
c-----------------------------------------------------------------------
      end
c-----end of fdreduce-----------------------------------------------------
c-----------------------------------------------------------------------
      subroutine fdaddbc(nx,ny,nz,a,ja,ia,iau,rhs,al,h)
      integer nx, ny, nz, ia(nx*ny*nz), ja(7*nx*ny*nz), iau(nx*ny*nz)
      real*8  h, al(6), a(7*nx*ny*nz), rhs(nx*ny*nz)
c-----------------------------------------------------------------------
c This subroutine will add the boundary condition to the linear system
c consutructed without considering the boundary conditions
c
c The Boundary condition is specified in the following form:
c           du
c     alpha -- + beta u = gamma
c           dn
c Alpha is stored in array AL.  The six side of the boundary appares
c in AL in the following order: left(west), right(east), bottom(south),
c top(north), front, back(rear). (see also the illustration in gen57pt)
c Beta and gamma appears as the functions, betfun and gamfun.
c They have the following prototype
c
c real*8 function xxxfun(x, y, z)
c real*8 x, y, z
c
c where x, y, z are vales in the range of [0, 1][0, (ny-1)*h]
c [0, (nz-1)*h]
c
c At the corners or boundary lines, the boundary conditions are applied
c in the follow order:
c 1) if one side is Dirichlet boundary condition, the Dirichlet boundary
c    condition is used;
c 2) if more than one sides are Dirichlet, the Direichlet condition
c    specified for X direction boundary will overwrite the one specified
c    for Y direction boundary which in turn has priority over Z
c     direction boundaries.
c 3) when all sides are non-Dirichlet, the average values are used.
c-----------------------------------------------------------------------
c     some constants
c
      real*8   half,zero,one,two
      parameter(half=0.5D0,zero=0.0D0,one=1.0D0,two=2.0D0)
c
c     local variables
c
      character*2 side
      integer  i,j,k,kx,ky,kz,node,nbr,ly,uy,lx,ux
      real*8   coeff, ctr, hhalf, x, y, z
      real*8   afun, bfun, cfun, dfun, efun, ffun, gfun, hfun
      external afun, bfun, cfun, dfun, efun, ffun, gfun, hfun
      real*8   betfun, gamfun
      integer  lctcsr
      external lctcsr, betfun, gamfun
c
      hhalf = half * h
      kx = 1
      ky = nx
      kz = nx*ny
c
c     In 3-D case, we need to go through all 6 faces one by one. If
c     the actual dimension is lower, test on ny is performed first.
c     If ny is less or equals to 1, then the value of nz is not
c     checked.
c-----
c     the surface on the left (west) side
c     Concentrate on the contribution from the derivatives related to x,
c     The terms with derivative of x was assumed to be:
c
c     a(3/2,j,k)*[u(2,j,k)-u(1,j,k)] + a(1/2,j,k)*[u(0,j,k)-u(1,j,k)] +
c     h*d(1,j,k)*[u(2,j,k)-u(0,j,k)]/2
c
c     But they actually are:
c
c     2*{a(3/2,j,k)*[u(2,j,k)-u(1,j,k)] -
c     h*a(1,j,k)*[beta*u(1,j,k)-gamma]/alpha]} +
c     h*h*d(1,j,k)*[beta*u(1,j,k)-gamma]/alpha
c
c     Therefore, in terms of local stencil the right neighbor of a node
c     should be changed to 2*a(3/2,j,k),
c     The matrix never contains the left neighbor on this border, nothing
c     needs to be done about it.
c     The following terms should be added to the center stencil:
c     -a(3/2,j,k) + a(1/2,j,k) + [h*d(1,j,k)-2*a(1,j,k)]*h*beta/alpha
c
c     And these terms should be added to the corresponding right-hand side
c     [h*d(1,j,k)-2*a(1,j,k)]*h*gamma/alpha
c
c     Obviously, the formula do not apply for the Dirichlet Boundary
c     Condition, where alpha will be zero. In that case, we simply set
c     all the elements in the corresponding row to zero(0), then let
c     the diagonal element be beta, and the right-hand side be gamma.
c     Thus the value of u at that point will be set. Later on point
c     like this will be removed from the matrix, since they are of
c     know value before solving the system.(not done in this subroutine)
c
      x = zero
      side = 'x1'
      do 20 k = 1, nz
         z = (k-1)*h
         do 21 j = 1, ny
            y = (j-1)*h
            node = 1+(j-1)*ky+(k-1)*kz
c
c     check to see if it's Dirichlet Boundary condition here
c
            if (al(1) .eq. zero) then
               call clrow(node, a, ja, ia)
               a(iau(node)) = betfun(side,x,y,z)
               rhs(node) = gamfun(side,x,y,z)
            else
c
c     compute the terms formulated above to modify the matrix.
c
c     the right neighbor is stroed in nbr'th posiiton in the a
               nbr = lctcsr(node, node+kx, ja, ia)
c
               coeff = two*afun(x,y,z)
               ctr = (h*dfun(x,y,z) - coeff)*h/al(1)
               rhs(node) = rhs(node) + ctr * gamfun(side,x,y,z)
               ctr = afun(x-hhalf,y,z) + ctr * betfun(side,x,y,z)
               coeff = afun(x+hhalf,y,z)
               a(iau(node)) = a(iau(node)) - coeff + ctr
               a(nbr) = two*coeff
            end if
 21      continue
 20   continue
c
c     the right (east) side boudary, similarly, the contirbution from
c     the terms containing the derivatives of x were assumed to be
c
c     a(nx+1/2,j,k)*[u(nx+1,j,k)-u(nx,j,k)] +
c     a(nx-1/2,j,k)*[u(nx-1,j,k)-u(nx,j,k)] +
c     d(nx,j,k)*[u(nx+1,j,k)-u(nx-1,j,k)]*h/2
c
c     Actualy they are:
c
c     2*{h*a(nx,j,k)*[gamma-beta*u(nx,j,k)]/alpha +
c     a(nx-1/2,j,k)*[u(nx-1,j,k)-u(nx,j,k)]} +
c     h*h*d(nx,j,k)*[gamma-beta*u(nx,j,k)]/alpha
c
c     The left stencil has to be set to 2*a(nx-1/2,j,k)
c
c     The following terms have to be added to the center stencil:
c
c     -a(nx-1/2,j,k)+a(nx+1/2,j,k)-[2*a(nx,j,k)+h*d(nx,j,k)]*beta/alpha
c
c     The following terms have to be added to the right-hand side:
c
c     -[2*a(nx,j,k)+h*d(nx,j,k)]*h*gamma/alpha
c
      x = one
      side = 'x2'
      do 22 k = 1, nz
         z = (k-1)*h
         do 23 j = 1, ny
            y = (j-1)*h
            node = (k-1)*kz + j*ky
c
            if (al(2) .eq. zero) then
               call clrow(node, a, ja, ia)
               a(iau(node)) = betfun(side,x,y,z)
               rhs(node) = gamfun(side,x,y,z)
            else
               nbr = lctcsr(node, node-kx, ja, ia)
c
               coeff = two*afun(x,y,z)
               ctr = (coeff + h*dfun(x,y,z))*h/al(2)
               rhs(node) = rhs(node) - ctr * gamfun(side,x,y,z)
               ctr = afun(x+hhalf,y,z) - ctr * betfun(side,x,y,z)
               coeff = afun(x-hhalf,y,z)
               a(iau(node)) = a(iau(node)) - coeff + ctr
               a(nbr) = two*coeff
            end if
 23      continue
 22   continue
c
c     If only one dimension, return now
c
      if (ny .le. 1) return
c
c     the bottom (south) side suface, This similar to the situation
c     with the left side, except all the function and realted variation
c     should be on the y.
c
c     These two block if statment here is to resolve the possible conflict
c     of assign the boundary value differently by different side of the
c     Dirichlet Boundary Conditions. They ensure that the edges that have
c     be assigned a specific value will not be reassigned.
c
      if (al(1) .eq. zero) then
         lx = 2
      else
         lx = 1
      end if
      if (al(2) .eq. zero) then
         ux = nx-1
      else
         ux = nx
      end if
      y = zero
      side = 'y1'
      do 24 k = 1, nz
         z = (k-1)*h
         do 25 i = lx, ux
            x = (i-1)*h
            node = i + (k-1)*kz
c
            if (al(3) .eq. zero) then
               call clrow(node, a, ja, ia)
               a(iau(node)) = betfun(side,x,y,z)
               rhs(node) = gamfun(side,x,y,z)
            else
               nbr = lctcsr(node, node+ky, ja, ia)
c
               coeff = two*bfun(x,y,z)
               ctr = (h*efun(x,y,z) - coeff)*h/al(3)
               rhs(node) = rhs(node) + ctr * gamfun(side,x,y,z)
               ctr = bfun(x,y-hhalf,z) + ctr * betfun(side,x,y,z)
               coeff = bfun(x,y+hhalf,z)
               a(iau(node)) = a(iau(node)) - coeff + ctr
               a(nbr) = two*coeff
            end if
 25      continue
 24   continue
c
c     The top (north) side, similar to the right side
c
      y = (ny-1) * h
      side = 'y2'
      do 26 k = 1, nz
         z = (k-1)*h
         do 27 i = lx, ux
            x = (i-1)*h
            node = (k-1)*kz+(ny-1)*ky + i
c
            if (al(4) .eq. zero) then
               call clrow(node, a, ja, ia)
               a(iau(node)) = betfun(side,x,y,z)
               rhs(node) = gamfun(side,x,y,z)
            else
               nbr = lctcsr(node, node-ky, ja, ia)
c
               coeff = two*bfun(x,y,z)
               ctr = (coeff + h*efun(x,y,z))*h/al(4)
               rhs(node) = rhs(node) - ctr * gamfun(side,x,y,z)
               ctr = bfun(x,y+hhalf,z) - ctr * betfun(side,x,y,z)
               coeff = bfun(x,y-hhalf,z)
               a(iau(node)) = a(iau(node)) - coeff + ctr
               a(nbr) = two*coeff
            end if
 27      continue
 26   continue
c
c     If only has two dimesion to work on, return now
c
      if (nz .le. 1) return
c
c     The front side boundary
c
c     If the edges of the surface has been decided by Dirichlet Boundary
c     Condition, then leave them alone.
c
      if (al(3) .eq. zero) then
         ly = 2
      else
         ly = 1
      end if
      if (al(4) .eq. zero) then
         uy = ny-1
      else
         uy = ny
      end if
c
      z = zero
      side = 'z1'
      do 28 j = ly, uy
         y = (j-1)*h
         do 29 i = lx, ux
            x = (i-1)*h
            node = i + (j-1)*ky
c
            if (al(5) .eq. zero) then
               call clrow(node, a, ja, ia)
               a(iau(node)) = betfun(side,x,y,z)
               rhs(node) = gamfun(side,x,y,z)
            else
               nbr = lctcsr(node, node+kz, ja, ia)
c
               coeff = two*cfun(x,y,z)
               ctr = (h*ffun(x,y,z) - coeff)*h/al(5)
               rhs(node) = rhs(node) + ctr * gamfun(side,x,y,z)
               ctr = cfun(x,y,z-hhalf) + ctr * betfun(side,x,y,z)
               coeff = cfun(x,y,z+hhalf)
               a(iau(node)) = a(iau(node)) - coeff + ctr
               a(nbr) = two*coeff
            end if
 29      continue
 28   continue
c
c     Similiarly for the top side of the boundary suface
c
      z = (nz - 1) * h
      side = 'z2'
      do 30 j = ly, uy
         y = (j-1)*h
         do 31 i = lx, ux
            x = (i-1)*h
            node = (nz-1)*kz + (j-1)*ky + i
c
            if (al(6) .eq. zero) then
               call clrow(node, a, ja, ia)
               a(iau(node)) = betfun(side,x,y,z)
               rhs(node) = gamfun(side,x,y,z)
            else
               nbr = lctcsr(node, node-kz, ja, ia)
c
               coeff = two*cfun(x,y,z)
               ctr = (coeff + h*ffun(x,y,z))*h/al(6)
               rhs(node) = rhs(node) - ctr * gamfun(side,x,y,z)
               ctr = cfun(x,y,z+hhalf) - ctr * betfun(side,x,y,z)
               coeff = cfun(x,y,z-hhalf)
               a(iau(node)) = a(iau(node)) - coeff + ctr
               a(nbr) = two*coeff
            end if
 31      continue
 30   continue
c
c     all set
c
      return
c-----------------------------------------------------------------------
      end
c-----end of fdaddbc----------------------------------------------------
c-----------------------------------------------------------------------
      subroutine clrow(i, a, ja, ia)
      integer i, ja(*), ia(*), k
      real *8 a(*)
c-----------------------------------------------------------------------
c     clear the row i to all zero, but still keep the structure of the
c     CSR matrix
c-----------------------------------------------------------------------
      do 10 k = ia(i), ia(i+1)-1
         a(k) = 0.0D0
 10   continue
c
      return
c-----end of clrow------------------------------------------------------
      end
c-----------------------------------------------------------------------
      function lctcsr(i,j,ja,ia)
      integer lctcsr, i, j, ja(*), ia(*), k
c-----------------------------------------------------------------------
c     locate the position of a matrix element in a CSR format
c     returns -1 if the desired element is zero
c-----------------------------------------------------------------------
      lctcsr = -1
      k = ia(i)
 10   if (k .lt. ia(i+1) .and. (lctcsr .eq. -1)) then
         if (ja(k) .eq. j) lctcsr = k
         k = k + 1
         goto 10
      end if
c
      return
c-----------------------------------------------------------------------
      end
c-----end of lctcsr-----------------------------------------------------


