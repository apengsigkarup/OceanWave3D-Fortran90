      program markov 
c-----------------------------------------------------------------------
c
c program to generate a Markov chain matrix (to test eigenvalue routines
c or algorithms for singular systems (in which case use I-A )) 
c the matrix models  simple random walk on a triangular grid.
c see additional comments in subroutine.
c -----
c just compile this segment and link to the rest of sparskit 
c (uses subroutine prtmt from MATGEN)
c will create a matrix in the HARWELL/BOEING format and put it in
c the file markov.mat
c
c-----------------------------------------------------------------------
      parameter (nmax=5000, nzmax= 4*nmax) 
      real*8 a(nzmax) 
      integer ja(nzmax), ia(nmax+1)
c 
      character title*72,key*8,type*3              
      open (unit=11,file='markov.mat') 
c
c read - in grid size - will not accept too large grids.
c
      write (6,'(17hEnter grid-size: ,$)') 
      read *, m
      if (m*(m+1) .gt. 2*nmax ) then
         print *, ' m too large - unable to produce matrix '
         stop
      endif
c 
c call generator. 
c
      call markgen (m, n, a, ja, ia)
c-----------------------------------------------------------------------
      title=' Test matrix from SPARSKIT - markov chain model           '
      key = 'randwk01'
      type = 'rua'
      iout = 11
      job = 2
      ifmt = 10
      call prtmt (n, n, a, ja, ia, x,'NN',title,
     *     key, type,ifmt, job, iout) 
      stop
      end
c     
      subroutine markgen (m, n, a, ja, ia) 
c-----------------------------------------------------------------------
c matrix generator for a markov model of a random walk on a triang. grid
c-----------------------------------------------------------------------
c this subroutine generates a test matrix that models a random
c walk on a triangular grid. This test example was used by 
c G. W. Stewart ["{SRRIT} - a FORTRAN subroutine to calculate the 
c dominant invariant subspaces of a real matrix", 
c Tech. report. TR-514, University of Maryland (1978).] and in a few
c papers on eigenvalue problems by Y. Saad [see e.g. LAA, vol. 34,
c pp. 269-295 (1980) ]. These matrices provide reasonably easy 
c test problems for eigenvalue algorithms. The transpose of the 
c matrix  is stochastic and so it is known that one is an exact 
c eigenvalue. One seeks the eigenvector of the transpose associated 
c with the eigenvalue unity. The problem is to calculate the 
c steady state probability distribution of the system, which is
c the eigevector associated with the eigenvalue one and scaled in
c such a way that the sum all the components is equal to one.
c-----------------------------------------------------------------------
c parameters
c------------
c on entry :
c----------
c m     = integer. number of points in each direction. 
c
c on return:
c----------
c n     = integer. The dimension of the matrix. (In fact n is known 
c         to be equal to (m(m+1))/2      ) 
c a, 
c ja,
c ia    = the matrix stored in CSR format.
c
c-----------------------------------------------------------------------
c Notes: 1) the code will actually compute the transpose of the 
c stochastic matrix that contains the transition probibilities.
c        2) It should also be possible to have a matrix generator
c with an additional parameter (basically redefining `half' below 
c to be another parameter and changing the rest accordingly, but 
c this is not as simple as it sounds). This is not likely to provide 
c any more interesting matrices. 
c-----------------------------------------------------------------------
      real*8 a(*), cst, pd, pu, half 
      integer ja(*), ia(*) 
c----------------------------------------------------------------------- 
      data half/0.5d0/
c
      cst = half/real(m-1)
c     
c     --- ix counts the grid point (natural ordering used), i.e.,
c     --- the row number of the matrix.
c     
      ix = 0
      jax = 1
      ia(1) = jax 
c
c     sweep y coordinates
c
      do 20 i=1,m
         jmax = m-i+1
c
c     sweep x coordinates
c
         do 10 j=1,jmax
            ix = ix + 1
            if (j .eq. jmax) goto 2
            pd = cst*real(i+j-1) 
c            
c     north
c
            a(jax) = pd
            if (i.eq. 1) a(jax) = a(jax)+pd
            ja(jax) =  ix + 1
            jax = jax+1
c     east
            a(jax) = pd
            if (j .eq. 1) a(jax) = a(jax)+pd
            ja(jax) = ix + jmax
            jax = jax+1
c     south
 2          pu = half - cst*real(i+j-3) 
            if ( j .gt. 1) then 
               a(jax) = pu
               ja(jax) = ix-1 
               jax = jax+1
            endif 
c     west
            if ( i .gt. 1) then
               a(jax) = pu
               ja(jax) = ix - jmax - 1 
               jax = jax+1
            endif 
            ia(ix+1) = jax 
 10      continue
 20   continue
      n = ix 
      return
      end
