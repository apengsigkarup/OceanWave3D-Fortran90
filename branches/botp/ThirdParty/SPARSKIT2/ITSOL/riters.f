      program riters
c-----------------------------------------------------------------------
c test program for iters -- the basic iterative solvers
c
c     this program generates a sparse matrix using
c     GEN57PT and then solves a linear system with an 
c     artificial rhs (the solution is a vector of (1,1,...,1)^T).
c-----------------------------------------------------------------------
c      implicit none
c      implicit real*8 (a-h,o-z)
      integer nmax, nzmax, maxits,lwk
      parameter (nmax=5000,nzmax=100000,maxits=60,lwk=nmax*40)
      integer ia(nmax),ja(nzmax),jau(nzmax),ju(nzmax),iw(nmax*3)
      integer ipar(16),nx,ny,nz,i,lfil,nwk,nrow,ierr
      real*8  a(nzmax),sol(nmax),rhs(nmax),au(nzmax),wk(nmax*40)
      real*8  xran(nmax), fpar(16), al(nmax)
      real*8  gammax,gammay,alpha,tol
      external gen57pt,cg,bcg,dbcg,bcgstab,tfqmr,gmres,fgmres,dqgmres
      external cgnr, fom, runrc, ilut
c     
      common /func/ gammax, gammay, alpha
c-----------------------------------------------------------------------  
c pde to be discretized is :
c---------------------------
c
c -Lap u + gammax exp (xy)delx u + gammay exp (-xy) dely u +alpha u
c
c where Lap = 2-D laplacean, delx = part. der. wrt x,
c dely = part. der. wrt y.
c gammax, gammay, and alpha are passed via the commun func.
c 
c-----------------------------------------------------------------------  
c
c data for PDE:
c
      nx = 6
      ny = 6
      nz = 1
      alpha = 0.0
      gammax = 0.0
      gammay = 0.0
c
c     set the parameters for the iterative solvers
c
      ipar(2) = 2
      ipar(3) = 1
      ipar(4) = lwk
      ipar(5) = 10 
      ipar(6) = maxits
      fpar(1) = 1.0D-5
      fpar(2) = 1.0D-10
c--------------------------------------------------------------
c call GEN57PT to generate matrix in compressed sparse row format
c
c     al(1:6) are used to store part of the boundary conditions
c     (see documentation on GEN57PT.)
c--------------------------------------------------------------
      al(1) = 0.0
      al(2) = 0.0
      al(3) = 0.0
      al(4) = 0.0
      al(5) = 0.0
      al(6) = 0.0
      nrow = nx * ny * nz
      call gen57pt(nx,ny,nz,al,0,nrow,a,ja,ia,ju,rhs)
      print *, 'RITERS: generated a finite difference matrix'
      print *, '        grid size = ', nx, ' X ', ny, ' X ', nz
      print *, '        matrix size = ', nrow
c
c     set-up the preconditioner ILUT(15, 1E-4)  ! new definition of lfil
c
      lfil = 3
      tol = 1.0D-4
      nwk = nzmax
      call ilut (nrow,a,ja,ia,lfil,tol,au,jau,ju,nwk,
     *     wk,iw,ierr)
      ipar(2) = 2
c
c     generate a linear system with known solution
c
      do i = 1, nrow
         sol(i) = 1.0D0
         xran(i) = 0.d0
      end do
      call amux(nrow, sol, rhs, a, ja, ia)
      print *, ' '
      print *, '	*** CG ***'
      call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
     +     cg)
      print *, ' '
      print *, '	*** BCG ***'
      call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
     +     bcg)
      print *, ' '
      print *, '	*** DBCG ***'
      call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
     +     dbcg)
      print *, ' '
      print *, '	*** CGNR ***'
      call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
     +     cgnr)
      print *, ' '
      print *, '	*** BCGSTAB ***'
      call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
     +     bcgstab)
      print *, ' '
      print *, '	*** TFQMR ***'
      call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
     +     tfqmr)
      print *, ' '
      print *, '	*** FOM ***'
      call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
     +     fom)
      print *, ' '
      print *, '	*** GMRES ***'
      call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
     +     gmres)
      print *, ' '
      print *, '	*** FGMRES ***'
      call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
     +     fgmres)
      print *, ' '
      print *, '	*** DQGMRES ***'
      call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju,
     +     dqgmres)
      stop
      end
c-----end-of-main
c-----------------------------------------------------------------------
