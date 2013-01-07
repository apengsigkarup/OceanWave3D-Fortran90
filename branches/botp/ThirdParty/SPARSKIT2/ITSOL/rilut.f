      program rilut
c-----------------------------------------------------------------------
c     test program for ilut preconditioned gmres.
c     this program generates a sparse matrix using
c     matgen and then solves a linear system with an 
c     artificial rhs.
c-----------------------------------------------------------------------
      implicit none 
c 
      integer nmax, nzmax 
      parameter (nmax=5000,nzmax=100000)
      integer ia(nmax),ja(nzmax),jau(nzmax),ju(nzmax),iw(nmax*3),
     &     iperm(nmax*2),ipar(16), levs(nzmax) 
      real*8  a(nzmax),x(nmax),y(nmax),au(nzmax),vv(nmax,20),
     *     xran(nmax),rhs(nmax),al(nmax),fpar(16)
c     real t(2), t1, etime
c
      integer nx,ny,nz,n,j,k,ierr,meth,lfil,nwk,im,maxits,iout
      real*8 tol,permtol,eps,alph,gammax,gammay,alpha 
      external gmres
c     
      common /func/ gammax, gammay, alpha
c-----------------------------------------------------------------------
c     pde to be discretized is :
c---------------------------
c     
c     -Lap u + gammax exp (xy)delx u + gammay exp (-xy) dely u +alpha u
c     
c     where Lap = 2-D laplacean, delx = part. der. wrt x,
c     dely = part. der. wrt y.
c     gammax, gammay, and alpha are passed via the commun func.
c     
c-----------------------------------------------------------------------
c     
c     data for PDE:
c     
      nx = 30 
      ny = 30 
      nz = 1
      alpha = -50.0
      gammax = 10.0
      gammay = 10.0
c     
c     data for preconditioner
c     
      nwk = nzmax 
c     
c     data for GMRES
c     
      im   = 10
      eps  = 1.0D-07
      maxits = 100 
      iout = 6
      permtol = 1.0
      ipar(2) = 2
      ipar(3) = 2
      ipar(4) = 20*nmax
      ipar(5) = im
      ipar(6) = maxits
      fpar(1) = eps
      fpar(2) = 2.22D-16
c     
c     same initial guess for gmres 
c     
c--------------------------------------------------------------
c     call gen57 to generate matrix in compressed sparse row format
c--------------------------------------------------------------
c     
c     define part of the boundary condition here
c     
      al(1) = 0.0
      al(2) = 1.0
      al(3) = 0.0
      al(4) = 0.0
      al(5) = 0.0
      al(6) = 0.0
      call gen57pt(nx,ny,nz,al,0,n,a,ja,ia,ju,rhs) 
c
c     zero initial guess to the iterative solvers
c
      do j=1, n
         xran(j) = 0.d0
      enddo
      print *, 'RILUT:  generated a finite difference matrix'
      print *, '        grid size = ', nx, ' X ', ny, ' X ', nz
      print *, '        matrix size = ', n
c--------------------------------------------------------------
c     gnerate right han side = A * (1,1,1,...,1)**T
c--------------------------------------------------------------
      do  k=1,n
         x(k) = 1.0
      enddo
      call amux(n, x, y, a, ja, ia)
c--------------------------------------------------------------
c     test all different methods available:
c     ILU0, MILU0, ILUT and with different values of tol and lfil
c     ( from cheaper to more expensive preconditioners)
c     The more accurate the preconditioner the fewer iterations 
c     are required in pgmres, in general. 
c     
c--------------------------------------------------------------
      do 200 meth = 1, 15 
         goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) meth
 1       continue
         write (iout,*) ' +++++ ILU(0) Preconditioner ++++ '
c        t1 = etime(t)    
         call ilu0 (n, a, ja, ia, au, jau, ju, iw, ierr)
c        t1 = etime(t) - t1 
         goto 100
 2       continue
         write (iout,*) ' +++++ MILU(0) Preconditioner ++++ '
c        t1 = etime(t)    
         call milu0 (n, a, ja, ia, au, jau, ju, iw, ierr)
c        t1 = etime(t) - t1 
         goto 100
 3       continue
         write (iout,*) ' +++++ ILUT Preconditioner ++++ '
         write (iout,*) ' +++++ tol = 0.0001, lfil=5 ++++ '
         tol  = 0.0001
         lfil = 5
c     
c        t1 = etime(t)       
         call ilut (n,a,ja,ia,lfil,tol,au,jau,ju,nwk,vv,iw,ierr)
c         t1 = etime(t) - t1 
         goto 100
 4       continue
         write (iout,*) ' +++++ ILUT Preconditioner ++++ '
         write (iout,*) ' +++++ tol = 0.0001, lfil=10 ++++ '
         tol = 0.0001
         lfil = 10
c     
c        t1 = etime(t)       
         call ilut (n,a,ja,ia,lfil,tol,au,jau,ju,nwk,vv,iw,ierr)

c        t1 = etime(t) - t1 
         goto 100
 5       continue
         write (iout,*) ' +++++ ILUT Preconditioner ++++ '
         write (iout,*) ' +++++ tol = .0001, lfil=15 ++++ '
         tol = 0.0001 
         lfil = 15 
c     
c        t1 = etime(t)    
         call ilut (n,a,ja,ia,lfil,tol,au,jau,ju,nwk,vv,iw,ierr)
c        t1 = etime(t) - t1 
         goto 100
 6       continue
         write (iout,*) ' +++++ ILUTP Preconditioner ++++ '
         write (iout,*) ' +++++ tol = 0.0001, lfil=5 ++++ '
         tol  = 0.0001
         lfil = 5
c     
c        t1 = etime(t)
         call ilutp(n,a,ja,ia,lfil,tol,permtol,n,au,jau,ju,nwk,
     *        vv,iw,iperm,ierr)
c        t1 = etime(t) - t1 
         goto 100
 7       continue
         write (iout,*) ' +++++ ILUTP Preconditioner ++++ '
         write (iout,*) ' +++++ tol = 0.0001, lfil=10 ++++ '
         tol = 0.0001
         lfil = 10
c     
c        t1 = etime(t)
         call ilutp(n,a,ja,ia,lfil,tol,permtol,n,au,jau,ju,nwk,
     *        vv,iw,iperm,ierr)
c        t1 = etime(t) - t1 
         goto 100
c-----------------------------------------------------------------------------
 8       continue
         write (iout,*) ' +++++ ILUTP Preconditioner ++++ '
         write (iout,*) ' +++++ tol = .0001, lfil=15 ++++ '
         tol = 0.0001
         lfil = 15
c     
c        t1 = etime(t)
         call ilutp(n,a,ja,ia,lfil,tol,permtol,n,au,jau,ju,nwk,
     *        vv,iw,iperm,ierr)
c        t1 = etime(t) - t1 
         goto 100 
 9      continue
         write (iout,*) ' +++++ ILUK Preconditioner ++++ '
         write (iout,*) ' +++++       lfil=0        ++++ '
         lfil = 0 
c        t1 = etime(t)
         call iluk(n,a,ja,ia,lfil,au,jau,ju,levs,nwk,vv,iw,ierr)
         print *, ' nnz for a =', ia(n+1) - ia(1) 
         print *, ' nnz for ilu =', jau(n+1) -jau(1) + n 
c        t1 = etime(t) - t1 
         goto 100 
c
 10      continue
         write (iout,*) ' +++++ ILUK Preconditioner ++++ '
         write (iout,*) ' +++++       lfil=1        ++++ '
         lfil = 1 
c        t1 = etime(t)
         call iluk(n,a,ja,ia,lfil,au,jau,ju,levs,nwk,vv,iw,ierr)
         print *, ' nnz for a =', ia(n+1) - ia(1) 
         print *, ' nnz for ilu =', jau(n+1) -jau(1) + n 
c        t1 = etime(t) - t1 
         goto 100 
c
 11      continue
         write (iout,*) ' +++++ ILUK Preconditioner ++++ '
         write (iout,*) ' +++++       lfil=3        ++++ '
         lfil = 3 
c        t1 = etime(t)
         call iluk(n,a,ja,ia,lfil,au,jau,ju,levs,nwk,vv,iw,ierr)
         print *, ' nnz for a =', ia(n+1) - ia(1) 
         print *, ' nnz for ilu =', jau(n+1) -jau(1) + n 
c        t1 = etime(t) - t1 
         goto 100 
c
 12      continue
         write (iout,*) ' +++++ ILUK Preconditioner ++++ '
         write (iout,*) ' +++++       lfil=6        ++++ '
         lfil = 6
c        t1 = etime(t)
         call iluk(n,a,ja,ia,lfil,au,jau,ju,levs,nwk,vv,iw,ierr)
         print *, ' nnz for a =', ia(n+1) - ia(1) 
         print *, ' nnz for ilu =', jau(n+1) -jau(1) + n 
c        t1 = etime(t) - t1 
         goto 100 
c
 13      continue
c----------------------------------------------------------------------- 
        write (iout,*) ' +++++ ILUD Preconditioner ++++ '
        write (iout,*) ' +++++ tol=0.075, alpha=0.0 ++++ '
        tol = 0.075 
        alph= 0.0 
c        t1 = etime(t)
        call ilud(n,a,ja,ia,alph,tol,au,jau,ju,nwk,vv,iw,ierr)
c
        print *, ' nnz for a =', ia(n+1) - ia(1) 
        print *, ' nnz for ilu =', jau(n+1) -jau(1) + n 
c        t1 = etime(t) - t1 
        goto 100 
c
 14      continue
        write (iout,*) ' +++++ ILUD Preconditioner ++++ '
        write (iout,*) ' +++++ tol=0.075, alpha=1.0 ++++ '
        tol = 0.075 
        alph=1.0 
c        t1 = etime(t)
        call ilud(n,a,ja,ia,alph,tol,au,jau,ju,nwk,vv,iw,ierr)
        print *, ' nnz for a =', ia(n+1) - ia(1) 
        print *, ' nnz for ilu =', jau(n+1) -jau(1) + n 
c        t1 = etime(t) - t1 
        goto 100 
 
 15      continue
        write (iout,*) ' +++++ ILUD Preconditioner ++++ '
        write (iout,*) ' +++++ tol=0.01, alpha=1.0 ++++ '
        tol = 0.01
c        t1 = etime(t)
        call ilud(n,a,ja,ia,alph,tol,au,jau,ju,nwk,vv,iw,ierr)
        print *, ' nnz for a =', ia(n+1) - ia(1) 
        print *, ' nnz for ilu =', jau(n+1) -jau(1) + n 
c        t1 = etime(t) - t1 
c        goto 100 
c
 100     continue
c     
c     check that return was succesful
c     
c        print *, ' ILU factorization time ', t1 
         print *, ' Precon set-up returned with ierr ', ierr
         if (ierr .ne. 0) goto 200
c--------------------------------------------------------------
c     call GMRES
c--------------------------------------------------------------
         call runrc(n,y,x,ipar,fpar,vv,xran,a,ja,ia,au,jau,ju,
     +        gmres)
         print *, 'GMRES return status = ', ipar(1)
         write (iout,*) '     ' 
 200  continue
c     
c----------------------------------------------------------
c     
      write (iout,*) ' **** SOLUTION ****  '
      write(iout, 111) (x(k),k=1,n) 
 111  format (5d15.5)
      stop
      end

