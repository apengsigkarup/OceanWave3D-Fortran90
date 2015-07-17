      program runall 
c-----------------------------------------------------------------------
c     not present: CG and FGMRES 
c-----------------------------------------------------------------------
c     program for running all examples -- for book tables.
c-----------------------------------------------------------------------
      parameter (nmax=10000,np1=nmax+1,nt2=2*nmax,nzmax = 100000,
     *     lw=31*nmax,nzumax = 20*nzmax) 
c-----------------------------------------------------------------------
      implicit none 
c
      integer ia(nt2),ja(nzmax),jlu(nzumax),levs(nzumax),ju(nmax+1) 
      real*8 wk(lw),rhs(nmax),a(nzmax),alu(nzumax)
      integer jw(nt2), iperm(nt2) 
      real*8 w(nt2),sol(nt2) 
      integer n, mbloc 
      character guesol*2, title*72, key*8, type*3 
c-----------------------------------------------------------------------
      integer iout,lfil,j,k,maxits,i,ncol,numat,mat,iwk,nnz,ipo,ipi,iin,
     *   job,nrhs,im,outf,ipre,its,ierr,ipre1,ipre2,lf,it,itp,
     *     icode,imat,kk,lfinc, niter  
      character*50 filnam
      integer iters(20) 
c
c     up to 10 matrices 
c
      real*8 droptol, dropinc, t, permtol, tol, eps, rand, alph, dnrm2 
c
c-----------------------------------------------------------------------
      data iin/4/, iout/8/ 
c-----------------------------------------------------------------------
      iwk = nzumax 
c
c read in matrix pathname from file. == then read matrix itself 
c
      open (2, file ='matfile')
      read (2,*)  numat
      i=20
c-----------------------------------------------------------------------
      open (iin,file='inputs')
      read (iin,*) alph 
      read (iin,*) im, maxits 
      read (iin,*) eps 
      read (iin,*) ipre1, ipre2, niter 
      read (iin,*) tol, dropinc, permtol 
      read (iin,*) lf, lfinc 
      outf = 30
c     
c     BIG MATRIX LOOP
c
c-----------------------------------------------------------------------
c     LOOP THROUGH MATRICES 
c-----------------------------------------------------------------------
      do 300 mat =1, numat
         read (2,'(a50)') filnam 
         imat = 20 
         open (imat,file =filnam) 
         write (7,*) filnam 
         job = 3
         nrhs = nmax 
         call readmt (nmax,nzmax,job,imat,a,ja,ia, rhs, nrhs,
     *        guesol,n,ncol,nnz,title,key,type,ierr)
c-----------------------------------------------------------------------
         rewind(imat)
         write (7,*) ' -- read -- returned with  ierr = ', ierr 
         write (7,*) ' matrix ', filnam 
         write (7,*) ' *** n  =  ', n, ' nnz ', nnz   
         if (ierr .ne. 0) stop ' not able to read *** ' 
c     
c     sort matrix 
c     
         call csort(n,a,ja,ia,jlu,.true.) 
         job = 0
         call roscal(n,job,1,a,ja,ia,wk,a,ja,ia,ierr)
         call coscal(n,job,1,a,ja,ia,wk,a,ja,ia,ierr)
c-----------------------------------------------------------------------
c     LOOP THROUGH PRECONDITONERS
c-----------------------------------------------------------------------
         droptol = tol 
         lfil = lf  
         kk = 0 
         do 500 ipre = ipre1, ipre2 
c----------------------------------------------------------------------- 
c     ipre = 1 --> ILUD
c     ipre = 2 --> ILUT
c     ipre = 3 --> ILUDP
c     ipre = 4 --> ILUTP
c     ipre = 5 --> ILUK 
c-----------------------------------------------------------------------
            lfil = lf
            droptol = tol 
            do 350 itp =1, niter 
               do k=1,n
                  sol(k) = 1.0 
               enddo
c     
               call amux (n, sol, rhs, a, ja, ia)
c     
               goto (1,2,3,4) ipre 
c     
 1             continue 
               call ilud (n,a,ja,ia,alph,droptol,alu,jlu,ju,iwk,
     *              w,jw,ierr)
               if (ierr .ne. 0) write (7,*) ' ilud ierr = ', ierr 
               call stb_test(n,rhs,alu,jlu,ju,t)
               write (7,*) ' ILUD-drpt =', droptol, ' stbtest = ', t
               write (7,*) '  --- fill-in =', jlu(n+1)-jlu(1) 
c     
               goto 100
 2             continue
c     
               call ilutn (n,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,
     *              ierr)
               if (ierr .ne. 0) write (7,*) ' ILUT ierr = ', ierr 
               call stb_test(n,rhs,alu,jlu,ju,t)
               write (7,*) ' ILUT-drpt =', droptol, 
     *              ' lfil = ', lfil, ' stbtest = ', t 
               write (7,*) '  --- fill-in =', jlu(n+1)-jlu(1) 
c     
               goto 100
c
c     pivoting codes
c 
 3             continue
               mbloc = n 
               call iludp (n,a,ja,ia,alph,droptol,permtol,mbloc,alu,
     *              jlu,ju,iwk,w,jw,iperm,ierr)
               if (ierr .ne. 0) write (7,*) ' ilud ierr = ', ierr 
               call stb_test(n,rhs,alu,jlu,ju,t)
               write (7,*) ' ILUD-drpt =', droptol, ' stbtest = ', t
               write (7,*) '  --- fill-in =', jlu(n+1)-jlu(1) 
               goto 100 
 4             continue
               mbloc = n
               call ilutpn (n,a,ja,ia,lfil,droptol,permtol,mbloc,alu,
     *              jlu,ju,iwk,w,jw,iperm,ierr) 
               if (ierr .ne. 0) write (7,*) ' milutp ierr = ', ierr 
c     
               call stb_test(n,rhs,alu,jlu,ju,t)
               write (7,*) ' ILUTP-drpt =', droptol, 
     *              ' lfil = ', lfil, ' stbtest = ', t 
               write (7,*) '  --- fill-in =', jlu(n+1)-jlu(1) 

               goto 100
 5             continue
c     
               call ilukn(n,a,ja,ia,lfil,alu,jlu,ju,levs,iwk,w,jw,ierr)
               if (ierr .ne. 0) write (7,*) ' ILUT ierr = ', ierr 
               call stb_test(n,rhs,alu,jlu,ju,t)
               write (7,*) ' ILUK -- lfil = ', lfil,' stbtest = ', t 
               write (7,*) '  --- fill-in =', jlu(n+1)-jlu(1) 
c               goto 100
c     ----------- 
c     big loop --
c     -----------
 100        continue 
            
            do 11 k=1,n
               sol(k) = 1.0 
 11         continue
            
            call amux (n, sol, rhs, a, ja, ia)
            do 12 j=1, n
               sol(j) = rand(j) 
 12         continue       
            its = 0 
c-----------------------------------------------------------------------
 10         continue 
c     
            call fgmr (n,im,rhs,sol,it,wk,ipo,ipi,eps,maxits,0,icode)
c-----------------------------------------------------------------------
            if (icode.eq.1) then
               call amux(n, wk(ipo), wk(ipi), a, ja, ia)
               its =  its+1 
               goto 10
            else if (icode.eq.3 .or. icode.eq.5) then
c     
c     output the residuals here ...
c     
               if (ipre .eq. 0) then 
c     NOPRE
                  call dcopy(n, wk(ipo),1,wk(ipi),1)
               else
                  call lusol(n,wk(ipo),wk(ipi),alu,jlu,ju)
               endif
               goto 10
            endif
c     
c     done *** 
c     
            call amux(n,sol,wk,a,ja,ia)
c     
            do i = 1, n
               wk(n+i) = sol(i) -1.0D0
               wk(i) = wk(i) - rhs(i)
            end do
c     
            kk = kk+1 
            iters(kk) = its
c     res(kk) = dnrm2(n,wk,1) 
c     err(kk) = dnrm2(n,wk(n+1),1) 
            write (iout, *) '# actual residual norm ', dnrm2(n,wk,1)
            write (iout, *) '# error norm ', dnrm2(n,wk(1+n),1)
c     
            droptol = droptol*dropinc
            lfil = lfil + lfinc 
c
c     pivoting codes
c 
            if (ipre .ge. 3)  then
               do k=1, ia(n+1)-1 
                  ja(k) = iperm(ja(k)) 
               enddo
            endif
 350        continue
c     
 500     continue
         call entline(outf,filnam,iters,kk) 
         close (imat) 
 300  continue 
      stop
c-------------end-of-main-program-ilut_solve_new------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      function distdot(n,x,ix,y,iy)
      integer n, ix, iy
      real*8 distdot, x(*), y(*), ddot
      external ddot
      distdot = ddot(n,x,ix,y,iy)
      return
c-----end-of-distdot
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine entline(outf,mat,its,kk) 
      implicit none 
      integer outf, kk,its(kk),k
      character mat*70
c      real*8 err(kk), res(kk) 
c----------------------------------------------------------------------- 
      write(outf,100)  mat(19:29),(its(k),k=1,8)
 100  format
     * (a,' & ', 8(i4,' & '),'\\\\ \\hline') 
      return
      end 
c-----------------------------------------------------------------------
      subroutine stb_test(n,sol,alu,jlu,ju,tmax)      
      implicit none
      integer n, jlu(*),ju(*) 
      real*8 sol(n),alu(*),tmax
c----------------------------------------------------------------------- 
      real*8 max, t, abs
      integer j
c
      call lusol(n,sol,sol,alu,jlu,ju)
c      
      tmax = 0.0 
      do j=1, n
         t = abs(sol(j)) 
         tmax = max(t,tmax) 
      enddo
      return
      end 
