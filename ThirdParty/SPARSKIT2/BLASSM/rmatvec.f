      program rmatvec
      parameter (nmax=10000, nzmax=80000) 
      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------------
c This test program tests all the subroutines in matvec.
c it generates matrices and transforms them in appropriate formats
c and then call the appropriate routines. 
c-----------------------------------------------------------------------
      integer ia1(nmax), ia2(nmax), ja1(nzmax), ja2(nzmax), 
     *     jad(nzmax), iwk1(nmax), iwk2(nmax), idim(11), ioff(10) 
      real*8 a1(nzmax), a2(nzmax), 
     *     x(nmax), y(nmax), y0(nmax), y1(nmax), stencil(100) 
c common used only to generate nonsymmetric matrices      
c      common /gam/ gamma, gamma1, cvar
      data idim /4, 10, 15, 40, 50, 60, 70, 80, 90, 100, 200 /
      data iout /6/
c     
c     initialize common gam 
c     
c      gamma = 0.5
c      gamma1 = 1.0
c      cvar = 1.0
c-----------------------------------------------------------------------
c     ii loop corresponds to size of problem 
c-----------------------------------------------------------------------
      do 100 ii = 1, 3 
         write (iout,*) '---------------- ii ',ii,'--------------------' 
         nfree = 1
         nx = idim(ii) 
         ny = nx
c-----------------------------------------------------------------------
c     jj loop corresponds to 2-D and 3-D problems.
c-----------------------------------------------------------------------
         do 150 jj=1, 2
         write (iout,*) '     ----------- jj ',jj,' -------------' 
            nz = 1
            if (jj .eq. 2) nz = 10
c     
c     call matrix generation routine --
c     (strange to use block version to generate 1 x 1 blocks...)
c     
            call gen57bl (nx,ny,nz,1,1,n,a1,ja1,ia1,ia2,stencil)
c     
c     initialize x
c     
            do 1 j=1, n
               x(j) = real(j)
 1          continue
c     
c initial call to get `` exact '' answer in y0
c
            call amux(n,x,y0, a1, ja1, ia1) 
c-----------------------------------------------------------------------
c TESTING AMUXE
c-----------------------------------------------------------------------
c     
c     convert to itpack format -----
c     
            call csrell (n,a1,ja1,ia1,7,a2,jad,n,ndiag,ierr)
            call amuxe (n, x, y, n, ndiag, a2,jad)
            call errpr (n, y, y0,iout,'amuxe ') 
c-----------------------------------------------------------------------
c TESTING AMUXD
c-----------------------------------------------------------------------
c     
c     convert to diagonal format
c
            idiag = 7
            call csrdia (n, idiag,10,a1, ja1, ia1, nmax, a2,
     *           ioff, a2, ja2, ia2, jad) 
            call amuxd (n,x,y,a2,nmax,idiag,ioff) 
            call errpr (n, y, y0,iout,'amuxd ')
c-----------------------------------------------------------------------
c TESTING ATMUX
c-----------------------------------------------------------------------
c
c    convert to csc format (transpose)
c    
            call csrcsc (n,1,1,a1,ja1,ia1,a2,ja2,ia2)
            call atmux (n, x, y, a2, ja2, ia2)
            call errpr (n, y, y0,iout,'atmux ')
c-----------------------------------------------------------------------
c TESTING AMUXJ
c-----------------------------------------------------------------------
c     
c     convert to jagged diagonal format
c     
            call csrjad (n,a1,ja1,ia1, jdiag, jad, a2, ja2, ia2) 
            call amuxj (n, x, y, jdiag, a2, ja2, ia2) 
            call dvperm (n, y, jad) 
            call errpr (n, y, y0,iout,'amuxj ') 
c
c convert back

            call jadcsr (n, jdiag, a2, ja2, ia2, jad, a1, ja1, ia1)
            call amux (n, x, y, a1, ja1, ia1)
            call errpr (n, y, y0,iout,'jadcsr')
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c       triangular systems solutions
c-----------------------------------------------------------------------
c TESTING LDSOL
c-----------------------------------------------------------------------
            call getl (n, a1, ja1, ia1, a2, ja2, ia2)
            call amux (n,x,y0, a2, ja2, ia2)
            call atmux(n,x,y1, a2, ja2, ia2)
            call csrmsr (n, a2, ja2, ia2, a2, ja2, y, iwk2) 
           do 2 k=1,n
               a2(k) = 1.0d0/ a2(k)
 2        continue
            call ldsol (n, y, y0, a2, ja2)
            call errpr (n, x, y, iout,'ldsol ') 
c-----------------------------------------------------------------------
c TESTING LDSOLL
c----------------------------------------------------------------------- 
            call levels (n, ja2, ja2, nlev, jad, iwk1, iwk2) 
            call ldsoll (n, y, y0, a2, ja2, nlev, jad, iwk1) 
            call errpr (n, x, y, iout,'ldsoll') 
c-----------------------------------------------------------------------
c TESTING UDSOLC
c-----------------------------------------------------------------------
c here we take advantage of the fact that the MSR format for U
c is the MSC format for L
c 
            call udsolc (n, y, y1, a2, ja2)
            call errpr (n, x, y, iout,'udsolc') 
c-----------------------------------------------------------------------
c TESTING LSOL 
c-----------------------------------------------------------------------
c here we exploit the fact that with MSR format a, ja, ja is actually
c the correct data structure for the strict lower triangular part of
c the CSR format. First rescale matrix. 
c
            scal = 0.1
            do 3 k=ja2(1), ja2(n+1)-1
               a2(k)=a2(k)*scal
 3          continue
            call amux(n, x, y0, a2, ja2, ja2) 
            do 4 j=1,n
               y0(j) = x(j) + y0(j) 
 4          continue
            call lsol (n, y, y0, a2, ja2, ja2) 
            call errpr (n, x, y, iout,'lsol  ')             
c-----------------------------------------------------------------------
c TESTING UDSOL
c-----------------------------------------------------------------------
            call getu (n, a1, ja1, ia1, a2, ja2, ia2)
            call amux (n,x,y0, a2, ja2, ia2)
            call atmux(n,x,y1, a2, ja2, ia2)
            call csrmsr (n, a2, ja2, ia2, a2, ja2, y, jad) 
            do 5 k=1,n
               a2(k) = 1.0d0/ a2(k)
 5          continue
            call udsol (n, y, y0, a2, ja2)
            call errpr (n, x, y, iout,'udsol ') 
c-----------------------------------------------------------------------
c TESTING LDSOLC
c-----------------------------------------------------------------------
c here we take advantage of the fact that the MSR format for L
c is the MSC format for U 
c 
            call ldsolc (n, y, y1, a2, ja2)
            call errpr (n, x, y, iout,'ldsolc') 
c-----------------------------------------------------------------------
c TESTING USOL 
c-----------------------------------------------------------------------
c here we exploit the fact that with MSR format a, ja, ja is actually
c the correct data structure for the strict lower triangular part of
c the CSR format. First rescale matrix. 
c
            scal = 0.1
            do 6 k=ja2(1), ja2(n+1)-1
               a2(k)=a2(k)*scal
 6          continue
            call amux(n, x, y1, a2, ja2, ja2) 
            do 7 j=1,n
               y1(j) = x(j) + y1(j) 
 7          continue
            call usol (n, y, y1, a2, ja2, ja2) 
            call errpr (n, x, y, iout,'usol  ')             
c-----------------------------------------------------------------------
c   --  END --
c----------------------------------------------------------------------- 
 150     continue
 100  continue
      stop 
c---------------end-of-main---------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
        subroutine errpr (n, y, y1,iout,msg) 
        real*8 y(*), y1(*), t, sqrt
        character*6 msg
        t = 0.0d0
        do 1 k=1,n
           t = t+(y(k)-y1(k))**2
 1      continue
        t = sqrt(t)
        write (iout,*) ' 2-norm of difference in ',msg,' =', t
        return
        end
       
