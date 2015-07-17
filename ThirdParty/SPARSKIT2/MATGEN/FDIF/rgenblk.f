      program bfivept
c-----------------------------------------------------------------------
c main program for generating BLOCK 5 point and 7-point matrices in the
c Harwell-Boeing format.  Creates a file with containing a 
c harwell-boeing matrix. 
c
c max block size = 5
c max number of grid points = 8000  = ( nx * ny * nz .le. 8000) 
c matrix dimension =  (nx*ny*nz* Block-size**2) .le. 8000 * 25= 200,000
c
c typical session:
c Enter nx, ny, nz : 10 10 1 
c enter block-size : 4
c enter filename for matrix: test.mat 
c output matrix in data file : test.mat 
c
c nz =1 will create a 2-D problem
c-----------------------------------------------------------------------
      parameter (nxmax = 20, nmx = nxmax*nxmax*nxmax, ntot=nmx*25) 
      integer ia(ntot),ja(ntot),iau(ntot), iao(ntot),jao(ntot)
      real*8 stencil(7,100), a(ntot), ao(ntot) 
      character title*72,key*8,type*3, matfile*50, guesol*2
c-----------------------------------------------------------------------
      write (6,*)  '  '
      write(6,'(22hEnter  nx, ny, nz   : ,$)') 
      read (5,*) nx, ny, nz 
      write(6,'(22hnfree (Block size)  : ,$)') 
      read (5,*) nfree 
     
      write(6,'(22hFilename for matrix : ,$)') 

      read(5,'(a50)') matfile
      open (unit=7,file=matfile) 
c
      write (6,*) ' output in data file : ', matfile

c------------------------------------------------------
      na = nfree*nfree
c
      call gen57bl (nx,ny,nz,nfree,na,n,a,ja,ia,iau,stencil)
c------------------------------------------------------

      print *, ' n=', n, ' nfree ', nfree, ' na =', na

      call bsrcsr(1,n,nfree, na, a, ja, ia, ao, jao, iao)
      n = n * nfree ! Apr. 21, 1995

      guesol='NN'

      title =
     *     ' BLOCK 5-POINT TEST MATRIX FROM SPARSKIT               '
      type  = 'RUA'
      key   = 'BLOCK5PT'
C              12345678 
      ifmt = 15
      job = 2
      iout = 7
      call prtmt (n,n,ao,jao,iao,rhs,guesol,title,key,type,
     1     ifmt,job,iout) 
      print *, ' output in data file : ', matfile
c
      stop
      end
