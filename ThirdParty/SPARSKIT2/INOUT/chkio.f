        program chkio                                               
c------------------------------------------------------------------c
c test suite for Part I : I/O routines.                            c
c tests the following : gen5pt.f, prtmt, readmt, amd pltmt.        c
c 1) generates a 100 x 100 5pt matrix,                             c 
c 2) prints it with a given format in file 'first.mat'             c
c 3) reads the matrix from 'first.mat' using readmat               c 
c 4) prints it again in file 'second.mat' in  a different format   c 
c 5) makes 4 pic files to show the different options of pltmt.     c
c    these are in job0.pic, job01.pic, job10.pic, job11.pic        c
c                          coded by Y. Saad, RIACS, 08/31/1989.    c
c------------------------------------------------------------------c
      parameter (nxmax = 20, nmx = nxmax*nxmax)
      implicit real*8 (a-h,o-z)
      integer ia(nmx),ja(7*nmx),iau(nmx) 
      real*8 a(7*nmx),rhs(3*nmx),al(6)
      character title*72, key*8, type*3, guesol*2 
c----- open statements ---------------- 
      open (unit=7,file='first.mat')
      open (unit=8,file='second.mat')
      open (unit=20,file='job00.pic')
      open (unit=21,file='job01.pic')
      open (unit=22,file='job10.pic')
      open (unit=23,file='job11.pic')
c
c---- dimension of grid 
c
      nx = 10
      ny = 10
      nz = 1
      al(1) = 1.0D0
      al(2) = 0.0D0
      al(3) = 2.3D1
      al(4) = 0.4D0
      al(5) = 0.0D0
      al(6) = 8.2D-2
c
c---- generate grid problem.
c
      call gen57pt (nx,ny,nz,al,0,n,a,ja,ia,iau,rhs)
c
c---- create the Harwell-Boeing matrix. Start by defining title, 
c     and type. them define format and print it.
c
      write (title,9) nx, ny
 9    format('Five-point matrix on a square region',
     *  ' using a ',I2,' by ',I2,' grid *SPARSKIT*')
      key = 'Fivept10'
      type= 'RSA'
      ifmt = 5
      job = 3
      guesol = 'GX'
c
c define a right hand side of ones, an initial guess of two's
c and an exact solution of three's.
c 
      do 2 k=1, 3*n
	 rhs(k) = real( 1 +  (k-1)/n ) 
 2    continue
c
      call prtmt (n,n,a,ja,ia,rhs,guesol,title,key,type,
     1	                  ifmt,job,7)
c---- read it again in same matrix a, ja, ia  
       nmax = nmx
       nzmax = 7*nmx
      do 3 k=1, 3*n
	 rhs(k) = 0.0
 3    continue
       job = 3
c
       rewind 7
c
       nrhs = 3*n   
c
       call readmt (nmax,nzmax,job,7,a,ja,ia,rhs,nrhs,guesol,
     1             nrow,ncol,nnz,title,key,type,ierr)
	 print *,  ' ierr = ', ierr, ' nrhs ' , nrhs
c
c matrix read.  print it again in a different format
c
      ifmt = 102
      ncol = nrow
      job = 3 
c
      call prtmt (nrow,ncol,a,ja,ia,rhs,guesol,title,key,type,
     1	                  ifmt,job,8)
c
c---- print four pic files 
c
      mode = 0
      do 10 i=1, 2
	 do 11 j=1, 2
	 job = (i-1)*10 +j-1
	 iout = 20+(i-1)*2+j-1
      call pltmt (nrow,ncol,mode,ja,ia,title,key,type,job,iout)
 11   continue
 10   continue
c-------- 
      stop
      end
