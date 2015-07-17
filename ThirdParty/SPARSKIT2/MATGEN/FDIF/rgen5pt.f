      program fivept
c-----------------------------------------------------------------------
c main program for generating 5 point and 7-point matrices in the
c Harwell-Boeing format.  Creates a file with containing a 
c harwell-boeing matrix. typical session:
c user answer are after the colon
c Enter nx, ny, nz  : 10 10 1
c Filename for matrix: test.mat
c output matrix in data file : test.mat
c
c nz = 1 will create a 2-D problem
c
c-----------------------------------------------------------------------
      integer nmx, nxmax
      parameter (nxmax = 50, nmx = nxmax*nxmax)
c      implicit none 
      integer ia(nmx),ja(7*nmx),iau(nmx)
      real*8 a(7*nmx),rhs(nmx),al(6)
      character title*72, key*8, type*3, matfile*50, guesol*2
c-----------------------------------------------------------------------    
      integer nx, ny, nz, iout, n, ifmt, job
      write (6,*)  '  '
      write(6,'(22hEnter  nx, ny, nz   : ,$)') 
      read (5,*) nx, ny, nz 
      write(6,'(22hFilename for matrix : ,$)') 
      read(5,'(a50)') matfile
      open (unit=7,file=matfile)
c
c     boundary condition is partly specified here
c
c      al(1) = 1.0D0
c      al(2) = 0.0D0
c      al(3) = 2.3D1
c      al(4) = 0.4D0
c      al(5) = 0.0D0
c      al(6) = 8.2D-2
       al(1) = 0.0D0
       al(2) = 0.0D0
       al(3) = 0.0D1
       al(4) = 0.0D0
       al(5) = 0.0D0
       al(6) = 0.0D0 
c
      call gen57pt (nx,ny,nz,al,0,n,a,ja,ia,iau,rhs)
      iout = 7
c
c     write out the matrix
c
      guesol='NN'
      title = 
     *     ' 5-POINT TEST MATRIX FROM SPARSKIT                    '
c          '123456789012345678901234567890123456789012345678901234567890
      type  = 'RUA' 
      key ='SC5POINT'
C           12345678 
      ifmt = 15
      job = 2
c   upper part only??  
c      call getu (n, a, ja, ia, a, ja, ia) 
      call prtmt (n,n,a,ja,ia,rhs,guesol,title,key,type,
     1     ifmt,job,iout) 
      write (6,*) ' output matrix in data file : ', matfile
c     
      stop
      end
      
