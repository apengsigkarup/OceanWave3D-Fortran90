      program zlatev
c-----------------------------------------------------------------------
c
c test suite for zlatev matrices. generates three matrices and 
c writes them in three different files in Harwell-Boeing format.
c      zlatev1.mat  produced from matrf2 
c      zlatev2.mat  produced from  dcn 
c      zlatev3.mat  produced from  ecn
c
c-----------------------------------------------------------------------
      parameter (nmax = 1000, nzmax=20*nmax)
      implicit real*8 (a-h,o-z)
      integer ia(nzmax), ja(nzmax), iwk(nmax)
      real*8 a(nzmax)
      character title*72, key*3,type*8, guesol*2
c
      open (unit=7,file='zlatev1.mat')
      open (unit=8,file='zlatev2.mat')
      open (unit=9,file='zlatev3.mat')
c
      m = 100
      n = m
      ic = n/2 
      index = 10
      alpha = 5.0
      nn = nzmax
c
c call matrf2
c 
      call matrf2(m,n,ic,index,alpha,nn,nz,a,ia,ja,ierr) 
      job = 1
c      do 110 i = 1, nz
c         print *, ia(i), ja(i), a(i)
c 110  continue
      call coicsr(n, nz, job, a, ja, ia, iwk) 
c-----
      title = ' 1st matrix from zlatev examples                '
      type  = 'RUA'
      key   = ' ZLATEV1'
      iout  = 7
      guesol='NN'
c
      ifmt = 3
      job = 2
c
c write result in H-B format. 
c 
c Replaces prtmt with smms in order to print matrix in format for 
c SMMS instead.
c      call smms (n,1,n,0,a,ja,ia,iout)
      call prtmt (n,n,a,ja,ia,rhs,guesol,title,type,key,
     1	                  ifmt,job,iout) 

c-------- second type of matrices dcn matrices ---------------
       n = 200 
       nn = nzmax 
       ic = 20
c-------------------------------------------------------
c matrix of the type e(c,n)
c-------------------------------------------------------
       call dcn(a,ia,ja,n,ne,ic,nn,ierr)
c---------------------------------------------------
      call coicsr(n, ne, job, a, ja, ia, iwk) 
      title = ' 2nd matrix from zlatev examples                '
      iout = iout+1
      guesol='NN'
      type  = 'RUA'
      key = ' ZLATEV2'
c
      ifmt = 3
      job = 2
c
c write result in second file 
c
c Replaced prtmt with smms in order to print matrix in format for 
c SMMS instead.
c      call smms (n,1,n,0,a,ja,ia,iout)
      call prtmt (n,n,a,ja,ia,rhs,guesol,title,type,key,
     1	                  ifmt,job,iout) 
c-------------------------------------------------------
c matrix of the type e(c,n)
c-------------------------------------------------------
       n = 200 
       ic = 20
       nn = nzmax 
c
c call ecn
c 
      call ecn(n,ic,ne,ia,ja,a,nn,ierr)
      call coicsr(n, ne, job, a, ja, ia, iwk) 
      title = ' 3nd matrix from zlatev examples                '
      guesol='NN'
      type  = 'RUA'
      key = ' ZLATEV3'
      iout = iout+1
c
      ifmt = 3
      job = 2
c
c write resulting matrix in third file
c
c Replaced prtmt with smms in order to print matrix in format for 
c SMMS instead.
c      call smms (n,1,n,0,a,ja,ia,iout)
      call prtmt (n,n,a,ja,ia,rhs,guesol,title,type,key,
     1	                  ifmt,job,iout) 
      stop
      end
