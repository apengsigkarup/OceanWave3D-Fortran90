        program chkfmt 
c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c test suite for the unary routines.                                   c
c tests some of the routines in the module unary. Still needs to tests c
c many other routines.                                                 c
c Last update: May 2, 1994.
c----------------------------------------------------------------------c 
      parameter (nxmax = 10, nmx = nxmax*nxmax, nnzmax=10*nmx)
      implicit real*8 (a-h,o-z)
      integer ia(nmx+1),ja(nnzmax),ia1(nnzmax),ja1(nnzmax),
     *         iwk(nmx*2+1), perm(16), qperm(16) ,
     *         iwork(nnzmax*2) 
      real*8 a(nnzmax),a1(nnzmax),dns(20,20),rhs(nnzmax),al(6)
      data ndns/20/
      data qperm /1, 3, 6, 8, 9, 11, 14, 16, 2, 4, 5, 7, 10, 12, 13, 15/
c
c define correct permutation
c 
      do 1 i=1, 16
         perm(qperm(i)) = i
 1    continue
c----- open statements ---------------- 
      open (unit=7,file='unary.mat') 
c
c---- dimension of grid 
c
c     generate a 16 by 16 matrix (after eliminating boundary)
      nx = 6
      ny = 6
      nz = 1
      al(1) = 0.D0
      al(2) = 0.D0
      al(3) = 0.D0
      al(4) = 0.D0
      al(5) = 0.D0
      al(6) = 0.D0
c
c---- generate grid problem.
c
      call gen57pt (nx,ny,nz,al,0,n,a,ja,ia,iwk,rhs)
c
c---- write out the matrix
c
	iout = 7
	nnz = ia(n+1)-1
	write (iout,*) '-----------------------------------------'
        write (iout,*) '  +++  initial matrix in CSR format +++ '
	write (iout,*) '-----------------------------------------'
        call dump(1,n,.true.,a,ja,ia,iout)
c
c call csrdns
c
	 call csrdns(n,n,a,ja,ia,dns,ndns,ierr)
c
c write it out as a dense matrix. 
c
	 write (iout,*) '-----------------------------------------'
         write (iout,*) '  +++ initial matrix in DENSE format+++ '
	 write (iout,*) '-----------------------------------------'
         call dmpdns(n, n, ndns, dns, iout)
c
c red black ordering 
c
         job = 1
c
        call dperm (n,a,ja,ia,a1,ja1,ia1,perm,perm,job)
c
	nnz = ia(n+1)-1
	write (iout,*) '-----------------------------------------'
        write (iout,*) '  +++ red-black matrix in CSR format +++ '
	write (iout,*) '-----------------------------------------'
        call dump(1,n,.true.,a1,ja1,ia1,iout)
c
c sort matrix
c 
        call csort (n,a1,ja1,ia1,iwork,.true.) 
	nnz = ia(n+1)-1
	write (iout,*) '-----------------------------------------'
        write (iout,*) '  +++     matrix after sorting    +++ '
	write (iout,*) '-----------------------------------------'
        call dump(1,n,.true.,a1,ja1,ia1,iout)
c
c
c convert into dense format
c
         call csrdns(n, n, a1,ja1,ia1,dns,ndns,ierr) 
	 write (iout,*) '-----------------------------------------'
         write (iout,*) '  +++ red-black matrix in DENSE format+++ '
	 write (iout,*) '-----------------------------------------'
         call dmpdns(n,n, ndns, dns, iout)
         stop
         end
c----------------------------------------------------------------------- 
         subroutine dmpdns(nrow, ncol, ndns, dns, iout) 
         integer nrow, ncol, ndns, iout
         real*8 dns(ndns,*) 
c-----------------------------------------------------------------------
c this subroutine prints out a dense matrix in a simple format.
c the zero elements of the matrix are omitted. The format for the
c nonzero elements is f4.1, i.e., very little precision is provided.
c-----------------------------------------------------------------------
c on entry
c --------
c nrow = row dimension of matrix
c ncol = column dimension of matrix
c ndns = first dimension of array dns.
c dns  = double dimensional array of size n x n containing the matrix 
c iout = logical unit where to write matrix
c
c on return
c ---------
c matrix will be printed out on unit output iout.
c------------------------------------------------------------------------ 
c         local variables
         integer j, j1, j2, last, i
         character*80 fmt 
c
c prints out a dense matrix -- without the zeros. 
c
	 write (iout,'(4x,16i4)') (j,j=1,ncol) 
         fmt(1:5) = '    |'
         j1 = 6
         do 1 j=1, ncol
            j2 = j1+4
            fmt(j1:j2) = '----'
            j1 = j2
 1       continue 
         last = j1 
         fmt(last:last) = '|'
            write (iout,*) fmt
c
c undo loop 1 --- 
c
         j1 = 6
         do 2 j=1,ncol
            j2 = j1+4
            fmt(j1:j2) = '   '
            j1 = j2
 2          continue            
c           
         do 4  i=1, nrow 
         j1 = 6 
         write (fmt,101) i
 101     format(1h ,i2,2h |)
         do 3 j=1, ncol 
            j2= j1+4
            if (dns(i,j) .ne. 0.0) then
               write (fmt(j1:j2),102) dns(i,j)
 102           format(f4.1) 
            endif
           j1 = j2
 3       continue
           fmt(last:last) = '|'
           write (iout,*) fmt
 4       continue
         fmt(1:5) = '    |'
         j1 = 6
         do 5 j=1, ncol 
            j2 = j1+4
            fmt(j1:j2) = '----'
            j1 = j2
 5       continue 
         fmt(last:last) = '|'
            write (iout,*) fmt
         return
         end
