        program chkfmt 
c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c test suite for the Formats routines.                                 c
c tests all of the routines in the module formats.                     c
c----------------------------------------------------------------------c 
c Note: the comments may not have been updated. 
c
c Here is the sequence of what is done by this program.            
c 1) call gen57bl togenerate a block matrix associated with a simple
c     5-point matrix on a 4 x 2 grid (2-D) with 2 degrees of freedom
c     per grid point.  Thus N = 16. This is produced in BSR format.
c    the pattern of the reduced  matrix is written in csr.mat. 
c
c 2) the block format is translated into a compressed sparse row
c    format by bsrcsr. The result is dumped in file csr.mat
c    matrix is converted back to bsr format. block pattern shown again
c
c 3) the matrix is translated in dense format by csrdns.
c    result a 16 x 16 matrix is written in unit dns.mat.
c    This is a good file to look at to see what the matrix is
c    and to compare results of other formats with.
c 4) the dense matrix obtained in 3) is reconverted back to
c    csr format using dnscsr. Result appended to file csr.mat
c 5) The matrix obtained in 4) is converted in coordinate format
c    and the resulting matrix is written in file coo.mat
c 6) the result is converted back to csr format. matrix 
c    appended to csr.mat.
c 7) result of 6) is converted to symmetric sparse row storage
c    (ssr) and the result is appended to csr.mat
c 8) result of 7) converted back to csr format and result is
c    appended to csr.mat
c 9) matrix resulting from 8) is converted to modified sparse
c    row format using csrmsr and result is written in msr.mat.
c10) the resulting matrix is converted back to csrformat and
c    result is appended to csr.mat
c11) result is converted to ellpack-itpack format with 
c    csrell and result is printed in itp.mat
c12) result is converted back to csr format and appended to csr.mat 
c12) result converted to csc format (transposition) using csrcsc
c    which should produce the same matrix here. result appended 
c    to csr.mat. A second call to csrcsc is made on resulting
c    matrix.
c13) the subroutine csrdia is used to extract two diagonals
c    (offsets -1 and 0) and then all the diagonals of matrix.
c    results in dia.mat
c14) diacsr is then called to convert the diagonally stored matrix
c    back to csr format. result appended to csr.mat
c15) result is converted to band format (bnd) by calling
c    csrbnd. result dumped to bnd.mat
c16) result is converted back to csr format and appended to csr.mat
c17) result sorted by a call to csrcsc and then converted to
c    block format (csrbsr) and then back to csr format again.
c    result appedned to csr.mat. 
c18) matrix converted to symmetric skyline format. result appended
c    to file band.mat 
c19) matrix converted back to csr format and result appended to
c    csr.mat.
c20) result converted to jad format. result output in jad.mat
c21) result concverted back to csr fromat. appended to csr.mat
c------------------------------------------------------------------
      implicit none 
      integer ndns, nxmax, nmx, nnzmax
      parameter (nxmax = 10, nmx = nxmax*nxmax, nnzmax=10*nmx)
      integer ia(nmx+1),ja(nnzmax),ia1(nnzmax),ja1(nnzmax),
     *         iwk(nmx*2+1),ioff(20),iwk2(nmx+1)
      real*8 stencil(7,100),a(nnzmax),a1(nnzmax),dns(20,20),wk(nmx)
c
      integer k,nx,ny,nz,nfree,na,nr,n,iout,nnz,i,ierr,
     *     maxcol, ndiag,job,lowd,len,imod,j,kend,k1,k2,
     *     kstart,mu,ml,idiag,idiag0,kdiag 
c-----------------------------------------------------------------------
      data ndns/20/
c----- open statements ---------------- 
      open (unit=7,file='csr.mat')
      open (unit=8,file='dns.mat')
      open (unit=9,file='coo.mat')
      open (unit=10,file='msr.mat')
      open (unit=11,file='itp.mat')
      open (unit=12,file='dia.mat')
      open (unit=13,file='bnd.mat')
      open (unit=14,file='jad.mat')
c
c---- dimension of grid 
c
      nx = 4
      ny = 2
      nz = 1
      nfree = 2
c
c---- generate grid problem.
c
c      na = nx*ny*nz*5
      na = nfree*nfree 
      call gen57bl (nx,ny,nz,nfree,na,nr,a1,ja1,ia1,iwk,stencil)
c      nr = n / nfree 
      n = nr*nfree 
c
c---- dump the reduced matrix 
c 
      iout = 7
      nnz = ia(n+1)-1
      write (iout,*) '-----------------------------------------'
      write (iout,*) '  +++ Pattern of block matrix  (CSR) +++ '
      write (iout,*) '-----------------------------------------'
c
      call dump(1,nr,.false.,a1,ja1,ia1,7)
c
c---- convert to CSR format and dump result 
c
      call bsrcsr (1,nr,nfree,na,a1,ja1,ia1,a,ja,ia)
      iout = 7
      nnz = ia(n+1)-1
      write (iout,*) '-----------------------------------------'
      write (iout,*) '  +++  initial matrix in CSR format +++ '
      write (iout,*) '-----------------------------------------'
      call dump(1,n,.true.,a,ja,ia,7)
c
c---- convert back to BSR format and dump pattern again.
c
      call csrbsr (1,n,nfree,na,a,ja,ia,a1,ja1,ia1,iwk2,ierr) 
      write (iout,*) '-----------------------------------------'
      write (iout,*) '  +++ Pattern of block matrix  (CSR) +++ '
      write (iout,*) '-----------------------------------------'
      call dump(1,nr,.false.,a1,ja1,ia1,7)
c
c----- convert to BSR format. 
c
      call bsrcsr (1,nr,nfree,na,a1,ja1,ia1,a,ja,ia)
      iout = 7
      nnz = ia(n+1)-1
      write (iout,*) '-----------------------------------------'
      write (iout,*) ' +++  matrix after BSRCSR conversion +++ '
      write (iout,*) '-----------------------------------------'
      call dump(1,n,.true.,a,ja,ia,7)
c
c----- convert to dense format. 
c
      call csrdns(n,n,a,ja,ia,dns,ndns,ierr)
c
      iout = iout+1 
      write (iout,*) '-----------------------------------------'
      write (iout,*) '  +++  initial matrix in DENSE format+++ '
      write (iout,*) '-----------------------------------------'
      write (iout,'(4x,16i4)') (j,j=1,n)
      write (iout,'(3x,65(1h-))') 
      do 3  i=1,n
         write (8,102) i,(dns(i,j), j=1,n)
 102     format(1h ,i2,1h|,16f4.1)
 3    continue
c     
c----- convert back to sparse format.
c
      call dnscsr(n,n,nnzmax,dns,ndns,a1,ja1,ia1,ierr)
      write (7,*) '-----------------------------------------'
      write (7,*) '  +++ matrix after conversion from dnscsr +++ '
      write (7,*) '-----------------------------------------'
      if (ierr .ne. 0) write (7,*) ' ***** ERROR FROM DNSCSR'
      if (ierr .ne. 0) write (7,*)  '     IERR = ', ierr
      call dump(1,n,.true.,a1,ja1,ia1,7)
c     
c     convert it to coordinate format.
c     
      call csrcoo(n,3,nnzmax,a,ja,ia,nnz,a1,ia1,ja1,ierr)
      iout = iout+ 1
      if (ierr .ne. 0) write (iout,*) ' ***** ERROR IN CSRCOO'
      if (ierr .ne. 0) write (iout,*)   '     IERR = ', ierr
      write (iout,*) '-----------------------------------------'
      write(iout,*) ' +++ Matrix in coordinate format +++ '
      write (iout,*) '-----------------------------------------'
      write(iout,103) (ia1(j),ja1(j),a1(j),j=1,nnz) 
 103  format (' i =', i3,'    j = ',i3,'     a(i,j) = ',f4.1) 
c     
c     convert it back again to csr format 
c     
      call coocsr(n,nnz,a1,ia1,ja1,a,ja,ia)
      
      write (7,*) '-----------------------------------------'
      write (7,*) '  +++ matrix after conversion from coocsr +++ '
      write (7,*) '-----------------------------------------'
      call dump(1,n,.true.,a,ja,ia,7)
c     
c     going to srs format
c     
      call csrssr(n,a,ja,ia,nnzmax,a1,ja1,ia1,ierr)
      write (7,*) '-----------------------------------------'
      write (7,*) '  +++ matrix after conversion to ssr format +++ '
      write (7,*) '      (lower part only stored in csr format)    '
      write (7,*) '-----------------------------------------'
      call dump(1,n,.true.,a1,ja1,ia1,7)
c     back to csr
      call ssrcsr (3,1,n,a1,ja1,ia1,nnzmax,a,ja,ia,iwk,iwk2,ierr)
      if (ierr .ne. 0) write(7,*) ' error in ssrcsr-IERR=',ierr
      write (7,*) '-----------------------------------------'
      write (7,*) '  +++ matrix after conversion from ssrcsr +++ '
      write (7,*) '-----------------------------------------'
      call dump(1,n,.true.,a,ja,ia,7)
c---- msr format
      iout = iout+1
      call csrmsr (n,a,ja,ia,a1,ja1,a1,ja1)
      write (iout,*) '-----------------------------------------'
      write (iout,*) '  +++ matrix in modified sparse row format +++'
      write (iout,*) '-----------------------------------------'
      write (iout,*) ' ** MAIN DIAGONAL '
      write (iout,'(16f4.1)') (a1(k),k=1,n)
      write (iout,*) ' ** POINTERS: '
      write (iout,'(17i4)') (ja1(k),k=1,n+1)
      write (iout,*) ' ** REMAINDER :'
      call dump(1,n,.true.,a1,ja1,ja1,iout)
c-------
      call msrcsr (n,a1,ja1,a,ja,ia,wk,iwk2)
      write (7,*) '-----------------------------------------'
      write (7,*) ' +++ matrix after conversion from msrcsr +++'
      write (7,*) '-----------------------------------------'
c     
      call dump(1,n,.true.,a,ja,ia,7)
c     
      maxcol = 13 
c     
      call csrell (n,a,ja,ia,maxcol,a1,ja1,n,ndiag,ierr)
      iout = iout+1 
      if (ierr .ne. 0) write (iout,*) ' ***** ERROR IN CSRELL'
      if (ierr .ne. 0) write (iout,*)   '     IERR = ', ierr
      write (iout,*) '-----------------------------------------'
      write (iout,*) '  +++ matrix in ELLPACK-ITPACK format +++ '
      write (iout,*) '-----------------------------------------'
      do 12 i=1,ndiag
         write (iout,*) ' Column number: ', i
         write (iout,104) (a1(n*(i-1)+k),k=1,n)
 104     format(9h COEF  = ,16f4.0)
         write (iout,105) (ja1(n*(i-1)+k),k=1,n)
 105     format (9h JCOEF = ,16i4) 
 12   continue
      call ellcsr (n,a1,ja1,n,ndiag,a,ja,ia,nnzmax,ierr)
      if (ierr .ne. 0) write (7,*) ' ***** ERROR IN ELLCSR'
      if (ierr .ne. 0) write (7,*)   '     IERR = ', ierr
      write (7,*) '-----------------------------------------'
      write (7,*) '  +++ matrix after conversion from ellcsr +++'
      write (7,*) '-----------------------------------------'
      call dump(1,n,.true.,a,ja,ia,7)
c     
      call csrcsc(n,1,1,a,ja,ia, a1,ja1,ia1)
      write (7,*) '-----------------------------------------'
      write (7,*) '  +++ matrix after conversion from csrcsc  +++ '
      write (7,*) '-----------------------------------------'
      call dump(1,n,.true.,a1,ja1,ia1,7)
      call csrcsc(n,1,1,a1,ja1,ia1, a,ja,ia)	
c     
c--------test 1 : get main diagonal and subdiagonal
c     get some info on diagonals
      call infdia(n,ja,ia,iwk,idiag0)
      job = 0
      ioff(1) = 0
      ioff(2) = -1
      idiag = 2
      call csrdia (n,idiag,job,a,ja,ia,ndns,dns,ioff,a1,ja1,ia1,iwk)
      iout = iout+1
      write (iout,*) '-----------------------------------------'
      write (iout,*) '  +++  diagonal format +++ '
      write (iout,*) '-----------------------------------------'
      write (iout,*) '  diagonals ioff = 0 and ioff = -1 '
      write (iout,*) ' number of diag.s returned from csrdia=',idiag
      do 13 kdiag = 1, idiag
	 write (iout,*) ' diagonal offset = ', ioff(kdiag)
         write (iout,'(16f4.1)') (dns(k,kdiag),k=1,n)
 13   continue
c     reverse conversion
      ndiag = ndns
      idiag = idiag0 
      job = 10
      call csrdia (n,idiag,job,a,ja,ia,ndns,dns,ioff,a1,ja1,ia1,iwk)
      write (iout,*) '-----------------------------------------'
      write (iout,*) '  +++  second test diagonal format +++ '
      write (iout,*) '         ** all diagonals of A  ** '
      write (iout,*) '-----------------------------------------'
      write (iout,*) ' number of diagonals on return from csrdia=',
     *     idiag
      do 131 kdiag = 1, idiag
	 write (iout,*) ' diagonal offset = ', ioff(kdiag)
         write (iout,'(16f4.1)') (dns(k,kdiag),k=1,n)
 131  continue
c     
c     reverse conversion
c     
      job = 0
      call diacsr (n,job,idiag,dns,ndns,ioff,a,ja,ia)
c--------
      write (7,*) '-----------------------------------------'
      write (7,*) '  +++ matrix after conversion from diacsr  +++ '
      write (7,*) '-----------------------------------------'
      call dump(1,n,.true.,a,ja,ia,7)
c     
c     checking the banded format
c     
      lowd = 0
      job = 1
      call csrbnd(n,a,ja,ia,job,dns,ndns,lowd,ml,mu,ierr)
      iout = iout+1
      if (ierr .ne. 0) write (iout,*) ' ***** ERROR IN CSRBND'
      if (ierr .ne. 0) write (iout,*)   '     IERR = ', ierr
      write (iout,*) '-----------------------------------------'
      write (iout,*) '       +++  banded  format +++ '
      write (iout,*) ' bandwidth values found ml=',ml,'  mu=',mu
      write (iout,*) '-----------------------------------------' 
      write (iout,'(4x,16i4)') (j,j=1,n)
      write (iout,'(3x,65(1h-))') 
      do 14  i=1, lowd
         write (iout,102) i, (dns(i,j), j=1,n)
 14   continue
c     
c     convert back to a, ja, ia format.
c
      len = nnzmax 
c--------
      call bndcsr(n,dns,ndns,lowd,ml,mu,a,ja,ia,len,ierr)
      write (7,*) ' IERR IN BNDCSR = ', ierr
      write (7,*) '-----------------------------------------'
      write (7,*) '  +++ matrix after conversion from bndcsr +++'
      write (7,*) '-----------------------------------------'
      call dump(1,n,.true.,a,ja,ia,7)
c     
c     make sure it is sorted
c     
      call csrcsc(n,1,1,a,ja,ia, a1,ja1,ia1)	
c     
c     checking skyline format.
c     
      imod = 1
      call csrssk (n,imod,a1,ja1,ia1,a,ia,nnzmax,ierr)
c     
      if (ierr .ne. 0) write (iout,*)   '     IERR = ', ierr
      write (iout,*) '-----------------------------------------'
      write (iout,*) '    +++ Sym. Skyline format +++ '
      write (iout,*) '-----------------------------------------' 
      write (iout,'(3x,65(1h-))') 
c---------------------
c     create column values.
c---------------------
      do 15  i=1, n
         kend = ia(i+1)-1 
         kstart = ia(i)  
         do k=kstart,kend             
            ja(k) =  i-(kend-k)
         enddo
 15   continue
c     
      call dump(1,n,.true.,a,ja,ia,iout)
c     
c     back to ssr format..
c     
      call sskssr (n,imod,a,ia,a1,ja1,ia1,nnzmax,ierr)
      write (7,*) '-----------------------------------------'
      write (7,*) '  +++ matrix after conversion from sskcsr +++'
      write (7,*) '-----------------------------------------'
      call dump(1,n,.true.,a1,ja1,ia1,7)
c     
c     checking jad format -----
c     
c     first go back to the csr format ----
c     
      call ssrcsr (3,1,n,a1,ja1,ia1,nnzmax,a,ja,ia,iwk,iwk2,ierr)
c     
      call csrjad (n, a, ja, ia, ndiag, iwk, a1, ja1, ia1) 
c     
      iout = iout+1 
      write (iout,*) '-----------------------------------------'
      write (iout,*) '   +++   matrix in JAD format +++ '
      write (iout,*) '-----------------------------------------'
c     
c     permutation array
c     
      write (iout,*) ' ** PERMUTATION ARRAY '
      write (iout,'(17i4)') (iwk(k),k=1,n)
c     ------ diagonals 
      do 16 i=1,ndiag
         write (iout,*) ' J-diagonal number: ', i
         k1 = ia1(i)
         k2 = ia1(i+1)-1 
         write (iout,104) (a1(k),k=k1,k2) 
         write (iout,105) (ja1(k),k=k1,k2) 
 16   continue
c     
c     back to csr format..
c     
      call jadcsr (n, ndiag, a1, ja1, ia1, iwk, a, ja, ia) 
c     
      write (7,*) '-----------------------------------------'
      write (7,*) '  +++ matrix after conversion from jadcsr +++'
      write (7,*) '-----------------------------------------'
      call dump(1,n,.true.,a,ja,ia,7)
c-----------------------------------------------------------------------
c     checking the linked list format 
c-----------------------------------------------------------------------
c     
      nnz = ia(n+1) - ia(1) 
      call csrlnk (n, a, ja, ia, iwk) 
c     
c     print links in file 7 (no need for another file) 
c     
      iout = 7 
      write (iout,*) '-----------------------------------------'
      write (iout,*) '   +++   matrix in LNK format +++ '
      write (iout,*) '-----------------------------------------'
c     
c     permutation array
c     
      write (iout,*) ' LINK ARRAY '
      write (iout,*) ' ---------- ' 
      write (iout,'(17i4)') (iwk(k),k=1,nnz) 
c     
c     back to csr format..
c     
      call lnkcsr (n, a, ja, ia, iwk, a1, ja1, ia1) 
c     
      write (7,*) '-----------------------------------------'
      write (7,*) '  +++ matrix after conversion from lnkcsr +++'
      write (7,*) '-----------------------------------------'
      call dump(1,n,.true.,a,ja,ia,7)
c------------------------------------------------------------
      stop
c-----------------------------------------------------------------------
c-----------end-of-chkfmt1---------------------------------------------- 
      end
