       subroutine dinfo1(n,iout,a,ja,ia,valued,
     *		          title,key,type,ao,jao,iao)
	implicit real*8 (a-h,o-z)
        real*8 a(*),ao(*)
        integer ja(*),ia(n+1),jao(*),iao(n+1),nzdiag 
        character title*72,key*8,type*3 
        logical valued 
c----------------------------------------------------------------------c
c  SPARSKIT:  ELEMENTARY INFORMATION ROUTINE.                          c
c----------------------------------------------------------------------c
c info1 obtains a number of statistics on a sparse matrix and writes   c
c it into the output unit iout. The matrix is assumed                  c
c to be stored in the compressed sparse COLUMN format sparse a, ja, ia c
c----------------------------------------------------------------------c
c Modified Nov 1, 1989. 1) Assumes A is stored in column               
c format. 2) Takes symmetry into account, i.e., handles Harwell-Boeing
c            matrices correctly. 
c          ***  (Because of the recent modification the words row and 
c            column may be mixed-up at occasions... to be checked...
c
c bug-fix July 25: 'upper' 'lower' mixed up in formats 108-107.
c
c On entry :
c-----------
c n	= integer. column dimension of matrix	
c iout  = integer. unit number where the information it to be output.	
c a	= real array containing the nonzero elements of the matrix
c	  the elements are stored by columns in order 
c	  (i.e. column i comes before column i+1, but the elements
c         within each column can be disordered).
c ja	= integer array containing the row indices of elements in a
c ia	= integer array containing of length n+1 containing the 
c         pointers to the beginning of the columns in arrays a and ja.
c	  It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
c 
c valued= logical equal to .true. if values are provided and .false.
c         if only the pattern of the matrix is provided. (in that
c         case a(*) and ao(*) are dummy arrays.
c
c title = a 72-character title describing the matrix
c         NOTE: The first character in title is ignored (it is often
c         a one).
c
c key   = an 8-character key for the matrix
c type  = a 3-character string to describe the type of the matrix.
c         see harwell/Boeing documentation for more details on the 
c         above three parameters.
c 
c on return
c---------- 
c 1) elementary statistics on the matrix is written on output unit 
c    iout. See below for detailed explanation of typical output. 
c 2) the entries of a, ja, ia are sorted.
c
c---------- 
c 
c ao	= real*8 array of length nnz used as work array.
c jao	= integer work array of length max(2*n+1,nnz) 
c iao   = integer work array of length n+1
c
c Note  : title, key, type are the same paramaters as those
c         used for Harwell-Bowing matrices.
c 
c-----------------------------------------------------------------------
c Output description:
c--------------------
c *** The following info needs to be updated.
c
c + A header containing the Title, key, type of the matrix and, if values
c   are not provided a message to that effect.
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c    * SYMMETRIC STRUCTURE MEDIEVAL RUSSIAN TOWNS                       
c    *                    Key = RUSSIANT , Type = SSA                   
c    * No values provided - Information of pattern only                 
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c +  dimension n, number of nonzero elements nnz, average number of
c    nonzero elements per column, standard deviation for this average.
c +  if the matrix is upper or lower triangular a message to that effect
c    is printed. Also the number of nonzeros in the strict upper 
c    (lower) parts and the main diagonal are printed.
c +  weight of longest column. This is the largest number of nonzero 
c    elements in a column encountered. Similarly for weight of 
c    largest/smallest row.
c +  lower dandwidth as defined by
c          ml = max ( i-j, / all  a(i,j).ne. 0 ) 
c +  upper bandwidth as defined by
c          mu = max ( j-i, / all  a(i,j).ne. 0 ) 
c    NOTE that ml or mu can be negative. ml .lt. 0 would mean
c    that A is confined to the strict upper part above the diagonal 
c    number -ml. Similarly for mu.
c
c +  maximun bandwidth as defined by
c    Max (  Max [ j ; a(i,j) .ne. 0 ] - Min [ j ; a(i,j) .ne. 0 ] )
c     i  
c +  average bandwidth = average over all columns of the widths each column.
c
c +  If there are zero columns /or rows a message is printed
c    giving the number of such columns/rows.
c    
c +  matching elements in A and transp(A) :this counts the number of
c    positions (i,j) such that if a(i,j) .ne. 0 then a(j,i) .ne. 0.
c    if this number is equal to nnz then the matrix is symmetric.
c +  Relative symmetry match : this is the ratio of the previous integer
c    over nnz. If this ratio is equal to one then the matrix has a
c    symmetric structure. 
c
c +  average distance of a given element from the diagonal, standard dev.
c    the distance of a(i,j) is defined as iabs(j-i). 
c
c +  Frobenius norm of A
c    Frobenius norm of 0.5*(A + transp(A))
c    Frobenius norm of 0.5*(A - transp(A))
c    these numbers provide information on the degree of symmetry
c    of the matrix. If the norm of the nonsymmetric part is
c    zero then the matrix is symmetric.
c
c + 90% of matrix is in the band of width k, means that
c   by moving away and in a symmetric manner from the main
c   diagonal you would have to include exactly k diagonals 
c   (k is always odd), in order to include 90% of the nonzero 
c   elements of A.  The same thing is then for 80%.
c 
c + The total number of nonvoid diagonals, i.e., among the
c   2n-1 diagonals of the matrix which have at least one nonxero
c   element.
c
c +  Most important diagonals. The code selects a number of k
c    (k .le. 10) diagonals that are the most important ones, i.e.
c    that have the largest number of nonzero elements. Any diagonal
c    that has fewer than 1% of the nonzero elements of A is dropped.
c    the numbers printed are the offsets with respect to the 
c    main diagonal, going from left tp right. 
c    Thus 0 means the main diagonal -1 means the subdiagonal, and 
c    +10 means the 10th upper diagonal.
c +  The accumulated percentages in the next line represent the 
c    percentage of the nonzero elements represented by the diagonals
c    up the current one put together. 
c    Thus:
c    *  The 10 most important diagonals are (offsets)    :             *
c    *     0     1     2    24    21     4    23    22    20    19     *
c    *  The accumulated percentages they represent are   :             *
c    *  40.4  68.1  77.7  80.9  84.0  86.2  87.2  88.3  89.4  90.4     *
c    *-----------------------------------------------------------------*
c    shows the offsets of the most important  diagonals and
c    40.4 represent ratio of the number of nonzero elements in the
c    diagonal zero (main diagonal) over the total number of nonzero
c    elements. the second number indicates that the diagonal 0 and the 
c    diagonal 1 together hold 68.1% of the matrix, etc..
c
c +  Block structure:
c    if the matrix has a block structure then the block size is found
c    and printed. Otherwise the info1 will say that the matrix
c    does not have a block structure. Note that block struture has
c    a very specific meaning here. the matrix has a block structure
c    if it consists of square blocks that are dense. even if there
c    are zero elements in the blocks  they should be represented 
c    otherwise it would be possible to determine the block size.
c
c-----------------------------------------------------------------------
	real*8 dcount(20),amx
	integer ioff(20)
	character*61 tmpst
	logical sym
c-----------------------------------------------------------------------
	data ipar1 /1/
	write (iout,99) 
        write (iout,97) title(2:72), key, type
 97     format(2x,' * ',a71,' *'/,
     *         2x,' *',20x,'Key = ',a8,' , Type = ',a3,25x,' *')
	if (.not. valued) write (iout,98)
 98    format(2x,' * No values provided - Information on pattern only',
     *   23x,' *')
c---------------------------------------------------------------------
	nnz = ia(n+1)-ia(1)
	sym = ((type(2:2) .eq. 'S') .or. (type(2:2) .eq. 'Z')
     *    .or. (type(2:2) .eq. 's') .or. (type(2:2) .eq. 'z')) 
c
        write (iout, 99)
        write(iout, 100) n, nnz
	job = 0
	if (valued) job = 1
	ipos = 1
        call csrcsc(n, job, ipos, a, ja, ia, ao, jao, iao)
        call csrcsc(n, job, ipos, ao, jao, iao, a, ja, ia)
c-------------------------------------------------------------------
c computing max bandwith, max number of nonzero elements per column
c min nonzero elements per column/row, row/column diagonal dominance 
c occurences, average distance of an element from diagonal, number of 
c elemnts in lower and upper parts, ...
c------------------------------------------------------------------
c    jao will be modified later, so we call skyline here
          call skyline(n,sym,ja,ia,jao,iao,nsky)
          call nonz_lud(n,ja,ia,nlower, nupper, ndiag) 
          call avnz_col(n,ja,ia,iao, ndiag, av, st)
c------ write out info ----------------------------------------------
	  if (sym)  nupper = nlower 
          write(iout, 101) av, st
          if (nlower .eq. 0 ) write(iout, 105)
 1        if (nupper .eq. 0) write(iout, 106)
          write(iout, 107) nlower
          write(iout, 108) nupper
          write(iout, 109) ndiag
c
          call nonz(n,sym, ja, ia, iao, nzmaxc, nzminc,
     *                  nzmaxr, nzminr, nzcol, nzrow)
          write(iout, 1020) nzmaxc, nzminc
c
	  if (.not. sym) write(iout, 1021) nzmaxr, nzminr
c
	  if (nzcol .ne. 0) write(iout,116) nzcol
	  if (nzrow .ne. 0) write(iout,115) nzrow
c     
          call diag_domi(n,sym,valued,a, ja,ia,ao, jao, iao,
     *                           ddomc, ddomr)
c----------------------------------------------------------------------- 
c symmetry and near symmetry - Frobenius  norms
c-----------------------------------------------------------------------
          call frobnorm(n,sym,a,ja,ia,Fnorm)   
          call ansym(n,sym,a,ja,ia,ao,jao,iao,imatch,av,fas,fan)
          call distaij(n,nnz,sym,ja,ia,dist, std)
          amx = 0.0d0
          do 40 k=1, nnz
            amx = max(amx, abs(a(k)) )
 40       continue
          write (iout,103) imatch, av, dist, std
          write(iout,96) 
          if (valued) then
             write(iout,104) Fnorm, fas, fan, amx, ddomr, ddomc
             write (iout,96)
          endif
c-----------------------------------------------------------------------
c--------------------bandedness- main diagonals ----------------------- -
c-----------------------------------------------------------------------
	n2 = n+n-1
	do 8 i=1, n2
             jao(i) = 0
 8        continue
          do 9 i=1, n
             k1 = ia(i)
             k2 = ia(i+1) -1
             do 91 k=k1, k2
                j = ja(k)
                jao(n+i-j) = jao(n+i-j) +1
 91          continue 
 9        continue
c
          call bandwidth(n,ja, ia, ml, mu, iband, bndav)
c
c     write bandwidth information .
c     
          write(iout,117)  ml, mu, iband, bndav
c     
          write(iout,1175) nsky
c     
c         call percentage_matrix(n,nnz,ja,ia,jao,90,jb2)
c         call percentage_matrix(n,nnz,ja,ia,jao,80,jb1)
          nrow = n
          ncol = n
          call distdiag(nrow,ncol,ja,ia,jao)
          call bandpart(n,ja,ia,jao,90,jb2)
          call bandpart(n,ja,ia,jao,80,jb1)
          write (iout,112) 2*jb2+1, 2*jb1+1
c-----------------------------------------------------------------
          nzdiag = 0
          n2 = n+n-1
          do 42 i=1, n2
             if (jao(i) .ne. 0) nzdiag=nzdiag+1
 42       continue
          call n_imp_diag(n,nnz,jao,ipar1, ndiag,ioff,dcount)
          write (iout,118) nzdiag
          write (tmpst,'(10i6)') (ioff(j),j=1,ndiag)
          write (iout,110) ndiag,tmpst
          write (tmpst,'(10f6.1)')(dcount(j), j=1,ndiag)
          write (iout,111) tmpst
          write (iout, 96)
c     jump to next page -- optional //
c     write (iout,'(1h1)') 
c-----------------------------------------------------------------------
c     determine block size if matrix is a block matrix..
c-----------------------------------------------------------------------
          call blkfnd(n, ja, ia, nblk)
          if (nblk .le. 1) then 
             write(iout,113)  
          else 
             write(iout,114) nblk
          endif 
          write (iout,96)
c     
c---------- done. Next define all the formats -------------------------- 
c
 99    format (2x,38(2h *))
 96    format (6x,' *',65(1h-),'*')
c-----------------------------------------------------------------------
 100   format(
     * 6x,' *  Dimension N                                      = ',
     * i10,'  *'/
     * 6x,' *  Number of nonzero elements                       = ',
     * i10,'  *')
 101   format(
     * 6x,' *  Average number of nonzero elements/Column        = ',
     * f10.4,'  *'/
     * 6x,' *  Standard deviation for above average             = ',
     * f10.4,'  *')
c-----------------------------------------------------------------------
 1020       format(
     * 6x,' *  Weight of longest column                         = ',
     * i10,'  *'/
     * 6x,' *  Weight of shortest column                        = ',
     * i10,'  *')
 1021       format(
     * 6x,' *  Weight of longest row                            = ',
     * i10,'  *'/
     * 6x,' *  Weight of shortest row                           = ',
     * i10,'  *')
 117        format(
     * 6x,' *  Lower bandwidth  (max: i-j, a(i,j) .ne. 0)       = ',
     * i10,'  *'/
     * 6x,' *  Upper bandwidth  (max: j-i, a(i,j) .ne. 0)       = ',
     * i10,'  *'/
     * 6x,' *  Maximum Bandwidth                                = ',
     * i10,'  *'/
     * 6x,' *  Average Bandwidth                                = ',
     * e10.3,'  *')
 1175       format(
     * 6x,' *  Number of nonzeros in skyline storage            = ',
     * i10,'  *')
 103   format(
     * 6x,' *  Matching elements in symmetry                    = ',
     * i10,'  *'/
     * 6x,' *  Relative Symmetry Match (symmetry=1)             = ',
     * f10.4,'  *'/
     * 6x,' *  Average distance of a(i,j)  from diag.           = ',
     * e10.3,'  *'/
     * 6x,' *  Standard deviation for above average             = ',
     * e10.3,'  *') 
 104   format(
     * 6x,' *  Frobenius norm of A                              = ',
     * e10.3,'  *'/
     * 6x,' *  Frobenius norm of symmetric part                 = ',
     * e10.3,'  *'/
     * 6x,' *  Frobenius norm of nonsymmetric part              = ',
     * e10.3,'  *'/
     * 6x,' *  Maximum element in A                             = ',
     * e10.3,'  *'/
     * 6x,' *  Percentage of weakly diagonally dominant rows    = ',
     * e10.3,'  *'/
     * 6x,' *  Percentage of weakly diagonally dominant columns = ',
     * e10.3,'  *')
 105        format(
     * 6x,' *  The matrix is lower triangular ...       ',21x,' *')
 106        format(
     * 6x,' *  The matrix is upper triangular ...       ',21x,' *')
 107        format(
     * 6x,' *  Nonzero elements in strict lower part            = ',
     * i10,'  *')
 108       format(
     * 6x,' *  Nonzero elements in strict upper part            = ',
     * i10,'  *')
 109       format(
     * 6x,' *  Nonzero elements in main diagonal                = ',
     * i10,'  *')
 110   format(6x,' *  The ', i2, ' most important',
     *	   ' diagonals are (offsets)    : ',10x,'  *',/,  
     * 6x,' *',a61,3x,' *')
 111   format(6x,' *  The accumulated percentages they represent are ',
     * '  : ', 10x,'  *',/,
     * 6x,' *',a61,3x,' *')
c 111	format( 
c     * 6x,' *  They constitute the following % of A             = ',
c     * f8.1,' %  *') 
 112	format( 
     * 6x,' *  90% of matrix is in the band of width            = ',
     * i10,'  *',/,
     * 6x,' *  80% of matrix is in the band of width            = ',
     * i10,'  *')
 113 	format( 
     * 6x,' *  The matrix does not have a block structure ',19x,
     *    ' *')
 114 	format( 
     * 6x,' *  Block structure found with block size            = ',
     * i10,'  *')
 115 	format( 
     * 6x,' *  There are zero rows. Number of such rows         = ', 
     * i10,'  *')
 116 	format( 
     * 6x,' *  There are zero columns. Number of such columns   = ', 
     * i10,'  *')
 118 	format( 
     * 6x,' *  The total number of nonvoid diagonals is         = ', 
     * i10,'  *')
c-------------------------- end of dinfo --------------------------
        return
        end
