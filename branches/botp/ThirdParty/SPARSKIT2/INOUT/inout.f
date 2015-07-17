c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
C                        INPUT-OUTPUT MODULE                           c
c----------------------------------------------------------------------c
c contents:                                                            c
c----------                                                            c
c  readmt : reads matrices in the Boeing/Harwell format.               c
c  prtmt  : prints matrices in the Boeing/Harwell format.              c
c  dump   : outputs matrix rows in a simple format (debugging purposes)c 
c  pspltm : generates a post-script plot of the non-zero pattern of A  c
c  pltmt  : produces a 'pic' file for plotting a sparse matrix         c
c  smms   : write the matrx in a format used in SMMS package           c
c  readsm : reads matrics in coordinate format (as in SMMS package)    c
c  readsk : reads matrices in CSR format (simplified H/B formate).     c
c  skit   : writes matrics to a file, format same as above.            c
c  prtunf : writes matrics (in CSR format) unformatted                 c
c  readunf: reads unformatted data of matrics (in CSR format)          c
c----------------------------------------------------------------------c
	subroutine readmt (nmax,nzmax,job,iounit,a,ja,ia,rhs,nrhs,
     *                     guesol,nrow,ncol,nnz,title,key,type,ierr)
c-----------------------------------------------------------------------
c this subroutine reads  a boeing/harwell matrix. handles right hand 
c sides in full format only (no sparse right hand sides).
c Also the matrix must be in assembled forms.
c Author: Youcef Saad - Date: Sept., 1989
c         updated Oct 31, 1989.
c-----------------------------------------------------------------------
c on entry:
c---------
c nmax 	 =  max column dimension  allowed for matrix. The array ia should 
c	    be of length at least ncol+1 (see below) if job.gt.0
c nzmax	 = max number of nonzeros elements allowed. the arrays a, 
c          and ja should be of length equal to nnz (see below) if these
c          arrays are to be read (see job).
c          
c job	 = integer to indicate what is to be read. (note: job is an
c          input and output parameter, it can be modified on return)
c          job = 0    read the values of ncol, nrow, nnz, title, key,
c                     type and return. matrix is not read and arrays
c                     a, ja, ia, rhs are not touched.
c          job = 1    read srtucture only, i.e., the arrays ja and ia.
c          job = 2    read matrix including values, i.e., a, ja, ia
c          job = 3    read matrix and right hand sides: a,ja,ia,rhs.
c		      rhs may contain initial guesses and exact 
c                     solutions appended to the actual right hand sides.
c		      this will be indicated by the output parameter
c                     guesol [see below]. 
c                     
c nrhs   = integer. nrhs is an input as well as ouput parameter.
c          at input nrhs contains the total length of the array rhs.
c          See also ierr and nrhs in output parameters.
c
c iounit = logical unit number where to read the matrix from.
c
c on return:
c---------- 
c job    = on return job may be modified to the highest job it could
c          do: if job=2 on entry but no matrix values are available it
c          is reset to job=1 on return. Similarly of job=3 but no rhs 
c          is provided then it is rest to job=2 or job=1 depending on 
c          whether or not matrix values are provided.
c          Note that no error message is triggered (i.e. ierr = 0 
c          on return in these cases. It is therefore important to
c          compare the values of job on entry and return ).
c
c a	 = the a matrix in the a, ia, ja (column) storage format
c ja 	 = row number of element a(i,j) in array a.
c ia     = pointer  array. ia(i) points to the beginning of column i.
c
c rhs    = real array of size nrow + 1 if available (see job)
c
c nrhs   = integer containing the number of right-hand sides found
c          each right hand side may be accompanied with an intial guess
c          and also the exact solution.
c
c guesol = a 2-character string indicating whether an initial guess 
c          (1-st character) and / or the exact solution (2-nd
c          character) is provided with the right hand side.
c	   if the first character of guesol is 'G' it means that an
c          an intial guess is provided for each right-hand side.
c          These are appended to the right hand-sides in the array rhs.
c	   if the second character of guesol is 'X' it means that an
c          exact solution is provided for each right-hand side.
c          These are  appended to the right hand-sides 
c          and the initial guesses (if any) in the array rhs.
c
c nrow   = number of rows in matrix
c ncol	 = number of columns in matrix 
c nnz	 = number of nonzero elements in A. This info is returned
c          even if there is not enough space in a, ja, ia, in order
c          to determine the minimum storage needed. 
c
c title  = character*72 = title of matrix test ( character a*72). 
c key    = character*8  = key of matrix 
c type   = charatcer*3  = type of matrix.
c          for meaning of title, key and type refer to documentation 
c          Harwell/Boeing matrices.
c
c ierr   = integer used for error messages 
c         * ierr  =  0 means that  the matrix has been read normally. 
c         * ierr  =  1 means that  the array matrix could not be read 
c         because ncol+1 .gt. nmax
c         * ierr  =  2 means that  the array matrix could not be read 
c         because nnz .gt. nzmax 
c         * ierr  =  3 means that  the array matrix could not be read 
c         because both (ncol+1 .gt. nmax) and  (nnz .gt. nzmax )
c         * ierr  =  4 means that  the right hand side (s) initial 
c         guesse (s) and exact solution (s)   could  not be
c         read because they are stored in sparse format (not handled
c         by this routine ...) 
c         * ierr  =  5 means that the right-hand-sides, initial guesses
c         and exact solutions could not be read because the length of 
c         rhs as specified by the input value of nrhs is not 
c         sufficient to store them. The rest of the matrix may have
c         been read normally.
c 
c Notes:
c-------
c 1) The file inout must be open (and possibly rewound if necessary)
c    prior to calling readmt.
c 2) Refer to the documentation on the Harwell-Boeing formats
c    for details on the format assumed by readmt.
c    We summarize the format here for convenience.
c  
c    a) all lines in inout are assumed to be 80 character long.
c    b) the file consists of a header followed by the block of the 
c       column start pointers followed by the block of the
c       row indices, followed by the block of the real values and
c       finally the numerical values of the right-hand-side if a 
c       right hand side is supplied. 
c    c) the file starts by a header which contains four lines if no
c       right hand side is supplied and five lines otherwise.
c       * first line contains the title (72 characters long) followed by
c         the 8-character identifier (name of the matrix, called key)
c        [ A72,A8 ]
c       * second line contains the number of lines for each
c         of the following data blocks (4 of them) and the total number 
c         of lines excluding the header.
c        [5i4]
c       * the third line contains a three character string identifying
c         the type of matrices as they are referenced in the Harwell
c         Boeing documentation [e.g., rua, rsa,..] and the number of
c         rows, columns, nonzero entries.
c         [A3,11X,4I14]
c       * The fourth line contains the variable fortran format
c         for the following data blocks.
c         [2A16,2A20] 
c       * The fifth line is only present if right-hand-sides are 
c         supplied. It consists of three one character-strings containing
c         the storage format for the right-hand-sides 
c         ('F'= full,'M'=sparse=same as matrix), an initial guess 
c         indicator ('G' for yes), an exact solution indicator 
c         ('X' for yes), followed by the number of right-hand-sides
c         and then the number of row indices. 
c         [A3,11X,2I14] 
c     d) The three following blocks follow the header as described 
c        above.
c     e) In case the right hand-side are in sparse formats then 
c        the fourth block uses the same storage format as for the matrix
c        to describe the NRHS right hand sides provided, with a column
c        being replaced by a right hand side.
c-----------------------------------------------------------------------
      character title*72, key*8, type*3, ptrfmt*16, indfmt*16,
     1       valfmt*20, rhsfmt*20, rhstyp*3, guesol*2
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     1     nnz, neltvl, nrhs, nmax, nzmax, nrwindx
      integer ia (nmax+1), ja (nzmax) 
      real*8 a(nzmax), rhs(*) 
c-----------------------------------------------------------------------
      ierr = 0
      lenrhs = nrhs
c
      read (iounit,10) title, key, totcrd, ptrcrd, indcrd, valcrd, 
     1     rhscrd, type, nrow, ncol, nnz, neltvl, ptrfmt, indfmt, 
     2     valfmt, rhsfmt
 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
c
      if (rhscrd .gt. 0) read (iounit,11) rhstyp, nrhs, nrwindx
 11   format (a3,11x,i14,i14)
c
c anything else to read ?
c
      if (job .le. 0) return
c     ---- check whether matrix is readable ------ 
      n = ncol
      if (ncol .gt. nmax) ierr = 1
      if (nnz .gt. nzmax) ierr = ierr + 2
      if (ierr .ne. 0) return
c     ---- read pointer and row numbers ---------- 
      read (iounit,ptrfmt) (ia (i), i = 1, n+1)
      read (iounit,indfmt) (ja (i), i = 1, nnz)
c     --- reading values of matrix if required....
      if (job .le. 1)  return
c     --- and if available ----------------------- 
      if (valcrd .le. 0) then
	 job = 1
	 return
      endif
      read (iounit,valfmt) (a(i), i = 1, nnz)
c     --- reading rhs if required ---------------- 
      if (job .le. 2)  return
c     --- and if available ----------------------- 
      if ( rhscrd .le. 0) then
	 job = 2
	 nrhs = 0
	 return
      endif
c     
c     --- read right-hand-side.-------------------- 
c     
      if (rhstyp(1:1) .eq. 'M') then 
         ierr = 4
         return
      endif
c
      guesol = rhstyp(2:3) 
c     
      nvec = 1 
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') nvec=nvec+1
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') nvec=nvec+1
c     
      len = nrhs*nrow 
c     
      if (len*nvec .gt. lenrhs) then
         ierr = 5
         return
      endif
c
c read right-hand-sides
c
      next = 1
      iend = len
      read(iounit,rhsfmt) (rhs(i), i = next, iend)
c
c read initial guesses if available
c
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') then
         next = next+len
         iend = iend+ len
         read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c     
c read exact solutions if available
c
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') then
         next = next+len
         iend = iend+ len
         read(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c     
      return
c--------- end of readmt -----------------------------------------------
c----------------------------------------------------------------------- 
      end
c-----------------------------------------------------------------------
      subroutine prtmt (nrow,ncol,a,ja,ia,rhs,guesol,title,key,type,
     1     ifmt,job,iounit)
c-----------------------------------------------------------------------
c writes a matrix in Harwell-Boeing format into a file.
c assumes that the matrix is stored in COMPRESSED SPARSE COLUMN FORMAT.
c some limited functionality for right hand sides. 
c Author: Youcef Saad - Date: Sept., 1989 - updated Oct. 31, 1989 to
c cope with new format. 
c-----------------------------------------------------------------------
c on entry:
c---------
c nrow   = number of rows in matrix
c ncol	 = number of columns in matrix 
c a	 = real*8 array containing the values of the matrix stored 
c          columnwise
c ja 	 = integer array of the same length as a containing the column
c          indices of the corresponding matrix elements of array a.
c ia     = integer array of containing the pointers to the beginning of 
c	   the row in arrays a and ja.
c rhs    = real array  containing the right-hand-side (s) and optionally
c          the associated initial guesses and/or exact solutions
c          in this order. See also guesol for details. the vector rhs will
c          be used only if job .gt. 2 (see below). Only full storage for
c          the right hand sides is supported. 
c
c guesol = a 2-character string indicating whether an initial guess 
c          (1-st character) and / or the exact solution (2-nd)
c          character) is provided with the right hand side.
c	   if the first character of guesol is 'G' it means that an
c          an intial guess is provided for each right-hand sides. 
c          These are assumed to be appended to the right hand-sides in 
c          the array rhs.
c	   if the second character of guesol is 'X' it means that an
c          exact solution is provided for each right-hand side.
c          These are assumed to be appended to the right hand-sides 
c          and the initial guesses (if any) in the array rhs.
c
c title  = character*72 = title of matrix test ( character a*72 ).
c key    = character*8  = key of matrix 
c type   = charatcer*3  = type of matrix.
c
c ifmt	 = integer specifying the format chosen for the real values
c	   to be output (i.e., for a, and for rhs-guess-sol if 
c          applicable). The meaning of ifmt is as follows.
c	  * if (ifmt .lt. 100) then the D descriptor is used,
c           format Dd.m, in which the length (m) of the mantissa is 
c           precisely the integer ifmt (and d = ifmt+6)
c	  * if (ifmt .gt. 100) then prtmt will use the 
c           F- descriptor (format Fd.m) in which the length of the 
c           mantissa (m) is the integer mod(ifmt,100) and the length 
c           of the integer part is k=ifmt/100 (and d = k+m+2)
c	    Thus  ifmt= 4   means  D10.4  +.xxxxD+ee    while
c	          ifmt=104  means  F7.4   +x.xxxx
c	          ifmt=205  means  F9.5   +xx.xxxxx
c	    Note: formats for ja, and ia are internally computed.
c
c job	 = integer to indicate whether matrix values and
c	   a right-hand-side is available to be written
c          job = 1   write srtucture only, i.e., the arrays ja and ia.
c          job = 2   write matrix including values, i.e., a, ja, ia
c          job = 3   write matrix and one right hand side: a,ja,ia,rhs.
c	   job = nrhs+2 write matrix and nrhs successive right hand sides
c	   Note that there cannot be any right-hand-side if the matrix
c	   has no values. Also the initial guess and exact solutions when 
c          provided are for each right hand side. For example if nrhs=2 
c          and guesol='GX' there are 6 vectors to write.
c          
c
c iounit = logical unit number where to write the matrix into.
c
c on return:
c---------- 
c the matrix a, ja, ia will be written in output unit iounit
c in the Harwell-Boeing format. None of the inputs is modofied.
c  
c Notes: 1) This code attempts to pack as many elements as possible per
c        80-character line. 
c        2) this code attempts to avoid as much as possible to put
c        blanks in the formats that are written in the 4-line header
c	 (This is done for purely esthetical reasons since blanks
c        are ignored in format descriptors.)
c        3) sparse formats for right hand sides and guesses are not
c        supported.
c-----------------------------------------------------------------------
      character title*72,key*8,type*3,ptrfmt*16,indfmt*16,valfmt*20,
     1     guesol*2, rhstyp*3
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     1     nnz, nrhs, len, nperli, nrwindx
      integer ja(*), ia(*) 	
      real*8 a(*),rhs(*)
c--------------
c     compute pointer format
c--------------
       nnz    = ia(ncol+1) -1
       if (nnz .eq. 0) then
	   return
	endif
      len    = int ( alog10(0.1+real(nnz+1))) + 1
      nperli = 80/len
      ptrcrd = ncol/nperli + 1
      if (len .gt. 9) then
         assign 101 to ix
      else
         assign 100 to ix
      endif
      write (ptrfmt,ix) nperli,len
 100  format(1h(,i2,1HI,i1,1h) )
 101  format(1h(,i2,1HI,i2,1h) )
c----------------------------
c compute ROW index format
c----------------------------
      len    = int ( alog10(0.1+real(nrow) )) + 1
      nperli = min0(80/len,nnz)
      indcrd = (nnz-1)/nperli+1
      write (indfmt,100) nperli,len
c---------------
c compute values and rhs format (using the same for both)
c--------------- 
      valcrd	= 0
      rhscrd  = 0
c quit this part if no values provided.
      if (job .le. 1) goto 20
c     
      if (ifmt .ge. 100) then
         ihead = ifmt/100
         ifmt = ifmt-100*ihead
         len = ihead+ifmt+2
         nperli = 80/len
c     
         if (len .le. 9 ) then
            assign 102 to ix
         elseif (ifmt .le. 9) then
            assign 103 to ix
         else 
            assign 104 to ix
         endif
c     
         write(valfmt,ix) nperli,len,ifmt
 102     format(1h(,i2,1hF,i1,1h.,i1,1h) )
 103     format(1h(,i2,1hF,i2,1h.,i1,1h) )
 104     format(1h(,i2,1hF,i2,1h.,i2,1h) )
C
      else
         len = ifmt + 6
         nperli = 80/len
c     try to minimize the blanks in the format strings.
         if (nperli .le. 9) then
	    if (len .le. 9 ) then
	       assign 105 to ix
	    elseif (ifmt .le. 9) then
	       assign 106 to ix
	    else 
	       assign 107 to ix
	    endif
	 else 
	    if (len .le. 9 ) then
	       assign 108 to ix
	    elseif (ifmt .le. 9) then
	       assign 109 to ix
	    else 
               assign 110 to ix
            endif
         endif
c-----------
         write(valfmt,ix) nperli,len,ifmt
 105     format(1h(,i1,1hD,i1,1h.,i1,1h) )
 106     format(1h(,i1,1hD,i2,1h.,i1,1h) )
 107     format(1h(,i1,1hD,i2,1h.,i2,1h) )
 108     format(1h(,i2,1hD,i1,1h.,i1,1h) )
 109     format(1h(,i2,1hD,i2,1h.,i1,1h) )
 110     format(1h(,i2,1hD,i2,1h.,i2,1h) )
c     
      endif 	    
      valcrd = (nnz-1)/nperli+1
      nrhs   = job -2
      if (nrhs .ge. 1) then
         i = (nrhs*nrow-1)/nperli+1
         rhscrd = i
         if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g')
     +      rhscrd = rhscrd+i
         if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x')
     +      rhscrd = rhscrd+i
         rhstyp = 'F'//guesol
      endif 
 20   continue
c     
      totcrd = ptrcrd+indcrd+valcrd+rhscrd
c     write 4-line or five line header
      write(iounit,10) title,key,totcrd,ptrcrd,indcrd,valcrd,
     1     rhscrd,type,nrow,ncol,nnz,nrhs,ptrfmt,indfmt,valfmt,valfmt
c-----------------------------------------------------------------------
      nrwindx = 0
      if (nrhs .ge. 1) write (iounit,11) rhstyp, nrhs, nrwindx
 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
 11   format(A3,11x,i14,i14)
c     
      write(iounit,ptrfmt) (ia (i), i = 1, ncol+1)
      write(iounit,indfmt) (ja (i), i = 1, nnz)
      if (job .le. 1) return
      write(iounit,valfmt) (a(i), i = 1, nnz)
      if (job .le. 2) return 
      len = nrow*nrhs 
      next = 1
      iend = len
      write(iounit,valfmt) (rhs(i), i = next, iend)
c     
c     write initial guesses if available
c     
      if (guesol(1:1) .eq. 'G' .or. guesol(1:1) .eq. 'g') then
         next = next+len
         iend = iend+ len
         write(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c     
c     write exact solutions if available
c     
      if (guesol(2:2) .eq. 'X' .or. guesol(2:2) .eq. 'x') then
         next = next+len
         iend = iend+ len
         write(iounit,valfmt) (rhs(i), i = next, iend)
      endif
c     
      return
c----------end of prtmt ------------------------------------------------ 
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine dump (i1,i2,values,a,ja,ia,iout)
      integer i1, i2, ia(*), ja(*), iout
      real*8 a(*) 
      logical values 
c-----------------------------------------------------------------------
c outputs rows i1 through i2 of a sparse matrix stored in CSR format 
c (or columns i1 through i2 of a matrix stored in CSC format) in a file, 
c one (column) row at a time in a nice readable format. 
c This is a simple routine which is useful for debugging. 
c-----------------------------------------------------------------------
c on entry:
c---------
c i1    = first row (column) to print out
c i2    = last row (column) to print out 
c values= logical. indicates whether or not to print real values.
c         if value = .false. only the pattern will be output.
c a,
c ja, 
c ia    =  matrix in CSR format (or CSC format) 
c iout  = logical unit number for output.
c---------- 
c the output file iout will have written in it the rows or columns 
c of the matrix in one of two possible formats (depending on the max 
c number of elements per row. The values are output with only 
c two digits of accuracy (D9.2). )
c-----------------------------------------------------------------------
c     local variables
      integer maxr, i, k1, k2 
c
c select mode horizontal or vertical 
c
        maxr = 0
        do 1 i=i1, i2
           maxr = max0(maxr,ia(i+1)-ia(i))
 1      continue
        
        if (maxr .le. 8) then
c
c able to do one row acros line
c
        do 2 i=i1, i2
           write(iout,100) i
	   k1=ia(i)
	   k2 = ia(i+1)-1
	   write (iout,101) (ja(k),k=k1,k2)
	   if (values) write (iout,102) (a(k),k=k1,k2)
 2      continue
      else 
c
c unable to one row acros line. do three items at a time
c across a line 
         do 3 i=i1, i2
            if (values) then
               write(iout,200) i
            else
               write(iout,203) i               
            endif
            k1=ia(i)
            k2 = ia(i+1)-1
            if (values) then
               write (iout,201) (ja(k),a(k),k=k1,k2)
            else
                write (iout,202) (ja(k),k=k1,k2)
             endif
 3       continue
      endif 
c
c formats :
c
 100  format (1h ,34(1h-),' row',i6,1x,34(1h-) )
 101  format(' col:',8(i5,6h     : ))
 102  format(' val:',8(D9.2,2h :) )
 200  format (1h ,30(1h-),' row',i3,1x,30(1h-),/
     *     3('  columns :    values  * ') )
c-------------xiiiiiihhhhhhddddddddd-*-
 201  format(3(1h ,i6,6h   :  ,D9.2,3h * ) )
 202  format(6(1h ,i5,6h  *    ) ) 
 203  format (1h ,30(1h-),' row',i3,1x,30(1h-),/
     *     3('  column  :  column   *') )
      return
c----end-of-dump--------------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine pspltm(nrow,ncol,mode,ja,ia,title,ptitle,size,munt,
     *     nlines,lines,iunt)
c-----------------------------------------------------------------------
      integer nrow,ncol,ptitle,mode,iunt, ja(*), ia(*), lines(nlines) 
      real size
      character title*(*), munt*2 
c----------------------------------------------------------------------- 
c PSPLTM - PostScript PLoTer of a (sparse) Matrix
c This version by loris renggli (renggli@masg1.epfl.ch), Dec 1991
c and Youcef Saad 
c------
c Loris RENGGLI, Swiss Federal Institute of Technology, Math. Dept
c CH-1015 Lausanne (Switzerland)  -- e-mail:  renggli@masg1.epfl.ch
c Modified by Youcef Saad -- June 24, 1992 to add a few features:
c separation lines + acceptance of MSR format.
c-----------------------------------------------------------------------
c input arguments description :
c
c nrow   = number of rows in matrix
c
c ncol   = number of columns in matrix 
c
c mode   = integer indicating whether the matrix is stored in 
c           CSR mode (mode=0) or CSC mode (mode=1) or MSR mode (mode=2) 
c
c ja     = column indices of nonzero elements when matrix is
c          stored rowise. Row indices if stores column-wise.
c ia     = integer array of containing the pointers to the 
c          beginning of the columns in arrays a, ja.
c
c title  = character*(*). a title of arbitrary length to be printed 
c          as a caption to the figure. Can be a blank character if no
c          caption is desired.
c
c ptitle = position of title; 0 under the drawing, else above
c
c size   = size of the drawing  
c
c munt   = units used for size : 'cm' or 'in'
c
c nlines = number of separation lines to draw for showing a partionning
c          of the matrix. enter zero if no partition lines are wanted.
c
c lines  = integer array of length nlines containing the coordinates of 
c          the desired partition lines . The partitioning is symmetric: 
c          a horizontal line across the matrix will be drawn in 
c          between rows lines(i) and lines(i)+1 for i=1, 2, ..., nlines
c          an a vertical line will be similarly drawn between columns
c          lines(i) and lines(i)+1 for i=1,2,...,nlines 
c
c iunt   = logical unit number where to write the matrix into.
c----------------------------------------------------------------------- 
c additional note: use of 'cm' assumes european format for paper size
c (21cm wide) and use of 'in' assumes american format (8.5in wide).
c The correct centering of the figure depends on the proper choice. Y.S.
c-----------------------------------------------------------------------
c external 
      integer LENSTR
      external LENSTR
c local variables ---------------------------------------------------
      integer n,nr,nc,maxdim,istart,ilast,ii,k,ltit
      real lrmrgn,botmrgn,xtit,ytit,ytitof,fnstit,siz
      real xl,xr, yb,yt, scfct,u2dot,frlw,delt,paperx,conv,xx,yy
      logical square 
c change square to .true. if you prefer a square frame around
c a rectangular matrix
      data haf /0.5/, zero/0.0/, conv/2.54/,square/.false./
c-----------------------------------------------------------------------
      siz = size
      nr = nrow
      nc = ncol
      n = nc
      if (mode .eq. 0) n = nr
c      nnz = ia(n+1) - ia(1) 
      maxdim = max(nrow, ncol)
      m = 1 + maxdim
      nc = nc+1
      nr = nr+1
c
c units (cm or in) to dot conversion factor and paper size
c 
      if (munt.eq.'cm' .or. munt.eq.'CM') then
         u2dot = 72.0/conv
        paperx = 21.0
      else
        u2dot = 72.0
        paperx = 8.5*conv
        siz = siz*conv
      end if
c
c left and right margins (drawing is centered)
c 
      lrmrgn = (paperx-siz)/2.0
c
c bottom margin : 2 cm
c
      botmrgn = 2.0
c scaling factor
      scfct = siz*u2dot/m
c matrix frame line witdh
      frlw = 0.25
c font size for title (cm)
      fnstit = 0.5
      ltit = LENSTR(title)
c position of title : centered horizontally
c                     at 1.0 cm vertically over the drawing
      ytitof = 1.0
      xtit = paperx/2.0
      ytit = botmrgn+siz*nr/m + ytitof
c almost exact bounding box
      xl = lrmrgn*u2dot - scfct*frlw/2
      xr = (lrmrgn+siz)*u2dot + scfct*frlw/2
      yb = botmrgn*u2dot - scfct*frlw/2
      yt = (botmrgn+siz*nr/m)*u2dot + scfct*frlw/2
      if (ltit.gt.0) then
        yt = yt + (ytitof+fnstit*0.70)*u2dot
      end if
c add some room to bounding box
      delt = 10.0
      xl = xl-delt
      xr = xr+delt
      yb = yb-delt
      yt = yt+delt
c
c correction for title under the drawing
      if (ptitle.eq.0 .and. ltit.gt.0) then
      ytit = botmrgn + fnstit*0.3
      botmrgn = botmrgn + ytitof + fnstit*0.7
      end if
c begin of output
c
      write(iunt,10) '%!'
      write(iunt,10) '%%Creator: PSPLTM routine'
      write(iunt,12) '%%BoundingBox:',xl,yb,xr,yt
      write(iunt,10) '%%EndComments'
      write(iunt,10) '/cm {72 mul 2.54 div} def'
      write(iunt,10) '/mc {72 div 2.54 mul} def'
      write(iunt,10) '/pnum { 72 div 2.54 mul 20 string'
      write(iunt,10) 'cvs print ( ) print} def'
      write(iunt,10)
     1  '/Cshow {dup stringwidth pop -2 div 0 rmoveto show} def'
c
c we leave margins etc. in cm so it is easy to modify them if
c needed by editing the output file
      write(iunt,10) 'gsave'
      if (ltit.gt.0) then
      write(iunt,*) '/Helvetica findfont ',fnstit,
     &             ' cm scalefont setfont '
      write(iunt,*) xtit,' cm ',ytit,' cm moveto '
      write(iunt,'(3A)') '(',title(1:ltit),') Cshow'
      end if
      write(iunt,*) lrmrgn,' cm ',botmrgn,' cm translate'
      write(iunt,*) siz,' cm ',m,' div dup scale '
c------- 
c draw a frame around the matrix
      write(iunt,*) frlw,' setlinewidth'
      write(iunt,10) 'newpath'
      write(iunt,11) 0, 0, ' moveto'
      if (square) then
      write(iunt,11) m,0,' lineto'
      write(iunt,11) m, m, ' lineto'
      write(iunt,11) 0,m,' lineto'
      else
      write(iunt,11) nc,0,' lineto'
      write(iunt,11) nc,nr,' lineto'
      write(iunt,11) 0,nr,' lineto'
      end if
      write(iunt,10) 'closepath stroke'
c
c     drawing the separation lines 
c 
      write(iunt,*)  ' 0.2 setlinewidth'
      do 22 kol=1, nlines 
         isep = lines(kol) 
c
c     horizontal lines 
c
         yy =  real(nrow-isep) + haf 
         xx = real(ncol+1) 
         write(iunt,13) zero, yy, ' moveto '
         write(iunt,13)  xx, yy, ' lineto stroke '
c
c vertical lines 
c
         xx = real(isep) + haf 
         yy = real(nrow+1)  
         write(iunt,13) xx, zero,' moveto '
         write(iunt,13) xx, yy, ' lineto stroke '             
 22     continue
c 
c----------- plotting loop ---------------------------------------------
c
      write(iunt,10) '1 1 translate'
      write(iunt,10) '0.8 setlinewidth'
      write(iunt,10) '/p {moveto 0 -.40 rmoveto '
      write(iunt,10) '           0  .80 rlineto stroke} def'
c     
      do 1 ii=1, n
        istart = ia(ii)
        ilast  = ia(ii+1)-1 
        if (mode .eq. 1) then
          do 2 k=istart, ilast
            write(iunt,11) ii-1, nrow-ja(k), ' p'
 2        continue 
        else
          do 3 k=istart, ilast
            write(iunt,11) ja(k)-1, nrow-ii, ' p'
 3        continue          
c add diagonal element if MSR mode.
          if (mode .eq. 2) 
     *         write(iunt,11) ii-1, nrow-ii, ' p' 
c
        endif
 1    continue
c-----------------------------------------------------------------------
      write(iunt,10) 'showpage'
      return
c
 10   format (A)
 11   format (2(I6,1x),A)
 12   format (A,4(1x,F9.2))
 13   format (2(F9.2,1x),A)
c-----------------------------------------------------------------------
      end
c
      integer function lenstr(s)
c-----------------------------------------------------------------------
c return length of the string S
c-----------------------------------------------------------------------
      character*(*) s
      integer len
      intrinsic len
      integer n
c----------------------------------------------------------------------- 
      n = len(s)
10    continue
        if (s(n:n).eq.' ') then
          n = n-1
          if (n.gt.0) go to 10
        end if
      lenstr = n
c
      return
c--------end-of-pspltm--------------------------------------------------
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine pltmt (nrow,ncol,mode,ja,ia,title,key,type,
     1     job, iounit)
c-----------------------------------------------------------------------
c this subroutine creates a 'pic' file for plotting the pattern of 
c a sparse matrix stored in general sparse format. it is not intended
c to be a means of plotting large matrices (it is very inefficient).
c It is however useful for small matrices and can be used for example
c for inserting matrix plots in a text. The size of the plot can be
c 7in x 7in or 5 in x 5in .. There is also an option for writing a 
c 3-line header in troff (see description of parameter job).
c Author: Youcef Saad - Date: Sept., 1989
c See SPARSKIT/UNSUPP/ for a version of this to produce a post-script
c file. 
c-----------------------------------------------------------------------
c nrow   = number of rows in matrix
c
c ncol	 = number of columns in matrix 
c
c mode   = integer indicating whether the matrix is stored
c          row-wise (mode = 0) or column-wise (mode=1)
c
c ja     = column indices of nonzero elements when matrix is
c	   stored rowise. Row indices if stores column-wise.
c ia     = integer array of containing the pointers to the 
c	   beginning of the columns in arrays a, ja.
c
c title  = character*71 = title of matrix test ( character a*71 ).
c key    = character*8  = key of matrix 
c type   = character*3  = type of matrix. 
c
c job    = this integer parameter allows to set a few minor 
c          options. First it tells pltmt whether or not to
c          reduce the plot. The standard size of 7in is then
c          replaced by a 5in plot. It also tells pltmt whether or
c          not to append to the pic file a few 'troff' lines that 
c          produce a centered caption includingg the title, key and 
c          types as well as the size and number of nonzero elements.
c          job =  0 : do not reduce and do not make caption.
c          job =  1 : reduce and do not make caption.
c          job = 10 : do not reduce and make caption
c          job = 11 : reduce and make caption.
c          (i.e. trailing digit for reduction, leading digit for caption)
c
c iounit = logical unit number where to write the matrix into.
c
c-----------------------------------------------------------------------
c example of usage . 
c-----------------
c In the fortran code:
c  a) read a Harwell/Boeing matrix
c          call readmt (.....)
c	   iout = 13
c  b) generate pic file:
c          call  pltmt (nrow,ncol,mode,ja,ia,title,key,type,iout)
c	   stop
c ---------
c Then in a unix environment plot the matrix by the command
c
c	pic FOR013.DAT | troff -me | lpr -Ppsx
c
c-----------------------------------------------------------------------
c notes: 1) Plots square as well as rectangular matrices.
c            (however not as much tested with rectangular matrices.)
c	  2) the dot-size is adapted according to the size of the
c            matrix.
c	  3) This is not meant at all as a way of plotting large
c            matrices. The pic file generaled will have one line for
c            each nonzero element. It is  only meant for use in
c	     such things as document poreparations etc..
c         4) The caption written will print the 71 character long
c            title. This may not be centered correctly if the
c            title has trailing blanks (a problem with Troff).
c            if you want the title centered then you can center
c            the string in title before calling pltmt. 
c       
c-----------------------------------------------------------------------
      integer ja(*), ia(*)
      character key*8,title*72,type*3
      real x, y
c-------
      n = ncol
      if (mode .eq. 0) n = nrow
      nnz = ia(n+1) - ia(1) 
      maxdim = max0 (nrow, ncol)
      xnrow = real(nrow)
      ptsize = 0.08
      hscale = (7.0 -2.0*ptsize)/real(maxdim-1) 
      vscale = hscale 
      xwid  = ptsize + real(ncol-1)*hscale + ptsize
      xht   = ptsize + real(nrow-1)*vscale + ptsize
      xshift = (7.0-xwid)/2.0
      yshift = (7.0-xht)/2.0 
c------
      if (mod(job,10) .eq. 1) then
         write (iounit,88)
      else
         write (iounit,89)
      endif
 88   format('.PS 5in',/,'.po 1.8i')
 89   format('.PS',/,'.po 0.7i')
      write(iounit,90) 
 90   format('box invisible wid 7.0 ht 7.0 with .sw at (0.0,0.0) ') 
      write(iounit,91) xwid, xht, xshift, yshift
 91   format('box wid ',f5.2,' ht ',f5.2,
     *     ' with .sw at (',f5.2,',',f5.2,')' )
c     
c     shift points slightly to account for size of dot , etc..
c     
      tiny = 0.03
      if (mod(job,10) .eq. 1) tiny = 0.05
      xshift = xshift + ptsize - tiny
      yshift = yshift + ptsize + tiny
c     
c-----------------------------------------------------------------------
c     
      ips = 8
      if (maxdim .le. 500) ips = 10
      if (maxdim .le. 300) ips = 12
      if (maxdim .le. 100) ips = 16
      if (maxdim .lt. 50) ips = 24
      write(iounit,92) ips
 92   format ('.ps ',i2)
c     
c-----------plottingloop --------------------------------------------- 
c     
      do 1 ii=1, n
         istart = ia(ii)
         ilast  = ia(ii+1)-1 
         if (mode .ne. 0) then
            x = real(ii-1)
            do 2 k=istart, ilast
               y = xnrow-real(ja(k))
               write(iounit,128) xshift+x*hscale, yshift+y*vscale
 2          continue 
         else
            y = xnrow - real(ii)
            do 3 k=istart, ilast
               x = real(ja(k)-1)
               write(iounit,128) xshift+x*hscale, yshift+y*vscale
 3          continue	    
         endif
 1    continue
c-----------------------------------------------------------------------
 128  format(7h"." at ,f6.3,1h,,f6.3,8h ljust  )
      write (iounit, 129)
 129  format('.PE')
c     quit if caption not desired. 
      if ( (job/10) .ne. 1)  return
c     
      write(iounit,127) key, type, title
      write(iounit,130) nrow,ncol,nnz
 127  format('.sp 4'/'.ll 7i'/'.ps 12'/'.po 0.7i'/'.ce 3'/,
     *     'Matrix:  ',a8,',  Type:  ',a3,/,a72)
 130  format('Dimension: ',i4,' x ',i4,',  Nonzero elements: ',i5)
      return
c----------------end-of-pltmt ------------------------------------------
c----------------------------------------------------------------------- 
      end
c-----------------------------------------------------------------------
      subroutine smms (n,first,last,mode,a,ja,ia,iout)
      integer ia(*), ja(*), n, first, last, mode, iout
      real*8 a(*)
c-----------------------------------------------------------------------
c writes a matrix in Coordinate (SMMS) format -- 
c-----------------------------------------------------------------------
c on entry:
c---------
c n     = integer = size of matrix -- number of rows (columns if matrix
c         is stored columnwise) 
c first  = first row (column) to be output. This routine will output 
c          rows (colums) first to last. 
c last   = last row (column) to be output. 
c mode   = integer giving some information about the storage of the 
c          matrix. A 3-digit decimal number. 'htu' 
c         * u = 0 means that matrix is stored row-wise 
c         * u = 1 means that matrix is stored column-wise 
c         * t = 0 indicates that the matrix is stored in CSR format 
c         * t = 1 indicates that the matrix is stored in MSR format. 
c         * h = ... to be added. 
c a,
c ja,
c ia    =  matrix in CSR or MSR format (see mode) 
c iout  = output unit number.
c
c on return:
c----------
c the output file iout will have written in it the matrix in smms
c (coordinate format)
c
c-----------------------------------------------------------------------
        logical msr, csc 
c
c determine mode ( msr or csr )
c
        msr = .false.
        csc = .false. 
	if (mod(mode,10) .eq. 1) csc = .true.
        if ( (mode/10) .eq. 1) msr = .true. 
        
        write (iout,*) n
      do 2 i=first, last 
         k1=ia(i)
         k2 = ia(i+1)-1
c     write (iout,*) ' row ', i 
         if (msr) write(iout,'(2i6,e22.14)')  i, i, a(i) 
         do 10 k=k1, k2
            if (csc) then
               write(iout,'(2i6,e22.14)')  ja(k), i, a(k)
            else
               write(iout,'(2i6,e22.14)')  i, ja(k), a(k)
            endif 
 10      continue
 2    continue
c----end-of-smms--------------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine readsm (nmax,nzmax,n,nnz,ia,ja,a,iout,ierr)
      integer nmax, nzmax, row, n, iout, i, j, k, ierr
      integer ia(nmax+1), ja(nzmax)
      real*8  a(nzmax), x
c-----------------------------------------------------------------------
c     read a matrix in coordinate format as is used in the SMMS
c     package (F. Alvarado), i.e. the row is in ascending order.
c     Outputs the matrix in CSR format.
c-----------------------------------------------------------------------
c coded by Kesheng Wu on Oct 21, 1991 with the supervision of Y. Saad
c-----------------------------------------------------------------------
c on entry:
c---------
c nmax  = the maximum size of array
c nzmax = the maximum number of nonzeros
c iout  = the I/O unit that has the data file
c
c on return:
c---------- 
c n     = integer = size of matrix
c nnz   = number of non-zero entries in the matrix
c a,
c ja, 
c ia    = matrix in CSR format
c ierr  = error code,
c         0 -- subroutine end with intended job done
c         1 -- error in I/O unit iout
c         2 -- end-of-file reached while reading n, i.e. a empty data file
c         3 -- n non-positive or too large
c         4 -- nnz is zero or larger than nzmax
c         5 -- data file is not orgnized in the order of ascending
c              row indices
c
c in case of errors:
c   n will be set to zero (0). In case the data file has more than nzmax
c   number of entries, the first nzmax entries will be read, and are not 
c   cleared on return. The total number of entry is determined.
c   Ierr is set.
c-----------------------------------------------------------------------
c
      rewind(iout)
      nnz = 0
      ia(1) = 1
      row = 1
c
      read (iout,*, err=1000, end=1010) n
      if ((n.le.0) .or. (n.gt.nmax)) goto 1020
c
 10   nnz = nnz + 1
      read (iout, *, err=1000, end=100) i, j, x

c     set the pointers when needed
      if (i.gt.row) then
         do 20 k = row+1, i
            ia(k) = nnz
 20      continue
         row = i
      else if (i.lt.row) then
         goto 1040
      endif

      ja(nnz) = j
      a (nnz) = x

      if (nnz.lt.nzmax) then
         goto 10
      else
         goto 1030
      endif

c     normal return -- end of file reached
 100  ia(row+1) = nnz
      nnz = nnz - 1
      if (nnz.eq.0) goto 1030
c
c     everything seems to be OK.
c
      ierr = 0
      return
c
c     error handling code
c
c     error in reading data entries
c
 1000 ierr = 1
      goto 2000
c
c     empty file
c
 1010 ierr  = 2
      goto 2000
c
c     problem with n
c
 1020 ierr = 3
      goto 2000
c
c     problem with nnz
c
 1030 ierr = 4
c
c     try to determine the real number of entries, in case needed
c
      if (nnz.ge.nzmax) then
 200     read(iout, *, err=210, end=210) i, j, x
         nnz = nnz + 1
         goto 200
 210     continue
      endif
      goto 2000
c
c     data entries not ordered
c
 1040 ierr = 5
 2000 n = 0
      return
c----end-of-readsm------------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine readsk (nmax,nzmax,n,nnz,a,ja,ia,iounit,ierr)
      integer nmax, nzmax, iounit, n, nnz, i, ierr
      integer ia(nmax+1), ja(nzmax) 
      real*8 a(nzmax)
c-----------------------------------------------------------------------
c Reads matrix in Compressed Saprse Row format. The data is supposed to
c appear in the following order -- n, ia, ja, a
c Only square matrices accepted. Format has following features
c (1) each number is separated by at least one space (or end-of-line), 
c (2) each array starts with a new line.
c-----------------------------------------------------------------------
c coded by Kesheng Wu on Oct 21, 1991 with supervision of Y. Saad
c-----------------------------------------------------------------------
c on entry:
c---------
c nmax 	 = max column dimension  allowed for matrix.
c nzmax	 = max number of nonzeros elements allowed. the arrays a, 
c          and ja should be of length equal to nnz (see below).
c iounit = logical unit number where to read the matrix from.
c
c on return:
c---------- 
c ia,
c ja,
c a      = matrx in CSR format
c n      = number of rows(columns) in matrix
c nnz	 = number of nonzero elements in A. This info is returned
c          even if there is not enough space in a, ja, ia, in order
c          to determine the minimum storage needed.
c ierr   = error code,
c          0 : OK;
c          1 : error when try to read the specified I/O unit.
c          2 : end-of-file reached during reading of data file.
c          3 : array size in data file is negtive or larger than nmax;
c          4 : nunmer of nonzeros in data file is negtive or larger than nzmax
c in case of errors:
c---------
c     n is set to 0 (zero), at the same time ierr is set.
c-----------------------------------------------------------------------
c     
c     read the size of the matrix
c
      rewind(iounit)
      read (iounit, *, err=1000, end=1010) n
      if ((n.le.0).or.(n.gt.nmax)) goto 1020
c     
c     read the pointer array ia(*)
c     
      read (iounit, *, err=1000, end=1010) (ia(i), i=1, n+1)
c     
c     Number of None-Zeros
c     
      nnz = ia(n+1) - 1
      if ((nnz.le.0).or.(nnz.gt.nzmax)) goto 1030
c     
c     read the column indices array
c     
      read (iounit, *, err=1000, end=1010) (ja(i), i=1, nnz)
c     
c     read the matrix elements
c     
      read (iounit, *, err=1000, end=1010) (a(i), i=1, nnz)
c     
c     normal return
c
      ierr = 0
      return
c     
c     error handling code
c     
c     error in reading I/O unit
 1000 ierr = 1
      goto 2000
c
c     EOF reached in reading
 1010 ierr =2
      goto 2000
c
c     n non-positive or too large
 1020 ierr = 3
      n = 0
      goto 2000
c
c     NNZ non-positive or too large
 1030 ierr = 4
c     
c     the real return statement
c     
 2000 n = 0
      return
c---------end of readsk ------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine skit (n, a, ja, ia, ifmt, iounit, ierr)
c-----------------------------------------------------------------------
c     Writes a matrix in Compressed Sparse Row format to an I/O unit.
c     It tryes to pack as many number as possible into lines of less than
c     80 characters. Space is inserted in between numbers for separation
c     to avoid carrying a header in the data file. This can be viewed
c     as a simplified Harwell-Boeing format. 
c-----------------------------------------------------------------------
c Modified from subroutine prtmt written by Y. Saad
c-----------------------------------------------------------------------
c on entry:
c---------
c n      = number of rows(columns) in matrix
c a      = real*8 array containing the values of the matrix stored 
c          columnwise
c ja     = integer array of the same length as a containing the column
c          indices of the corresponding matrix elements of array a.
c ia     = integer array of containing the pointers to the beginning of 
c          the row in arrays a and ja.
c ifmt   = integer specifying the format chosen for the real values
c          to be output (i.e., for a, and for rhs-guess-sol if 
c          applicable). The meaning of ifmt is as follows.
c          * if (ifmt .lt. 100) then the D descriptor is used,
c          format Dd.m, in which the length (m) of the mantissa is 
c          precisely the integer ifmt (and d = ifmt+6)
c          * if (ifmt .gt. 100) then prtmt will use the 
c          F- descriptor (format Fd.m) in which the length of the 
c          mantissa (m) is the integer mod(ifmt,100) and the length 
c          of the integer part is k=ifmt/100 (and d = k+m+2)
c          Thus  ifmt= 4   means  D10.4  +.xxxxD+ee    while
c          ifmt=104  means  F7.4   +x.xxxx
c          ifmt=205  means  F9.5   +xx.xxxxx
c          Note: formats for ja, and ia are internally computed.
c     
c iounit = logical unit number where to write the matrix into.
c     
c on return:
c----------
c ierr   = error code, 0 for normal 1 for error in writing to iounit.
c
c on error:
c--------
c     If error is encontacted when writing the matrix, the whole matrix
c     is written to the standard output.
c     ierr is set to 1.
c-----------------------------------------------------------------------
      character ptrfmt*16,indfmt*16,valfmt*20
      integer iounit, n, ifmt, len, nperli, nnz, i, ihead
      integer ja(*), ia(*), ierr
      real*8 a(*)
c--------------
c     compute pointer format
c--------------
      nnz    = ia(n+1)
      len    = int ( alog10(0.1+real(nnz))) + 2
      nnz    = nnz - 1
      nperli = 80/len
      
      print *, ' skit entries:', n, nnz, len, nperli

      if (len .gt. 9) then
         assign 101 to ix
      else
         assign 100 to ix
      endif
      write (ptrfmt,ix) nperli,len
 100  format(1h(,i2,1HI,i1,1h) )
 101  format(1h(,i2,1HI,i2,1h) )
c----------------------------
c     compute ROW index format
c----------------------------
      len    = int ( alog10(0.1+real(n) )) + 2
      nperli = min0(80/len,nnz)
      write (indfmt,100) nperli,len
c---------------------------
c     compute value format
c---------------------------
      if (ifmt .ge. 100) then
         ihead = ifmt/100
         ifmt = ifmt-100*ihead
         len = ihead+ifmt+3
         nperli = 80/len
c     
         if (len .le. 9 ) then
            assign 102 to ix
         elseif (ifmt .le. 9) then
            assign 103 to ix
         else 
            assign 104 to ix
         endif
c     
         write(valfmt,ix) nperli,len,ifmt
 102     format(1h(,i2,1hF,i1,1h.,i1,1h) )
 103     format(1h(,i2,1hF,i2,1h.,i1,1h) )
 104     format(1h(,i2,1hF,i2,1h.,i2,1h) )
C     
      else
         len = ifmt + 7
         nperli = 80/len
c     try to minimize the blanks in the format strings.
         if (nperli .le. 9) then
	    if (len .le. 9 ) then
	       assign 105 to ix
	    elseif (ifmt .le. 9) then
	       assign 106 to ix
	    else 
	       assign 107 to ix
	    endif
	 else 
	    if (len .le. 9 ) then
	       assign 108 to ix
	    elseif (ifmt .le. 9) then
	       assign 109 to ix
	    else 
               assign 110 to ix
            endif
         endif
c-----------
         write(valfmt,ix) nperli,len,ifmt
 105     format(1h(,i1,1hD,i1,1h.,i1,1h) )
 106     format(1h(,i1,1hD,i2,1h.,i1,1h) )
 107     format(1h(,i1,1hD,i2,1h.,i2,1h) )
 108     format(1h(,i2,1hD,i1,1h.,i1,1h) )
 109     format(1h(,i2,1hD,i2,1h.,i1,1h) )
 110     format(1h(,i2,1hD,i2,1h.,i2,1h) )
c     
      endif 	    
c     
c     output the data
c     
      write(iounit, *) n
      write(iounit,ptrfmt,err=1000) (ia(i), i = 1, n+1)
      write(iounit,indfmt,err=1000) (ja(i), i = 1, nnz)
      write(iounit,valfmt,err=1000) ( a(i), i = 1, nnz)
c
c     done, if no trouble is encounted in writing data
c
      ierr = 0
      return
c     
c     if can't write the data to the I/O unit specified, should be able to
c     write everything to standard output (unit 6)
c     
 1000 write(0, *) 'Error, Can''t write data to sepcified unit',iounit
      write(0, *) 'Write the matrix into standard output instead!'
      ierr = 1
      write(6,*) n
      write(6,ptrfmt) (ia(i), i=1, n+1)
      write(6,indfmt) (ja(i), i=1, nnz)
      write(6,valfmt) ( a(i), i=1, nnz)
      return
c----------end of skit ------------------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine prtunf(n, a, ja, ia, iout, ierr)
c-----------------------------------------------------------------------
c This subroutine dumps the arrays used for storing sparse compressed row
c format in machine code, i.e. unformatted using standard FORTRAN term.
c-----------------------------------------------------------------------
c First coded by Kesheng Wu on Oct 21, 1991 under the instruction of
c Prof. Y. Saad
c-----------------------------------------------------------------------
c On entry:
c     n: the size of the matrix (matrix is n X n)
c    ia: integer array stores the stariting position of each row.
c    ja: integer array stores the column indices of each entry.
c     a: the non-zero entries of the matrix.
c  iout: the unit number opened for storing the matrix.
c On return:
c  ierr: a error, 0 if everything's OK, else 1 if error in writing data.
c On error:
c  set ierr to 1.
c  No redirection is made, since direct the machine code to the standard
c output may cause unpridictable consequences.
c-----------------------------------------------------------------------
      integer iout, n, nnz, ierr, ia(*), ja(*)
      real*8  a(*)
      nnz = ia(n+1)-ia(1) 
c
      write(unit=iout, err=1000)  n
      write(unit=iout, err=1000) (ia(k),k=1,n+1) 
      if (nnz .gt. 0) then
         write(unit=iout, err=1000) (ja(k),k=1,nnz)
         write(unit=iout, err=1000) ( a(k),k=1,nnz)
      endif 
c
      ierr = 0
      return
c
 1000 ierr = 1
      return
      end
c---------end of prtunf ------------------------------------------------
c
c-----------------------------------------------------------------------
      subroutine readunf(nmax,nzmax,n,nnz,a,ja,ia,iounit,ierr)
c-----------------------------------------------------------------------
c This subroutine reads a matix store in machine code (FORTRAN
c unformatted form). The matrix is in CSR format.
c-----------------------------------------------------------------------
c First coded by Kesheng Wu on Oct 21, 1991 under the instruction of
c Prof. Y. Saad
c-----------------------------------------------------------------------
c On entry:
c    nmax: the maximum value of matrix size.
c   nzmax: the maximum number of non-zero entries.
c  iounit: the I/O unit that opened for reading.
c On return:
c       n: the actual size of array.
c     nnz: the actual number of non-zero entries.
c ia,ja,a: the matrix in CSR format.
c    ierr: a error code, it's same as that used in reaadsk
c          0 -- OK
c          1 -- error in reading iounit
c          2 -- end-of-file reached while reading data file
c          3 -- n is non-positive or too large
c          4 -- nnz is non-positive or too large
c On error:
c     return with n set to 0 (zero). nnz is kept if it's set already,
c     in case one want to use it to determine the size of array needed
c     to hold the data.
c-----------------------------------------------------------------------
c
      integer nmax, nzmax, n, iounit, nnz, k
      integer ia(nmax+1), ja(nzmax)
      real*8  a(nzmax)
c
      rewind iounit
c      
      read (unit=iounit, err=1000, end=1010) n
      if ((n.le.0) .or. (n.gt.nmax)) goto 1020
c
      read(unit=iounit, err=1000, end=1010) (ia(k),k=1,n+1)
c
      nnz = ia(n+1) - 1
      if ((nnz.le.0) .or. (nnz.gt.nzmax)) goto 1030
c
      read(unit=iounit, err=1000, end=1010) (ja(k),k=1,nnz)
      read(unit=iounit, err=1000, end=1010) (a(k),k=1,nnz)
c
c     everything seems to be OK.
c
      ierr = 0
      return
c
c     error handling
c
 1000 ierr = 1
      goto 2000
 1010 ierr = 2
      goto 2000
 1020 ierr = 3
      goto 2000
 1030 ierr = 4
 2000 n = 0
      return
      end
c---------end of readunf ----------------------------------------------
