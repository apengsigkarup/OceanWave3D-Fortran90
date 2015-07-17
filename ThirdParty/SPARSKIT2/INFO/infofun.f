c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c                  INFORMATION ROUTINES. INFO MODULE                   c
c----------------------------------------------------------------------c
c bandwidth :  Computes  ml     = lower_bandwidth(A)                   c
c                        mu     = upper_bandwidth(A)                   c
c                        iband  = max_bandwidth(A)                     c
c                        bndav  = average_bandwidth(A)                 c
c nonz      :  Computes  nzmaxc = max_column_length(A)                 c
c                        nzminc = min_column_length(A)                 c
c                        nzmaxr = max_row_length(A)                    c
c                        nzminr = min_row_length(A)                    c
c                        nzcol  = zero_column_number(A)                c
c                        nzrow  = zero_row_number(A)                   c
c diag_domi :  Computes  ddomc  = diag_domi_column_percentage(A)       c
c                        ddomr  = diag_domi_row_percentage(A)          c
c frobnorm  :  Computes  Fnorm  = Frobenius_norm(A)                    c
c ansym     :  Computes  fas    = sym_part_Frobenius_norm(A)           c
c                        fan    = nonsym_part_Frobenius_norm(A)        c
c                        imatch = matching_elements_number(A)          c
c                        av     = relative_sym_match(A)                c
c distaij   :  Computes  dist   = average_dist_of_a(i,j)(A)            c
c                        std    = standard_deviation(A)                c
c skyline   :  Computes  nsky   = nonzero_number_in_skyline(A)         c
c distdiag  :  Computes  dist   = element_number_in_eachdiag(A)        c
c bandpart  :  Computes  band   = bandwidth_width(A)                   c
c n_imp_diag:  Computes  ndiag  = important_diag_number(A)             c
c nonz_lud  :  Computes  nlower = nonzero_number_of_lower_part(A)      c
c                        nupper = nonzero_number_of_upper_part(A)      c
c                        ndiag  = nonzero_number_of_maindiag(A)        c
c avnz_col  :  Computes  av     = average_nonzero_number_in_column(A)  c
c                        st     = standard_deviation(A)                c
c vbrinfo   :  Print info on matrix in variable block row format       c
c----------------------------------------------------------------------c
      subroutine bandwidth(n,ja, ia, ml, mu, iband, bndav)
      implicit none 
      integer n, ml, mu, iband
      integer ja(*), ia(n+1)
      real*8 bndav
c-----------------------------------------------------------------------
c this routine computes the lower, upper, maximum, and average 
c bandwidths.     revised -- July 12, 2001  -- bug fix -- YS. 
c-----------------------------------------------------------------------
c On Entry:
c----------
c n     = integer. column dimension of matrix
c a     = real array containing the nonzero elements of the matrix
c         the elements are stored by columns in order
c         (i.e. column i comes before column i+1, but the elements
c         within each column can be disordered).
c ja    = integer array containing the row indices of elements in a
c ia    = integer array containing of length n+1 containing the
c         pointers to the beginning of the columns in arrays a and ja.
c         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
c
c on return
c----------
c ml    = lower bandwidth as defined by
c        ml = max(i-j | all  a(i,j).ne. 0)
c mu    = upper bandwidth as defined by
c        mu = max ( j-i | all  a(i,j).ne. 0 )
c iband =  maximum bandwidth as defined by
c         iband = Max (  Max [ j | a(i,j) .ne. 0 ] - 
c                        Min [ j | a(i,j) .ne. 0 ] )
c bndav = Average Bandwidth          
c-----------------------------------------------------------------------
c     locals
      integer max 
      integer j0, j1, jminc, jmaxc, i
c-----------------------------------------------------------------------
      ml = -n
      mu = -n
      bndav = 0.0d0
      iband = 0 
      do 10 i=1,n
         j0 = ia(i)
         j1 = ia(i+1) - 1
         jminc = ja(j0)
         jmaxc = ja(j1)
         ml = max(ml,i-jminc)
         mu = max(mu,jmaxc-i)
         iband = max(iband,jmaxc-jminc+1)
         bndav = bndav+real( jmaxc-jminc+1)
 10   continue
      bndav = bndav/real(n)
      return
c-----end-of-bandwidth--------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine nonz(n,sym, ja, ia, iao, nzmaxc, nzminc, 
     *     nzmaxr, nzminr, nzcol, nzrow)
      implicit none 
      integer n , nzmaxc, nzminc, nzmaxr, nzminr, nzcol, nzrow
      integer ja(*), ia(n+1), iao(n+1)
      logical sym
c----------------------------------------------------------------------
c     this routine computes maximum numbers of nonzero elements 
c     per column/row, minimum numbers of nonzero elements per column/row, 
c     and  numbers of zero columns/rows.
c----------------------------------------------------------------------
c     On Entry:
c----------
c     n     = integer column dimension of matrix
c     ja    = integer array containing the row indices of elements in a
c     ia    = integer array containing of length n+1 containing the
c     pointers to the beginning of the columns in arrays a and ja.
c     It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
c     iao   = similar array for the transpose of the matrix. 
c     sym   = logical variable indicating whether or not the matrix is 
c     stored in symmetric mode.
c     on return
c----------
c     nzmaxc = max length of columns
c     nzminc = min length of columns 
c     nzmaxr = max length of rows
c nzminr = min length of rows
c     nzcol  = number of zero columns
c     nzrow = number of zero rows
c-----------------------------------------------------------------------
      integer  i, j0, j0r, j1r, indiag, k, j1, lenc, lenr 
c
      nzmaxc = 0
      nzminc = n
      nzmaxr = 0
      nzminr = n
      nzcol  = 0
      nzrow  = 0
c-----------------------------------------------------------------------
      do 10 i = 1, n
         j0 = ia(i)
         j1 = ia(i+1) 
         j0r = iao(i)
         j1r = iao(i+1)
         indiag = 0
         do 20 k=j0, j1-1 
            if (ja(k) .eq. i) indiag = 1
 20      continue
c         
         lenc = j1-j0
         lenr = j1r-j0r
c         
         if (sym) lenc = lenc + lenr - indiag
         if (lenc .le. 0) nzcol = nzcol +1
         nzmaxc = max0(nzmaxc,lenc)
         nzminc = min0(nzminc,lenc)
         if (lenr .le. 0) nzrow = nzrow+1
         nzmaxr = max0(nzmaxr,lenr)
         nzminr = min0(nzminr,lenr)
 10   continue
      return
      end
c-----------------------------------------------------------------------
      subroutine diag_domi(n,sym,valued,a, ja,ia,ao,jao, iao, 
     *     ddomc, ddomr)
      implicit none 
      real*8 a(*), ao(*), ddomc, ddomr
      integer n, ja(*), ia(n+1), jao(*), iao(n+1)
      logical sym, valued
c-----------------------------------------------------------------
c     this routine computes the percentage of weakly diagonally 
c     dominant rows/columns
c-----------------------------------------------------------------
c     on entry:
c     ---------
c     n     = integer column dimension of matrix
c a     = real array containing the nonzero elements of the matrix
c     the elements are stored by columns in order
c     (i.e. column i comes before column i+1, but the elements
c     within each column can be disordered).
c     ja    = integer array containing the row indices of elements in a
c ia    = integer array containing of length n+1 containing the
c     pointers to the beginning of the columns in arrays a and ja.
c     It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
c     ao    = real array containing the nonzero elements of the matrix
c     the elements are stored by rows in order
c     (i.e. row i comes before row i+1, but the elements
c     within each row can be disordered).
c ao,jao, iao, 
c     structure for transpose of a 
c     sym   = logical variable indicating whether or not the matrix is
c     symmetric.
c     valued= logical equal to .true. if values are provided and .false.
c         if only the pattern of the matrix is provided. (in that
c     case a(*) and ao(*) are dummy arrays.
c     
c     ON RETURN
c----------
c     ddomc = percentage of weakly diagonally dominant columns
c     ddomr = percentage of weakly diagonally dominant rows
c-------------------------------------------------------------------
c     locals
      integer i, j0, j1, k, j 
      real*8 aii, dsumr, dsumc 
c     number of diagonally dominant columns
c     real arithmetic used to avoid problems.. YS. 03/27/01 
      ddomc = 0.0  
c     number of diagonally dominant rows
      ddomr = 0.0 
      do 10 i = 1, n
         j0 = ia(i)
         j1 = ia(i+1) - 1
         if (valued) then
            aii = 0.0d0
            dsumc = 0.0d0
            do 20 k=j0,j1
               j = ja(k) 
               if (j .eq. i) then
                  aii = abs(a(k))
               else
                  dsumc = dsumc + abs(a(k))
               endif
 20         continue
            dsumr = 0.0d0
            if (.not. sym) then
               do 30 k=iao(i), iao(i+1)-1
                  if (jao(k) .ne. i) dsumr = dsumr+abs(ao(k))
 30            continue 
            else
               dsumr = dsumc
            endif
            if (dsumc .le. aii) ddomc = ddomc + 1.0
            if (dsumr .le. aii) ddomr = ddomr + 1.0
         endif
 10   continue
      ddomr = ddomr / real(n)
      ddomc = ddomc / real(n)
      return
c-----------------------------------------------------------------------
c--------end-of-diag_moni-----------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine frobnorm(n,sym,a,ja,ia,Fnorm)
      implicit none 
      integer n 
      real*8 a(*),Fnorm
      integer ja(*),ia(n+1)
      logical sym
c--------------------------------------------------------------------------
c     this routine computes the Frobenius norm of A.
c--------------------------------------------------------------------------
c     on entry:
c-----------
c n      = integer colum dimension of matrix
c a      = real array containing the nonzero elements of the matrix
c          the elements are stored by columns in order
c          (i.e. column i comes before column i+1, but the elements
c          within each column can be disordered).
c ja     = integer array containing the row indices of elements in a.
c ia     = integer array containing of length n+1 containing the
c          pointers to the beginning of the columns in arrays a and 
c          ja. It is assumed that ia(*)= 1 and ia(n+1)  = nnz +1.
c sym    = logical variable indicating whether or not the matrix is
c          symmetric.
c
c on return
c-----------
c Fnorm  = Frobenius norm of A.
c--------------------------------------------------------------------------
      real*8 Fdiag
      integer i, k 
      Fdiag = 0.0
      Fnorm = 0.0
      do i =1,n
         do k = ia(i), ia(i+1)-1
            if (ja(k) .eq. i) then
               Fdiag = Fdiag + a(k)**2
            else
               Fnorm = Fnorm + a(k)**2
            endif
         enddo 
      enddo 
      if (sym) then
         Fnorm = 2*Fnorm +Fdiag
      else
        Fnorm = Fnorm + Fdiag
      endif
      Fnorm = sqrt(Fnorm)
      return
      end
c----------------------------------------------------------------------
      subroutine ansym(n,sym,a,ja,ia,ao,jao,iao,imatch,
     *     av,fas,fan)
c---------------------------------------------------------------------
c     this routine computes the Frobenius norm of the symmetric and
c     non-symmetric parts of A, computes number of matching elements
c     in symmetry and relative symmetry match. 
c---------------------------------------------------------------------
c on entry:
c----------
c n   = integer column dimension of matrix
c a   = real array containing the nonzero elements of the matrix
c       the elements are stored by columns in order
c       (i.e. column i comes before column i+1, but the elements
c       within each column can be disordered).
c ja  = integer array containing the row indices of elements in a
c ia  = integer array containing of length n+1 containing the
c       pointers to the beginning of the columns in arrays a and ja.
c       It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
c sym = logical variable indicating whether or not the matrix is
c       symmetric.
c on return
c----------
c fas   = Frobenius norm of symmetric part
c fan   = Frobenius norm of non-symmetric part
c imatch = number of matching elements in symmetry
c av     = relative symmetry match (symmetry = 1)
c ao,jao,iao = transpose of A just as a, ja, ia contains 
c              information of A.
c-----------------------------------------------------------------------
      implicit real*8 (a-h, o-z)
      real*8 a(*),ao(*),fas,fan,av, Fnorm, st
      integer n, ja(*), ia(n+1), jao(*), iao(n+1),imatch
      logical sym
c-----------------------------------------------------------------------
      nnz    = ia(n+1)-ia(1)
      call csrcsc(n,1,1,a,ja,ia,ao,jao,iao)
      if (sym) goto 7
      st     = 0.0d0
      fas    = 0.0d0
      fan    = 0.0d0
      imatch = 0
      do 6 i=1,n
         k1 = ia(i)
         k2 = iao(i)
         k1max = ia(i+1) - 1
         k2max = iao(i+1) - 1
c     
 5       if (k1 .gt. k1max .or. k2 .gt. k2max) goto 6
c     
         j1 = ja(k1)
         j2 = jao(k2)
         if (j1 .ne. j2 ) goto 51
         fas = fas + (a(k1)+ao(k2))**2
         fan = fan + (a(k1)-ao(k2))**2
         st  = st + a(k1)**2
         imatch = imatch + 1
 51      k1 = k1+1
         k2 = k2+1
         if (j1 .lt. j2)  k2 = k2 - 1
         if (j1 .gt. j2)  k1 = k1 - 1
         goto 5
 6    continue
      fas = 0.25D0 * fas
      fan = 0.25D0 * fan
 7    call frobnorm(n,sym,ao,jao,iao,Fnorm)
      if (sym) then
         imatch = nnz
         fas = Fnorm
         fan = 0.0d0
      else
         if (imatch.eq.nnz) then
            st = 0.0D0
         else
            st = 0.5D0 * (Fnorm**2 - st)
            if (st.lt.0.0D0) st = 0.0D0
         endif
         fas = sqrt(fas + st)
         fan = sqrt(fan + st)
      endif
      av = real(imatch)/real(nnz)
      return
      end
c------end-of-ansym-----------------------------------------------------
c-----------------------------------------------------------------------
      subroutine distaij(n,nnz,sym,ja,ia,dist, std)
      implicit real*8 (a-h, o-z)
      real*8  dist, std
      integer ja(*), ia(n+1)
c-----------------------------------------------------------------------
c     this routine computes the average distance of a(i,j) from diag and
c     standard deviation  for this average.
c-----------------------------------------------------------------------
c On entry :
c-----------
c n     = integer. column dimension of matrix
c nnz   = number of nonzero elements of matrix
c ja    = integer array containing the row indices of elements in a
c ia    = integer array containing of length n+1 containing the
c         pointers to the beginning of the columns in arrays a and ja.
c         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
c sym   = logical variable indicating whether or not the matrix is
c         symmetric.
c on return
c----------
c dist  = average distance of a(i,j) from diag.
c std   = standard deviation for above average.
c-----------------------------------------------------------------------
c
c distance of an element from diagonal.
c
      dist   = 0.0
      std = 0.0
      do 3 i=1,n
         j0 = ia(i)
         j1 = ia(i+1) - 1
         do 31 k=j0, j1
            j=ja(k)
            dist = dist + real(iabs(j-i) )
 31      continue
 3    continue
      dist = dist/real(nnz)
      do 6 i = 1, n 
         do 61 k=ia(i), ia(i+1) - 1
            std=std+(dist-real(iabs(ja(k)-i)))**2
 61      continue
 6    continue
      std = sqrt(std/ real(nnz))
      return
      end
c-----------------------------------------------------------------------
      subroutine skyline(n,sym,ja,ia,jao,iao,nsky)
      implicit real*8 (a-h, o-z)
      integer n, ja(*), ia(n+1), jao(*), iao(n+1)
      integer nskyl, nskyu, nsky
      logical sym
c-------------------------------------------------------------------
c this routine computes the number of nonzeros in the skyline storage.
c-------------------------------------------------------------------
c
c On entry :
c-----------
c n     = integer. column dimension of matrix
c ja    = integer array containing the row indices of elements in a
c ia    = integer array containing of length n+1 containing the
c         pointers to the beginning of the columns in arrays a and ja.
c         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
c iao   = integer array containing of length n+1 containing the
c         pointers to the beginning of the rows in arrays ao and jao.
c         It is assumed that iao(*) = 1 and iao(n+1) = nzz+1.
c jao   = integer array containing the column indices of elements in ao.
c sym   = logical variable indicating whether or not the matrix is
c         symmetric.
c on return
c----------
c nsky  = number of nonzeros in skyline storage
c-------------------------------------------------------------------
c
c nskyu = skyline storage for upper part
      nskyu = 0
c nskyl = skyline storage for lower part
      nskyl = 0
      do 10 i=1,n
         j0 = ia(i)
         j0r = iao(i)
         
         jminc = ja(j0)
         jminr = jao(j0r)
         if (sym) jminc = jminr
         
         nskyl = nskyl + i-jminr + 1
         nskyu = nskyu + i-jminc + 1
         
 10   continue
      nsky = nskyl+nskyu-n
      if (sym)  nsky = nskyl
      return
      end
c-----------------------------------------------------------------------
      subroutine distdiag(nrow,ncol,ja,ia,dist)
      implicit real*8 (a-h, o-z)
      integer nrow,ncol,ja(*), ia(nrow+1),dist(*)
c----------------------------------------------------------------------
c this routine computes the numbers of elements in each diagonal. 
c----------------------------------------------------------------------
c On entry :
c-----------
c nrow  = integer. row dimension of matrix
c ncol  = integer. column dimension of matrix
c ja    = integer array containing the row indices of elements in a
c ia    = integer array containing of length n+1 containing the
c         pointers to the beginning of the columns in arrays a and ja.
c         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
c on return
c----------
c dist  = integer array containing the numbers of elements in each of 
c         the nrow+ncol-1 diagonals of A. dist(k) contains the 
c         number of elements in diagonal '-nrow+k'.  k ranges from 
c         1 to (nrow+ncol-1).
c----------------------------------------------------------------------
      nnz  = ia(nrow+1)-ia(1)
      n2   = nrow+ncol-1
      do 8 i=1, n2
         dist(i) = 0
 8    continue
      do 9 i=1, nrow
         k1 = ia(i)
         k2 = ia(i+1) -1
         do 91 k=k1, k2
            j = ja(k)
            dist(nrow+j-i) = dist(nrow+j-i) +1
 91      continue
 9    continue
      return
      end
c-----------------------------------------------------------------------
      subroutine bandpart(n,ja,ia,dist,nper,band)
      implicit real*8 (a-h, o-z)
      integer n,ja(*), ia(n+1),dist(*)
      integer nper,band
c-------------------------------------------------------------------------
c this routine computes the bandwidth of the banded matrix, which contains
c 'nper' percent of the original matrix.
c-------------------------------------------------------------------------
c On entry :
c-----------
c n     = integer. column dimension of matrix
c ja    = integer array containing the row indices of elements in a
c ia    = integer array containing of length n+1 containing the
c         pointers to the beginning of the columns in arrays a and ja.
c         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
c dist  = integer array containing the numbers of elements in the 
c         matrix with different distance of row indices and column 
c         indices.
c nper  = percentage of matrix  within the bandwidth
c on return
c----------
c band  = the width of the bandwidth
c----------------------------------------------------------------------
      nnz  = ia(n+1)-ia(1)
      iacc  = dist(n)
      band  = 0
      j     = 0
 10   j     = j+1
      iacc  = iacc + dist(n+j) +dist(n-j) 
      if (iacc*100 .le. nnz*nper) then 
         band = band +1 
         goto 10         
      endif 
      return 
      end 
c-----------------------------------------------------------------------      
      subroutine  n_imp_diag(n,nnz,dist, ipar1,ndiag,ioff,dcount)
      implicit real*8 (a-h, o-z)
      real*8  dcount(*)
      integer n,nnz, dist(*), ndiag, ioff(*), ipar1
c-----------------------------------------------------------------------
c     this routine computes the most important diagonals.
c-----------------------------------------------------------------------
c
c On entry :
c-----------
c n     = integer. column dimension of matrix
c nnz   = number of nonzero elements of matrix
c dist  = integer array containing the numbers of elements in the 
c         matrix with different distance of row indices and column 
c         indices. ipar1 = percentage of nonzero  elements of A that 
c         a diagonal should have in order to be an important diagonal 
c on return
c----------
c ndiag = number of the most important diagonals
c ioff  = the offsets with respect to the main diagonal
c dcount= the accumulated percentages
c-----------------------------------------------------------------------
      n2 = n+n-1
      ndiag = 10
      ndiag = min0(n2,ndiag)
      itot  = 0
      ii = 0
      idiag = 0
c     sort diagonals by decreasing order of weights.
 40   jmax = 0
      i    = 1
      do 41 k=1, n2
         j = dist(k)
         if (j .lt. jmax) goto 41
         i = k
         jmax = j
 41   continue
c     permute ----
c     save offsets and accumulated count if diagonal is acceptable
c     (if it has at least ipar1*nnz/100 nonzero elements)
c     quite if no more acceptable diagonals --
c     
      if (jmax*100 .lt. ipar1*nnz) goto 4
      ii = ii+1
      ioff(ii) = i-n
      dist(i)   = - jmax
      itot = itot + jmax
      dcount(ii) = real(100*itot)/real(nnz)
      if (ii .lt. ndiag) goto 40
 4    continue
      ndiag = ii
      return
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
      subroutine nonz_lud(n,ja,ia,nlower, nupper, ndiag)
      implicit real*8 (a-h, o-z)
      integer n,  ja(*), ia(n+1)
      integer nlower, nupper, ndiag
c-----------------------------------------------------------------------
c this routine computes the number of nonzero elements in strict lower
c part, strict upper part, and main diagonal.
c-----------------------------------------------------------------------
c
c On entry :
c-----------
c n     = integer. column dimension of matrix
c ja    = integer array containing the row indices of elements in a
c ia    = integer array containing of length n+1 containing the
c         pointers to the beginning of the columns in arrays a and ja.
c         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
c on return
c----------
c nlower= number of nonzero elements in strict lower part
c nupper= number of nonzero elements in strict upper part
c ndiag = number of nonzero elements in main diagonal
c------------------------------------------------------------------- 
c
c number of nonzero elements in upper part
c
      nupper = 0
      ndiag = 0
      
      do 3 i=1,n
c     indiag = nonzero diagonal element indicator
         do 31 k=ia(i), ia(i+1)-1
            j=ja(k)
            if (j .lt. i) nupper = nupper+1
            if (j .eq. i) ndiag = ndiag + 1 
 31      continue
 3    continue
      nlower = ia(n+1)-1-nupper-ndiag
      return
      end
c-----------------------------------------------------------------------
      subroutine avnz_col(n,ja,ia,iao, ndiag, av, st)
      implicit real*8 (a-h, o-z)
      real*8 av, st
      integer n,  ja(*), ia(n+1), iao(n+1)
c---------------------------------------------------------------------
c     this routine computes average number of nonzero elements/column and
c     standard deviation for this average
c---------------------------------------------------------------------
c
c On entry :
c-----------
c n     = integer. column dimension of matrix
c ja    = integer array containing the row indices of elements in a
c ia    = integer array containing of length n+1 containing the
c         pointers to the beginning of the columns in arrays a and ja.
c         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
c ndiag = number of the most important diagonals
c On return
c----------
c av    = average number of nonzero elements/column
c st    = standard deviation for this average
c Notes
c---------
c standard deviation will not be correct for symmetric storage. 
c----------------------------------------------------------------------
c     standard deviatioan for the average
      st = 0.0d0
c     average and standard deviation
c     
      av = real(ia(n+1)-1)/real(n)
c     
c     will be corrected later.
c     
      do 3 i=1,n
         j0 = ia(i)
         j1 = ia(i+1) - 1
c     indiag = nonzero diagonal element indicator
         do 31 k=j0, j1
            j=ja(k)
 31      continue
         lenc = j1+1-j0
         st = st + (real(lenc) - av)**2
 3    continue
c     
      st = sqrt( st / real(n) )
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine vbrinfo(nr, nc, kvstr, kvstc, ia, ja, ka, iwk, iout)
c-----------------------------------------------------------------------
      integer nr, nc, kvstr(nr+1), kvstc(nc+1), ia(nr+1), ja(*), ka(*)
      integer iwk(nr+nc+2+nr), iout
c-----------------------------------------------------------------------
c     Print info on matrix in variable block row format.
c-----------------------------------------------------------------------
c     On entry:
c--------------
c     nr,nc   = matrix block row and block column dimension
c     kvstr   = first row number for each block row
c     kvstc   = first column number for each block column
c     ia,ja,ka,a = input matrix in VBR format
c     iout    = unit number for printed output
c
c     On return:
c---------------
c     Printed output to unit number specified in iout.  If a non-square
c     matrix is provided, the analysis will be performed on the block
c     rows, otherwise a row/column conformal partitioning will be used.
c
c     Work space:
c----------------
c     iwk(1:nb+1) = conformal block partitioning
c        (nb is unknown at start but is no more than nr+nc)
c     iwk(nb+2:nb+2+nr) = frequency of each blocksize
c     The workspace is not assumed to be initialized to zero, nor is it
c     left that way.
c
c-----------------------------------------------------------------------
c-----local variables
      integer n, nb, nnz, nnzb, i, j, neq, max, num
      character*101 tmpst
      integer bsiz(10), freq(10)
c-----------------------------------------------------------------------
      n = kvstr(nr+1)-1
      nnz = ka(ia(nr+1)) - ka(1)
      nnzb = ia(nr+1) - ia(1)
      write (iout, 96)
      write (iout, 100) n, nnz, real(nnz)/real(n)
      write (iout, 101) nr, nnzb, real(nnzb)/real(nr)
c-----if non-square matrix, do analysis on block rows,
c     else do analysis on conformal partitioning
      if (kvstr(nr+1) .ne. kvstc(nc+1)) then
         write (iout, 103)
         do i = 1, nr+1
            iwk(i) = kvstr(i)
         enddo
         nb = nr
      else
         call kvstmerge(nr, kvstr, nc, kvstc, nb, iwk)
         if ((nr .ne. nc) .or. (nc .ne. nb)) write (iout, 104) nb
      endif
c-----accumulate frequencies of each blocksize
      max = 1
      iwk(1+nb+2) = 0
      do i = 1, nb
         neq = iwk(i+1) - iwk(i)
         if (neq .gt. max) then
            do j = max+1, neq
               iwk(j+nb+2) = 0
            enddo
            max = neq
         endif
         iwk(neq+nb+2) = iwk(neq+nb+2) + 1
      enddo
c-----store largest 10 of these blocksizes
      num = 0
      do i = max, 1, -1
         if ((iwk(i+nb+2) .ne. 0) .and. (num .lt. 10)) then
            num = num + 1
            bsiz(num) = i
            freq(num) = iwk(i+nb+2)
         endif
      enddo
c-----print information about blocksizes
      write (iout, 109) num
      write (tmpst,'(10i6)') (bsiz(j),j=1,num)
      write (iout,110) num,tmpst
      write (tmpst,'(10i6)') (freq(j),j=1,num)
      write (iout,111) tmpst
      write (iout, 96)
c-----------------------------------------------------------------------
 99    format (2x,38(2h *))
 96    format (6x,' *',65(1h-),'*')
 100   format(
     * 6x,' *  Number of rows                                   = ',
     * i10,'  *'/
     * 6x,' *  Number of nonzero elements                       = ',
     * i10,'  *'/
     * 6x,' *  Average number of nonzero elements/Row           = ',
     * f10.4,'  *')
 101   format(
     * 6x,' *  Number of block rows                             = ',
     * i10,'  *'/
     * 6x,' *  Number of nonzero blocks                         = ',
     * i10,'  *'/
     * 6x,' *  Average number of nonzero blocks/Block row       = ',
     * f10.4,'  *')
 103   format(
     * 6x,' *  Non-square matrix.                                 ',
     *     '            *'/
     * 6x,' *  Performing analysis on block rows.                 ',
     *     '            *')
 104   format(
     * 6x,' *  Non row-column conformal partitioning supplied.    ',
     *     '            *'/
     * 6x,' *  Using conformal partitioning.  Number of bl rows = ',
     * i10,'  *')
 109   format(
     * 6x,' *  Number of different blocksizes                   = ',
     * i10,'  *')
 110   format(6x,' *  The ', i2, ' largest dimension nodes',
     *     ' have dimension    : ',10x,'  *',/,
     * 6x,' *',a61,3x,' *')
 111   format(6x,' *  The frequency of nodes these ',
     *     'dimensions are      : ',10x,'  *',/,
     * 6x,' *',a61,3x,' *')
c---------------------------------
      return
      end
c-----------------------------------------------------------------------
c-----------------------end-of-vbrinfo----------------------------------
