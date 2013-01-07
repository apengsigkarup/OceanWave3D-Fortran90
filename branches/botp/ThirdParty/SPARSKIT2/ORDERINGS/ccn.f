c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c          REORDERING ROUTINES -- STRONGLY CONNECTED COMPONENTS        c 
c----------------------------------------------------------------------c
c     Contributed by:
C     Laura C. Dutto - email: dutto@cerca.umontreal.ca
c                      July 1992 - Update: March 1994
C-----------------------------------------------------------------------
c     CONTENTS:
c     --------
c     blccnx : Driver routine to reduce the structure of a  matrix 
c              to its strongly connected components.
c     cconex : Main routine to compute the strongly connected components
c              of a (block diagonal) matrix.
c     anccnx : We put in ICCNEX the vertices marked in the component MCCNEX.
c     newcnx : We put in ICCNEX the vertices marked in the component
c              MCCNEX. We modify also the vector KPW.
c     blccn1 : Parallel computation of the connected components of a
c              matrix. The parallel loop is performed only if the matrix
c              has a block diagonal structure.
c     ccnicopy:We copy an integer vector into anothoer.
c     compos : We calculate the composition between two permutation
c              vectors.
c     invlpw : We calculate the inverse of a permutation vector.
c     numini : We initialize a vector to the identity.
c     tbzero : We initialize to ZERO an integer vector.
c     iplusa : Given two integers IALPHA and IBETA, for an integer vector 
c              IA we calculate IA(i) = ialpha + ibeta * ia(i)
C
c----------------------------------------------------------------------c
      subroutine BLCCNX(n, nbloc, nblcmx, nsbloc, job, lpw, amat, ja,
     *                  ia, iout, ier, izs, nw)
C-----------------------------------------------------------------------
c     
c     This routine determines if the matrix given by the structure
c     IA et JA is irreductible. If not, it orders the unknowns such
c     that all the consecutive unknowns in KPW between NSBLOC(i-1)+1
c     and NSBLOC(i) belong to the ith component of the matrix.
c     The numerical values of the matrix are in AMAT. They are modified
c     only if JOB = 1 and if we have more than one connected component.
c
c     On entry:
c     --------
c     n      = row and column dimension of the matrix
c     nblcmx = maximum number of connected components allowed. The size
c              of NSBLOC is nblcmx + 1 (in fact, it starts at 0).
c     job    = integer indicating the work to be done:
c              job = 1  if the permutation LPW is modified, we
c                       permute not only the structure of the matrix
c                       but also its numerical values.
c              job.ne.1 if the permutation LPW is modified, we permute 
c                       the structure of the matrix ignoring real values.
c     iout   = impression parameter. If 0 < iout < 100, we print
c              comments and error messages on unit IOUT.
c     nw     = length of the work vector IZS.
c
c     Input / output:
c     --------------
c     nbloc  = number of connected components of the matrix. If the
c              matrix is not irreductible, nbloc > 1. We allow
c              nbloc > 1 on entry; in this case we calculate the
c              number of connected components in each previous one.
c     nsbloc = integer array of length NBLOC + 1 containing the pointers
c              to the first node of each component on the old (input)
c              and on the new (output) ordering.
c     lpw    = integer array of length N corresponding to the
c              permutation of the unknowns. We allow LPW to be a vector 
c              different from the identity on input.
c     amat   = real*8 values of the matrix given by the structure IA, JA.
c     ja     = integer array of length NNZERO (= IA(N+1)-IA(1)) corresponding
c              to the column indices of nonzero elements of the matrix, stored
c              rowwise. It is modified only if the matrix has more
c              than one connected component.
c     ia     = integer array of length N+1 corresponding to the
c              pointer to the beginning of each row in JA (compressed
c              sparse row storage). It is modified only if
c              the matrix has more than one connected component.
c
c     On return:
c     ----------
c     ier    = integer. Error message. Normal return ier = 0.
c
c     Work space:
c     ----------
c     izs    = integer vector of length NW
c
C-----------------------------------------------------------------------
C     Laura C. Dutto - email: dutto@cerca.umontreal.ca
c                      July 1992 - Update: March 1994
C-----------------------------------------------------------------------
      integer izs(nw), lpw(n), nsbloc(0:nblcmx), ia(n+1), ja(*)
      real*8 amat(*)
      logical impr
      character*6 chsubr
C-----------------------------------------------------------------------
      ier    = 0
      impr   = iout.gt.0.and.iout.le.99
      ntb    = ia(n+1) - 1
      mxccex = max(nblcmx,20)
c.....The matrix AMAT is a real*8 vector
      ireal  = 2
c
c.....MXPTBL: maximal number of vertices by block
      mxptbl = 0
      do ibloc = 1, nbloc
         mxptbl = max( mxptbl, nsbloc(ibloc) - nsbloc(ibloc-1))
      enddo
c
      long1 = nbloc * mxptbl
      long2 = nbloc * (mxccex+1)
c.....Dynamic allocation of memory
      iend   = 1
      iiend  = iend
      ilpw   = iiend 
      ikpw   = ilpw   + n
      ilccnx = ikpw   + long1
      imark  = ilccnx + long2
      iend   = imark  + n
      if(iend .gt. nw) go to 220
c
      nbloc0 = nbloc
      chsubr = 'BLCCN1'
c.....We determine if the matrix has more than NBLOC0 connected components.
      call BLCCN1(n, nbloc, nblcmx, nsbloc, izs(ilpw), izs(ikpw), ia,
     *            ja, izs(imark), mxccex, izs(ilccnx), mxptbl, iout,
     *            ier)
      if(ier.ne.0) go to 210
c
      if(nbloc .gt. nbloc0) then
c..........The matrix has more than NBLOC0 conneted components. So, we
c..........modify the vectors IA and JA to take account of the new permutation.
           nfree = iend - ikpw
           call tbzero(izs(ikpw), nfree)
           iiat  = ikpw
           ijat  = iiat + n + 1
           iamat = ijat + ntb
           iend  = iamat
           if(job .eq. 1) iend = iamat + ireal * ntb
           if(iend .gt. nw) go to 220
c
c..........We copy IA and JA on IAT and JAT respectively
           call ccnicopy(n+1, ia, izs(iiat))
           call ccnicopy(ntb, ja, izs(ijat))
           if(job .eq. 1) call dcopy(ntb, amat, 1, izs(iamat), 1)
           call dperm(n, izs(iamat), izs(ijat), izs(iiat), amat,
     *                ja, ia, izs(ilpw), izs(ilpw), job)
           ipos = 1
c..........We sort columns inside JA.
           call csrcsc(n, job, ipos, amat, ja, ia, izs(iamat),
     *                 izs(ijat), izs(iiat))
           call csrcsc(n, job, ipos, izs(iamat), izs(ijat), izs(iiat),
     *                 amat, ja, ia)
      endif
c.....We modify the ordering of unknowns in LPW
      call compos(n, lpw, izs(ilpw))
c
 120  nfree = iend - iiend
      call tbzero(izs(iiend), nfree)
      iend = iiend
      return
c
 210  IF(IMPR) WRITE(IOUT,310) chsubr,ier
      go to 120
 220  IF(IMPR) WRITE(IOUT,320) nw, iend
      if(ier.eq.0) ier = -1
      go to 120
c
 310  FORMAT(' ***BLCCNX*** ERROR IN ',a6,'. IER = ',i8)
 320  FORMAT(' ***BLCCNX*** THERE IS NOT ENOUGH MEMORY IN THE WORK',
     1       ' VECTOR.'/13X,' ALLOWED MEMORY = ',I10,'  - NEEDED',
     2       ' MEMORY = ',I10)
      end 
c **********************************************************************
      subroutine CCONEX(n, icol0, mxccnx, lccnex, kpw, ia, ja, mark,
     *                  iout, ier)
C-----------------------------------------------------------------------
c
c     This routine determines if the matrix given by the structure
c     IA and JA is irreductible. If not, it orders the unknowns such
c     that all the consecutive unknowns in KPW between LCCNEX(i-1)+1
c     and LCCNEX(i) belong to the ith component of the matrix.
c     The structure of the matrix could be nonsymmetric.
c     The diagonal vertices (if any) will belong to the last connected
c     component (convention).
c
c     On entry:
c     --------
c     n      = row and column dimension of the matrix
c     icol0  = the columns of the matrix are between ICOL0+1 and ICOL0+N
c     iout   = impression parameter. If 0 < IOUT < 100, we print
c              comments and error messages on unit IOUT.
c     ia     = integer array of length N+1 corresponding to the
c              pointer to the beginning of each row in JA (compressed
c              sparse row storage).
c     ja     = integer array of length NNZERO (= IA(N+1)-IA(1))
c              corresponding to the column indices of nonzero elements
c              of the matrix, stored rowwise.
c
c     Input/Output:
c     ------------
c     mxccnx = maximum number of connected components allowed on input,
c              and number of connected components of the matrix, on output.
c
c     On return:
c     ----------
c     lccnex = integer array of length MXCCNX + 1 containing the pointers
c              to the first node of each component, in the vector KPW.
c     kpw    = integer array of length N corresponding to the
c              inverse of permutation vector.
c     ier    = integer. Error message. Normal return ier = 0.
c
c     Work space:
c     ----------
c     mark   = integer vector of length N
c
C-----------------------------------------------------------------------
C     Laura C. Dutto - email: dutto@cerca.umontreal.ca
c                      July 1992 - Update: March 1994
C-----------------------------------------------------------------------
      dimension ia(n+1), lccnex(0:mxccnx), kpw(n), ja(*), mark(n)
      logical impr
C-----------------------------------------------------------------------
      ier    = 0
      ipos = ia(1) - 1 
      impr   = iout.gt.0.and.iout.le.99
c
      nccnex = 0
c.....We initialize MARK to zero. At the end of the algorithm, it would
c.....indicate the number of connected component associated with the vertex.
c.....The number (-1) indicates that the row associated with this vertex
c.....is a diagonal row. This value could be modified because we accept
c.....a non symmetric matrix. All the diagonal vertices will be put in
c.....the same connected component.
      call tbzero(mark, n)
c
 5    do i = 1,n
         if(mark(i) .eq. 0) then
              ideb = i
              go to 15
         endif
      enddo
      go to 35
c
 15   if( ia(ideb+1) - ia(ideb) .eq. 1) then
c..........The row is a diagonal row.
           mark(ideb) = -1
           go to 5
      endif
      iccnex = nccnex + 1
      if(iccnex .gt. mxccnx) go to 220
      index  = 0
      newind = 0
      jref   = 0
      mark(ideb) = iccnex
      index  = index + 1
      kpw(index) = ideb
c
 20   jref = jref + 1
      ideb = kpw(jref)

      do 30 ir = ia(ideb)-ipos, ia(ideb+1)-ipos-1
         j = ja(ir) - icol0
         mccnex = mark(j)
         if(mccnex .le. 0) then
              index = index + 1
              kpw(index) = j
              mark(j) = iccnex
         else if( mccnex .eq. iccnex) then
              go to 30
         else if( mccnex .gt. iccnex) then
c.............We realize that the connected component MCCNX is,
c.............in fact, included in this one. We modify MARK and KPW.
              call NEWCNX(n, mccnex, iccnex, index, kpw, mark)
              if(mccnex .eq. nccnex) nccnex = nccnex - 1
         else
c.............We realize that the previously marked vertices belong,
c.............in fact, to the connected component ICCNX. We modify MARK.
              call ANCCNX(n, iccnex, mccnex, mark, nwindx)
              iccnex = mccnex
              newind = newind + nwindx
         endif
 30   continue
      if(jref .lt. index) go to 20
c
c.....We have finished with this connected component.
      index = index + newind
      if(iccnex .eq. nccnex+1) nccnex = nccnex + 1
      go to 5
c.......................................................................
c
c     We have partitioned the graph in its connected components!
c
c.......................................................................
 35   continue
c
c.....All the vertices have been already marked. Before modifying KPW
c.....(if necessary), we put the diagonal vertex (if any) in the last
c.....connected component.
      call tbzero(lccnex(1), nccnex)
c
      idiag = 0
      do i = 1, n
         iccnex = mark(i)
         if(iccnex .eq. -1) then
              idiag = idiag + 1
              if(idiag .eq. 1) then
                   nccnex = nccnex + 1
                   if(nccnex .gt. mxccnx) go to 220
                   if(impr) write(iout,340)
              endif
              mark(i) = nccnex
         else
              lccnex(iccnex) = lccnex(iccnex) + 1
         endif
      enddo
      if(idiag .ge. 1) lccnex(nccnex) = idiag
c
      if(nccnex .eq. 1) then
         lccnex(nccnex) = n
         go to 40
      endif
c
      iccnex = 1
 8    if(iccnex .gt. nccnex) go to 12
      if(lccnex(iccnex) .le. 0) then
           do i = 1, n
              if(mark(i) .ge. iccnex) mark(i) = mark(i) - 1
           enddo
           nccnex = nccnex - 1
           do mccnex = iccnex, nccnex
              lccnex(mccnex) = lccnex(mccnex + 1)
           enddo
      else
           iccnex = iccnex + 1
      endif
      go to 8
c
 12   index = 0
      do iccnex = 1, nccnex
         noeicc = lccnex(iccnex)
         lccnex(iccnex) = index
         index  = index + noeicc
      enddo
      if(index .ne. n) go to 210
c
c.....We define correctly KPW
      do i = 1,n
         iccnex = mark(i)
         index = lccnex(iccnex) + 1
         kpw(index) = i
         lccnex(iccnex) = index
      enddo
c
 40   mxccnx = nccnex
      lccnex(0) = nccnex
      if(nccnex .eq. 1) call numini(n, kpw)
      return
c
 210  if(impr) write(iout,310) index,n
      go to 235
 220  if(impr) write(iout,320) nccnex, mxccnx
      go to 235
 235  ier = -1
      return
c
 310  format(' ***CCONEX*** ERROR TRYING TO DETERMINE THE NUMBER',
     *       ' OF CONNECTED COMPONENTS.'/13X,' NUMBER OF MARKED',
     *       ' VERTICES =',i7,3x,'TOTAL NUMBER OF VERTICES =',I7)
 320  format(' ***CCONEX*** THE ALLOWED NUMBER OF CONNECTED COMPONENTS',
     *       ' IS NOT ENOUGH.'/13X,' NECESSARY NUMBER = ',I4,
     *       5x,' ALLOWED NUMBER = ',I4)
 323  format(' ***CCONEX*** ERROR IN ',A6,'. IER = ',I8)
 340  format(/' ***CCONEX*** THE LAST CONNECTED COMPONENT WILL',
     *       ' HAVE THE DIAGONAL VERTICES.')
      end
c **********************************************************************
      subroutine ANCCNX(n, mccnex, iccnex, mark, ncount)
C-----------------------------------------------------------------------
c     
c     We put in ICCNEX the vertices marked in the component MCCNEX.
C
C-----------------------------------------------------------------------
c     include "NSIMPLIC"
      dimension mark(n)
C-----------------------------------------------------------------------
C     Laura C. Dutto - email: dutto@cerca.umontreal.ca - December 1993
C-----------------------------------------------------------------------
      ncount = 0
      do i = 1, n
         if( mark(i) .eq. mccnex) then
              mark(i) = iccnex
              ncount  = ncount + 1
         endif
      enddo
c
      return
      end
c **********************************************************************
      subroutine NEWCNX(n, mccnex, iccnex, index, kpw, mark)
C-----------------------------------------------------------------------
c     
c     We put in ICCNEX the vertices marked in the component MCCNEX. We
c     modify also the vector KPW.
C
C-----------------------------------------------------------------------
c     include "NSIMPLIC"
      dimension kpw(*), mark(n)
C-----------------------------------------------------------------------
C     Laura C. Dutto - email: dutto@cerca.umontreal.ca - December 1993
C-----------------------------------------------------------------------
      do i = 1, n
         if( mark(i) .eq. mccnex) then
              mark(i) = iccnex
              index = index + 1
              kpw(index) = i
         endif
      enddo
c
      return
      end
c **********************************************************************
      subroutine BLCCN1(n, nbloc, nblcmx, nsbloc, lpw, kpw, ia, ja,
     *                  mark, mxccex, lccnex, mxptbl, iout, ier)
C-----------------------------------------------------------------------
c     
c     This routine determines if the matrix given by the structure
c     IA et JA is irreductible. If not, it orders the unknowns such
c     that all the consecutive unknowns in KPW between NSBLOC(i-1)+1
c     and NSBLOC(i) belong to the ith component of the matrix.
c
c     On entry:
c     --------
c     n      = row and column dimension of the matrix
c     nblcmx = The size of NSBLOC is nblcmx + 1 (in fact, it starts at 0).
c     ia     = integer array of length N+1 corresponding to the
c              pointer to the beginning of each row in JA (compressed
c              sparse row storage).
c     ja     = integer array of length NNZERO (= IA(N+1)-IA(1)) corresponding
c              to the column indices of nonzero elements of the matrix,
c              stored rowwise.
c     mxccex = maximum number of connected components allowed by block.
c     mxptbl = maximum number of points (or unknowns) in each connected
c              component (mxptbl .le. n).
c     iout   = impression parameter. If 0 < iout < 100, we print
c              comments and error messages on unit IOUT.
c
c     Input/Output:
c     ------------
c     nbloc  = number of connected components of the matrix. If the
c              matrix is not irreductible, nbloc > 1. We allow
c              nbloc > 1 on entry; in this case we calculate the
c              number of connected components in each previous one.
c     nsbloc = integer array of length NBLOC + 1 containing the pointers
c              to the first node of each component on the new ordering.
c              Normally, on entry you put: NBLOC = 1, NSBLOC(0) = 0,
c              NSBLOC(NBLOC) = N.
c
c     On return:
c     ----------
c     lpw    = integer array of length N corresponding to the
c              permutation vector (the row i goes to lpw(i)).
c     ier    = integer. Error message. Normal return ier = 0.
c
c     Work space:
c     ----------
c     kpw    = integer vector of length MXPTBL*NBLOC necessary for parallel
c              computation.
c     mark   = integer vector of length N
c     lccnex = integer vector of length (MXCCEX+1)*NBLOC necessary for parallel
c              computation.
c
C-----------------------------------------------------------------------
C     Laura C. Dutto - e-mail: dutto@cerca.umontreal.ca
c                      Juillet 1992. Update: March 1994
C-----------------------------------------------------------------------
      dimension lpw(n), kpw(mxptbl*nbloc), ia(n+1), ja(*),
     *          lccnex((mxccex+1)*nbloc), nsbloc(0:nbloc), mark(n)
      logical impr
      character chsubr*6
C-----------------------------------------------------------------------
      ier  = 0
      impr = iout.gt.0.and.iout.le.99
      isor = 0
c
      chsubr = 'CCONEX'
      newblc = 0
C$DOACROSS if(nbloc.gt.1), LOCAL(ibloc, ik0, ik1, ins0, ntb0,
C$&   nccnex, ilccnx, info, kpibl), REDUCTION(ier, newblc)
      do 100 ibloc = 1,nbloc
         ik0 = nsbloc(ibloc - 1)
         ik1 = nsbloc(ibloc)
         ntb0 = ia(ik0+1)
         if(ia(ik1+1) - ntb0 .le. 1) go to 100
         ntb0 = ntb0 - 1
         ins0 = ik1 - ik0
c........We need more memory place for KPW1 because of parallel computation
         kpibl = (ibloc-1) * mxptbl
         call numini( ins0, kpw(kpibl+1))
         nccnex = mxccex
         ilccnx = (mxccex+1) * (ibloc-1) + 1
c.......................................................................
c
c        Call to the main routine: CCONEX
c
c.......................................................................
         call cconex(ins0, ik0, nccnex, lccnex(ilccnx), kpw(kpibl+1), 
     *               ia(ik0+1), ja(ntb0+1), mark(ik0+1), isor, info)
         ier = ier + info
         if(info .ne. 0 .or. nccnex .lt. 1) go to 100
c
c........We add the new connected components on NEWBLC
         newblc = newblc + nccnex
c........We define LPW different from the identity only if there are more 
c........than one connected component in this block
         if(nccnex .eq. 1) then
              call numini(ins0, lpw(ik0+1))
         else
              call invlpw(ins0, kpw(kpibl+1), lpw(ik0+1))
         endif
         call iplusa(ins0, ik0, 1, lpw(ik0+1))
 100  continue
c
      if(ier .ne. 0) go to 218
      if(newblc .eq. nbloc) go to 120
      if(newblc .gt. nblcmx) go to 230
c
c.....We modify the number of blocks to indicate the number of connected 
c.....components in the matrix.
      newblc = 0
      nsfin = 0
CDIR$ NEXT SCALAR
      do ibloc = 1, nbloc
         ilccnx = (mxccex+1) * (ibloc-1) + 1
         nccnex = lccnex(ilccnx)
         if(nccnex .gt. 1 .and. impr) write(iout,420) ibloc,nccnex
         lcc0 = 0
CDIR$ NEXT SCALAR
         do icc = 1,nccnex
            newblc = newblc + 1
            nsb = lccnex(ilccnx+icc)
c...........Be careful! In LCCNEX we have the cumulated number of vertices 
            nsbloc(newblc) = nsfin + nsb
            if(nccnex .gt. 1 .and. impr) write(iout,425) icc,nsb-lcc0
            lcc0 = nsb
         enddo
         nsfin = nsfin + nsb
      enddo
      nbloc = newblc
c
 120  return
c
 218  if(impr) write(iout,318) chsubr,ier
      go to 120
 230  if(impr) write(iout,330) newblc,nblcmx
      if(ier.eq.0) ier = -1
      go to 120
c
 318  format(' ***BLCCN1*** ERROR IN ',a6,'. IER = ',i8)
 330  format(' ***BLCCN1*** THE MEMORY SPACE ALLOWED FOR NSBLOC IS',
     *       ' NOT ENOUGH.'/13X,' NUMBER (NECESSARY) OF CONNECTED',
     *       ' COMPONENTS = ',I5/13X,' MAXIMAL NUMBER OF BLOCKS',14x,
     *       '= ',i5)
 420  FORMAT(' *** The block ',i3,' has ',i3,' strongly connected',
     *       ' components. The number of vertices by component is:')
 425  format(5x,'Component No.',i3,' - Number of vertices = ',i6)
      end
C***********************************************************************
      SUBROUTINE CCNICOPY(N,IX,IY)
C.......................................................................
C     We copy the vector IX on the vector IY
C.......................................................................
      DIMENSION IX(n),IY(n)
C.......................................................................
      IF(N.LE.0) RETURN
C$DOACROSS if(n .gt. 250), local(i)
      DO 10 I = 1,N
         IY(I) = IX(I)
   10 CONTINUE
C
      RETURN
      END
c***********************************************************************
      SUBROUTINE COMPOS(n, lpw0, lpw1)
C-----------------------------------------------------------------------
c
c     We take account of the original order of unknowns. We put the
c     final result on LPW0.
c
C-----------------------------------------------------------------------
      DIMENSION lpw0(n), lpw1(n)
C-----------------------------------------------------------------------
c     Laura C. Dutto - Mars 1994
C-----------------------------------------------------------------------
C$DOACROSS if(n .gt. 250), local(i0)
      do i0 = 1, n
         lpw0(i0) = lpw1(lpw0(i0))
      enddo
c
      return
      end
C **********************************************************************
      SUBROUTINE INVLPW(n, lpw, kpw)
c.......................................................................
c 
c     KPW is the inverse of LPW
c
c.......................................................................
      dimension lpw(n), kpw(n)
c.......................................................................
c     Laura C. Dutto - Novembre 1993
c.......................................................................
C$DOACROSS if(n .gt. 200), local(i0, i1)
      do i0 = 1, n
         i1 = lpw(i0)
         kpw(i1) = i0
      enddo
c
      return
      end
C **********************************************************************
      subroutine NUMINI(n, lpw)
c.......................................................................
      dimension lpw(n)
c.......................................................................
c
c     The vector LPW is initialized as the identity.
c
c.......................................................................
c     Laura C. Dutto - Novembre 1993
c.......................................................................
C$DOACROSS if(n .gt. 250), local(i)
      do i=1,n
         lpw(i) = i
      enddo
c
      return
      end
C***********************************************************************
      SUBROUTINE TBZERO(M,NMOT)
C.......................................................................
C     We initialize to ZERO an integer vector of length NMOT.
C.......................................................................
      DIMENSION M(NMOT)
C.......................................................................
      IF(NMOT.le.0) return
C$DOACROSS if(nmot.gt.500), LOCAL(i)
      DO 1 I=1,NMOT
           M(I)=0
 1         CONTINUE
      RETURN
      END
C **********************************************************************
      SUBROUTINE IPLUSA (n, nalpha, nbeta, ia)
c.......................................................................
C
c     We add NALPHA to each element of NBETA * IA:
c
c            ia(i) = nalpha + nbeta * ia(i)
c
c.......................................................................
      integer ia(n)
c.......................................................................
c     Laura C. Dutto - February 1994
c.......................................................................
      if(n .le. 0) return
c
      nmax = 500
      if(nalpha .eq. 0) then
           if(nbeta .eq. 1) return
           if(nbeta .eq. -1) then
C$DOACROSS if(n .gt. nmax), local (i)
                do i = 1, n
                   ia(i) = - ia(i)
                enddo
           else
C$DOACROSS if(n .gt. nmax/2), local (i)
                do i = 1, n
                   ia(i) = nbeta * ia(i)
                enddo
           endif
           return
      endif
      if(nbeta .eq. 0) then
C$DOACROSS if(n .gt. nmax), local (i)
           do i = 1, n
              ia(i) = nalpha
           enddo
           return
      endif
      if(nbeta .eq. -1) then
C$DOACROSS if(n .gt. nmax/2), local (i)
           do i = 1, n
              ia(i) = nalpha - ia(i)
           enddo
      else if(nbeta .eq. 1) then
C$DOACROSS if(n .gt. nmax/2), local (i)
           do i = 1, n
              ia(i) = nalpha + ia(i)
           enddo
      else
C$DOACROSS if(n .gt. nmax/3), local (i)
           do i = 1, n
              ia(i) = nalpha + nbeta * ia(i)
           enddo
      endif
c
      return
      end
