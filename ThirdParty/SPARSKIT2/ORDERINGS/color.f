c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c          REORDERING ROUTINES -- COLORING BASED ROUTINES              c 
c----------------------------------------------------------------------c
c contents:                                                            c
c----------                                                            c
c multic  : greedy algorithm for multicoloring                         c
c indset0 : greedy algorithm for independent set ordering              c
c indset1 : independent set ordering using minimal degree traversal    c
c indset2 : independent set ordering with local minimization           c
c indset3 : independent set ordering by vertex cover algorithm         c
c HeapSort, FixHeap, HeapInsert, interchange, MoveBack, FiHeapM,       c
c           FixHeapM, HeapInsertM,indsetr,rndperm, are utility         c
c           routines for sorting, generating random permutations, etc. c
c----------------------------------------------------------------------c
      subroutine multic (n,ja,ia,ncol,kolrs,il,iord,maxcol,ierr) 
      integer n, ja(*),ia(n+1),kolrs(n),iord(n),il(maxcol+1),ierr 
c-----------------------------------------------------------------------
c     multicoloring ordering -- greedy algorithm -- 
c     determines the coloring permutation and sets up
c     corresponding data structures for it.
c-----------------------------------------------------------------------
c on entry
c -------- 
c n     = row and column dimention of matrix
c ja    = column indices of nonzero elements of matrix, stored rowwise.
c ia    = pointer to beginning of each row in ja.
c maxcol= maximum number of colors allowed -- the size of il is
c         maxcol+1 at least. Note: the number of colors does not
c         exceed the maximum degree of each node +1.
c iord  = en entry iord gives the order of traversal of the nodes
c         in the multicoloring algorithm. If there is no preference 
c         then set iord(j)=j for j=1,...,n
c
c on return
c --------- 
c ncol  = number of colours found 
c kolrs = integer array containing the color number assigned to each node 
c il    = integer array containing the pointers to the
c         beginning of each color set. In the permuted matrix
c         the rows /columns il(kol) to il(kol+1)-1 have the same color.
c iord  = permutation array corresponding to the multicolor ordering.
c         row number i will become row nbumber iord(i) in permuted 
c         matrix. (iord = destination permutation array).
c ierr  = integer. Error message. normal return ierr = 0. If ierr .eq.1
c         then the array il was overfilled. 
c 
c-----------------------------------------------------------------------
c     
      integer kol, i, j, k, maxcol, mycol 
c     
      ierr = 0
      do 1 j=1, n
         kolrs(j) = 0
 1    continue 
      do 11 j=1, maxcol
         il(j) = 0
 11   continue
c     
      ncol = 0
c     
c     scan all nodes 
c     
      do 4 ii=1, n
         i = iord(ii) 
c     
c     look at adjacent nodes to determine colors already assigned
c     
         mcol = 0
         do 2 k=ia(i), ia(i+1)-1
            j = ja(k)
            icol = kolrs(j)
            if (icol .ne. 0) then
               mcol = max(mcol,icol) 
c     
c     il used as temporary to record already assigned colors.
c     
               il(icol) = 1 
            endif
 2       continue
c     
c     taken colors determined. scan il until a slot opens up.
c     
         mycol = 1
 3       if (il(mycol) .eq. 1) then
            mycol = mycol+1 
            if (mycol .gt. maxcol) goto 99
            if (mycol .le. mcol) goto 3
         endif
c     
c     reset il to zero for next nodes
c     
         do 35 j=1, mcol
            il(j) = 0
 35      continue
c     
c     assign color and update number of colors so far
c     
         kolrs(i) = mycol
         ncol = max(ncol,mycol)
 4    continue
c     
c     every node has now been colored. Count nodes of each color
c     
      do 6 j=1, n
         kol = kolrs(j)+1
         il(kol) = il(kol)+1
 6    continue
c     
c     set pointers il
c     
      il(1) = 1
      do 7 j=1, ncol
         il(j+1) = il(j)+il(j+1)
 7    continue
c     
c     set iord
c     
      do 8 j=1, n
         kol = kolrs(j) 
         iord(j) = il(kol)
         il(kol) = il(kol)+1
 8    continue
c     
c     shift il back 
c     
      do 9 j=ncol,1,-1
         il(j+1) = il(j)
 9    continue
      il(1) = 1
c     
      return
 99   ierr = 1
      return
c----end-of-multic------------------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------

      subroutine indset0 (n,ja,ia,nset,iord,riord,sym,iptr) 
      integer n, nset, ja(*),ia(*),riord(*),iord(*) 
      logical sym
c---------------------------------------------------------------------- 
c greedy algorithm for independent set ordering
c----------------------------------------------------------------------
c parameters:
c ----------
c n      = row dimension of matrix
c ja, ia = matrix pattern in CRS format
c nset   = (output) number of elements in the independent set
c iord   = permutation array corresponding to the independent set 
c          ordering. Row number i will become row number iord(i) in 
c          permuted matrix.
c riord  = reverse permutation array. Row number i in the permutated 
c          matrix is row number riord(i) in original matrix. 
c----------------------------------------------------------------------
c notes: works for CSR, MSR, and CSC formats but assumes that the
c matrix has a symmetric structure. 
c---------------------------------------------------------------------- 
c local variables
c 
      integer j, k1, k2, nod, k, mat
      do 1 j=1, n
         iord(j) = 0
 1    continue
      nummat = 1
      if (.not. sym) nummat = 2
c     
c     iord used as a marker
c     
      nset = 0
      do 12  nod=1, n
         if (iord(nod) .ne. 0) goto 12 
         nset = nset+1
         iord(nod) = 1
c     
c     visit all neighbors of current nod 
c     
         ipos = 0
         do 45 mat=1, nummat
            do 4 k=ia(ipos+nod), ia(ipos+nod+1)-1 
               j = ja(k)
               if (j .ne. nod) iord(j) = 2
 4          continue
            ipos = iptr-1
 45      continue
 12   continue
c     
c     get permutation
c     
      k1 = 0 
      k2 = nset 
      do 6 j=1,n
         if (iord(j) .eq. 1) then
            k1 = k1+1
            k = k1
         else
            k2 = k2+1
            k = k2 
         endif
         riord(k) = j         
         iord(j) = k
 6    continue
      return
c---------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------
      subroutine indset1 (n,ja,ia,nset,iord,riord,iw,sym,iptr) 
      integer n, nset, iptr, ja(*),ia(*),riord(*),iord(*),iw(*)
      logical sym
c---------------------------------------------------------------------- 
c greedy algorithm for independent set ordering -- with intial 
c order of traversal given by that of min degree. 
c----------------------------------------------------------------------
c parameters:
c ----------
c n      = row dimension of matrix
c ja, ia = matrix pattern in CRS format
c nset   = (output) number of elements in the independent set
c iord   = permutation array corresponding to the independent set 
c          ordering. Row number i will become row number iord(i) in 
c          permuted matrix.
c riord  = reverse permutation array. Row number i in the permutated 
c          matrix is row number riord(i) in original matrix. 
c----------------------------------------------------------------------
c notes: works for CSR, MSR, and CSC formats but assumes that the
c matrix has a symmetric structure. 
c---------------------------------------------------------------------- 
c local variables
      integer j,k1,k2,nummat,nod,k,ipos 
c
c     nummat is the number of matrices to loop through (A in symmetric
c     pattern case (nummat=1) or A,and transp(A) otherwise (mummat=2) 
c
      if (sym) then
         nummat = 1
      else
         nummat = 2
      endif
      iptrm1 = iptr-1
c     
c     initialize arrays
c
      do 1 j=1,n
         iord(j) = j 
         riord(j) = j 
         iw(j) = 0 
 1    continue
c     
c     initialize degrees of all nodes
c     
      ipos = 0
      do 100 imat =1,nummat
         do 15 j=1,n
            iw(j) = iw(j) + ia(ipos+j+1)-ia(ipos+j) 
 15      continue
         ipos = iptrm1
 100  continue 
c     
c     call heapsort -- sorts nodes in increasing degree. 
c     
      call HeapSort (iw,iord,riord,n,n) 
c     
c     weights no longer needed -- use iw to store order of traversal.
c     
      do 16 j=1, n
         iw(n-j+1) = iord(j)
         iord(j) = 0
 16   continue
c     
c     iord used as a marker
c     
      nset = 0
      do 12  ii = 1, n
         nod = iw(ii) 
         if (iord(nod) .ne. 0) goto 12 
         nset = nset+1
         iord(nod) = 1
c     
c     visit all neighbors of current nod 
c     
         ipos = 0
         do 45 mat=1, nummat
            do 4 k=ia(ipos+nod), ia(ipos+nod+1)-1 
               j = ja(k)
               if (j .ne. nod) iord(j) = 2
 4          continue
            ipos = iptrm1 
 45      continue
 12   continue
c     
c     get permutation
c     
      k1 = 0 
      k2 = nset 
      do 6 j=1,n
         if (iord(j) .eq. 1) then
            k1 = k1+1
            k = k1
         else
            k2 = k2+1
            k = k2 
         endif
         riord(k) = j         
         iord(j) = k
 6    continue
      return
c----------------------------------------------------------------------
      end
c----------------------------------------------------------------------
      subroutine indset2(n,ja,ia,nset,iord,riord,iw,sym,iptr) 
      integer n,nset,iptr,ja(*),ia(*),riord(n),iord(n),iw(n) 
      logical sym 
c---------------------------------------------------------------------- 
c greedy algorithm for independent set ordering -- local minimization
c using heap strategy -- 
c----------------------------------------------------------------------
c This version for BOTH unsymmetric and symmetric patterns
c----------------------------------------------------------------------
c on entry
c -------- 
c n     = row and column dimension of matrix
c ja    = column indices of nonzero elements of matrix,stored rowwise.
c ia    = pointer to beginning of each row in ja.
c sym   = logical indicating whether the matrix has a symmetric pattern.
c         If not the transpose must also be provided -- appended to the
c         ja, ia structure -- see description of iptr next.
c iptr  = in case the matrix has an unsymmetric pattern,the transpose
c         is assumed to be stored in the same arrays ia,ja. iptr is the
c         location in ia of the pointer to the first row of transp(A).
c         more generally, ia(iptr),...,ia(iptr+n) are the pointers to 
c         the beginnings of rows 1, 2, ...., n+1 (row n+1 is fictitious)
c         of the transpose of A in the array ja. For example,when using 
c         the msr format,one can write:
c          iptr = ja(n+1)
c          ipos = iptr+n+2                ! get the transpose of A:
c          call csrcsc (n,0,ipos,a,ja,ja,a,ja,ja(iptr))    ! and then:
c          call indset(n,ja,ja,nset,iord,riord,iwk,.false.,iptr) 
c
c iw    = work space of length n.
c
c on return: 
c---------- 
c nset  = integer. The number of unknowns in the independent set. 
c iord  = permutation array corresponding to the new ordering. The 
c         first nset unknowns correspond to the independent set.
c riord = reverse permutation array.  
c----------------------------------------------------------------------
c local variables --
c 
      integer j,k1,k2,nummat,nod,k,ipos,i,last,lastlast,jold,jnew,
     *     jo,jn 
c
c     nummat is the number of matrices to loop through (A in symmetric
c     pattern case (nummat=1) or A,and transp(A) otherwise (mummat=2) 
c
      if (sym) then
         nummat = 1
      else
         nummat = 2
      endif
      iptrm1 = iptr-1
c
c     initialize arrays
c
      do 1 j=1,n
         iord(j) = j
         riord(j) = j
         iw(j) = 0 
 1    continue
c
c     initialize degrees of all nodes
c
      ipos = 0
      do 100 imat =1,nummat
         do 15 j=1,n
            iw(j) = iw(j) + ia(ipos+j+1)-ia(ipos+j) 
 15      continue
 100     ipos = iptrm1
c     
c start by constructing a heap
c 
      do 2 i=n/2,1,-1 
         j = i
         call FixHeap (iw,iord,riord,j,j,n) 
 2    continue
c     
c main loop -- remove nodes one by one. 
c 
      last = n
      nset = 0
 3    continue
      lastlast = last
      nod = iord(1) 
c     
c     move first element to end
c 
      call moveback (iw,iord,riord,last) 
      last = last -1 
      nset = nset + 1 
c
c     scan all neighbors of accepted node -- move them to back -- 
c     
      ipos = 0
      do 101 imat =1,nummat      
         do 5 k=ia(ipos+nod),ia(ipos+nod+1)-1
            jold = ja(k)
            jnew = riord(jold)
            if (jold .eq. nod .or. jnew .gt. last) goto 5 
            iw(jnew) = -1
            call HeapInsert (iw,iord,riord,jnew,ichild,jnew) 
            call moveback (iw,iord,riord,last) 
            last = last -1 
 5       continue 
         ipos = iptrm1
 101  continue
c
c update the degree of each edge
c 
         do 6 k=last+1,lastlast-1
            jold = iord(k) 
c     
c     scan the neighbors of current node
c     
         ipos = 0
         do 102 imat =1,nummat
            do 61 i=ia(ipos+jold),ia(ipos+jold+1)-1 
               jo = ja(i) 
               jn = riord(jo) 
c
c     consider this node only if it has not been moved          
c     
               if (jn .gt. last) goto 61
c     update degree of this neighbor 
               iw(jn) = iw(jn)-1
c     and fix the heap accordingly
               call HeapInsert (iw,iord,riord,jn,ichild,jn)
 61         continue
            ipos = iptrm1
 102     continue
 6    continue
c
c     stopping test -- end main "while"loop 
c
      if (last .gt. 1) goto 3
      nset = nset + last
c
c     rescan all nodes one more time to determine the permutations 
c
      k1 = 0 
      k2 = nset
      do 7 j=n,1,-1
         if (iw(j) .ge. 0) then
            k1 = k1+1
            k = k1
         else
            k2 = k2+1
            k = k2 
         endif
         riord(k) = iord(j) 
 7    continue
      do j=1,n
         iord(riord(j)) = j
      enddo
      return
c----------------------------------------------------------------------
      end
c----------------------------------------------------------------------
      subroutine indset3(n,ja,ia,nset,iord,riord,iw,sym,iptr) 
      integer n,nset,iptr,ja(*),ia(*),riord(n),iord(n),iw(n) 
      logical sym 
c---------------------------------------------------------------------- 
c greedy algorithm for independent set ordering -- local minimization
c using heap strategy -- VERTEX COVER ALGORITHM --  
c ASSUMES MSR FORMAT (no diagonal element) -- ADD A SWITCH FOR CSR -- 
c----------------------------------------------------------------------
c This version for BOTH unsymmetric and symmetric patterns
c----------------------------------------------------------------------
c on entry
c -------- 
c n     = row and column dimension of matrix
c ja    = column indices of nonzero elements of matrix,stored rowwise.
c ia    = pointer to beginning of each row in ja.
c sym   = logical indicating whether the matrix has a symmetric pattern.
c         If not the transpose must also be provided -- appended to the
c         ja, ia structure -- see description of iptr next.
c iptr  = in case the matrix has an unsymmetric pattern,the transpose
c         is assumed to be stored in the same arrays ia,ja. iptr is the
c         location in ia of the pointer to the first row of transp(A).
c         more generally, ia(iptr),...,ia(iptr+n) are the pointers to 
c         the beginnings of rows 1, 2, ...., n+1 (row n+1 is fictitious)
c         of the transpose of A in the array ja. For example,when using 
c         the msr format,one can write:
c          iptr = ja(n+1)
c          ipos = iptr+n+2                ! get the transpose of A:
c          call csrcsc (n,0,ipos,a,ja,ja,a,ja,ja(iptr))    ! and then:
c          call indset(n,ja,ja,nset,iord,riord,iwk,.false.,iptr) 
c
c iw    = work space of length n.
c
c on return: 
c---------- 
c nset  = integer. The number of unknowns in the independent set. 
c iord  = permutation array corresponding to the new ordering. The 
c         first nset unknowns correspond to the independent set.
c riord = reverse permutation array.  
c----------------------------------------------------------------------
c local variables --
c 
      integer j,nummat,nod,k,ipos,i,lastnset,jold,jnew
c
c     nummat is the number of matrices to loop through (A in symmetric
c     pattern case (nummat=1) or A,and transp(A) otherwise (mummat=2) 
c
      if (sym) then
         nummat = 1
      else
         nummat = 2
      endif
      iptrm1 = iptr-1
c
c     initialize arrays
c
      do 1 j=1,n
         riord(j) = j
         iord(j) = j
         iw(j) = 0 
 1    continue
c
c     initialize degrees of all nodes
c
      nnz = 0
      ipos = 0
      do 100 imat =1,nummat
         do 15 j=1,n
            ideg = ia(ipos+j+1)-ia(ipos+j) 
            iw(j) = iw(j) + ideg 
            nnz = nnz + ideg
 15      continue
 100     ipos = iptrm1
c
c     number of edges
c     
         if (sym) then nnz = 2*nnz
c     
c start by constructing a Max heap
c 
      do 2 i=n/2,1,-1 
         j = i
         call FixHeapM (iw,riord,iord,j,j,n) 
 2    continue
      nset = n
c----------------------------------------------------------------------     
c main loop -- remove nodes one by one. 
c----------------------------------------------------------------------
 3    continue
      lastnset = nset
      nod = riord(1) 
c     
c     move first element to end
c 
      call movebackM (iw,riord,iord,nset) 
      nnz = nnz - iw(nset) 
      nset = nset -1 
c
c     scan all neighbors of accepted node -- 
c     
      ipos = 0
      do 101 imat =1,nummat      
         do 5 k=ia(ipos+nod),ia(ipos+nod+1)-1
            jold = ja(k)
            jnew = iord(jold)
            if (jold .eq. nod .or. jnew .gt. nset) goto 5 
            iw(jnew) = iw(jnew) - 1
            nnz = nnz-1 
            call FixHeapM (iw,riord,iord,jnew,jnew,nset) 
 5       continue 
         ipos = iptrm1
 101  continue
c      
      if (nnz .gt. 0) goto 3
      return
c----------------------------------------------------------------------- 
      end
c-----------------------------------------------------------------------
      subroutine HeapSort (a,ind,rind,n,ncut)
      integer a(*),ind(n),rind(n),n, ncut 
c----------------------------------------------------------------------  
c integer version -- min heap sorts decreasinly. 
c---------------------------------------------------------------------- 
c sorts inger keys in array a increasingly and permutes the companion
c array ind rind accrodingly. 
c n    = size of array
c ncut = integer indicating when to cut the process.the process is
c        stopped after ncut outer steps of the heap-sort algorithm.
c        The first ncut values are sorted and they are the smallest
c        ncut values of the array.
c----------------------------------------------------------------------
c local variables 
c
      integer i,last, j,jlast
c   
c    Heap sort algorithm ---
c    
c    build heap 
       do 1 i=n/2,1,-1 
          j = i
          call FixHeap (a,ind,rind,j,j,n) 
 1    continue
c     
c   done -- now remove keys one by one 
c
      jlast = max(2,n-ncut+1)
      do 2 last=n,jlast,-1 
         call moveback (a,ind,rind,last) 
 2    continue
      return
      end 
c----------------------------------------------------------------------
      subroutine FixHeap (a,ind,rind,jkey,vacant,last)
      integer a(*),ind(*),rind(*),jkey,vacant,last
c----------------------------------------------------------------------
c     inserts a key (key and companion index) at the vacant position 
c     in a (min) heap -
c arguments
c     a(1:last)    = real array
c     ind(1:last)  = integer array -- permutation of initial data
c     rind(1:last) = integer array -- reverse permutation 
c     jkey         = position of key to be inserted. a(jkey) 
c                    will be inserted into the heap
c     vacant       = vacant where a key is to be inserted
c     last         = number of elements in the heap.
c----------------------------------------------------------------------
c local variables 
c
      integer child,lchild,rchild,xkey
      xkey = a(jkey) 
      ikey = ind(jkey) 
      lchild = 2*vacant
 1    continue 
      rchild = lchild+1  
      child = lchild  
      if (rchild .le. last .and. a(rchild) .lt. a(child))
     *     child = rchild 
      if (xkey .le. a(child) .or. child .gt. last) goto 2 
      a(vacant) = a(child) 
      ind(vacant) = ind(child) 
      rind(ind(vacant)) = vacant 
      vacant = child 
      lchild = 2*vacant 
      if (lchild .le.  last) goto 1
 2    continue 
      a(vacant) = xkey 
      ind(vacant) = ikey 
      rind(ikey) = vacant 
      return
c----------------------------------------------------------------------      
      end
c----------------------------------------------------------------------      
      subroutine HeapInsert (a,ind,rind,jkey,child,node) 
      integer a(*),ind(*),rind(*),jkey,child,node
c----------------------------------------------------------------------
c inserts a key to a heap from `node'. Checks values up
c only -- i.e.,assumes that the subtree (if any) whose root
c is node is such that the keys are all inferior to those
c to ge inserted. 
c 
c child is where the key ended up.
c---------------------------------------------------------------------- 
c---- local variables 
      integer parent,xkey,ikey
      xkey = a(jkey) 
      ikey = ind(jkey) 
c      node = node + 1 
      a(node) = xkey 
      ind(node) = ikey 
      rind(ikey) = node 
      if (node .le. 1) return
      child=node 
 1    parent = child/2  
      if (a(parent) .le. a(child)) goto 2
      call interchange(a,ind,rind,child,parent) 
      child = parent 
      if (child .gt. 1) goto 1 
 2    continue
      return
      end 
c-----------------------------------------------------------------------
      subroutine interchange (a,ind,rind,i,j)
      integer a(*),ind(*),rind(*),i,j 
      integer tmp,itmp
      tmp = a(i)  
      itmp = ind(i) 
c
      a(i) = a(j) 
      ind(i) = ind(j) 
c
      a(j) = tmp
      ind(j) = itmp 
      rind(ind(j)) = j
      rind(ind(i)) = i
c
      return
      end 
c----------------------------------------------------------------------
      subroutine moveback (a,ind,rind,last) 
      integer a(*),ind(*),rind(*),last 
c moves the front key to the back and inserts the last
c one back in from the top -- 
c 
c local variables 
c
      integer vacant,xmin 
c
         vacant = 1 
         xmin = a(vacant)
         imin = ind(vacant) 
         call FixHeap(a,ind,rind,last,vacant,last-1) 
         a(last) = xmin 
        ind(last) = imin 
         rind(ind(last)) = last
c
         return
         end 
c----------------------------------------------------------------------
      subroutine FixHeapM (a,ind,rind,jkey,vacant,last)
      integer a(*),ind(*),rind(*),jkey,vacant,last
c----      
c     inserts a key (key and companion index) at the vacant position 
c     in a heap -  THIS IS A MAX HEAP VERSION
c arguments
c     a(1:last)    = real array
c     ind(1:last)  = integer array -- permutation of initial data
c     rind(1:last) = integer array -- reverse permutation 
c     jkey         = position of key to be inserted. a(jkey) 
c                    will be inserted into the heap
c     vacant       = vacant where a key is to be inserted
c     last         = number of elements in the heap.
c----     
c local variables 
c
      integer child,lchild,rchild,xkey
      xkey = a(jkey) 
      ikey = ind(jkey) 
      lchild = 2*vacant
 1    continue 
      rchild = lchild+1  
      child = lchild  
      if (rchild .le. last .and. a(rchild) .gt. a(child))
     *     child = rchild 
      if (xkey .ge. a(child) .or. child .gt. last) goto 2
      a(vacant) = a(child) 
      ind(vacant) = ind(child) 
      rind(ind(vacant)) = vacant 
      vacant = child 
      lchild = 2*vacant 
      if (lchild .le.  last) goto 1
 2    continue 
      a(vacant) = xkey 
      ind(vacant) = ikey 
      rind(ikey) = vacant 
      return
      end
c      
      subroutine HeapInsertM (a,ind,rind,jkey,child,node) 
      integer a(*),ind(*),rind(*),jkey,child,node
c----------------------------------------------------------------------
c inserts a key to a heap from `node'. Checks values up
c only -- i.e.,assumes that the subtree (if any) whose root
c is node is such that the keys are all inferior to those
c to ge inserted. 
c 
c child is where the key ended up.
c---------------------------------------------------------------------- 
c---- local variables 
      integer parent,xkey,ikey
      xkey = a(jkey) 
      ikey = ind(jkey) 
c      node = node + 1 
      a(node) = xkey 
      ind(node) = ikey 
      rind(ikey) = node 
      if (node .le. 1) return
      child=node 
 1    parent = child/2  
      if (a(parent) .ge. a(child)) goto 2
      call interchange(a,ind,rind,child,parent) 
      child = parent 
      if (child .gt. 1) goto 1 
 2    continue
      return
      end 
c----------------------------------------------------------------------
      subroutine movebackM (a,ind,rind,last) 
      integer a(*),ind(*),rind(*),last 
c----------------------------------------------------------------------
c moves the front key to the back and inserts the last
c one back in from the top --  MAX HEAP VERSION 
c----------------------------------------------------------------------
c 
c local variables 
c
      integer vacant,xmin 
c
      vacant = 1 
      xmin = a(vacant)
      imin = ind(vacant) 
      call FixHeapM(a,ind,rind,last,vacant,last-1) 
      a(last) = xmin 
      ind(last) = imin 
      rind(ind(last)) = last
c----------------------------------------------------------------------
      return
      end 
c----------------------------------------------------------------------
      subroutine indsetr (n,ja,ia,nset,iord,riord,sym,iptr) 
      integer n, nset, ja(*),ia(*),riord(*),iord(*) 
      logical sym
c---------------------------------------------------------------------- 
c greedy algorithm for independent set ordering -- RANDOM TRAVERSAL -- 
c----------------------------------------------------------------------
c parameters:
c ----------
c n      = row dimension of matrix
c ja, ia = matrix pattern in CRS format
c nset   = (output) number of elements in the independent set
c iord   = permutation array corresponding to the independent set 
c          ordering. Row number i will become row number iord(i) in 
c          permuted matrix.
c riord  = reverse permutation array. Row number i in the permutated 
c          matrix is row number riord(i) in original matrix. 
c----------------------------------------------------------------------
c notes: works for CSR, MSR, and CSC formats but assumes that the
c matrix has a symmetric structure. 
c---------------------------------------------------------------------- 
c local variables
c 
      integer j, k1, k2, nod, k, mat
      do 1 j=1, n
         iord(j) = 0
 1    continue
c
c generate random permutation
c 
      iseed = 0 
      call rndperm(n, riord, iseed)
      write (8,'(10i6)') (riord(j),j=1,n) 
c
      nummat = 1
      if (.not. sym) nummat = 2
c     
c iord used as a marker
c 
      nset = 0
      do 12  ii=1, n
         nod = riord(ii) 
         if (iord(nod) .ne. 0) goto 12 
         nset = nset+1
         iord(nod) = 1
c
c visit all neighbors of current nod 
c     
         ipos = 0
         do 45 mat=1, nummat
            do 4 k=ia(ipos+nod), ia(ipos+nod+1)-1 
               j = ja(k)
               if (j .ne. nod) iord(j) = 2
 4          continue
            ipos = iptr-1
 45      continue
 12   continue
c
c get permutation
c     
      k1 = 0 
      k2 = nset 
      do 6 j=1,n
         if (iord(j) .eq. 1) then
            k1 = k1+1
            k = k1
         else
            k2 = k2+1
            k = k2 
         endif
            riord(k) = j         
            iord(j) = k
 6    continue
      return
c----------------------------------------------------------------------
      end
c----------------------------------------------------------------------
      subroutine rndperm(n,iord,iseed) 
      integer n, iseed, iord(n) 
c----------------------------------------------------------------------
c this subroutine will generate a pseudo random permutation of the
c n integers 1,2, ...,n.
c iseed is the initial seed. any integer.
c----------------------------------------------------------------------
c local
c
      integer i, j, itmp 
c----------------------------------------------------------------------
      do j=1, n
         iord(j) = j
      enddo
c
      do i=1, n
         j = mod(irand(0),n) + 1
         itmp = iord(i) 
         iord(i) = iord(j) 
         iord(j) = itmp
      enddo
c----------------------------------------------------------------------
      return
c----------------------------------------------------------------------
      end 
