c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c               REORDERING ROUTINES -- LEVEL SET BASED ROUTINES        c
c----------------------------------------------------------------------c
c dblstr   : doubled stripe partitioner 
c rdis     : recursive dissection partitioner
c dse2way  : distributed site expansion usuing sites from dblstr 
c dse      : distributed site expansion usuing sites from rdis
c------------- utility routines ----------------------------------------- 
c BFS      : Breadth-First search traversal algorithm 
c add_lvst : routine to add a level -- used by BFS 
c stripes  : finds the level set structure
c stripes0 : finds a trivial one-way partitioning from level-sets 
c perphn   : finds a pseudo-peripheral node and performs a BFS from it.
c mapper4  : routine used by dse and dse2way to do center expansion
c get_domns: routine to find subdomaine from linked lists found by 
c            mapper4. 
c add_lk   : routine to add entry to linked list -- used by mapper4. 
c find_ctr : routine to locate an approximate center of a subgraph. 
c rversp   : routine to reverse a given permutation (e.g., for RCMK)
c maskdeg  : integer function to compute the `masked' of a node
c-----------------------------------------------------------------------
      subroutine dblstr(n,ja,ia,ip1,ip2,nfirst,riord,ndom,map,mapptr,
     *     mask,levels,iwk) 
      implicit none
      integer ndom,ja(*),ia(*),ip1,ip2,nfirst,riord(*),map(*),mapptr(*),
     *     mask(*),levels(*),iwk(*),nextdom
c-----------------------------------------------------------------------
c     this routine does a two-way partitioning of a graph using 
c     level sets recursively. First a coarse set is found by a
c     simple cuthill-mc Kee type algorithm. Them each of the large
c     domains is further partitioned into subsets using the same 
c     technique. The ip1 and ip2 parameters indicate the desired number 
c     number of partitions 'in each direction'. So the total number of
c     partitions on return ought to be equal (or close) to ip1*ip2 
c----------------------parameters----------------------------------------
c on entry: 
c---------
c n      = row dimension of matrix == number of vertices in graph
c ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data
c          structure)
c ip1    = integer indicating the number of large partitions ('number of
c          paritions in first direction') 
c ip2    = integer indicating the number of smaller partitions, per 
c          large partition, ('number of partitions in second direction') 
c nfirst = number of nodes in the first level that is input in riord 
c riord  = (also an ouput argument). on entry riord contains the labels  
c          of the nfirst nodes that constitute the first level.   
c on return:
c-----------
c ndom   = total number of partitions found 
c map    = list of nodes listed partition by partition from partition 1
c          to paritition ndom.
c mapptr = pointer array for map. All nodes from position 
c          k1=mapptr(idom),to position k2=mapptr(idom+1)-1 in map belong
c          to partition idom.
c work arrays:
c-------------
c mask   = array of length n, used to hold the partition number of each 
c          node for the first (large) partitioning. 
c          mask is also used as a marker of  visited nodes. 
c levels = integer array of length .le. n used to hold the pointer 
c          arrays for the various level structures obtained from BFS. 
c 
c-----------------------------------------------------------------------
      integer n, j,idom,kdom,jdom,maskval,k,nlev,init,ndp1,numnod
      maskval = 1 
      do j=1, n
         mask(j) = maskval 
      enddo
      iwk(1) = 0 
      call BFS(n,ja,ia,nfirst,iwk,mask,maskval,riord,levels,nlev)      
c
c     init = riord(1) 
c     call perphn (ja,ia,mask,maskval,init,nlev,riord,levels) 
      call stripes (nlev,riord,levels,ip1,map,mapptr,ndom)
c-----------------------------------------------------------------------
      if (ip2 .eq. 1) return      
      ndp1 = ndom+1
c     
c     pack info into array iwk 
c 
      do j = 1, ndom+1
         iwk(j) = ndp1+mapptr(j)  
      enddo
      do j=1, mapptr(ndom+1)-1
         iwk(ndp1+j) = map(j) 
      enddo
      do idom=1, ndom 
         j = iwk(idom) 
         numnod = iwk(idom+1) - iwk(idom) 
         init = iwk(j) 
         do k=j, iwk(idom+1)-1 
         enddo
      enddo

      do idom=1, ndom 
         do k=mapptr(idom),mapptr(idom+1)-1 
            mask(map(k)) = idom
         enddo
      enddo
      nextdom = 1 
c
c     jdom = counter for total number of (small) subdomains 
c 
      jdom = 1
      mapptr(jdom) = 1 
c----------------------------------------------------------------------- 
      do idom =1, ndom
         maskval = idom
         nfirst = 1
         numnod = iwk(idom+1) - iwk(idom) 
         j = iwk(idom) 
         init = iwk(j) 
         nextdom = mapptr(jdom) 
c  note:    old version uses iperm array 
         call perphn(numnod,ja,ia,init,mask,maskval,
     *        nlev,riord,levels)
c          
         call stripes (nlev,riord,levels,ip2,map(nextdom),
     *        mapptr(jdom),kdom)
c          
         mapptr(jdom) = nextdom
         do j = jdom,jdom+kdom-1
            mapptr(j+1) = nextdom + mapptr(j+1)-1
         enddo
         jdom = jdom + kdom
      enddo
c
      ndom = jdom - 1
      return
      end 
c-----------------------------------------------------------------------
      subroutine rdis(n,ja,ia,ndom,map,mapptr,mask,levels,size,iptr) 
      implicit none
      integer n,ja(*),ia(*),ndom,map(*),mapptr(*),mask(*),levels(*),
     *     size(ndom),iptr(ndom)
c-----------------------------------------------------------------------
c     recursive dissection algorithm for partitioning.
c     initial graph is cut in two - then each time, the largest set
c     is cut in two until we reach desired number of domains.
c----------------------------------------------------------------------- 
c     input
c     n, ja, ia = graph 
c     ndom      = desired number of subgraphs
c     output 
c     ------
c     map, mapptr  = pointer array data structure for domains. 
c             if k1 = mapptr(i), k2=mapptr(i+1)-1 then 
c             map(k1:k2) = points in domain number i
c    work arrays:
c    -------------
c    mask(1:n)    integer 
c    levels(1:n)  integer 
c    size(1:ndom) integer 
c    iptr(1:ndom) integer 
c----------------------------------------------------------------------- 
      integer idom,maskval,k,nlev,init,nextsiz,wantsiz,lev,ko,
     *     maxsiz,j,nextdom  
c-----------------------------------------------------------------------
      idom = 1
c-----------------------------------------------------------------------
c     size(i) = size of domnain  i
c     iptr(i)  = index of first element of domain i
c-----------------------------------------------------------------------
      size(idom) = n
      iptr(idom) = 1
      do j=1, n
         mask(j) = 1 
      enddo
c     
c     domain loop
c
 1    continue
c
c     select domain with largest size
c     
      maxsiz = 0 
      do j=1, idom
         if (size(j) .gt. maxsiz) then
            maxsiz = size(j)
            nextdom = j
         endif
      enddo
c
c     do a Prphn/ BFS on nextdom
c     
      maskval = nextdom
      init = iptr(nextdom) 
      call perphn(n,ja,ia,init,mask,maskval,nlev,map,levels) 
c
c     determine next subdomain
c
      nextsiz = 0
      wantsiz = maxsiz/2 
      idom = idom+1
      lev = nlev 
      do while (nextsiz .lt. wantsiz) 
         do k = levels(lev), levels(lev+1)-1
            mask(map(k)) = idom
        enddo
        nextsiz = nextsiz + levels(lev+1) - levels(lev) 
        lev = lev-1
      enddo
c
      size(nextdom) = size(nextdom) - nextsiz
      size(idom) = nextsiz
c
c     new initial point = last point of previous domain
c
       iptr(idom) = map(levels(nlev+1)-1) 
c       iptr(idom) = map(levels(lev)+1) 
c      iptr(idom) = 1 
c
c alternative 
c      lev = 1 
c      do while (nextsiz .lt. wantsiz) 
c         do k = levels(lev), levels(lev+1)-1
c            mask(map(k)) = idom
c         enddo
c         nextsiz = nextsiz + levels(lev+1) - levels(lev) 
c         lev = lev+1
c      enddo
c
c     set size of new domain and adjust previous one
c
c      size(idom) = nextsiz 
c      size(nextdom) = size(nextdom) - nextsiz 
c      iptr(idom) = iptr(nextdom) 
c      iptr(nextdom) = map(levels(lev))

      if (idom .lt. ndom) goto 1
c
c     domains found -- build data structure 
c     
      mapptr(1) = 1
      do idom=1, ndom 
         mapptr(idom+1) = mapptr(idom) + size(idom) 
      enddo
      do k=1, n
         idom = mask(k) 
         ko = mapptr(idom) 
         map(ko) = k
         mapptr(idom) = ko+1
      enddo
c
c     reset pointers
c     
      do j = ndom,1,-1
         mapptr(j+1) = mapptr(j) 
      enddo
      mapptr(1) = 1 
c
      return
      end 
c-----------------------------------------------------------------------
      subroutine dse2way(n,ja,ia,ip1,ip2,nfirst,riord,ndom,dom,idom,
     *     mask,jwk,link) 
c-----------------------------------------------------------------------
c     uses centers obtained from dblstr partition to get new partition
c----------------------------------------------------------------------- 
c     input: n, ja, ia   = matrix
c     nfirst = number of first points 
c     riord  = riord(1:nfirst) initial points 
c     output 
c     ndom   = number of domains
c     dom, idom = pointer array structure for domains. 
c     mask , jwk, link = work arrays,
c-----------------------------------------------------------------------
      implicit none 
      integer n, ja(*), ia(*), ip1, ip2, nfirst, riord(*), dom(*),
     *     idom(*), mask(*), jwk(*),ndom,link(*)  
c
c-----------------------------------------------------------------------
c     local variables
      integer i, mid,nsiz, maskval,init, outer, nouter, k
      call dblstr(n,ja,ia,ip1,ip2,nfirst,riord,ndom,dom,idom,mask,
     *     link,jwk)
c
      nouter = 3
c----------------------------------------------------------------------- 

      do outer =1, nouter 
c
c     set masks 
c
      do i=1, ndom
         do k=idom(i),idom(i+1)-1
            mask(dom(k)) = i
         enddo
      enddo
c
c     get centers 
c 
      do i =1, ndom
         nsiz = idom(i+1) - idom(i) 
         init = dom(idom(i))
         maskval = i 
c
c         use link for local riord -- jwk for other arrays -- 
c 
         call find_ctr(n,nsiz,ja,ia,init,mask,maskval,link, 
     *    jwk,mid,jwk(nsiz+1)) 
         riord(i) = mid 
      enddo
c
c     do level-set expansion from centers -- save previous diameter 
c 
      call mapper4(n,ja,ia,ndom,riord,jwk,mask,link) 
      call get_domns2(ndom,riord,link,jwk,dom,idom)
c----------------------------------------------------------------------- 
      enddo 
      return 
      end 
c----------------------------------------------------------------------- 
      subroutine dse(n,ja,ia,ndom,riord,dom,idom,mask,jwk,link) 
      implicit none 
      integer n, ja(*), ia(*), ndom, riord(*), dom(*),
     *     idom(*), mask(*), jwk(*),link(*)  
c-----------------------------------------------------------------------
c     uses centers produced from rdis to get a new partitioning -- 
c     see calling sequence in rdis.. 
c-----------------------------------------------------------------------
c     local variables
      integer i, mid, nsiz, maskval,init, outer, nouter, k 
c-----------------------------------------------------------------------
      nouter = 3
c 
      call rdis(n,ja,ia,ndom,dom,idom,mask,link,jwk,jwk(ndom+1)) 
c
c     initial points = 
c
      do outer =1, nouter 
c
c     set masks 
c
      do i=1, ndom
         do k=idom(i),idom(i+1)-1
            mask(dom(k)) = i
         enddo
      enddo
c
c     get centers 
c 
      do i =1, ndom
         nsiz = idom(i+1) - idom(i) 
         init = dom(idom(i))
         maskval = i 
c
c         use link for local riord -- jwk for other arrays -- 
c 

         call find_ctr(n,nsiz,ja,ia,init,mask,maskval,link, 
     *    jwk,mid,jwk(nsiz+1)) 
         riord(i) = mid 
      enddo
c
c     do level-set expansion from centers -- save previous diameter 
c 
      call mapper4(n,ja,ia,ndom,riord,jwk,mask,link) 
      call get_domns2(ndom,riord,link,jwk,dom,idom)
c----------------------------------------------------------------------- 
      enddo 
      return 
      end 
c----------------------------------------------------------------------- 
      subroutine BFS(n,ja,ia,nfirst,iperm,mask,maskval,riord,levels,
     *     nlev)
      implicit none 
      integer n,ja(*),ia(*),nfirst,iperm(n),mask(n),riord(*),levels(*),
     *     nlev,maskval 
c-----------------------------------------------------------------------
c finds the level-structure (breadth-first-search or CMK) ordering for a
c given sparse matrix. Uses add_lvst. Allows an set of nodes to be 
c the initial level (instead of just one node). 
c-------------------------parameters------------------------------------
c on entry:
c---------
c     n      = number of nodes in the graph 
c     ja, ia = pattern of matrix in CSR format (the ja,ia arrays of csr data
c     structure)
c     nfirst = number of nodes in the first level that is input in riord
c     iperm  = integer array indicating in which order to  traverse the graph
c     in order to generate all connected components. 
c     if iperm(1) .eq. 0 on entry then BFS will traverse the nodes
c     in the  order 1,2,...,n.
c     
c     riord  = (also an ouput argument). On entry riord contains the labels  
c     of the nfirst nodes that constitute the first level.      
c     
c     mask   = array used to indicate whether or not a node should be 
c     condidered in the graph. see maskval.
c     mask is also used as a marker of  visited nodes. 
c     
c     maskval= consider node i only when:  mask(i) .eq. maskval 
c     maskval must be .gt. 0. 
c     thus, to consider all nodes, take mask(1:n) = 1. 
c     maskval=1 (for example) 
c     
c     on return
c     ---------
c     mask   = on return mask is restored to its initial state. 
c     riord  = `reverse permutation array'. Contains the labels of the nodes
c     constituting all the levels found, from the first level to
c     the last. 
c     levels = pointer array for the level structure. If lev is a level
c     number, and k1=levels(lev),k2=levels(lev+1)-1, then
c     all the nodes of level number lev are:
c     riord(k1),riord(k1+1),...,riord(k2) 
c     nlev   = number of levels found
c-----------------------------------------------------------------------
c     
      integer j, ii, nod, istart, iend 
      logical permut
      permut = (iperm(1) .ne. 0) 
c     
c     start pointer structure to levels 
c     
      nlev   = 0 
c     
c     previous end
c     
      istart = 0 
      ii = 0
c     
c     current end 
c     
      iend = nfirst
c     
c     intialize masks to zero -- except nodes of first level -- 
c     
      do 12 j=1, nfirst 
         mask(riord(j)) = 0 
 12   continue
c-----------------------------------------------------------------------
 13   continue 
c     
 1    nlev = nlev+1
      levels(nlev) = istart + 1
      call add_lvst (istart,iend,nlev,riord,ja,ia,mask,maskval) 
      if (istart .lt. iend) goto 1
 2    ii = ii+1 
      if (ii .le. n) then
         nod = ii         
         if (permut) nod = iperm(nod)          
         if (mask(nod) .eq. maskval) then
c     
c     start a new level
c
            istart = iend
            iend = iend+1 
            riord(iend) = nod
            mask(nod) = 0
            goto 1
         else 
            goto 2
         endif
      endif
c----------------------------------------------------------------------- 
 3    levels(nlev+1) = iend+1 
      do j=1, iend
         mask(riord(j)) = maskval 
      enddo
c----------------------------------------------------------------------- 
      return
      end
c-----------------------------------------------------------------------
      subroutine add_lvst(istart,iend,nlev,riord,ja,ia,mask,maskval) 
      integer nlev, nod, riord(*), ja(*), ia(*), mask(*) 
c-------------------------------------------------------------
c     adds one level set to the previous sets.. 
c     span all nodes of previous mask
c-------------------------------------------------------------
      nod = iend
      do 25 ir = istart+1,iend 
         i = riord(ir)		
         do 24 k=ia(i),ia(i+1)-1
            j = ja(k)
            if (mask(j) .eq. maskval) then
               nod = nod+1 
               mask(j) = 0
               riord(nod) = j
            endif 
 24      continue
 25   continue
      istart = iend 
      iend   = nod 
      return
      end 
c----------------------------------------------------------------------- 
      subroutine stripes (nlev,riord,levels,ip,map,mapptr,ndom)
      implicit none
      integer nlev,riord(*),levels(nlev+1),ip,map(*),
     *    mapptr(*), ndom
c-----------------------------------------------------------------------
c    this is a post processor to BFS. stripes uses the output of BFS to 
c    find a decomposition of the adjacency graph by stripes. It fills 
c    the stripes level by level until a number of nodes .gt. ip is 
c    is reached. 
c---------------------------parameters-----------------------------------
c on entry: 
c --------
c nlev   = number of levels as found by BFS 
c riord  = reverse permutation array produced by BFS -- 
c levels = pointer array for the level structure as computed by BFS. If 
c          lev is a level number, and k1=levels(lev),k2=levels(lev+1)-1, 
c          then all the nodes of level number lev are:
c                      riord(k1),riord(k1+1),...,riord(k2) 
c  ip    = number of desired partitions (subdomains) of about equal size.
c 
c on return
c ---------
c ndom     = number of subgraphs (subdomains) found 
c map      = node per processor list. The nodes are listed contiguously
c            from proc 1 to nproc = mpx*mpy. 
c mapptr   = pointer array for array map. list for proc. i starts at 
c            mapptr(i) and ends at mapptr(i+1)-1 in array map.
c-----------------------------------------------------------------------
c local variables. 
c
      integer ib,ktr,ilev,k,nsiz,psiz 
      ndom = 1 
      ib = 1
c to add: if (ip .le. 1) then ...
      nsiz = levels(nlev+1) - levels(1) 
      psiz = (nsiz-ib)/max(1,(ip - ndom + 1)) + 1 
      mapptr(ndom) = ib 
      ktr = 0 
      do 10 ilev = 1, nlev
c
c     add all nodes of this level to domain
c     
         do 3 k=levels(ilev), levels(ilev+1)-1
            map(ib) = riord(k)
            ib = ib+1
            ktr = ktr + 1 
            if (ktr .ge. psiz  .or. k .ge. nsiz) then 
               ndom = ndom + 1
               mapptr(ndom) = ib 
               psiz = (nsiz-ib)/max(1,(ip - ndom + 1)) + 1 
               ktr = 0
            endif
c
 3       continue
 10   continue
      ndom = ndom-1
      return 
      end
c-----------------------------------------------------------------------
      subroutine stripes0 (ip,nlev,il,ndom,iptr)
      integer ip, nlev, il(*), ndom, iptr(*)
c-----------------------------------------------------------------------
c     This routine is a simple level-set partitioner. It scans
c     the level-sets as produced by BFS from one to nlev.
c     each time the number of nodes in the accumulated set of
c     levels traversed exceeds the parameter ip, this set defines 
c     a new subgraph. 
c-------------------------parameter-list---------------------------------
c on entry:
c --------
c ip     = desired number of nodes per subgraph.
c nlev   = number of levels found  as output by BFS
c il     = integer array containing the pointer array for
c          the level data structure as output by BFS. 
c          thus il(lev+1) - il(lev) = the number of 
c          nodes that constitute the level numbe lev.
c on return
c ---------
c ndom   = number of sungraphs found
c iptr   = pointer array for the sugraph data structure. 
c          thus, iptr(idom) points to the first level that 
c          consistutes the subgraph number idom, in the 
c          level data structure. 
c-----------------------------------------------------------------------
      ktr = 0
      iband = 1 
      iptr(iband) = 1 
c-----------------------------------------------------------------------

      do 10 ilev = 1, nlev
         ktr = ktr + il(ilev+1) - il(ilev)
         if (ktr .gt. ip) then
            iband = iband+1 
            iptr(iband) = ilev+1
            ktr = 0
         endif
c
 10   continue
c-----------returning --------------------
      iptr(iband) = nlev + 1 
      ndom = iband-1
      return
c-----------------------------------------------------------------------
c-----end-of-stripes0--------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      integer function maskdeg (ja,ia,nod,mask,maskval) 
      implicit none 
      integer ja(*),ia(*),nod,mask(*),maskval
c-----------------------------------------------------------------------
      integer deg, k 
      deg = 0 
      do k =ia(nod),ia(nod+1)-1
         if (mask(ja(k)) .eq. maskval) deg = deg+1 
      enddo
      maskdeg = deg 
      return
      end 
c-----------------------------------------------------------------------
      subroutine perphn(n,ja,ia,init,mask,maskval,nlev,riord,levels) 
      implicit none
      integer n,ja(*),ia(*),init,mask(*),maskval,
     *     nlev,riord(*),levels(*)
c-----------------------------------------------------------------------
c     finds a peripheral node and does a BFS search from it. 
c-----------------------------------------------------------------------
c     see routine  dblstr for description of parameters
c input:
c-------
c ja, ia  = list pointer array for the adjacency graph
c mask    = array used for masking nodes -- see maskval
c maskval = value to be checked against for determing whether or
c           not a node is masked. If mask(k) .ne. maskval then
c           node k is not considered.
c init    = init node in the pseudo-peripheral node algorithm.
c
c output:
c-------
c init    = actual pseudo-peripherial node found.
c nlev    = number of levels in the final BFS traversal.
c riord   =
c levels  =
c----------------------------------------------------------------------- 
      integer j,nlevp,deg,nfirst,mindeg,nod,maskdeg
      integer iperm(1) 
      nlevp = 0 
 1    continue
      riord(1) = init
      nfirst = 1 
      iperm(1) = 0
c
      call BFS(n,ja,ia,nfirst,iperm,mask,maskval,riord,levels,nlev)
      if (nlev .gt. nlevp) then 
         mindeg = n+1 
         do j=levels(nlev),levels(nlev+1)-1
            nod = riord(j) 
            deg = maskdeg(ja,ia,nod,mask,maskval)
            if (deg .lt. mindeg) then
               init = nod
               mindeg = deg
            endif 
         enddo
         nlevp = nlev 
         goto 1 
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine mapper4 (n,ja,ia,ndom,nodes,levst,marker,link)
      implicit none 
      integer n,ndom,ja(*),ia(*),marker(n),levst(2*ndom), 
     *     nodes(*),link(*) 
c-----------------------------------------------------------------------
c     finds domains given ndom centers -- by doing a level set expansion 
c-----------------------------------------------------------------------
c     on entry:
c     ---------
c     n      = dimension of matrix 
c     ja, ia = adajacency list of matrix (CSR format without values) -- 
c     ndom   = number of subdomains (nr output by coarsen)  
c     nodes  = array of size at least n. On input the first ndom entries
c              of nodes should contain the labels of the centers of the 
c              ndom domains from which to do the expansions. 
c     
c     on return 
c     --------- 
c     link  = linked list array for the ndom domains. 
c     nodes = contains the list of nodes of the domain corresponding to
c             link. (nodes(i) and link(i) are related to the same node). 
c    
c     levst = levst(j) points to beginning of subdomain j in link.
c
c     work arrays:
c     ----------- 
c     levst : work array of length 2*ndom -- contains the beginning and
c     end of  current level in link. 
c     beginning of last level in link for each processor.
c     also ends in levst(ndom+i) 
c     marker : work array of length n. 
c
c     Notes on implementation: 
c     ----------------------- 
c     for j .le. ndom link(j) is <0  and indicates the end of the
c     linked list. The most recent element added to the linked
c     list is added at the end of the list (traversal=backward) 
c     For  j .le. ndom, the value of -link(j) is the size of 
c     subdomain j. 
c
c-----------------------------------------------------------------------
c     local variables 
      integer mindom,j,lkend,nod,nodprev,idom,next,i,kk,ii,ilast,nstuck,
     *     isiz, nsize 
c     
c     initilaize nodes and link arrays
c
      do 10 j=1, n
         marker(j) = 0
 10   continue
c     
      do 11 j=1, ndom
         link(j) = -1 
         marker(nodes(j)) = j
         levst(j) = j
         levst(ndom+j) = j 
 11   continue
c
c     ii = next untouched node for restarting new connected component. 
c 
      ii = 0
c     
      lkend = ndom
      nod   = ndom  
      nstuck = 0 
c-----------------------------------------------------------------------
 100  continue 
      idom = mindom(n,ndom,link) 
c-----------------------------------------------------------------------
c     begin level-set loop 
c-----------------------------------------------------------------------
 3    nodprev = nod 
      ilast = levst(ndom+idom) 
      levst(ndom+idom) = lkend       
      next = levst(idom) 
c     
c     linked list traversal loop
c 
      isiz = 0 
      nsize = link(idom) 
 1    i = nodes(next) 
      isiz = isiz + 1 
c     
c     adjacency list traversal loop 
c     
      do 2 kk=ia(i), ia(i+1)-1
         j = ja(kk) 
         if (marker(j) .eq. 0) then 
            call add_lk(j,nod,idom,ndom,lkend,levst,link,nodes,marker) 
         endif
 2    continue
c     
c     if last element of the previous level not reached continue
c     
      if (next .gt. ilast) then
         next = link(next) 
         if (next .gt. 0) goto 1
      endif
c-----------------------------------------------------------------------
c     end level-set traversal --  
c-----------------------------------------------------------------------
      if (nodprev .eq. nod) then
c     
c     link(idom) >0 indicates that set is stuck  --  
c
         link(idom) = -link(idom) 
         nstuck = nstuck+1
      endif
c     
      if (nstuck .lt. ndom) goto 100 
c
c     reset sizes -- 
c 
      do j=1, ndom
         if (link(j) .gt. 0) link(j) = -link(j)
      enddo
c
      if (nod .eq. n) return 
c
c     stuck. add first non assigned point to smallest domain
c     
 20   ii = ii+1
      if (ii .le. n) then
         if (marker(ii) .eq. 0) then 
            idom = 0 
            isiz = n+1 
            do 30 kk=ia(ii), ia(ii+1)-1
               i = marker(ja(kk)) 
               if (i .ne. 0) then 
                  nsize = abs(link(i)) 
                  if (nsize .lt. isiz) then 
                     isiz = nsize
                     idom = i 
                  endif
               endif
 30            continue
c
c     if no neighboring domain select smallest one 
c
               if (idom .eq. 0) idom = mindom(n,ndom,link) 
c     
c     add ii to sudomain idom at end of linked list  
c     
            call add_lk(ii,nod,idom,ndom,lkend,levst,link,nodes,marker) 
            goto 3 
         else
            goto 20 
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine get_domns2(ndom,nodes,link,levst,riord,iptr)
      implicit none 
      integer ndom,nodes(*),link(*),levst(*),riord(*),iptr(*) 
c-----------------------------------------------------------------------
c     constructs the subdomains from its linked list data structure
c-----------------------------------------------------------------------
c     input:
c     ndom  = number of subdomains
c     nodes = sequence of nodes are produced by mapper4. 
c     link  = link list array as produced by mapper4.
c     on return:
c----------
c     riord = contains the nodes in each subdomain in succession.
c     iptr  = pointer in riord for beginnning of each subdomain.
c     Thus subdomain number i consists of nodes 
c     riord(k1),riord(k1)+1,...,riord(k2) 
c     where k1 = iptr(i), k2= iptr(i+1)-1
c     
c-----------------------------------------------------------------------
c     local variables 
      integer nod, j, next, ii 
      nod = 1
      iptr(1) = nod 
      do 21 j=1, ndom 
         next = levst(j)
 22      ii = nodes(next)
         riord(nod) = ii 
         nod = nod+1 
         next = link(next) 
         if (next .gt.  0) goto 22
         iptr(j+1) = nod 
 21   continue
c
      return
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      function mindom(n, ndom, link) 
      implicit none 
      integer mindom, n, ndom, link(n) 
c-----------------------------------------------------------------------
c     returns  the domain with smallest size
c----------------------------------------------------------------------- 
c      locals
c
      integer i, nsize, isiz 
c
      isiz = n+1 
      do 10 i=1, ndom
         nsize = - link(i) 
         if (nsize .lt. 0) goto 10 
         if (nsize .lt. isiz) then 
            isiz = nsize
            mindom = i
         endif
 10   continue
      return
      end 
c-----------------------------------------------------------------------
      subroutine add_lk(new,nod,idom,ndom,lkend,levst,link,nodes,marker) 
      implicit none
      integer new,nod,idom,ndom,lkend,levst(*),link(*),nodes(*),
     *     marker(*) 
c----------------------------------------------------------------------- 
c     inserts new element to linked list from the tail.
c----------------------------------------------------------------------- 
c     adds one entry (new) to linked list and ipdates everything.
c     new  = node to be added
c     nod  = current number of marked nodes
c     idom = domain to which new is to be added
c     ndom = total number of domains
c     lkend= location of end of structure (link and nodes)
c     levst= pointer array for link, nodes
c     link = link array 
c     nodes= nodes array -- 
c     marker = marker array == if marker(k) =0 then node k is not
c              assigned yet. 
c----------------------------------------------------------------------- 
c      locals
c     
      integer ktop  
      lkend = lkend + 1
      nodes(lkend) = new
      nod = nod+1 
      marker(new) = idom 
      ktop = levst(idom) 
      link(lkend) = ktop 
      link(idom) = link(idom)-1 
      levst(idom) = lkend 
      return
c-----------------------------------------------------------------------
c-------end-of-add_lk--------------------------------------------------- 
      end 
c----------------------------------------------------------------------- 
      subroutine find_ctr(n,nsiz,ja,ia,init,mask,maskval,riord,
     *     levels,center,iwk) 
      implicit none
      integer n,nsiz,ja(*),ia(*),init,mask(*),maskval,riord(*),
     *     levels(*),center,iwk(*) 
c-----------------------------------------------------------------------
c     finds a center point of a subgraph -- 
c-----------------------------------------------------------------------
c     n, ja, ia = graph
c     nsiz = size of current domain.
c     init = initial node in search
c     mask
c     maskval 
c-----------------------------------------------------------------------
c     local variables 
      integer midlev, nlev,newmask, k, kr, kl, init0, nlev0  
      call perphn(n,ja,ia,init,mask,maskval,nlev,riord,levels)
c-----------------------------------------------------------------------
c     midlevel = level which cuts domain into 2 roughly equal-size 
c     regions 
c
      midlev = 1
      k = 0 
 1    continue
      k = k + levels(midlev+1)-levels(midlev) 
      if (k*2 .lt. nsiz) then
         midlev = midlev+1
         goto 1 
      endif
c-----------------------------------------------------------------------
      newmask = n+maskval
c     
c     assign temporary masks to mid-level elements
c     
      do k=levels(midlev),levels(midlev+1)-1
         mask(riord(k)) = newmask
      enddo
c     
c     find pseudo-periph node for mid-level `line'
c
      kr = 1
      kl = kr + nsiz 
      init0 = riord(levels(midlev))
      call perphn(n,ja,ia,init0,mask,newmask,nlev0,iwk(kr),iwk(kl)) 
c-----------------------------------------------------------------------
c     restore  mask to initial state 
c-----------------------------------------------------------------------
      do k=levels(midlev),levels(midlev+1)-1
         mask(riord(k)) = maskval 
      enddo
c-----------------------------------------------------------------------
c     define center 
c-----------------------------------------------------------------------  
      midlev = 1 + (nlev0-1)/2
      k = iwk(kl+midlev-1)
      center = iwk(k) 
c----------------------------------------------------------------------- 
      return 
      end 
c-----------------------------------------------------------------------
      subroutine rversp (n, riord)
      integer n, riord(n)
c-----------------------------------------------------------------------
c     this routine does an in-place reversing of the permutation array
c     riord --
c-----------------------------------------------------------------------
      integer j, k
      do 26 j=1,n/2
         k = riord(j)
         riord(j) = riord(n-j+1)
         riord(n-j+1) = k
 26   continue
      return
      end
c-----------------------------------------------------------------------
