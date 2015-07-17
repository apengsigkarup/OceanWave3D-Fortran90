c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c      MATRIX GENERATION ROUTINES - FINITE ELEMENT MATRICES            c
c----------------------------------------------------------------------c
c contents:                                                            c
c----------                                                            c
c genfea       : generates finite element matrices in assembled form   c
c genfea_wbc   : generates finite element matrices in assembled form   c
c                without applying the boundary conditions              c
c genfeu       : generates finite element matrices in unassembled form c
c genfeu_wbc   : generates finite element matrices in unassembled form c
c                without applying the boundary conditions              c
c genfeu_lstif : generates finite element matrices in unassembled form c
c                using the lstif problem appearing in elmtlib2.f       c
c assmb1       : assembles an unassembled matrix (produced by genfeu)  c
c----------------------------------------------------------------------c
      subroutine genfea (nx,nelx,node,job,x,y,ijk,nodcode,fs,nint,
     *     a,ja,ia,f,iwk,jwk,ierr,xyk)
c-----------------------------------------------------------------------
c this subroutine generates a finite element matrix in assembled form.
c the matrix is assembled in compressed sparse row format. See genfeu
c for matrices in unassembled form. The user must provide the grid, 
c (coordinates x, y and connectivity matrix ijk) as well as some 
c information on the nodes (nodcode) and the material properties 
c (the function K(x,y) above) in the form of a subroutine xyk. 
c----------------------------------------------------------------------
c
c on entry:
c ---------
c 
c nx	    = integer . the number of nodes in the grid . 
c nelx	    = integer . the number of elements in the grid.
c node      = integer = the number of nodes per element (should be
c             set to three in this version). also the first dimension
c             of ijk
c job	    = integer. If job=0, it is assumed that there is no heat
c             source (i.e. fs = 0) and the right hand side
c             produced will therefore be a zero vector.
c             If job = 1 on entry then the contributions from the
c             heat source in each element are taken into account.
c 
c x, y      = two real arrays containing the coordinates of the nodes.
c 
c ijk       =  an integer array containing the connectivity matrix.
c              ijk(i,nel), i=1,2,..node, is the list of the nodes
c              constituting the element nel, ans listed in 
c              counter clockwise order.
c
c nodcode   = an integer array containing the boundary information for
c             each node with the following meaning.
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner node. [This node and the
c             corresponmding element are discarded.]
c
c fs	    = real array of length nelx on entry containing the heat 
c             source for each element (job = 1 only) 
c
c xyk	    = subroutine defining the material properties at each 
c	      element. Form: 
c 	      call xyk(nel,xyke,x,y,ijk,node) with on return
c             xyke =  material constant matrices. 
c	      for each element nel, xyke(1,nel),xyke(2,nel) 
c             and xyke(3,nel) represent the constants
c             K11, K22, and K12 at that element.
c                           
c on return
c --------- 
c nint	    = integer. The number of active (nonboundary) nodes. Also 
c             equal to the dimension of the assembled matrix.
c
c a, ja, ia = assembled matrix in compressed sparse row format.
c
c f	    = real array containing the right hand for the linears 
c             system to solve.
c 
c ierr	    = integer. Error message. If (ierr .ne. 0) on return
c             it means that one of the elements has a negative or zero
c             area probably because of a bad ordering of the nodes 
c             (see ijk above). Use the subroutine chkelmt to reorder
c             the nodes properly if necessary.
c iwk, jwk  = two integer work arrays of length nx each.
c
c-----------------------------------------------------------------------
      real*8 a(*),x(*),y(*),f(*),fs(*)
      integer ijk(node,*), nodcode(*),ia(*),ja(*),iwk(*),jwk(*)
      external xyk, funb, func, fung 
c     
      ierr = 0 
c     
c     take into boundary conditions to remove boundary nodes.
c     
      call bound (nx,nelx,ijk,nodcode,node,nint,jwk,
     *     x,y,f,iwk)
c     
c     assemble the matrix
c     
       call assmbo (nx,nelx,node,ijk,nodcode,x,y,
     *     a,ja,ia,f,iwk,jwk,ierr,xyk, funb, func, fung)
c     
c     if applicable (job .eq. 1) get heat source function
c     
      indic = 1
      if (job .eq. 1) 
     *     call hsourc (indic,nx,nelx,node,x,y,ijk,fs,f) 
c     
c     call diric for Dirichlet conditions
c     
      call diric(nx,nint,a,ja,ia,f)
c     done
      return
c------end of genfea ---------------------------------------------------
c-----------------------------------------------------------------------
      end
      subroutine genfea_wbc (nx,nelx,node,job,x,y,ijk,nodcode,fs,
     *     a,ja,ia,f,iwk,jwk,ierr,xyk)
c-----------------------------------------------------------------------
c this subroutine generates a finite element matrix in assembled form.
c the matrix is assembled in compressed sparse row format. See genfeu
c for matrices in unassembled form. The user must provide the grid, 
c (coordinates x, y and connectivity matrix ijk) as well as some 
c information on the nodes (nodcode) and the material properties 
c (the function K(x,y) above) in the form of a subroutine xyk. 
c----------------------------------------------------------------------
c Irene Moulitsas, moulitsa@cs.umn.edu  : It does not apply boundary
c                 conditions; variable nint is eliminated
c----------------------------------------------------------------------
c
c on entry:
c ---------
c 
c nx	    = integer . the number of nodes in the grid . 
c nelx	    = integer . the number of elements in the grid.
c node      = integer = the number of nodes per element (should be
c             set to three in this version). also the first dimension
c             of ijk
c job	    = integer. If job=0, it is assumed that there is no heat
c             source (i.e. fs = 0) and the right hand side
c             produced will therefore be a zero vector.
c             If job = 1 on entry then the contributions from the
c             heat source in each element are taken into account.
c 
c x, y      = two real arrays containing the coordinates of the nodes.
c 
c ijk       =  an integer array containing the connectivity matrix.
c              ijk(i,nel), i=1,2,..node, is the list of the nodes
c              constituting the element nel, ans listed in 
c              counter clockwise order.
c
c nodcode   = an integer array containing the boundary information for
c             each node with the following meaning.
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner node. [This node and the
c             corresponmding element are discarded.]
c
c fs	    = real array of length nelx on entry containing the heat 
c             source for each element (job = 1 only) 
c
c xyk	    = subroutine defining the material properties at each 
c	      element. Form: 
c 	      call xyk(nel,xyke,x,y,ijk,node) with on return
c             xyke =  material constant matrices. 
c	      for each element nel, xyke(1,nel),xyke(2,nel) 
c             and xyke(3,nel) represent the constants
c             K11, K22, and K12 at that element.
c                           
c on return
c --------- 
c a, ja, ia = assembled matrix in compressed sparse row format.
c
c f	    = real array containing the right hand for the linears 
c             system to solve.
c 
c ierr	    = integer. Error message. If (ierr .ne. 0) on return
c             it means that one of the elements has a negative or zero
c             area probably because of a bad ordering of the nodes 
c             (see ijk above). Use the subroutine chkelmt to reorder
c             the nodes properly if necessary.
c iwk, jwk  = two integer work arrays of length nx each.
c
c-----------------------------------------------------------------------
      real*8 a(*),x(*),y(*),f(*),fs(*)
      integer ijk(node,*), nodcode(*),ia(*),ja(*),iwk(*),jwk(*)
      external xyk, funb, func, fung 
c     
      ierr = 0 
c     
c     assemble the matrix
c     
       call assmbo (nx,nelx,node,ijk,nodcode,x,y,
     *     a,ja,ia,f,iwk,jwk,ierr,xyk, funb, func, fung)
c     
c     if applicable (job .eq. 1) get heat source function
c     
      indic = 1
      if (job .eq. 1) 
     *     call hsourc (indic,nx,nelx,node,x,y,ijk,fs,f) 
c
c     done
      return
c------end of genfea_wbc -----------------------------------------------
c-----------------------------------------------------------------------
       end
c----------------------------------------------------------------------- 
      subroutine genfeu (nx,nelx,node,job,x,y,ijk,nodcode,fs,
     *     nint,a,na,f,iwk,jwk,ierr,xyk)
c-----------------------------------------------------------------------
c this subroutine generates finite element matrices for heat 
c condution problem 
c
c                  - Div ( K(x,y) Grad u ) = f
c                    u = 0 on boundary 
c 
c (with Dirichlet boundary conditions). The matrix is returned 
c in unassembled form. The user must provide the grid, 
c (coordinates x, y and connectivity matrix ijk) as well as some 
c information on the nodes (nodcode) and the material properties 
c (the function K(x,y) above) in the form of a subroutine xyk. 
c
c----------------------------------------------------------------------
c
c on entry:
c ---------
c 
c nx	    = integer . the number of nodes in the grid . 
c nelx	    = integer . the number of elements in the grid.
c node      = integer = the number of nodes per element (should be
c             set to three in this version). also the first dimension
c             of ijk
c job	    = integer. If job=0, it is assumed that there is no heat
c             source (i.e. fs = 0) and the right hand side
c             produced will therefore be a zero vector.
c             If job = 1 on entry then the contributions from the
c             heat source in each element are taken into account.
c 
c na	    = integer. The first dimension of the array a. 
c             a is declared as an array of dimension a(na,node,node).
c
c x, y      = two real arrays containing the coordinates of the nodes.
c 
c ijk       =  an integer array containing the connectivity matrix.
c              ijk(i,nel), i=1,2,..node, is the list of the nodes
c              constituting the element nel, ans listed in 
c              counter clockwise order.
c
c xyk	    = subroutine defining the material properties at each 
c	      element. Form: 
c 	      call xyk(nel,xyke,x,y,ijk,node) with on return
c             xyke =  material constant matrices. 
c	      for each element nel, xyke(1,nel),xyke(2,nel) 
c             and xyke(3,nel) represent the constants
c             K11, K22, and K12 at that element.
c                           
c nodcode   = an integer array containing the boundary information for
c             each node with the following meaning.
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner node. [This node and the
c             corresponmding element are discarded.]
c
c fs	    = real array of length nelx on entry containing the heat 
c             source for each element (job = 1 only) 
c                           
c on return
c --------- 
c nint	    = integer. The number of active (nonboundary) nodes. Also 
c             equal to the dimension of the assembled matrix.
c
c a         = matrix in unassembled form. a(nel,*,*) contains the 
c             element matrix for element nel.
c
c f	    = real array containing the right hand for the linears 
c             system to solve, in assembled form. 
c 
c ierr	    = integer. Error message. If (ierr .ne. 0) on return
c             it means that one of the elements has a negative or zero
c             area probably because of a bad ordering of the nodes 
c             (see ijk above). Use the subroutine chkelmt to reorder
c             the nodes properly if necessary.
c iwk, jwk  = two integer work arrays of length nx each.
c
c-----------------------------------------------------------------------
      real*8 a(na,node,node),x(*),y(*),f(*), fs(*)
      integer ijk(node,*), nodcode(*),iwk(*),jwk(*)
      external xyk
c     
      ierr = 0 
c     
c     take boundary conditions into account to move boundary nodes to
c     the end..
c     
      call bound (nx,nelx,ijk,nodcode,node,nint,jwk,
     *     x,y,f,iwk)
c     
c     assemble the matrix
c     
      call unassbl (a,na,f,nx,nelx,ijk,nodcode,
     *     node,x,y,ierr,xyk) 
c     
c     if applicable (job .eq. 1) get heat source function
c     
      indic = 0
      if (job .eq. 1) 
     *     call hsourc (indic,nx,nelx,node,x,y,ijk,fs,f) 
c     
c     done
c     
      return
      end
c----- end of genfeu ---------------------------------------------------- 
      subroutine genfeu_wbc (nx,nelx,node,job,x,y,ijk,nodcode,fs,
     *          a,na,f,iwk,jwk,ierr,xyk)
c-----------------------------------------------------------------------
c this subroutine generates finite element matrices for heat 
c condution problem 
c
c                  - Div ( K(x,y) Grad u ) = f
c                    u = 0 on boundary 
c 
c (with Dirichlet boundary conditions). The matrix is returned 
c in unassembled form. The user must provide the grid, 
c (coordinates x, y and connectivity matrix ijk) as well as some 
c information on the nodes (nodcode) and the material properties 
c (the function K(x,y) above) in the form of a subroutine xyk. 
c
c----------------------------------------------------------------------
c moulitsa@cs   : It does not apply boundary conditions
c                 variable nint is eliminated
c----------------------------------------------------------------------
c
c on entry:
c ---------
c 
c nx	    = integer . the number of nodes in the grid . 
c nelx	    = integer . the number of elements in the grid.
c node      = integer = the number of nodes per element (should be
c             set to three in this version). also the first dimension
c             of ijk
c job	    = integer. If job=0, it is assumed that there is no heat
c             source (i.e. fs = 0) and the right hand side
c             produced will therefore be a zero vector.
c             If job = 1 on entry then the contributions from the
c             heat source in each element are taken into account.
c 
c na	    = integer. The first dimension of the array a. 
c             a is declared as an array of dimension a(na,node,node).
c
c x, y      = two real arrays containing the coordinates of the nodes.
c 
c ijk       =  an integer array containing the connectivity matrix.
c              ijk(i,nel), i=1,2,..node, is the list of the nodes
c              constituting the element nel, ans listed in 
c              counter clockwise order.
c
c xyk	    = subroutine defining the material properties at each 
c	      element. Form: 
c 	      call xyk(nel,xyke,x,y,ijk,node) with on return
c             xyke =  material constant matrices. 
c	      for each element nel, xyke(1,nel),xyke(2,nel) 
c             and xyke(3,nel) represent the constants
c             K11, K22, and K12 at that element.
c                           
c nodcode   = an integer array containing the boundary information for
c             each node with the following meaning.
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner node. [This node and the
c             corresponmding element are discarded.]
c
c fs	    = real array of length nelx on entry containing the heat 
c             source for each element (job = 1 only) 
c                           
c on return
c --------- 
c a         = matrix in unassembled form. a(nel,*,*) contains the 
c             element matrix for element nel.
c
c f	    = real array containing the right hand for the linears 
c             system to solve, in assembled form. 
c 
c ierr	    = integer. Error message. If (ierr .ne. 0) on return
c             it means that one of the elements has a negative or zero
c             area probably because of a bad ordering of the nodes 
c             (see ijk above). Use the subroutine chkelmt to reorder
c             the nodes properly if necessary.
c iwk, jwk  = two integer work arrays of length nx each.
c
c-----------------------------------------------------------------------
      real*8 a(na,node,node),x(*),y(*),f(*), fs(*)
      integer ijk(node,*), nodcode(*),iwk(*),jwk(*)
      external xyk
c     
      ierr = 0 
c     
c     assemble the matrix
c     
      call unassbl (a,na,f,nx,nelx,ijk,nodcode,
     *     node,x,y,ierr,xyk) 
c     
c     if applicable (job .eq. 1) get heat source function
c     
      indic = 0
      if (job .eq. 1) 
     *     call hsourc (indic,nx,nelx,node,x,y,ijk,fs,f) 
c     
c     done
c     
      return
      end
c----- end of genfeu_wbc ----------------------------------------------- 
      subroutine genfeu_lstif (nx,nelx,node,job,x,y,ijk,nodcode,fs,
     *          a,na,f,iwk,jwk,ierr,xyk)
c-----------------------------------------------------------------------
c this subroutine generates finite element matrices using unassmbl_lstif.
c The matrix is returned in unassembled form.
c The user must provide the grid, coordinates x, y and connectivity matrix
c ijk) as well as some information on the nodes (nodcode) and the material
c properties (the function K(x,y) above) in the form of a subroutine xyk. 
c
c----------------------------------------------------------------------
c moulitsa@cs.umn.edu   : It does not apply boundary conditions
c                         variable nint is eliminated
c----------------------------------------------------------------------
c
c on entry:
c ---------
c 
c nx	    = integer . the number of nodes in the grid . 
c nelx	    = integer . the number of elements in the grid.
c node      = integer = the number of nodes per element (should be
c             set to three in this version). also the first dimension
c             of ijk
c job	    = integer. If job=0, it is assumed that there is no heat
c             source (i.e. fs = 0) and the right hand side
c             produced will therefore be a zero vector.
c             If job = 1 on entry then the contributions from the
c             heat source in each element are taken into account.
c 
c na	    = integer. The first dimension of the array a. 
c             a is declared as an array of dimension a(na,node,node).
c
c x, y      = two real arrays containing the coordinates of the nodes.
c 
c ijk       =  an integer array containing the connectivity matrix.
c              ijk(i,nel), i=1,2,..node, is the list of the nodes
c              constituting the element nel, ans listed in 
c              counter clockwise order.
c
c xyk	    = subroutine defining the material properties at each 
c	      element. Form: 
c 	      call xyk(nel,xyke,x,y,ijk,node) with on return
c             xyke =  material constant matrices. 
c	      for each element nel, xyke(1,nel),xyke(2,nel) 
c             and xyke(3,nel) represent the constants
c             K11, K22, and K12 at that element.
c                           
c nodcode   = an integer array containing the boundary information for
c             each node with the following meaning.
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner node. [This node and the
c             corresponmding element are discarded.]
c
c fs	    = real array of length nelx on entry containing the heat 
c             source for each element (job = 1 only) 
c                           
c on return
c --------- 
c a         = matrix in unassembled form. a(nel,*,*) contains the 
c             element matrix for element nel.
c
c f	    = real array containing the right hand for the linears 
c             system to solve, in assembled form. 
c 
c ierr	    = integer. Error message. If (ierr .ne. 0) on return
c             it means that one of the elements has a negative or zero
c             area probably because of a bad ordering of the nodes 
c             (see ijk above). Use the subroutine chkelmt to reorder
c             the nodes properly if necessary.
c iwk, jwk  = two integer work arrays of length nx each.
c
c-----------------------------------------------------------------------
      real*8 a(na,node,node),x(*),y(*),f(*), fs(*)
      integer ijk(node,*), nodcode(*),iwk(*),jwk(*)
      external xyk, funb, func, fung
c     
      ierr = 0 
c     
c     assemble the matrix
c     
      call unassbl_lstif (a,na,f,nx,nelx,ijk,nodcode,
     *     node,x,y,ierr,xyk,funb,func,fung) 
c     
c     if applicable (job .eq. 1) get heat source function
c     
      indic = 0
      if (job .eq. 1) 
     *     call hsourc (indic,nx,nelx,node,x,y,ijk,fs,f) 
c     
c     done
c     
      return
      end
c----- end of genfeu_lstif ---------------------------------------------
c-----------------------------------------------------------------------
      subroutine assmb1 (u,nu,a,ja,ia,fu,f,nx,nelx,ijk,nodcode,
     *     node,iwk,jwk)
c--------------------------------------------------------------
c u	 = unassembled matrix u(na,node,node)
c nu	 = 1-st dimension of u
c a,ja,ia= assembled matrix on output
c fu	 = unassembled right hand side
c f      = right hand side (global load vector) assembled
c nx     = number of nodes at input
c nelx	 = number of elements at input
c ijk	 = connectivity matrix: for node k, ijk(*,k) point to the
c          nodes of element k.
c node	 = total number of nodal points in each element
c
c nodcode= boundary information list for each node with the
c	   following meaning:
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner point (corner points
c
c x,y   = real*8 arrays containing the $x$ and $y$ coordinates 
c	  resp. of the nodes.
c         K11, K22, and K12 at that element.
c iwk,jwk = two integer work arrays.
c ierr	= error message integer . 
c	  ierr = 0 --> normal return
c	  ierr = 1 --> negative area encountered (due to bad 
c	           numbering of nodes of an element- see
c		   message printed in unit iout). not used..
c iout	= output unit (not used here).
c--------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      real*8 u(nu,node,node),a(*),fu(node,*),f(*)
      integer ja(*),ia(*),ijk(node,*),iwk(*),jwk(*),nodcode(*)
c     max number of nonzeros per row allowed  = 200
c--------------------------------------------------------------
c     initialize
c--------------------------------------------------------------
      do 100 i=1,nx 
         f(i) = 0.0d0 
 100  continue
c     
c     initialize  pointer arrays.
c     
      do 5 k=1,nx+1
         ia(k) = 1
         jwk(k) = 0
 5    continue
      do 6 k=1,nelx
         do 59 j=1,node
            knod = ijk(j,k)
            ia(knod) = ia(knod) + 1
 59      continue
 6    continue
c---------------------------------------------------
      do 7 k=1, nx
         if (nodcode(k) .ge.1 ) ia(k)=ia(k)+1
 7    continue
c     
      ksav = ia(1)
      ia(1) = 1
      do 101 j=2, nx+1
         ksavn = ia(j)
         ia(j) = ia(j-1) +  ksav
         iwk(j-1) = ia(j-1)-1
         ksav = ksavn
 101  continue
c-----------------
c     main loop
c-----------------
      do 102 nel=1, nelx
c     
c     get nodal points
c     
         do 120 ka=1, node
            ii = ijk(ka,nel)
            f(ii) = f(ii) + fu(ka,nel)
c     
c     unpack row into jwk1
c     
            irowst = ia(ii)
            ilast  = iwk(ii) 
            do 109 k=irowst,ilast 
               jwk(ja(k)) = k
 109        continue
c     
            do 108 kb = 1,node
c     
c     column number = jj
c     
               jj = ijk(kb,nel)
               k = jwk(jj)
               if (k .eq. 0) then
                  ilast = ilast+1
                  jwk(jj) = ilast
                  ja(ilast) = jj
                  a(ilast) = u(nel,ka,kb) 
               else 
                  a(k) = a(k) + u(nel,ka,kb)
               endif
 108        continue
c     refresh jwk
            do 119 k=irowst,ilast 
               jwk(ja(k)) = 0
 119        continue
            iwk(ii) = ilast 
 120     continue
c     
 102  continue
      return
c---------end-of-assmb1---------------------------------------------- 
      end
c--------------------------------------------------------------------
