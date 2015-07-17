      subroutine inmesh (nmesh,iin,nx,nelx,node,x,y,nodcode,ijk,iperm)
      implicit none 
      real*8 x(*),y(*)
      integer nmesh,iin,nx,nelx,node,nodcode(nx),ijk(node,nelx), 
     *     iperm(nx) 
c-----------------------------------------------------------------------
c this subroutine selects and initializes a mesh among a few
c choices. So far there are 9 initial meshes provided and the user can
c also enter his own mesh as a 10th option. 
c 
c on entry:
c--------- 
c nmesh	    = integer indicating the mesh chosen. nmesh=1,...,9
c             corresponds to one of the 9 examples supplied by
c             SPARSKIT. nmesh = 0 is a user supplied initial mesh.
c             see below for additional information for the format.
c iin       = integer containing the I/O unit number where to read
c             the data from in case nmesh = 1. A dummy integer 
c             otherwise.
c node      = integer = the number of nodes per element (should be
c             set to node=3 in this version). node is also the first 
c             dimension of the array ijk.
c
c on return
c ---------
c nx	    = integer . the number of nodes
c nelx	    = integer . the number of elements
c x, y      = two real arrays containing the coordinates of the nodes.
c nodcode   = an integer array containing the boundary information for
c             each node with the following meaning.
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner node. 
c                           
c ijk(node,*)= an integer array containing the connectivity matrix.
c
c-----------------------------------------------------------------------
c format for user supplied mesh (when nmesh = 7) 
c
c option nmesh = 0, is a user definied initial mesh. 
c--------- 
c format is as follows:
c line 1: two integers, the first containing the number of nodes 
c         the second the number of elements.
c line 2: to line nx+1:  node information.
c        enter the following, one line per node:
c     * the number of the node in the numbering chosen (integer 
c          taking the values 1 to nx), followed by,
c        * the coordinates of the nodes (2 reals)  followed by 
c	   the boundary information, an integer taking one of the 
c          values 0, 1, or 2,  with the meaning explained above.
c
c line nx+2 to nx+nelx+1: connectivity matrix
c       enter the following one line per element: 
c       * number of the element in the numbering chosen, followed by
c       * The three numbers of the nodes (according to the above numbering
c       of the nodes) that constitute the element, in a counter clock-wise
c       order (this is in fact not important since it is checked by the
c       subroutine chkelemt). 
c
c AN EXAMPLE: consisting on one single element (a triangle) 
c------------ 
c    3    1
c    1    0.0000    0.0000    2    
c    2    4.0000    0.0000    2     
c    3    0.0000    4.0000    2                  
c    1    1    2    3
c
c----------------------------------------------------------------------- 
c     local variables      
      integer i, j, ii
c
c     print *, ' ----- nmesh = ', nmesh 
      goto (10,1,2,3,4,5,6,7,8,9) nmesh+1 
 1    continue
      call fmesh1  (nx,nelx,node,x,y,nodcode,ijk)
      goto 18
 2    continue
      call fmesh2  (nx,nelx,node,x,y,nodcode,ijk)
      goto 18
 3    continue
      call fmesh3  (nx,nelx,node,x,y,nodcode,ijk)
      goto 18
 4    continue
      call fmesh4  (nx,nelx,node,x,y,nodcode,ijk)
      goto 18
 5    continue
      call fmesh5  (nx,nelx,node,x,y,nodcode,ijk)
      goto 18
 6    continue
      call fmesh6 (nx,nelx,node,x,y,nodcode,ijk)
      goto 18
 7    continue
      call fmesh7 (nx,nelx,node,x,y,nodcode,ijk,iperm) 
      goto 18
 8    continue
      call fmesh8 (nx,nelx,node,x,y,nodcode,ijk,iperm) 
      goto 18
 9    continue
      call fmesh9(nx,nelx,node,x,y,nodcode,ijk,iperm) 
      goto 18
 10   continue
c
c-------option 0 : reading mesh from IO unit iin.
c
      read (iin,*) nx, nelx
c     
      do 16 i=1,nx
         read(iin,*) ii,x(ii),y(ii),nodcode(ii)
 16   continue
      do 17 i=1,nelx
         read(iin,*) ii,(ijk(j,ii),j=1,node)
         if (ii. gt. nelx) nelx = ii
 17   continue
c-----------------------------------------------------------------------
 18    continue
c     and return
      return
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
      subroutine fmesh1 (nx,nelx,node,x,y,nodcode,ijk)
c--------------------------------------------------------------
c 
c initial mesh for a simple square with two elemnts
c      3             4
c       --------------
c       |          . |
c       |   2    .   |
c       |      .     |
c       |   .    1   |
c       | .          |
c       --------------
c      1              2
c--------------------------------------------------------------
c input parameters: node = first dimensoin of ijk (must be .ge. 3)
c output parameters:
c    nx    = number of nodes
c    nelx = number of elemnts
c    (x(1:nx), y(1:nx)) = coordinates of nodes
c    nodcode(1:nx) = integer code for each node with the 
c	    following meening:
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner point.
c   ijk(1:3,1:nelx) = connectivity matrix. for a given element
c	    number nel, ijk(k,nel), k=1,2,3 represent the nodes
c	    composing the element nel.
c--------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*),nodcode(*),ijk(node,*)
      real*8 x1(4),y1(4)
      integer ijk1(2),ijk2(2),ijk3(2)
c--------------------------------------------------------------
c coordinates of nodal points
c--------------------------------------------------------------
      data  x1/0.0, 1.0, 0.0, 1.0/
      data  y1/0.0, 0.0, 1.0, 1.0/
c
c------------------|--|
c elements         1  2
c------------------|--|
      data ijk1   /1, 1/
      data ijk2   /2, 4/
      data ijk3   /4, 3/
c
      nx = 4
c
      do 1 k=1, nx
         x(k) = x1(k)
         y(k) = y1(k)
         nodcode(k) = 1
 1    continue
c     
      nodcode(2) = 2
      nodcode(3) = 2
c
      nelx = 2
c
      do  2 k=1,nelx
         ijk(1,k) = ijk1(k) 
         ijk(2,k) = ijk2(k)
         ijk(3,k) = ijk3(k)
 2    continue
c
      return
      end
c-----------------------------------------------------------------------
      subroutine fmesh2 (nx,nelx,node,x,y,nodcode,ijk)
c---------------------------------------------------------------
c initial mesh for a simple D-shaped region with 4 elemnts
c       6
c       | .
c       |    .
c       |      .
c       |   4     .
c       |           .
c     4 -------------- 5
c       |          . |
c       |   3    .   |
c       |      .     |
c       |   .    2   |
c       | .          |
c       --------------
c       | 2         . 3
c       |         .
c       |   1   .
c       |     .
c       |  .
c       |.
c       1  
c--------------------------------------------------------------
c input parameters: node = first dimensoin of ijk (must be .ge. 3)
c output parameters:
c    nx    = number of nodes
c    nelx = number of elemnts
c    (x(1:nx), y(1:nx)) = coordinates of nodes
c    nodcode(1:nx) = integer code for each node with the 
c	    following meening:
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner point.
c   ijk(1:3,1:nelx) = connectivity matrix. for a given element
c	    number nel, ijk(k,nel), k=1,2,3 represent the nodes
c	    composing the element nel.
c
c--------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*),nodcode(*),ijk(node,*)
      real*8 x1(6),y1(6)
      integer ijk1(4),ijk2(4),ijk3(4)
c--------------------------------------------------------------
c coordinates of nodal points
c--------------------------------------------------------------
      data x1/0.0, 0.0, 1.0, 0.0, 1.0, 0.0/
      data y1/0.0, 1.0, 1.0, 2.0, 2.0, 3.0/
c
c------------------|--|--|--|
c elements         1  2  3  4
c------------------|--|--|--|
      data ijk1   /1, 2, 2, 4/
      data ijk2   /3, 3, 5, 5/
      data ijk3   /2, 5, 4, 6/
c
      nx = 6
c
      do 1 k=1, nx
         x(k) = x1(k)
         y(k) = y1(k)
         nodcode(k) = 1
 1    continue
c
      nelx = 4
c     
      do  2 k=1,nelx
         ijk(1,k) = ijk1(k) 
         ijk(2,k) = ijk2(k)
         ijk(3,k) = ijk3(k)
 2    continue
c
      return
      end
c-----------------------------------------------------------------------
      subroutine fmesh3 (nx,nelx,node,x,y,nodcode,ijk)
c---------------------------------------------------------------
c initial mesh for a C-shaped region composed of 10 elements --
c
c
c      10           11            12             
c       ---------------------------
c       |          . |          . |
c       |  7     .   |   9    .   |
c       |      .     |      .     |
c       |   .    8   |   .   10   |
c       | .          | .          |
c     7 ---------------------------
c       |          . |8           9
c       |   5    .   |
c       |      .     |
c       |   .    6   |
c     4 | .          |5           6
c       ---------------------------
c       |          . |          . |
c       |   1    .   |  3     .   |
c       |      .     |      .     |
c       |   .    2   |   .   4    |
c       | .          | .          |
c       ---------------------------
c      1             2            3
c
c--------------------------------------------------------------
c input parameters: node = first dimensoin of ijk (must be .ge. 3)
c    nx    = number of nodes
c    nelx = number of elemnts
c    (x(1:nx), y(1:nx)) = coordinates of nodes
c    nodcode(1:nx) = integer code for each node with the 
c	    following meening:
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner point.
c   ijk(1:3,1:nelx) = connectivity matrix. for a given element
c	    number nel, ijk(k,nel), k=1,2,3 represent the nodes
c	    composing the element nel.
c
c--------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*),nodcode(*),ijk(node,*)
      real*8 x1(12),y1(12)
      integer ijk1(10),ijk2(10),ijk3(10)
c--------------------------------------------------------------
c coordinates of nodal points
c--------------------------------------------------------------
      data x1/0.0,1.0,2.0,0.0,1.0,2.0,0.0,1.0,2.0,0.0,1.0,2.0/
      data y1/0.0,0.0,0.0,1.0,1.0,1.0,2.0,2.0,2.0,3.0,3.0,3.0/
c
c------------------|--|--|--|--|--|--|---|---|---|
c elements         1  2  3  4  5  6  7   8   9  10
c------------------|--|--|--|--|--|--|---|---|---|
      data ijk1   /1, 1, 2, 2, 4, 4, 7,  7,  8, 8/
      data ijk2   /5, 2, 6, 3, 8, 5, 11, 8, 12, 9/
      data ijk3   /4, 5, 5, 6, 7, 8, 10, 11,11, 12/
c
      nx = 12
c
      do 1 k=1, nx
         x(k) = x1(k)
         y(k) = y1(k)
         nodcode(k) = 1
 1    continue
c
      nodcode(3) = 2
      nodcode(10) = 2
      nodcode(9) = 2
c     
      nelx = 10
c     
      do  2 k=1,nelx
         ijk(1,k) = ijk1(k) 
         ijk(2,k) = ijk2(k)
         ijk(3,k) = ijk3(k)
 2    continue
c
      return
      end
c-----------------------------------------------------------------------
      subroutine fmesh4 (nx,nelx,node,x,y,nodcode,ijk)
c----------------------------------------------------------------------- 
c initial mesh for a C-shaped region composed of 10 elements --
c      10                   11 
c       +------------------+ .
c       | .                |    .
c       |    .       8     |       . 12
c       |        .         |  9   . |
c       |     7      .     |   .    |
c     7 |                . | .   10 |
c       -------------------+--------+ 9
c       |                 .| 8 
c       |     5       .    |
c       |         .        |
c       |    .       6     |
c       |.                 | 5      6
c    4  +------------------+--------+ 
c       |               .  | .   4  |
c       |    1       .     |    .   |
c       |        .         |  3    .| 3
c       |    .        2    |    .
c       | .                | .
c       -------------------- 
c       1                  2 
c--------------------------------------------------------------
c input parameters: node = first dimensoin of ijk (must be .ge. 3)
c    nx    = number of nodes
c    nelx = number of elemnts
c    (x(1:nx), y(1:nx)) = coordinates of nodes
c    nodcode(1:nx) = integer code for each node with the 
c	    following meening:
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner point.
c   ijk(1:3,1:nelx) = connectivity matrix. for a given element
c	    number nel, ijk(k,nel), k=1,2,3 represent the nodes
c	    composing the element nel.
c
c--------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*),nodcode(*),ijk(node,*)
      real*8 x1(12),y1(12)
      integer ijk1(10),ijk2(10),ijk3(10)
c--------------------------------------------------------------
c coordinates of nodal points
c--------------------------------------------------------------
      data x1/0.0,1.0,1.5,0.0,1.0,1.5,0.0,1.0,1.5,0.0,1.0,1.5/
      data y1/0.0,0.0,0.5,1.0,1.0,1.0,2.0,2.0,2.0,3.0,3.0,2.5/
c
c------------------|--|--|--|--|--|--|---|---|---|
c elements         1  2  3  4  5  6  7   8   9  10
c------------------|--|--|--|--|--|--|---|---|---|
      data ijk1   /1, 1, 2, 5, 4, 4, 7, 10,  8, 8/
      data ijk2   /5, 2, 3, 3, 8, 5, 8,  8, 12, 9/
      data ijk3   /4, 5, 5, 6, 7, 8, 10, 11,11, 12/
c
      nx = 12
c     
      do 1 k=1, nx
         x(k) = x1(k)
         y(k) = y1(k)
         nodcode(k) = 1
 1    continue
c
      nodcode(6) = 2
      nodcode(9) = 2
c
      nelx = 10
c     
      do  2 k=1,nelx
         ijk(1,k) = ijk1(k) 
         ijk(2,k) = ijk2(k)
         ijk(3,k) = ijk3(k)
 2    continue
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine fmesh5 (nx,nelx,node,x,y,nodcode,ijk)
c---------------------------------------------------------------
c     initial mesh for a whrench shaped region composed of 14 elements --
c     
c                                      13            15    
c                                        . ----------.           |-3
c                                      .   .   13  .   .         |
c                                   .   12   .   .  14    .      |
c 9        10        11       12  .            . 14        . 16  |
c ----------------------------------------------------------     |-2
c |       . |       . |       . |            . |                 |  
c | 1   .   |  3  .   |  5  .   |    7   .     |                 |
c |   .  2  |   .  4  |   .  6  |     .    8   |                 |
c |.        |.        |.        | .            |                 |
c -----------------------------------------------------------    |-1
c 1         2         3       4  .           6 .           . 8   |
c                                   .   9    .   .   11   .      |
c                                      .   .  10    .   .        |
c                                        .___________.           |-0
c                                       5             7
c
c 0---------1--------2----------3--------------4-------------5
c--------------------------------------------------------------
c input parameters: node = first dimensoin of ijk (must be .ge. 3)
c    nx    = number of nodes
c    nelx = number of elemnts
c    (x(1:nx), y(1:nx)) = coordinates of nodes
c    nodcode(1:nx) = integer code for each node with the 
c	    following meening:
c	nodcode(i) = 0 -->  node i is internal
c	nodcode(i) = 1 -->  node i is a boundary but not a corner point
c	nodcode(i) = 2 -->  node i is a corner point.
c   ijk(1:3,1:nelx) = connectivity matrix. for a given element
c	    number nel, ijk(k,nel), k=1,2,3 represent the nodes
c	    composing the element nel.
c
c--------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*),nodcode(*),ijk(node,*)
      real*8 x1(16),y1(16)
      integer ijk1(14),ijk2(14),ijk3(14)
c--------------------------------------------------------------
c     coordinates of nodal points
c--------------------------------------------------------------
      data x1/0.,1.,2.,3.,3.5,4.,4.5,5.,0.,1.,2.,3.,3.5,4.,4.5,5./
      data y1/1.,1.,1.,1.,0.,1.,0.,1.,2.,2.,2.,2.,3.,2.,3.,2./
c     
c------------------|--|--|--|--|--|--|---|---|---|--|---|---|---|
c elements         1  2  3  4  5  6  7   8   9  10  11  12  13  14
c------------------|--|--|--|--|--|--|---|---|---|--|---|---|---|
      data ijk1   /1, 1, 2, 2, 3, 3, 4,  4,  4,  5, 6, 12, 14, 14/
      data ijk2   /10,2,11, 3,12, 4,14,  6,  5,  7, 7, 14, 15, 16/
      data ijk3   /9,10,10,11,11,12,12, 14,  6,  6, 8, 13, 13, 15/
c
      nx = 16
c     
      do 1 k=1, nx
         x(k) = x1(k)
         y(k) = y1(k)
         nodcode(k) = 1
 1    continue
c     
      nodcode(9) = 2
      nodcode(8) = 2
      nodcode(16) = 2
c     
      nelx = 14
c     
      do  2 k=1,nelx
         ijk(1,k) = ijk1(k) 
         ijk(2,k) = ijk2(k)
         ijk(3,k) = ijk3(k)
 2    continue
c     
      return
      end
c----------------------------------------------------------------------- 
      subroutine fmesh6 (nx,nelx,node,x,y,nodcode,ijk)
c---------------------------------------------------------------
c this generates a finite element mesh for an ellipse-shaped 
c domain.
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*),nodcode(*),ijk(node,*), ijktr(200,3)
      integer nel(200) 
c--------------------------------------------------------------
c     coordinates of nodal points
c--------------------------------------------------------------
      nd = 8 
      nr = 3
c     
c     define axes of ellipse
c     
      a = 2.0
      b = 1.30 
c     
      nx = 1
      pi = 4.0* atan(1.0) 
      theta = 2.0 * pi / real(nd)  
      x(1) = 0.0
      y(1) = 0.0 
      delr = a / real(nr)
      nx = 0
      do i = 1, nr 
         ar = real(i)*delr
         br = ar*b / a
         do j=1, nd 
            nx = nx+1
            x(nx) = a +ar*cos(real(j)*theta)
            y(nx) = b +br*sin(real(j)*theta)
c            write (13,*) ' nod ', nx, ' x,y', x(nx), y(nx) 
            nodcode(nx) = 0
            if (i .eq. nr) nodcode(nx) = 1
         enddo
      enddo
c     
      nemax = 200
      call dlauny(x,y,nx,ijktr,nemax,nelx)
c     
c     print *, ' delauny -- nx, nelx ', nx, nelx
      do 3 j=1,nx
         nel(j) = 0
 3    continue
c     transpose ijktr into ijk and count the number of 
c     elemnts to which each node belongs
c     
      do 4 j=1, nelx
         do 41 k=1, node
	    i = ijktr(j,k)
	    ijk(k,j) = i
	    nel(i) = nel(i)+1
 41      continue
 4    continue
c     
c     take care of ordering within each element
c     
      call chkelmt (nx, x, y, nelx, ijk, node)
c     
      return
      end
c--------------------------------------------------------
      subroutine fmesh7 (nx,nelx,node,x,y,nodcode,ijk,iperm)
      implicit none 
      real*8 x(*),y(*) 
      integer nx,nelx,node,nodcode(nx),ijk(node,nelx),iperm(nx) 
c---------------------------------------------------------------
c     this generates a U-shaped domain with an elliptic inside.
c     then a Delauney triangulation is used to generate the mesh.
c     mesh needs to be post-processed -- see inmesh -- 
c---------------------------------------------------------------
      integer nr,nsec,i,k,nemax,j,nel(200),ijktr(200,3),nodexc
      real*8 a,b,x1, y1, x2, y2, xcntr,ycntr, rad,pi,delr,xnew,ynew,
     *     arx, ary, cos, sin, theta,excl   
c--------------------------------------------------------------
c     coordinates of nodal points
c--------------------------------------------------------------
      data x1/0.0/,y1/0.0/,x2/6.0/,y2/6.0/,nodexc/1/
      xcntr = x2/3.0
      ycntr = y2/2.0 
      rad = 1.8
      nsec = 20 
      nr = 3 
c
c     exclusion zone near the boundary 
c 
      excl = 0.02*x2 
c-----------------------------------------------------------------------
c     enter the four corner points. 
c-----------------------------------------------------------------------
      nx = 1 
      x(nx) = x1
      y(nx) = y1
      nodcode(nx) = 1
      nx = nx+1
      x(nx) = x2
      y(nx) = y1
      nodcode(nx) = 1
      nx = nx+1
      x(nx) = x2
      y(nx) = y2 
      nodcode(nx) = 1
      nx = nx+1
      x(nx) = x1
      y(nx) = y2
      nodcode(nx) = 1
c     
c     define axes of ellipse 
c     
      a = 2.0
      b = 1.30
c-----------------------------------------------------------------------
      pi = 4.0*atan(1.0) 
      delr = a / real(nr)
      do 2 i = 1, nsec 
         theta = 2.0 * real(i-1) * pi / real(nsec)  
         xnew = xcntr + rad*cos(theta)
         ynew = ycntr + rad*b*sin(theta)/a 
         if ((xnew .ge. x2) .or. (xnew .le. x1) .or. (ynew .ge. y2)
     *        .or. (ynew .le. y1)) goto 2 
         nx = nx+1 
         x(nx) = xnew 
         y(nx) = ynew 
         nodcode(nx) = nodexc 
         arx = delr*cos(theta)
         ary = delr*b*sin(theta)/a 
c     
c     while inside domain do: 
c     
 1       continue 
         xnew = x(nx) + arx 
         ynew = y(nx) + ary 
         if (xnew .ge. x2) then 
            x(nx) = x2 
            nodcode(nx) = 1
         else if (xnew .le. x1) then
            x(nx) = x1 
            nodcode(nx) = 1         
         else if (ynew .ge. y2) then 
            y(nx) = y2 
            nodcode(nx) = 1
         else if (ynew .le. y1) then
            y(nx) = y1 
            nodcode(nx) = 1         
         else
            nx = nx+1 
            x(nx) = xnew 
            y(nx) = ynew 
            nodcode(nx) = 0
            call clos2bdr(nx,xnew,ynew,x,y,x1,x2,y1,y2,excl,nodcode)
         endif 
c        write (13,*) ' nod ', nx, ' x,y', x(nx), y(nx)
c     *         ,' arx--ary ', arx, ary
         arx = arx*1.2 
         ary = ary*1.2
         if (nodcode(nx) .le. 0) goto 1 
 2    continue 
c     
      nemax = 200
      call dlauny(x,y,nx,ijktr,nemax,nelx)
c     
c     print *, ' delauny -- nx, nelx ', nx, nelx
      do 3 j=1,nx
         nel(j) = 0
 3    continue
c     
c     transpose ijktr into ijk and count the number of 
c     elemnts to which each node belongs
c     
      do 4 j=1, nelx
         do 41 k=1, node
	    i = ijktr(j,k)
	    ijk(k,j) = i
	    nel(i) = nel(i)+1
 41      continue
 4    continue
c
c     this mesh needs cleaning up --
c 
      call cleanel (nelx,ijk, node,nodcode,nodexc) 
      call cleannods(nx,x,y,nelx,ijk,node,nodcode,iperm)       
c
c     take care of ordering within each element
c
      call chkelmt (nx, x, y, nelx, ijk, node)
      return
      end
c-----------------------------------------------------------------------
      subroutine fmesh8 (nx,nelx,node,x,y,nodcode,ijk,iperm)
      implicit none 
      real*8 x(*),y(*) 
      integer nx,nelx,node,nodcode(nx),ijk(node,nelx),iperm(nx) 
c---------------------------------------------------------------
c     this generates a small rocket type shape inside a rectangle
c     then a Delauney triangulation is used to generate the mesh.
c     mesh needs to be post-processed -- see inmesh -- 
c---------------------------------------------------------------
      integer nr,nsec,i,k,nemax,j,nel(1500),ijktr(1500,3),nodexc 
      real*8 a,b,x1, y1, x2, y2, xcntr,ycntr, rad,pi,delr,xnew,ynew,
     *     arx, ary, cos, sin, theta,radi,excl 
c--------------------------------------------------------------
c     coordinates of corners + some additional data 
c--------------------------------------------------------------
      data x1/0.0/,y1/0.0/,x2/6.0/,y2/6.0/,nodexc/3/  
      xcntr = 4.0 
      ycntr = y2/2.0 
      rad = 0.6
c
c     exclusion zone near the boundary. 
c 
      excl = 0.02*x2 
      nsec = 30
      nr = 4
c-----------------------------------------------------------------------
c     enter the four corner points. 
c-----------------------------------------------------------------------
      nx = 1 
      x(nx) = x1
      y(nx) = y1
      nodcode(nx) = 1
      nx = nx+1
      x(nx) = x2
      y(nx) = y1
      nodcode(nx) = 1
      nx = nx+1
      x(nx) = x2
      y(nx) = y2 
      nodcode(nx) = 1
      nx = nx+1
      x(nx) = x1
      y(nx) = y2
      nodcode(nx) = 1
c     
c     define axes of ellipse /circle / object
c     
      a = 2.0
      b = 1.0
c-----------------------------------------------------------------------
      pi = 4.0*atan(1.0) 
      delr = 2.0*rad / real(nr)
      do 2 i = 1, nsec 
         theta = 2.0*real(i-1) * pi / real(nsec) 
         if (theta .gt. pi) theta = theta - 2.0*pi 
         radi=rad*(1.0+0.05*((pi/2.0)**2-theta**2)**2)
     *        /(1.0+0.05*(pi/2.0)**4) 
         arx = radi*cos(theta)/real(nr) 
         ary = radi*sin(theta)/real(nr) 
c      a hack!
c         arx = (abs(theta)+0.25)*cos(theta)/real(nr) 
c         ary = (abs(theta)+0.25)*sin(theta)/real(nr) 
c
         xnew = xcntr + radi*cos(theta)
         ynew = ycntr + radi*b*sin(theta)/a 
         if ((xnew .ge. x2) .or. (xnew .le. x1) .or. (ynew .ge. y2)
     *        .or. (ynew .le. y1)) goto 2 
         nx = nx+1 
         x(nx) = xnew 
         y(nx) = ynew 
         nodcode(nx) = nodexc 
c     
c     while inside domain do: 
c     
 1       continue 
         xnew = xnew  + arx 
         ynew = ynew  + ary 
         if (xnew .ge. x2) then 
            x(nx) = x2 
            nodcode(nx) = 1
         else if (xnew .le. x1) then
            x(nx) = x1 
            nodcode(nx) = 1         
         else if (ynew .ge. y2) then 
            y(nx) = y2 
            nodcode(nx) = 1
         else if (ynew .le. y1) then
            y(nx) = y1 
            nodcode(nx) = 1         
c
c     else we can add this as interior point
c 
         else
            nx = nx+1 
            x(nx) = xnew 
            y(nx) = ynew 
            nodcode(nx) = 0
c
c     do something if point is too close to boundary 
c 
            call clos2bdr(nx,xnew,ynew,x,y,x1,x2,y1,y2,excl,nodcode)
         endif
c            
         arx = arx*1.1
         ary = ary*1.1 
         if (nodcode(nx) .eq. 0) goto 1 
 2    continue 
c     
      nemax = 1500 
      call dlauny(x,y,nx,ijktr,nemax,nelx)
c     
      print *, ' delauney -- nx, nelx ', nx, nelx
      do 3 j=1,nx
         nel(j) = 0
 3    continue
c-----------------------------------------------------------------------     
c     transpose ijktr into ijk and count the number of 
c     elemnts to which each node belongs
c----------------------------------------------------------------------- 
      do 4 j=1, nelx
         do 41 k=1, node
	    i = ijktr(j,k)
	    ijk(k,j) = i
	    nel(i) = nel(i)+1
 41      continue
 4    continue
c
c     this mesh needs cleaning up --
c 
      call cleanel (nelx,ijk, node,nodcode,nodexc) 
      call cleannods(nx,x,y,nelx,ijk,node,nodcode,iperm)       
c
c     take care of ordering within each element
c
      call chkelmt (nx, x, y, nelx, ijk, node)
      return
      end
c-----------------------------------------------------------------------  
      subroutine fmesh9 (nx,nelx,node,x,y,nodcode,ijk,iperm)
      implicit none 
      real*8 x(*),y(*) 
      integer nx,nelx,node,nodcode(nx),ijk(node,nelx),iperm(nx) 
c---------------------------------------------------------------
c     this generates a U-shaped domain with an elliptic inside.
c     then a Delauney triangulation is used to generate the mesh.
c     mesh needs to be post-processed -- see inmesh -- 
c---------------------------------------------------------------
      integer nr,nsec,i,k,nemax,j,nel(1500),ijktr(1500,3),nodexc 
      real*8 x1, y1, x2, y2, xcntr,ycntr, rad,pi,delr,xnew,ynew,
     *     arx, ary, cos, sin, theta,excl   
c--------------------------------------------------------------
c     coordinates of nodal points
c--------------------------------------------------------------
      data x1/0.0/,y1/0.0/,x2/11.0/,y2/5.5/,nodexc/3/ 
      xcntr = 1.50
      ycntr = y2/2.0 
      rad = 0.6
      nsec = 30
      nr = 3 
c
c-----------------------------------------------------------------------
c     enter the four corner points. 
c-----------------------------------------------------------------------
      nx = 1 
      x(nx) = x1
      y(nx) = y1
      nodcode(nx) = 1
      nx = nx+1
      x(nx) = x2
      y(nx) = y1
      nodcode(nx) = 1
      nx = nx+1
      x(nx) = x2
      y(nx) = y2 
      nodcode(nx) = 1
      nx = nx+1
      x(nx) = x1
      y(nx) = y2
      nodcode(nx) = 1
c     
c     define axes of ellipse 
c     
c-----------------------------------------------------------------------
      pi = 4.0*atan(1.0) 
      delr = rad / real(nr) 
      do 2 i = 1, nsec 
         theta = 2.0 * real(i-1) * pi / real(nsec)  
         xnew = xcntr + rad*cos(theta)
         ynew = ycntr + rad*sin(theta) 
         if ((xnew .ge. x2) .or. (xnew .le. x1) .or. (ynew .ge. y2)
     *        .or. (ynew .le. y1)) goto 2 
         nx = nx+1 
         x(nx) = xnew 
         y(nx) = ynew 
         nodcode(nx) = nodexc 
         arx = delr*cos(theta)
         ary = delr*sin(theta)
c 
c     exclusion zone near the boundary 
c
c         excl = 0.1*delr 
         excl = 0.15*delr
c     
c     while inside domain do: 
c     
 1       continue 
         xnew = x(nx) + arx 
         ynew = y(nx) + ary 
         if (xnew .ge. x2) then 
            x(nx) = x2 
            nodcode(nx) = 1
         else if (xnew .le. x1) then
            x(nx) = x1 
            nodcode(nx) = 1         
         else if (ynew .ge. y2) then 
            y(nx) = y2 
            nodcode(nx) = 1
         else if (ynew .le. y1) then
            y(nx) = y1 
            nodcode(nx) = 1         
         else
            nx = nx+1 
            x(nx) = xnew 
            y(nx) = ynew 
            nodcode(nx) = 0
            call clos2bdr(nx,xnew,ynew,x,y,x1,x2,y1,y2,excl,nodcode)
         endif 
         arx = arx*1.1 
         ary = ary*1.1
         excl = excl*1.1 
         if (nodcode(nx) .le. 0) goto 1 
 2    continue 
c     
      nemax = 1500
      call dlauny(x,y,nx,ijktr,nemax,nelx)
c     
c     print *, ' delauny -- nx, nelx ', nx, nelx
      do 3 j=1,nx
         nel(j) = 0
 3    continue
c     
c     transpose ijktr into ijk and count the number of 
c     elemnts to which each node belongs
c     
      do 4 j=1, nelx
         do 41 k=1, node
	    i = ijktr(j,k)
	    ijk(k,j) = i
	    nel(i) = nel(i)+1
 41      continue
 4    continue
c
c     this mesh needs cleaning up --
c 
      call cleanel (nelx,ijk, node,nodcode,nodexc) 
      call cleannods(nx,x,y,nelx,ijk,node,nodcode,iperm)       
c
c     take care of ordering within each element
c
      call chkelmt (nx, x, y, nelx, ijk, node)
      return
      end
c-----------------------------------------------------------------------
      subroutine clos2bdr (nx,xnew,ynew,x,y,x1,x2,y1,y2,excl,nodcode) 
      implicit none 
      integer nx,nodcode(nx)
      real*8 x(nx),y(nx),xnew,ynew,x1,x2,y1,y2,excl 
c----------------------------------------------------------------------- 
c     takes care of case where a point generated is too close to the
c     boundary  -- in this case projects the previous point to the
c     rectangle boundary == that makes some exclusion criterion
c     violated... does a simple job.
c----------------------------------------------------------------------- 
      if (xnew .ge. x2-excl) then
         x(nx) = x2
         y(nx) = y(nx-1) 
         nodcode(nx) = 1 
      endif
      if (xnew .le. x1+excl) then
         x(nx) = x1 
         y(nx) = y(nx-1) 
         nodcode(nx) = 1 
      endif
      if (ynew .ge. y2-excl) then
         y(nx) = y2
         x(nx) = x(nx-1) 
         nodcode(nx) = 1 
      endif
      if (ynew .le. y1+excl) then
         y(nx) = y1 
         x(nx) = x(nx-1) 
         nodcode(nx) = 1 
      endif
c     
      return
      end
