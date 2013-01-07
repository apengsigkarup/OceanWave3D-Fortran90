c-----------------------------------------------------------------------
c     contains the functions needed for defining the PDE poroblems. 
c
c     first for the scalar 5-point and 7-point PDE 
c-----------------------------------------------------------------------
      function afun (x,y,z)
      real*8 afun, x,y,z 
      afun = -1.0d0
      return 
      end

      function bfun (x,y,z)
      real*8 bfun, x,y,z 
      bfun = -1.0d0
      return 
      end

      function cfun (x,y,z)
      real*8 cfun, x,y,z 
      cfun = -1.0d0
      return 
      end

      function dfun (x,y,z)
      real*8 dfun, x,y,z 
      data gamma /100.0/ 
c     dfun = gamma * exp( x * y )
      dfun = 10.d0
      return 
      end

      function efun (x,y,z)
      real*8 efun, x,y,z
      data gamma /100.0/ 
c     efun = gamma * exp( (- x) * y ) 
      efun = 0.d0
      return 
      end

      function ffun (x,y,z)
      real*8 ffun, x,y,z 
      ffun = 0.0
      return 
      end

      function gfun (x,y,z)
      real*8 gfun, x,y,z 
      gfun = 0.0 
      return 
      end

      function hfun(x, y, z)
      real*8 hfun, x, y, z
      hfun = 0.0
      return
      end

      function betfun(side, x, y, z)
      real*8 betfun, x, y, z
      character*2 side
      betfun = 1.0
      return
      end

      function gamfun(side, x, y, z)
      real*8 gamfun, x, y, z
      character*2 side
      if (side.eq.'x2') then
         gamfun = 5.0
      else if (side.eq.'y1') then
         gamfun = 2.0
      else if (side.eq.'y2') then
         gamfun = 7.0
      else
         gamfun = 0.0
      endif
      return
      end

c-----------------------------------------------------------------------
c     functions for the block PDE's 
c-----------------------------------------------------------------------
      subroutine afunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
            coeff((j-1)*nfree+i) = 0.0d0
 1       continue
         coeff((j-1)*nfree+j) = -1.0d0
 2    continue
      return 
      end

      subroutine bfunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
            coeff((j-1)*nfree+i) = 0.0d0
 1       continue
         coeff((j-1)*nfree+j) = -1.0d0
 2    continue
      return 
      end

      subroutine cfunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
            coeff((j-1)*nfree+i) = 0.0d0
 1       continue
         coeff((j-1)*nfree+j) = -1.0d0
 2    continue
      return 
      end

      subroutine dfunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
            coeff((j-1)*nfree+i) = 0.0d0
 1       continue
 2    continue
      return 
      end

      subroutine efunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
            coeff((j-1)*nfree+i) = 0.0d0
 1       continue
 2    continue
      return 
      end

      subroutine ffunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
            coeff((j-1)*nfree+i) = 0.0d0
 1       continue
 2    continue
      return 
      end

      subroutine gfunbl (nfree,x,y,z,coeff)
      real*8 x, y, z, coeff(100) 
      do 2 j=1, nfree
         do 1 i=1, nfree
            coeff((j-1)*nfree+i) = 0.0d0
 1       continue
 2    continue
      return 
      end
c-----------------------------------------------------------------------
c     The material property function xyk for the 
c     finite element problem 
c-----------------------------------------------------------------------
      subroutine xyk(nel,xyke,x,y,ijk,node)
      implicit real*8 (a-h,o-z)
      dimension xyke(2,2), x(*), y(*), ijk(node,*)
c     
c     this is the identity matrix.
c     
      xyke(1,1) = 1.0d0
      xyke(2,2) = 1.0d0
      xyke(1,2) = 0.0d0
      xyke(2,1) = 0.0d0

      return
      end
