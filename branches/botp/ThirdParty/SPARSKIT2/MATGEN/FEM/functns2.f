c-----------------------------------------------------------------------
c     contains the functions needed for defining the PDE poroblems. 
c
c     first for the scalar 5-point and 7-point PDE 
c-----------------------------------------------------------------------
      function afun (x,y,z)
      real*8 afun, x,y,z 
      afun = -1.0
      return 
      end

      function bfun (x,y,z)
      real*8 bfun, x,y,z 
      bfun = -1.0
      return 
      end

      function cfun (x,y,z)
      real*8 cfun, x,y,z 
      cfun = -1.0d0
      return 
      end

      function dfun (x,y,z)
      real*8 dfun, x,y,z 
      dfun = 10.d0
      return 
      end

      function efun (x,y,z)
      real*8 efun, x,y,z
      efun = 0.0d0
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
c        subroutine xyk(nel,xyke,x,y,ijk,node)
c        implicit real*8 (a-h,o-z)
c        dimension xyke(2,2), x(*), y(*), ijk(node,*)
cc
cc this is the identity matrix.
cc
c        xyke(1,1) = 1.0d0
c        xyke(2,2) = 1.0d0
c        xyke(1,2) = 0.0d0
c        xyke(2,1) = 0.0d0
cc
c        return
c        end
c
        subroutine xyk (xyke, x, y)
        implicit real*8(a-h,o-z)
        dimension xyke(2,2)

        xyke(1,1) = 1.
        xyke(1,1) = exp(x+y)
        xyke(1,2) = 0.
        xyke(2,1) = 0.
        xyke(2,2) = 1.
        xyke(2,2) = exp(x+y)
        return
        end

        function funb(x,y)
        implicit real*8(a-h,o-z)

        funb = 0.
        funb = 2.5
        funb = 2*x
        return
        end

        function func(x,y)
        implicit real*8(a-h,o-z)

        func = 0.
        func = -1.5
        func = -5*y
        return
        end

        function fung(x,y)
c  Right hand side corresponding to the exact solution of
c   u = exp(x+y)*x*(1.-x)*y*(1.-y)
c   (That exact solution is defined in the function exact)
        implicit real*8(a-h,o-z)

c        fung = 1. + x
c        fung = 2*y*(1.-y) + 2*x*(1.-x)
c        fung = 2*y*(1.-y) + 2*x*(1.-x) +2.5*(1-2*x)*y*(1-y) 
c     1        - 1.5*(1.-2*y)*x*(1-x)
        r = exp(x+y)
        fung = r*r*((x*x+3.*x)*y*(1.-y) +
     1                    (y*y+3.*y)*x*(1.-x)) 
     2          + r*(r-2.*x)*(x*x+x-1.)*y*(1.-y)
     3                + r*(r+5.*y)*(y*y+y-1.)*x*(1.-x)
        return
        end

         function exact(x,y)
C  Exact Solution.
         implicit real*8(a-h,o-z)

         exact = exp(x+y)*x*(1.-x)*y*(1.-y)

         return
         end
