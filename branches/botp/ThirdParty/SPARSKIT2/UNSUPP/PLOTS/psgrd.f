      subroutine psgrid (npts,ja,ia,xx,yy,title,ptitle,size,munt,iunt) 
c-----------------------------------------------------------------------
c     plots a symmetric graph defined by ja,ia and the coordinates
c     xx(*),yy(*) 
c----------------------------------------------------------------------- 
c npts   = number of points in mesh
c ja, ia = connectivity of mesh -- as given by pattern of sparse
c            matrix.
c xx, yy = cordinates of the points. 
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
c     iunt   = logical unit number where to write the matrix into.
c----------------------------------------------------------------------- 
      implicit none 
      integer npts,ptitle,iunt, ja(*), ia(*) 
      character title*(*), munt*2 
      real*8 xx(npts),yy(npts) 
      real size
c----------------------------------------------------------------------- 
c     local variables --------------------------------------------------
c----------------------------------------------------------------------- 
      integer nr,ii,k,ltit
      real xi,yi,xj,yj,lrmrgn,botmrgn,xtit,ytit,ytitof,fnstit,siz,
     *     xl,xr, yb,yt, scfct,u2dot,frlw,delt,paperx,conv,dimen,
     *     haf,zero, xdim, ydim
      real*8 xmin,xmax,ymin,ymax,max,min
      integer j
c-----------------------------------------------------------------------
      integer LENSTR
      external LENSTR
c----------------------------------------------------------------------- 
      data haf /0.5/, zero/0.0/, conv/2.54/
      siz = size
c
c     get max and min dimensions
c
      xmin = xx(1) 
      xmax = xmin 
      ymin = yy(1)
      ymax = ymin
      do j=2, npts
         xmax = max(xmax,xx(j))
         xmin = min(xmin,xx(j))
         ymax = max(ymax,yy(j))
         ymin = min(ymin,yy(j))
      enddo
c----------------------------------------------------------------------- 
c      n = npts
      nr = npts 
      xdim = xmax -xmin
      ydim = ymax -ymin
      dimen = max(xdim,ydim) 
c-----------------------------------------------------------------------
      print *, ' xmin', xmin, ' xmax', xmax 
      print *, ' ymin', ymin, ' ymax', ymax, ' dimen ', dimen
c-----------------------------------------------------------------------
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
      botmrgn = 2.0/dimen
c     scaling factor
      scfct = siz*u2dot/dimen
c     frame line witdh
      frlw = 0.25/dimen
c     font siz for title (cm)
      fnstit = 0.5/dimen
      ltit = LENSTR(title)
c
c     position of title : centered horizontally
c                     at 1.0 cm vertically over the drawing
      ytitof = 1.0/dimen
      xtit = paperx/2.0
      ytit = botmrgn+siz*nr/dimen + ytitof
c almost exact bounding box
      xl = lrmrgn*u2dot - scfct*frlw/2
      xr = (lrmrgn+siz)*u2dot + scfct*frlw/2
      yb = botmrgn*u2dot - scfct*frlw/2
      yt = (botmrgn+siz*ydim/dimen)*u2dot + scfct*frlw/2
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
c
      if (ptitle.eq.0 .and. ltit.gt.0) then
      ytit = botmrgn + fnstit*0.3
      botmrgn = botmrgn + ytitof + fnstit*0.7
      end if
c
c     begin output
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
c     we leave margins etc. in cm so it is easy to modify them if
c     needed by editing the output file
c
      write(iunt,10) 'gsave'
      if (ltit.gt.0) then
      write(iunt,*) '/Helvetica findfont',fnstit,' cm scalefont setfont'
      write(iunt,*) xtit,' cm',ytit,' cm moveto'
      write(iunt,'(3A)') '(',title(1:ltit),') Cshow'
      end if
c
      write(iunt,*) lrmrgn,' cm ',botmrgn,' cm translate'
      write(iunt,*) siz,' cm ',dimen,' div dup scale '
c     
c     draw a frame around the matrix  // REMOVED 
c
c       del = 0.005
c       del2 = del*2.0 
c       write(iunt,*) del, ' setlinewidth'
c       write(iunt,10) 'newpath'
c       write(iunt,11) -del2, -del2, ' moveto'
c       write(iunt,11) dimen+del2,-del2,' lineto'
c       write(iunt,11) dimen+del2, dimen+del2, ' lineto'
c       write(iunt,11) -del2,dimen+del2,' lineto'
c       write(iunt,10) 'closepath stroke'
c     
c----------- plotting loop ---------------------------------------------
       write(iunt,*)  ' 0.01 setlinewidth'
c
       do 1 ii=1, npts 
       
c     if (mask(ii) .eq. 0) goto 1
          xi = xx(ii) - xmin 
          yi = yy(ii) - ymin
c     write (iout+1,*) ' ******** ii pt', xi, yi, xmin, ymin 
          do 2 k=ia(ii),ia(ii+1)-1 
             j = ja(k) 
             if (j .le. ii) goto 2 
             xj = xx(j) - xmin 
             yj = yy(j) - ymin
c     write (iout+1,*) ' j pt -- j= ',j, 'pt=', xj, yj
c
c     draw a line from ii to j
c
             write(iunt,11) xi, yi, ' moveto '
             write(iunt,11) xj, yj, ' lineto'
             write(iunt,10) 'closepath stroke'
 2        continue 
 1     continue
c-----------------------------------------------------------------------
       write(iunt,10) 'showpage'
       return
c
 10   format (A)
 11   format (2F9.2,A)
 12   format (A,4F9.2)
 13   format (2F9.2,A)
c-----------------------------------------------------------------------
      end
