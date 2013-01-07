      subroutine texgrd(npts,ja,ia,xx,yy,munt,size,vsize,hsize,
     *     xleft,bot,job,title,ptitle,ijk,node,nelx,iunt)
c-----------------------------------------------------------------------
      integer npts,iunt,ptitle,ja(*),ia(*), ijk(node,*) 
      character title*(*), munt*2 
      real*8 xx(npts), yy(npts) 
c-----------------------------------------------------------------------
c     allows to have several grids in same picture by calling texgrd 
c     several times and exploiting job and different shifts.  
c-----------------------------------------------------------------------
c     input arguments description :
c     
c npts   = number of rows in matrix
c
c ncol   = number of columns in matrix 
c
c mode   = integer indicating whether the matrix is stored in 
c          CSR mode (mode=0) or CSC mode (mode=1) or MSR mode (mode=2) 
c
c ja     = column indices of nonzero elements when matrix is
c          stored rowise. Row indices if stores column-wise.
c ia     = integer array of containing the pointers to the 
c          beginning of the columns in arrays a, ja.
c
c munt   = units used for sizes : either 'cm' or 'in'
c
c size   = size of the matrix box in 'munt' units 
c
c vsize  = vertical size of the frame containing the picture
c          in 'munt' units 
c
c hsize  = horizontal size of the frame containing the picture 
c          in 'munt' units 
c
c xleft  =  position of left border of matrix in 'munt' units 
c
c bot    =  position of bottom border of matrix in 'munt' units 
c
c job    = job indicator for preamble and post process
c          can be viewed as a 2-digit number job = [job1,job2]      
c          where job1 = job /10 , job2 = job - 10*job1 = mod(job,10)
c          job2  relates to preamble/post processing:
c          job2 = 0: all preambles+end-document lines 
c          job2 = 1: preamble only
c          job2 = 2: end-document only
c          anything else:  no preamble or end-docuiment lines
c       Useful for plotting several matrices in same frame.
c
c          job1 relates to the way in which the nodes and elements must
c          be processed:
c          job1  relates to options for the plot. 
c          job1 = 0 : only a filled circle for the nodes, no labeling 
c          job1 = 1 : labels the nodes (in a circle) 
c          job1 = 2 : labels both nodes and elements. 
c          job1 = 3 : no circles, no labels 
c
c title  = character*(*). a title of arbitrary length to be printed 
c          as a caption to the matrix. Can be a blank character if no
c          caption is desired. Can be put on top or bottom of matrix
c          se ptitle.
c
c ptitle = position of title; 0 under the frame, else above
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
      real*8 xmin, xmax, ymin, ymax

      n = npts 
      siz = size
      job1 = job /10 
      job2 = job - 10*job1 
c
c     get max and min dimensions
c
      xmin = xx(1) 
      xmax = xmin 
      ymin = yy(1)
      ymax = ymin
c
      do j=2, npts
         xmax = max(xmax,xx(j))
         xmin = min(xmin,xx(j))
         ymax = max(ymax,yy(j))
         ymin = min(ymin,yy(j))
      enddo
c----------------------------------------------------------------------- 
      n = npts
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
      tdim = max(ydim,xdim) 
      unit0 = size/tdim 
      hsiz = hsize/unit0
      vsiz = vsize/unit0
      siz  = size/unit0 
c
c     size of little circle for each node -- cirr = in local units
c     cirabs in inches (or cm) -- rad = radius in locl units --
c 
      cirabs = 0.15
      if (job1 .le. 0) cirabs = 0.08
      if (job1 .eq. 3) cirabs = 0.0 
      cirr   = cirabs/unit0 
      rad = cirr/2.0
c
c     begin document generation
c 
      if (job2 .le. 1) then 
         write (iunt,*) ' \\documentstyle[epic,eepic,12pt]{article} '
         write (iunt,*) ' \\begin{document}'
         write (iunt,100) unit0, munt 
         write (iunt,99) hsiz, vsiz 
      else 
c     redeclare unitlength 
         write (iunt,100) unit0, munt 
      endif 
      if (job1 .le. 0) then 
         write (iunt,101) cirr
      else
         write (iunt,102) cirr
      endif 
 102  format('\\def\\cird{\\circle{',f5.2,'} }') 
 101   format('\\def\\cird{\\circle*{',f5.2,'} }') 
 99   format('\\begin{picture}(',f8.2,1h,,f8.2,1h)) 
 100  format(' \\setlength{\\unitlength}{',f5.3,a2,'}')
c
c     bottom margin between cir and title
c
      xs = xleft /  unit0 + (tdim - xdim)*0.5 + (hsiz - siz)*0.5 
      ys = bot / unit0 - ymin
c
      xmargin = 0.30/unit0
      if (munt .eq. 'cm' .or. munt .eq. 'CM') xmargin = xmargin*2.5
      xtit = xs + xmin 
      if (ptitle .eq. 0) then 
         ytit = ys + ymin - xmargin 
      else
         ytit = ys + ymax + xmargin 
      endif 
      ltit = LENSTR(title)
      write(iunt,111) xtit,ytit,xdim,xmargin,title(1:ltit) 
c      
 111  format ('\\put(',F6.2,',',F6.2,
     *     '){\\makebox(',F6.2,1h,,F6.2,'){',A,2h}}) 
c      
c     print all the circles if needed
c 
c ##### temporary for showing f
c      write (iunt,102) 0.01      
c      write (iunt,112)xs+xx(1), ys+yy(1) 
c-----------------------------------------------------------------------
      if (job1 .eq. 3) goto 230
      do 22 i=1, npts 
         x = xs + xx(i) 
         y = ys + yy(i) 
         write(iunt,112) x, y 
 112     format ('\\put(',F6.2,',',F6.2,'){\\cird}')
c         write(iunt,113) x-rad, y-rad, cirr, cirr, i 
c 113     format ('\\put(',F6.2,',',F6.2,
c     *     '){\\makebox(',F6.2,1h,,F6.2,'){\\scriptsize ',i4,2h}}) 
         if (job1 .ge. 1) then 
            write(iunt,113) x, y, i 
         endif
 113     format ('\\put(',F6.2,',',F6.2,
     *        '){\\makebox(0.0,0.0){\\scriptsize ',i4,2h}}) 
 22   continue
 230  continue 
c         
c     number the elements if needed
c
      if (job1 .eq. 2) then
         do 23 iel = 1, nelx 
            x = 0.0
            y = 0.0
            do j=1, node
               x = x+xx(ijk(j,iel)) 
               y = y+yy(ijk(j,iel)) 
            enddo
            x = xs + x / real(node) 
            y = ys + y / real(node) 
            write(iunt,113) x, y, iel 
 23      continue
      endif 
c
c     draw lines
c
      write (iunt,*) ' \\Thicklines '
      do 1 ii=1, npts 
         xi = xs+ xx(ii) 
         yi = ys+ yy(ii) 
         do 2 k=ia(ii),ia(ii+1)-1 
            j = ja(k) 
            if (j .le. ii) goto 2
            xj = xs + xx(j) 
            yj = ys + yy(j) 
            xspan = xj - xi
            yspan = yj - yi
            tlen = sqrt(xspan**2 + yspan**2) 
            if (abs(xspan) .gt. abs(yspan)) then
               ss = yspan / tlen
               cc = sqrt(abs(1.0 - ss**2))
               cc = sign(cc,xspan) 
            else
               cc = xspan / tlen 
               ss = sqrt(abs(1.0 - cc**2))
               ss = sign(ss,yspan) 
            endif
c            print *, ' ss -- cc ', ss, cc 
c            write(iunt,114)xi,yi,xj,yj
            write(iunt,114)xi+cc*rad,yi+ss*rad,xj-cc*rad,yj-ss*rad 
 114        format('\\drawline(',f6.2,1h,,f6.2,')(',f6.2,1h,,f6.2,1h))
c            tlen = tlen - 2.0*cirr
c            write(iunt,114) xi+cc*cirr,yi+ss*cirr,cc,ss,tlen 
c 114        format('\\put(',f6.2,1h,,f6.2,'){\\line(',
c     *           f6.2,1h,,f6.2,'){',f6.2,'}}') 
 2       continue 
 1    continue
c-----------------------------------------------------------------------
      if (job2 .eq. 0 .or. job2 .eq. 2) then 
         write (iunt,*) ' \\end{picture} '
         write (iunt,*) ' \\end{document} '
      endif
c
      return
      end
