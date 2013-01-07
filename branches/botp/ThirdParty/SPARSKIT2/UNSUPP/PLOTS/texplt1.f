      subroutine texplt(nrow,ncol,mode,ja,ia,munt,size,vsize,hsize,
     *     xleft,bot,job,title,ptitle,nlines,lines,iunt)
c-----------------------------------------------------------------------
      integer nrow,ncol,mode,iunt,ptitle,ja(*),ia(*),lines(nlines) 
      character title*(*), munt*2 
c-----------------------------------------------------------------------
c     allows to have several matrices in same picture by calling texplt
c     several times and exploiting job and different shifts.  
c-----------------------------------------------------------------------
c     input arguments description :
c     
c nrow   = number of rows in matrix
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
c          
c          can be thought of as a 2-digit number job = [job1,job2]
c          where job1 = job /10 , job2 = job - 10*job1 = mod(job,10)
c          job2 = 0: all preambles+end-document lines 
c          job2 = 1: preamble only
c          job2 = 2: end-document only
c          anything else for job2:  no preamble or end-docuiment lines
c       Useful for plotting several matrices in same frame.
c
c          job1 indicates what to put for a nonzero dot.
c          job1  relates to preamble/post processing:
c          job1 = 0 : a filled squate 
c          job1 = 1 : a filled circle 
c          job1 = 2 : the message $a_{ij}$ where i,j are the trow/column
c                     positions of the nonzero element. 
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
      data haf /0.5/, zero/0.0/, conv/2.54/
c-----------------------------------------------------------------------
      n = ncol 
      if (mode .eq. 0) n = nrow 
      job1 = job /10 
      job2 = job - 10*job1 
      maxdim = max(nrow, ncol)
      rwid = real(ncol-1)
      rht  = real(nrow-1)
      unit0 = size/real(maxdim) 
      hsiz = hsize/unit0
      vsiz = vsize/unit0
      siz  = size/unit0
c     
c     size of little box for each dot -- boxr = in local units
c     boxabs in inches (or cm)
c 
      boxr   = 0.6
      boxabs = unit0*boxr 
c
c     spaces between frame to nearest box
c
      space = 0.03/unit0+(1.0-boxr)/2.0 
c
c     begin document generation
c  for very first call better have \unitlength set first..
      if (job2 .le. 1) then 
         write (iunt,*) ' \\documentstyle[epic,12pt]{article} '
         write (iunt,*) ' \\begin{document}'
         write (iunt,100) unit0, munt 
         write (iunt,99) hsiz, vsiz 
      else 
c     redeclare unitlength 
         write (iunt,100) unit0, munt 
      endif 
c----- always redefine units 

      if (job1 .eq. 0) then
         write (iunt,101) boxabs, boxabs
      else 
         write (iunt,102) boxabs/unit0
      endif 
 100  format(' \\setlength{\\unitlength}{',f5.3,a2,'}') 
 99   format('\\begin{picture}(',f8.2,1h,,f8.2,1h)) 
 101  format('\\def\\boxd{\\vrule height',f7.4,'in width',f7.4,'in }') 
 102  format('\\def\\boxd{\\circle*{',f7.4,'}}')
c
c     draw a frame around the matrix
c     get shifts from real inches to local units
c
      xs = xleft/unit0 + (hsiz-siz)*0.5
      ys = bot/unit0
c     
      eps = 0.0 
      xmin = xs 
      xmax = xs +rwid + boxr + 2.0*space 
      ymin = ys 
      ymax = ys+rht + boxr + 2.0*space
c     
c     bottom margin between box and title
c
      xmargin = 0.30/unit0
      if (munt .eq. 'cm' .or. munt .eq. 'CM') xmargin = xmargin*2.5
      xtit = 0.5*(xmin+xmax) 
      xtit = xmin
      ytit = ymax 
      if (ptitle .eq. 0) ytit = ymin - xmargin 
      xdim = xmax-xmin
      ltit = LENSTR(title)
      write(iunt,111) xtit,ytit,xdim,xmargin,title(1:ltit) 
c      
 111  format ('\\put(',F6.2,',',F6.2,
     *     '){\\makebox(',F6.2,1h,F6.2,'){',A,2h}}) 
c      
      write(iunt,*)  ' \\thicklines'
      write (iunt,108) xmin,ymin,xmax,ymin,xmax,ymax,
     *     xmin,ymax,xmin,ymin
 108  format('\\drawline',1h(,f8.2,1h,,f8.2,1h), 
     *     1h(,f8.2,1h,,f8.2,1h), 1h(,f8.2,1h,,f8.2,1h), 
     *     1h(,f8.2,1h,,f8.2,1h), 1h(,f8.2,1h,,f8.2,1h)) 
c     
c     draw the separation lines 
c     
c      if (job1 .gt.0) then
c         xs = xs + 0.25
c         ys = ys + 0.25
c      endif
      write(iunt,*)  ' \\thinlines'
      do 22 kol=1, nlines 
         isep = lines(kol)
c     
c     horizontal lines 
c     
         yy =  ys + real(nrow-isep)  
         write(iunt,109) xmin, yy, xmax, yy
c     
c     vertical lines 
c     
         xx = xs+real(isep) 
         write(iunt,109) xx, ymin, xx, ymax 
 22   continue
c     
 109  format('\\drawline',
     *     1h(,f8.2,1h,,f8.2,1h), 1h(,f8.2,1h,,f8.2,1h)) 
      
c-----------plotting loop ---------------------------------------------
c     
c     add some space right of the frame and up from the bottom
c     
      xs = xs+space
      ys = ys+space
c-----------------------------------------------------------------------
      do 1 ii=1, n
         istart = ia(ii)
         ilast  = ia(ii+1)-1 
         if (mode .eq. 1) then
            do 2 k=istart, ilast
               if (job1 .le. 1) then               
                  write(iunt,12) xs+real(ii-1),ys+real(nrow-ja(k)) 
               else
                  write(iunt,13) xs+real(ii-1),ys+real(nrow-ja(k)),
     *                 ii,ja(k) 
               endif
 2          continue 
         else
            y = ys+real(nrow-ii)
            do 3 k=istart, ilast
               if (job1 .le. 1) then
                  write(iunt,12) xs+real(ja(k)-1), y 
               else 
                  write(iunt,13) xs+real(ja(k)-1), y, ii, ja(k) 
               endif
 3          continue          
c     add diagonal element if MSR mode.
            if (mode .eq. 2) 
     *           write(iunt,12) xs+real(ii-1), ys+real(nrow-ii) 
         endif
 1    continue
c-----------------------------------------------------------------------
 12   format ('\\put(',F6.2,',',F6.2,')','{\\boxd}') 
 13   format ('\\put(',F6.2,',',F6.2,')','{$a_{',i3,1h,,i3,'}$}')
c-----------------------------------------------------------------------
      if (job2 .eq. 0 .or. job2 .eq. 2) then 
         write (iunt,*) ' \\end{picture} '
         write (iunt,*) ' \\end{document} '
      endif
c
      return
      end
      integer function lenstr0(s)
c-----------------------------------------------------------------------
c return length of the string S
c-----------------------------------------------------------------------
      character*(*) s
      integer len
      intrinsic len
      integer n
c----------------------------------------------------------------------- 
      n = len(s)
10    continue
        if (s(n:n).eq.' ') then
          n = n-1
          if (n.gt.0) go to 10
        end if
      lenstr0 = n
c
      return
c-----------------------------------------------------------------------
      end
c----------------------------------------------------------------------- 
