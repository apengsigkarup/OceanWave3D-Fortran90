      program hb2ps 
c----------------------------------------------------------------------
c translates a harwell - boeing file into a post-script file. Usage:   
c                   hb2ps < HB_file > Postscript_file 
c where hb2ps is the executable generated from this program,
c HB_file is a file containing a matrix stored in Harwell-Boeing
c format and Postscript_file is a file to contain the post-script file. 
c---------------------------------------------------------------------- 
      parameter (nmax = 10000, nzmax = 100000)
      integer ia(nmax+1),ja(nzmax), idummy(1), ptitle
      real*8  a(1),rhs(1)
      real size
      character title*72, key*8, guesol*2, munt*2 
      data iin /5/, iout/6/, size/5.0/, nlines/0/, ptitle/0/,mode/0/
      data munt/'in'/
c-----------------------------------------------------------------------
      job = 1 
      nrhs = 0
c
c     read matrix in Harwell-Boeing format 
c
      call readmt (nmax,nzmax,job,iin,a,ja,ia, rhs, nrhs,
     *     guesol,nrow,ncol,nnz,title,key,type,ierr)
c
c     if   not readable return 
c
      if (ierr .ne. 0) then 
         write (iout,100) ierr
         stop
      endif
c
c     call post script generator
c     
      call pspltm(nrow,ncol,mode,ja,ia,title,ptitle,size,munt,
     *     nlines,idummy,iout) 
c
 100     format(' **ERROR: Unable to read matrix',/,
     *        ' Message returned fom readmt was ierr =',i3)
c-----------------------------------------------------------------------  
      stop
      end

