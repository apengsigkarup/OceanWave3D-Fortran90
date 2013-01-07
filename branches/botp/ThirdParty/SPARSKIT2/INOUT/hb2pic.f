        program hb2pic                                               
c------------------------------------------------------------------c
c
c reads a harwell-Boeing matrix and creates a pic file for pattern.
c 
c------------------------------------------------------------------c
        implicit real*8 (a-h,o-z)
        parameter (nmax = 5000, nzmax = 70000)
        integer ia(nmax+1), ja(nzmax)
	real*8  a(nzmax), rhs(1)
        character title*72, key*8, guesol*2 
        logical valued
c--------------
      data iin /5/, iout/6/
c--------------
      job = 2
      nrhs = 0
      call readmt (nmax,nzmax,job,iin,a,ja,ia, rhs, nrhs,
     *     guesol,nrow,ncol,nnz,title,key,type,ierr)
c---- if not readable return 
      if (ierr .ne. 0) then
         write (iout,100) ierr
 100     format(' **ERROR: Unable to read matrix',/,
     *        ' Message returned fom readmt was ierr =',i3)
         stop
      endif
      valued = (job .ge. 2)
c-------
      mode = 1
      iounit = 6
      job = 11
      call pltmt (nrow,ncol,mode,ja,ia,title,key,type,job,iout)
c----------------------------------------------------------------------- 
      stop
      end
