      program rsobel
      integer n, ia(1:200), ja(1:1000), ib(1:200), jb(1:1000)
      integer nrowc, ncolc
      integer ic(1:200), jc(1:1000), ierr
      real*8  a(1:1000), b(1:1000), c(1:1000)

      write (*, '(1x, 9hInput n:  ,$)')
      read *, n
      call sobel(n,nrowc,ncolc,c,jc,ic,a,ja,ia,b,jb,ib,1000,ierr)
      print *, 'ierr =', ierr
      print *, 'Nrow =', nrowc, '	Ncol =', ncolc
      call dump(1, nrowc, .true., c, jc, ic, 6)
      end

