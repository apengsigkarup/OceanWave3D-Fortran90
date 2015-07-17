       program matprod
c-----------------------------------------------------------------------
c      test program for some routines in BLASSM.f
c----------------------------------------------------------------------- 
c      Last update: May 2, 1994
c----------------------------------------------------------------------- 
       implicit real*8 (a-h,o-z)
       parameter (nxmax = 30,nmx = nxmax*nxmax,nzmax=7*nmx)	
       integer ia(nmx+1),ib(nmx+1),ic(nzmax),
     *	       ja(nzmax),jb(nzmax),jc(nzmax),iw(nmx)
c-----------------------------------------------------------------------
       real*8 a(nzmax),b(nzmax),c(nzmax),
     *	      x(nmx),y(nmx),y1(nmx),rhs(nmx), al(6)
       character title*71,key*8,type*3 

       nx = 20
       ny = 20
       nz = 1
       al(1) = 1.0D0
       al(2) = 0.0D0
       al(3) = 2.3D1
       al(4) = 0.4D0
       al(5) = 0.0D0
       al(6) = 8.2D-2
       iout = 8
c-----------------------------------------------------------------------
       call gen57pt (nx,ny,nz,al,0,n,a,ja,ia,iw,rhs)
       call gen57pt (ny,nx,nz,al,0,n,b,jb,ib,iw,rhs)
c
       s = 3.812
c        
       call aplsb1(n,n,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,ierr)
       if (ierr .ne. 0) print *,' ierr = ',ierr
c
c       call dump (1,n,.true.,c,jc,ic,9) 
c
       do 1 k=1,n
          x(k) = real(k)/real(n)
 1     continue
c
       call ope (n,x,y1,a,ja,ia)
       call ope (n,x,y,b,jb,ib)
       do 2 j=1, n
          y1(j) = s*y(j) + y1(j)
 2     continue
c     
       call ope (n,x,y,c,jc,ic)
c------------------------------------------------------
       write (6,*) ' ------------ checking APLSB --------------'
       call ydfnorm(n,y1,y,6)
c------------------------------------------------------
       
       type = '--------'
       title=' test matrix for blassm c = a+b '
       key = 'rua'
c     
       ifmt = 103
c     
       job = -1
c--------
       do 121 jj=1,2
	  write (9,*) 'DUMP A____________________________'
          call dump (1,n,.true.,a,ja,ia,9) 
	  write (9,*) 'DUMP B____________________________'
          call dump (1,n,.true.,b,jb,ib,9) 
          call apmbt(n,n,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
	  write (9,*) 'DUMP C____________________________'
          call dump (1,n,.true.,c,jc,ic,9) 
          if (ierr .ne. 0) print *,' ierr = ',ierr
          call ope (n,x,y1,a,ja,ia)
          call opet (n,x,y,b,jb,ib)
          s = real(job) 
          do 3 j=1, n
 3           y1(j) = y1(j) + s*y(j)
c     
             call ope (n,x,y,c,jc,ic)
c------------xs------------------------------------------
	 write (6,*) '  '
	 write (6,*) ' ------------ checking APMBT---------------'
	 write (6,*) ' ------------ with JOB = ',job,' -------------'
	 call ydfnorm(n,y1,y,6)
c------------------------------------------------------
	 job = job + 2 
 121	 continue
c
	 type = '--------'
	 title=' test matrix for blassm c = a+b^T '
c 
c
c
	s = 0.1232445
	call aplsbt(n,n,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
c
	if (ierr .ne. 0) print *,' ierr = ',ierr
         call ope (n,x,y1,a,ja,ia)
    	 call opet (n,x,y,b,jb,ib)
	 do 4 j=1, n
 4	      y1(j) = y1(j) + s*y(j)
c
 	 call ope (n,x,y,c,jc,ic)
c------------------------------------------------------
c------------------------------------------------------
	 write (6,*) '  '
	 write (6,*) ' ------------ checking APLSBT---------------'
	 call ydfnorm(n,y1,y,6)
c----------------------------------------------------------------------- 
c testing products
c-----------------------------------------------------------------------
        job = 1 
	call amub (n,n,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
c
	if (ierr .ne. 0) print *,' ierr = ',ierr
    	 call ope (n,x,y,b,jb,ib)
    	 call ope (n,y,y1,a,ja,ia)
c
 	 call ope (n,x,y,c,jc,ic)
c----------------------------------------------------------------------- 
	 write (6,*) '  '
	 write (6,*) ' ------------ checking AMUB  ---------------'
	 call ydfnorm(n,y1,y,6)
c
       stop
       end
c
c
    	subroutine ope (n,x,y,a,ja,ia) 
	implicit real*8 (a-h,o-z)
	real*8  x(1),y(1),a(1) 
	integer ia(*),ja(*)
c sparse matrix * vector multiplication
c	
	do 100 i=1,n
		k1 = ia(i)
		k2 = ia(i+1) -1
	y(i) = 0.0
	do 99 k=k1,k2
		y(i) = y(i) + a(k)*x(ja(k)) 
 99	continue
 100	continue
	return
	end
c
    	subroutine opet (n,x,y,a,ja,ia) 
	implicit real*8 (a-h,o-z)
	real*8  x(1),y(1),a(1) 
	integer ia(*),ja(*)
c sparse matrix * vector multiplication
c	
	do 1 j=1, n
 1	y(j) = 0.0d0 
c
	do 100 i=1,n
	do 99 k=ia(i), ia(i+1)-1 
		y(ja(k)) = y(ja(k)) + x(i)*a(k)
 99	continue
 100	continue
	return
	end
c
	subroutine ydfnorm(n,y1,y,iout)
	implicit real*8 (a-h,o-z)
	real*8 y(*),y1(*)
c
	 t = 0.0d0
	 do 21 k=1,n
	 t = t+(y(k)-y1(k))**2
 21	 continue
	 t = sqrt(t) 
       write(iout,*) '2-norm of error (exact answer-tested answer)=',t
c-----------------------------------------------------------------------
	return
	end
