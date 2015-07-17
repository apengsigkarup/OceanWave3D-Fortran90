      	program exptest
c-------------------------------------------------------------------
c
c Test program for exponential propagator using Arnoldi approach
c This main program is a very simple test using diagonal matrices
c (Krylov subspace methods are blind to the structure of the matrix
c except for symmetry). This provides a good way of testing the
c accuracy of the method as well as the error estimates.
c
c-------------------------------------------------------------------
        implicit real*8 (a-h,o-z)
        parameter (nmax = 400, ih0=60, ndmx=20,nzmax = 7*nmax)
        real*8 a(nzmax), u(ih0*nmax), w(nmax),w1(nmax),x(nmax),y(nmax)
	integer ioff(10) 
	data iout/6/, a0/0.0/, b0/1.0/, epsmac/1.d-10/, eps /1.d-10/
c
c set dimension of matrix
c
	n = 100
c--------------------------- define matrix -----------------------------
c A is a single diagonal matrix (ndiag = 1 and ioff(1) = 0 ) 
c-----------------------------------------------------------------------
        ndiag = 1
	ioff(1) = 0
c
c-------- entries in the diagonal are uniformly distributed. 
c
	h = 1.0d0 / real(n+1)
        do 1 j=1, n
	a(j) = real(j+1)* h 
 1	continue
c-------- 
        write (6,'(10hEnter tn: ,$)')
	read (5,*) tn 
c
        write (6,'(36hEpsilon (desired relative accuracy): ,$)')
	read (5,*)  eps 
c------- 
        write (6,'(36h m (= dimension of Krylov subspace): ,$)')
	read (5,*)  m 
c------- 
c define initial conditions: chosen so that solution = (1,1,1,1..1)^T
c-------
	do 2 j=1,n
	w(j) = dexp(a(j)*tn) 
 2	w1(j) = w(j) 
c
         call expprod (n, m, eps, tn, u, w, x, y, a, ioff, ndiag)
c
        print *, ' final answer '
        print *, (w(k),k=1,20) 
c        
	do 4 k=1,n
 4	w1(k) = dexp(-a(k)*tn) * w1(k) 
        print *, ' exact solution '
        print *, (w1(k),k=1,20) 
c
c---------- computing actual 2-norm of error ------------------
c
	t = 0.0d0
	do 47 k=1,n
 47	t = t+ (w1(k)-w(k))**2
        t = dsqrt(t / ddot(n, w,1,w,1) )
c
        write (6,*) ' final error', t
c--------------------------------------------------------------
	stop
	end
c------- 
	subroutine oped(n,x,y,diag,ioff,ndiag)
c======================================================
c this kernel performs a matrix by vector multiplication
c for a diagonally structured  matrix stored in diagonal
c format
c======================================================
	implicit real*8 (a-h,o-z)
	real*8 x(n), y(n), diag(n,ndiag)
	common nope, nmvec
	integer j, n, ioff(ndiag)
CDIR$ IVDEP
	do 1 j=1, n
		y(j) = 0.00
 1	continue	
c
	do 10 j=1,ndiag
	io = ioff(j)
	i1=max0(1,1-io)
	i2=min0(n,n-io)
CDIR$ IVDEP
	do 9 k=i1,i2
	y(k) = y(k)+diag(k,j)*x(k+io)
 9	continue
 10	continue
	nmvec  = nmvec + 1 
	nope = nope + 2*ndiag*n
	return
	end
c
