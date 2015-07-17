      subroutine phiprod (n, m, eps, tn, u, w, r, x, y, a, ioff, ndiag)
      real*8 eps, tn 
      real*8 a(n,ndiag), u(n,m+1), w(n), r(n), x(n), y(n)
      integer n, m, ndiag, ioff(ndiag)

c-----------------------------------------------------------------------
c this subroutine computes an approximation to the vector
c
c   	w(tn) = w(t0) + tn *  phi( - A * tn ) * (r - A w(t0))
c
c where phi(z) = (1-exp(z)) / z
c 
c i.e. solves   dw/dt = - A w + r in [t0,t0+ tn] (returns only w(t0+tn))
c 
c for matrices stored in diagonal (DIA) format.
c
c this routine constitutes an interface for the routine phipro for
c matrices stored in diagonal (DIA) format. The phipro routine uses
c reverse communication and as a result does not depend on any
c data structure of the matrix.

c-----------------------------------------------------------------------
c ARGUMENTS
c---------- 
c see phipro for meaning of parameters n, m, eps, tn, u, w, x, y.
c 
c a, ioff, and ndiag are the arguments of the matrix:
c
c a(n,ndiag) = a rectangular array with a(*,k) containing the diagonal 
c              offset by ioff(k) (negative or positive or zero), i.e., 
c              a(i,jdiag) contains the element A(i,i+ioff(jdiag)) in 
c              the usual dense storage scheme.
c 
c ioff	     = integer array containing the offsets  of the ndiag diagonals
c ndiag      = integer. the number of diagonals.
c 
c-----------------------------------------------------------------------
c local variables 
c 
      integer indic, ierr

      indic = 0
 101  continue
      call phipro (n, m, eps, tn, w, r, u, x, y, indic, ierr)
      if (indic .eq. 1) goto 102
c     
c     matrix vector-product for diagonal storage --
c     
      call oped(n, x, y, a, ioff, ndiag)
      goto 101
 102  continue
      return
      end
c----------end-of-phiprod----------------------------------------------- 
c-----------------------------------------------------------------------
      subroutine phipro (n, m, eps, tn, w, r, u, x, y, indic, ierr)
c     implicit  real*8 (a-h,o-z)
      integer n, m, indic, ierr
      real*8 eps, tn, w(n), r(n), u(n,m+1), x(n), y(n) 
c-----------------------------------------------------------------------
c
c this subroutine computes an approximation to the vector
c
c     w(tn) = w(t0) + tn *  phi( - A * tn ) * (r - A w(t0))
c     where phi(z) = (1-exp(z)) / z
c 
c     i.e. solves dw/dt=-Aw+r in [t0,t0+tn] (returns w(t0+tn))
c     t0 need not be known.
c 
c note that for w(t0)=0 the answer is    w=tn *phi(-tn * A) r
c in other words this allows to compute phi(A tn) v. 
c This code will work well only for cases where eigenvalues are
c real (or nearly real) and positive. It has also been coded to 
c work for cases where tn .lt. 0.0 (and A has real negative spectrum)
c 
c----------------------------------------------------------------------- 
c
c THIS IS A REVERSE COMMUNICATION IMPLEMENTATION. 
c------------------------------------------------- 
c USAGE: (see also comments on argument indic below).
c------ 
c
c      indic = 0
c 1    continue
c      call phipro (n, m, eps, tn, u, w, x, y, indic)
c      if (indic .eq. 1) goto 2 <-- indic .eq. 1 means phipro has finished
c      call matvec(n, x, y)     <--- user's matrix-vec. product
c                                    with x = input vector, and
c                                     y = result = A * x.
c      goto 1
c 2    continue
c      .....
c
c-----------------------------------------------------------------------
c
c en entry:
c---------- 
c n	= dimension of matrix
c
c m	= dimension of Krylov subspace (= degree of polynomial 
c         approximation to the exponential used. )
c
c eps   = scalar indicating the relative error tolerated for the result. 
c         the code will try to compute an answer such that 
c         norm2(exactanswer-approximation) / norm2(w) .le. eps 
c
c tn	= scalar by which to multiply matrix. (may be .lt. 0)
c         the code will compute a solution to dw/dt = -A w + r,
c         and overwrite the result w(tn) onto in w.
c
c w	= real array of length n. Initial condition for the ODE system
c         on input, result w(tn) on output (input and output argument) 
c
c r     = real array of length n. the constant term in the system 
c         dw/dt = -A w + r to be solved.
c
c u	= work array of size n*(m+1) (used to hold the Arnoldi basis )
c
c x, y  = two real work vectors of length n each. x and y are used to
c         carry the input and output vectors for the matrix-vector
c         products y=Ax in the reverse communication protocole.
c         see argument indic (return) below for details on their usage.
c
c indic = integer used as indicator for the reverse communication.
c         in the first call enter indic = 0.
c
c ierr  = error indicator. 
c         ierr = 1 means phipro was called with indic=1 (not allowed)
c         ierr = -1 means that the input is zero the solution has been 
c         unchanged.
c 
c on return:
c-----------
c 
c w     = contains the result w(tn)=w(t0)+tn*phi(-A*tn)*(r-Aw(t0))
c         when phipro has finished (as indicated by indic see below)
c
c indic = indicator for the reverse communication protocole.
c       * INDIC .eq. 1  means that phipro has finished and w contains the
c         result. 
c       * INDIC .gt. 1 means that phipro has not finished and that 
c         it is requesting another matrix vector product before
c         continuing. The user must compute Ax where A is the matrix
c         and x is the vector provided by phipro and return the 
c         result in y. Then phipro must be called again without
c         changing any other argument. typically this is best 
c         implemented in a loop with phipro being called as long
c         indic is returned with a value .ne. 1.
c 
c NOTES:  m should not exceed 60 in this version  (see mmax below)
c-----------------------------------------------------------------------
c local variables 
c 
      integer mmax
      parameter (mmax=60) 
      real*8 errst, tcur, told, dtl, beta, red, dabs, dble
      real*8 hh(mmax+2,mmax+1), z(mmax+1)
      complex*16   wkc(mmax+1) 
      integer ih, k, job
      logical verboz
      data verboz/.true./
      save
c-----------------------------------------------------------------------
c indic = 4  means  getting y=Ax needed in phipro 
c indic = 3  means  passing through only with result of y= Ax to phihes
c indic = 2  means phihes has finished its job
c indic = 1  means phipro has finished its job (real end)/
c-----------------------------------------------------------------------
      ierr = 0 
      if (indic .eq. 3) goto 101 
      if (indic .eq. 4) goto 11
      if (indic .eq. 1) then 
         ierr = 1
         return
      endif
c----- 
      ih = mmax 
      m  = min0(m,mmax) 
      tcur = 0.0d0
      dtl = tn - tcur
      job = -1
c-------------------- outer loop ----------------------------- 
 100  continue
      if(verboz) print *,'In PHIPRO, current time = ', tcur ,'---------'
c------------------------------------------------------------- 
c ---- call phionential propagator --------------------------- 
c------------------------------------------------------------- 
      told = tcur 
c
c     if (told + dtl .gt. tn) dtl = tn-told
c      construct initial vector for Arnoldi:  r - A w(old) 
c
      do 10 k=1, n
         x(k) = w(k)
 10   continue
      indic = 4
      return
 11   continue
      do 12 k=1, n
         u(k,1) = r(k) - y(k)
 12   continue
c
 101  continue
      call phihes (n,m,dtl,eps,u,job,z,wkc,beta,errst,hh,ih,x, y,indic,
     *            ierr) 
c-----------------------------------------------------------------------
      if (ierr .ne. 0) return
      if (indic .eq. 3) return
      tcur = told + dtl 
      if(verboz) print *, ' tcur now = ', tcur, ' dtl = ', dtl
c
c     relative error 
c      if(verboz) print *, ' beta', beta
      errst = errst / beta
c---------
      if ((errst .le. eps) .and. ( (errst .gt. eps/100.0) .or.
     *     (tcur .eq. tn))) goto 102
c
c     use approximation :  [ new err ] = fact**m  * [cur. error]
c 
      red =  (0.5*eps / errst)**(1.0d0 /dble(m) ) 
      dtl = dtl*red
      if (dabs(told+dtl) .gt. dabs(tn) )  dtl = tn-told
      if(verboz) print *, ' red =',red,' , reducing dt to: ', dtl
c-------
      job = 1 
      goto 101
c-------
 102  continue 
c
      call project(n, m, w, dtl, u, z) 
c never go beyond tcur
      job = 0
      dtl = dmin1(dtl, tn-tcur)
      if (dabs(tcur+dtl) .gt. dabs(tn)) dtl = tn-tcur 
      if (dabs(tcur) .lt. dabs(tn)) goto 100
      indic = 1
      return
      end
c----------end-of-phipro------------------------------------------------
c-----------------------------------------------------------------------
      subroutine phihes (n,m,dt,eps,u,job,z,wkc,beta,errst,hh,ih,
     *                   x, y, indic,ierr) 
c     implicit  real*8 (a-h,o-z)
      integer n, m, job, ih, indic, ierr
      real*8 hh(ih+2,m+1), u(n,m+1), z(m+1), x(n), y(n)
      complex*16 wkc(m+1) 
      real*8 dt, eps, beta, errst
c-----------------------------------------------------------------------
c this subroutine computes the Arnoldi basis Vm and the corresponding 
c coeffcient vector ym in the approximation 
c 
c        	w  ::= beta  Vm  ym 
c               where ym = phi(- Hm * dt) * e1
c
c to the vector phi(-A * dt) w where A is an arbitary matrix and 
c w is a given input vector. The phi function is defined   by
c               phi(z) = (1 - exp(z) ) / z  
c 
c In case job .lt.0 the arnoldi basis is recomputed. Otherwise the
c code assumes assumes that  u(*) contains an already computed 
c arnoldi basis and computes only the y-vector (which is stored in
c v(*)). Three different options are available through the argument job. 
c-----------------------------------------------------------------------
c on entry:
c---------- 
c n	= dimension of matrix
c
c m	= dimension of Krylov subspace (= degree of polynomial 
c         approximation to the phionential used. )
c
c dt	= scalar by which to multiply matrix. Can be viewed
c         as a time step. dt must be positive [to be fixed].
c
c eps   = scalar indicating the relative error tolerated for the result. 
c         the code will try to compute an answer such that 
c         norm2(exactanswer-approximation) / norm2(w) .le. eps 
c
c u	= work array of size n*(m+1) to contain the Arnoldi basis
c
c w	= real array of length n = input vector to  which phi(-A) is
c         to be applied. 
c
c job	= integer. job indicator. If job .lt.  0 then the Arnoldi
c         basis is recomputed. If job .gt. 0 then it is assumed
c         that the user wants to use a previously computed Krylov
c         subspace but a different dt. Thus the Arnoldi basis and
c         the Hessenberg matrix Hm are not recomputed. 
c	  In that case the user should not modify the values of beta
c         and the matrices hh and u(n,*) when recalling phipro. 
c         job = -1 : recompute basis and get an initial estimate for 
c                    time step dt to be used.
c         job = 0  : recompute basis and do not alter dt.
c         job = 1  : do not recompute arnoldi basis. 
c
c z     = real work array of  size (m+1)
c wkc   = complex*16 work array of size (m+1) 
c
c hh    = work array of size size at least (m+2)*(m+1)
c
c ih+2	= first dimension of hh as declared in the calling program.
c         ih must be .ge. m.
c 
c-----------------------------------------------------------------------
c on return:
c-----------
c w2	= resulting vector w2 = phi(-A *dt) * w
c beta  = real equal to the 2-norm of w. Needed if phipro will
c         be recalled with the same Krylov subspace and a different dt. 
c errst = rough estimates of the 2-norm of the error.
c hh	= work array of dimension at least (m+2) x (m+1) 
c
c----------------------------------------------------------------------- 
c local variables 
c 
      integer ndmax
      parameter (ndmax=20) 
      real*8 alp0, fnorm, t, rm, ddot, dabs, dsqrt, dsign,dble
      complex*16 alp(ndmax+1), rd(ndmax+1)
      integer i, j, k, ldg, i0, i1, m1 
      logical verboz
      data verboz/.true./
      save
c------use degree 14 chebyshev all the time -------------------------- 
      if (indic .eq. 3) goto 60
c
c------get partial fraction expansion of rational function -----------
c-----------------------------------------------------------------------
c chebyshev (14,14) 
c      ldg= 7
c      alp0 =  0.183216998528140087E-11
c      alp(1)=( 0.557503973136501826E+02,-0.204295038779771857E+03)
c      rd(1)=(-0.562314417475317895E+01, 0.119406921611247440E+01)
c      alp(2)=(-0.938666838877006739E+02, 0.912874896775456363E+02)
c      rd(2)=(-0.508934679728216110E+01, 0.358882439228376881E+01)
c      alp(3)=( 0.469965415550370835E+02,-0.116167609985818103E+02)
c      rd(3)=(-0.399337136365302569E+01, 0.600483209099604664E+01)
c      alp(4)=(-0.961424200626061065E+01,-0.264195613880262669E+01)
c      rd(4)=(-0.226978543095856366E+01, 0.846173881758693369E+01)
c      alp(5)=( 0.752722063978321642E+00, 0.670367365566377770E+00)
c      rd(5)=( 0.208756929753827868E+00, 0.109912615662209418E+02)
c      alp(6)=(-0.188781253158648576E-01,-0.343696176445802414E-01)
c      rd(6)=( 0.370327340957595652E+01, 0.136563731924991884E+02)
c      alp(7)=( 0.143086431411801849E-03, 0.287221133228814096E-03)
c      rd(7)=( 0.889777151877331107E+01, 0.166309842834712071E+02)
c-----------------------------------------------------------------------
c Pade of  degree =  (4,4) 
c
c        ldg= 2
c        alp(1)=(-0.132639894655051648E+03,-0.346517448171383875E+03)
c        rd(1)=(-0.579242120564063611E+01, 0.173446825786912484E+01)
c        alp(2)=( 0.926398946550511936E+02, 0.337809095284865179E+02)
c        rd(2)=(-0.420757879435933546E+01, 0.531483608371348736E+01)
c
c Pade of degree =  8
c
        ldg= 4
        alp(1)=( 0.293453004361944040E+05, 0.261671093076405813E+05)
        rd(1)=(-0.104096815812822569E+02, 0.523235030527069966E+01)
        alp(2)=(-0.212876889060526154E+05,-0.764943398790569044E+05)
        rd(2)=(-0.111757720865218743E+02, 0.173522889073929320E+01)
        alp(3)=(-0.853199767523084301E+04,-0.439758928252937039E+03)
        rd(3)=(-0.873657843439934822E+01, 0.882888500094418304E+01)
        alp(4)=( 0.330386145089576530E+03,-0.438315990671386316E+03)
        rd(4)=(-0.567796789779646360E+01, 0.127078225972105656E+02)
c
      do 102 k=1, ldg
         alp(k) = - alp(k) / rd(k)
 102  continue
      alp0 = 0.0d0
c     
c     if job .gt. 0 skip arnoldi process:
c     
      if (job .gt. 0) goto 2
c------normalize vector u and put in first column of u -- 
      beta = dsqrt(ddot(n,u,1,u,1))
c-----------------------------------------------------------------------
      if(verboz) print *, ' In PHIHES, beta ', beta 
      if (beta .eq. 0.0d0) then
         ierr = -1
         indic = 1
         return
      endif 
c
      t = 1.0d0/beta 
      do 25 j=1, n
         u(j,1) = u(j,1)*t 
 25   continue
c------------------Arnoldi loop ----------------------------------------- 
c      fnorm = 0.0d0
      i1 = 1
 58   i = i1
      i1 = i + 1
      do 59 k=1, n
         x(k) = u(k,i)
 59   continue
      indic = 3
      return
 60   continue
      do 61 k=1, n
         u(k,i1) = y(k) 
 61   continue
      i0 =1
c switch  for Lanczos version 
c     i0 = max0(1, i-1)
      call mgsr (n, i0, i1, u, hh(1,i))
      fnorm = fnorm + ddot(i1, hh(1,i),1, hh(1,i),1)
      if (hh(i1,i) .eq. 0.0) m = i
      if  (i .lt. m) goto 58
c--------------done with arnoldi loop ---------------------------------
      rm = dble(m) 
      fnorm = dsqrt( fnorm / rm )
c------- put beta*e1 into z ------------------------------------------- 
      m1 = m+1 
      do 4 i=1,m1
         hh(i,m1) = 0.0
 4    continue
c     
c     compute initial dt when  job .lt. 1 
c
      if (job .ge. 0) goto 2
c
      t = 2.0*eps 
      do 41 k=1, m
         t = 2.0*t*dble(k+1)/rm
 41   continue
c
      t = rm* (t**(1.0d0/rm) ) / fnorm 
      if(verboz) print *, ' t, dt = ', t, dt
      t = dmin1(dabs(dt),t) 
      dt = dsign(t, dt) 
c---------------------- get the vector phi(Hm)e_1 + estimate ----------- 
 2    continue 
      z(1) = beta 
      do 3 k=2, m1
         z(k) = 0.0d0 
 3    continue
c-------get  : phi(H) * beta*e1
      call hes(ldg,m1,hh,ih,dt,z,rd,alp,alp0,wkc)
c-------error estimate 
      errst = dabs(z(m1))
      if(verboz) print *, ' error estimate =', errst 
c----------------------------------------------------------------------- 
      indic = 2
      return
      end
c-----------------------------------------------------------------------
      subroutine mgsr (n, i0, i1, ss, r)
c     implicit  real*8 (a-h,o-z)
      integer n, i0, i1
      real*8 ss(n,i1), r(i1)
c-----------------------------------------------------------------------
c modified gram - schmidt  with  partial  reortho. the vector ss(*,i1) is
c orthogonalized against the first i vectors  of ss  (which  are  already
c orthogonal).  the coefficients of the orthogonalization are returned in
c the array r
c------------------------------------------------------------------------
c local variables 
c 
      integer i, j, k, it
      real*8 hinorm, tet, ddot, t, dsqrt
      data  tet/10.0d0/

      do 53 j=1, i1
         r(j) = 0.0d0
 53   continue
      i = i1-1
      it = 0
 54   hinorm = 0.0d0
      it = it +1
      if (i .eq. 0) goto 56
c     
      do 55 j=i0, i
         t = ddot(n, ss(1,j),1,ss(1,i1),1)
         hinorm = hinorm + t**2
         r(j) = r(j) + t
         call daxpy(n,-t,ss(1,j),1,ss(1,i1),1)
 55   continue
      t = ddot(n, ss(1,i1), 1, ss(1,i1), 1)
 56   continue
c     
c     test for reorthogonalization see daniel et. al.
c     two reorthogonalization allowed ---
c     
      if (t*tet .le. hinorm .and. it .lt. 2) goto 54
      t =dsqrt(t)
      r(i1)= t
      if (t .eq. 0.0d0) return
      t = 1.0d0/t
      do 57  k=1,n
         ss(k,i1) = ss(k,i1)*t
 57   continue
      return
      end
c----------end-of-mgsr--------------------------------------------------
c-----------------------------------------------------------------------
      subroutine project(n, m, w, t, u, v) 
      integer n, m
      real*8 u(n,m), v(m), w(n), t, scal
c
c     computes the vector w = w + t * u * v
c
c local variables 
c
      integer j, k

      do 100 j=1,m
         scal = t*v(j) 
         do 99 k=1,n
            w(k) = w(k) + scal*u(k,j) 
 99      continue
 100  continue
      return
      end
c-----------------------------------------------------------------------
      subroutine hes (ndg,m1,hh,ih,dt,y,root,coef,coef0,w2)
c     implicit  real*8 (a-h,o-z)
      integer ndg, m1, ih 
      real*8 hh(ih+2,m1), y(m1)
      complex*16   coef(ndg), root(ndg), w2(m1)
      real*8 dt, coef0
c--------------------------------------------------------------------  
c computes phi ( H dt) * y    (1)
c where H = Hessenberg matrix (hh)
c y	  = arbitrary vector.
c ----------------------------
c ndg	= number of poles as determined by getrat
c m1    = dimension of hessenberg matrix
c hh	= hessenberg matrix (real)
c ih+2	= first dimension of hh
c dt	= scaling factor used for hh (see (1)) 
c y	= real vector. on return phi(H dt ) y is computed
c         and overwritten on y.
c root  = poles of the rational approximation to phi as
c         computed by getrat
c coef, 
c coef0 = coefficients of partial fraction phiansion 
c         
c  phi(t) ~ coef0 +  sum     Real [   coef(i) / (t - root(i)  ]
c                  i=1,ndg  
c
c valid for real t.
c coef0 is real, coef(*) is a complex array.
c
c--------------------------------------------------------------------  
c local variables 
c
      integer m1max
      parameter (m1max=70) 
      complex*16   hloc(m1max+1,m1max), t, zpiv, dcmplx
      real*8 yloc(m1max), dble
      integer i, j, ii
c     
c      if (m1 .gt. m1max) print *, ' *** ERROR : In HES, M+1 TOO LARGE'
c     
c     loop associated with the poles.
c     
      do 10 j=1,m1
         yloc(j) = y(j)
         y(j)    = y(j)*coef0
 10   continue
c     
      do 8 ii = 1, ndg
c     
c     copy Hessenberg matrix into temporary
c     
         do 2 j=1, m1
            do 1 i=1, j+1
               hloc(i,j) = dcmplx( dt*hh(i,j) )
 1          continue
            hloc(j,j) = hloc(j,j) - root(ii) 
            w2(j)     = dcmplx(yloc(j)) 
 2       continue 
c
c forward solve 
c 
         do 4 i=2,m1
            zpiv  = hloc(i,i-1) / hloc(i-1,i-1)
            do 3 j=i,m1
               hloc(i,j) = hloc(i,j) - zpiv*hloc(i-1,j)
 3          continue
            w2(i)     = w2(i) - zpiv*w2(i-1)
 4       continue 
c     
c     backward solve
c     
         do 6 i=m1,1,-1 
            t=w2(i)
            do 5 j=i+1,m1
               t = t-hloc(i,j)*w2(j)
 5          continue
            w2(i) = t/hloc(i,i)
 6       continue
c     
c     accumulate result in y.
c     
         do 7 i=1,m1
            y(i) = y(i) + dble ( coef(ii) * w2(i) ) 
 7       continue
 8    continue
      return
      end
c----------end-of-hes---------------------------------------------------
c-----------------------------------------------------------------------
      subroutine daxpy(n,t,x,indx,y,indy)
      integer n, indx, indy
      real*8 x(n), y(n), t
c-------------------------------------------------------------------
c does the following operation
c y <--- y + t * x ,   (replace by the blas routine daxpy )
c indx and indy are supposed to be one here
c-------------------------------------------------------------------
       integer k

       do 1 k=1,n
          y(k) = y(k) + x(k)*t
1      continue
       return
       end
c----------end-of-daxpy-------------------------------------------------
c----------------------------------------------------------------------- 
       function ddot(n,x,ix,y,iy)
       integer n, ix, iy
       real*8 ddot, x(n), y(n)
c-------------------------------------------------------------------
c computes the inner product t=(x,y) -- replace by blas routine ddot
c-------------------------------------------------------------------
       integer j
       real*8 t
        
       t = 0.0d0
       do 10 j=1,n
          t = t + x(j)*y(j)
10     continue
       ddot=t
       return
       end
c----------end-of-ddot--------------------------------------------------
c----------------------------------------------------------------------- 

