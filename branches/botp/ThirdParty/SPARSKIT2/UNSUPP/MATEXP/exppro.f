      subroutine expprod (n, m, eps, tn, u, w, x, y, a, ioff, ndiag)
      real*8 eps, tn 
      real*8 a(n,ndiag), u(n,m+1), w(n), x(n), y(n)
      integer n, m, ndiag, ioff(ndiag)

c-----------------------------------------------------------------------
c this subroutine computes an approximation to the vector
c
c        	w :=  exp( - A * tn ) * w
c
c for matrices stored in diagonal (DIA) format.
c
c this routine constitutes an interface for the routine exppro for
c matrices stored in diagonal (DIA) format.
c-----------------------------------------------------------------------
c ARGUMENTS
c---------- 
c see exppro for meaning of parameters n, m, eps, tn, u, w, x, y.
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
      call exppro (n, m, eps, tn, u, w, x, y, indic, ierr) 
      if (indic .eq. 1) goto 102
c     
c     matrix vector-product for diagonal storage --
c     
      call oped(n, x, y, a, ioff, ndiag)
      goto 101
 102  continue
      return
      end
c----------end-of-expprod-----------------------------------------------
c-----------------------------------------------------------------------
      subroutine exppro (n, m, eps, tn, u, w, x, y, indic, ierr)
c     implicit  real*8 (a-h,o-z)
      integer n, m, indic, ierr
      real*8 eps, tn, u(n,m+1), w(n), x(n), y(n) 
c-----------------------------------------------------------------------
c     
c this subroutine computes an approximation to the vector
c
c        	w :=  exp( - A * tn ) * w
c
c where A is an arbitary matrix and w is a given input vector 
c uses a dynamic estimation of internal time advancement (dt) 
c----------------------------------------------------------------------- 
c THIS IS A REVERSE COMMUNICATION IMPLEMENTATION. 
c------------------------------------------------- 
c USAGE: (see also comments on indic below).
c------ 
c
c      indic = 0
c 1    continue
c      call exppro (n, m, eps, tn, u, w, x, y, indic)
c      if (indic .eq. 1) goto 2 <-- indic .eq. 1 means job is finished
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
c         the code will compute an approximation to exp(- tn * A) w
c         and overwrite the result onto w.
c
c u	= work array of size n*(m+1) (used to hold the Arnoldi basis )
c
c w	= real array of length n = input vector to  which exp(-A) is
c         to be applied. this is also an output argument 
c
c x, y  = two real work vectors of length at least  n each. 
c         see indic for usage.
c
c indic = integer used as indicator for the reverse communication.
c         in the first call enter indic = 0. See below for more.
c
c on return:
c-----------
c 
c w     = contains the resulting vector exp(-A * tn ) * w when 
c         exppro has finished (see indic) 
c
c indic = indicator for the reverse communication protocole.
c       * INDIC .eq. 1  means that exppro has finished and w contains the
c         result. 
c       * INDIC .gt. 1 ,  means that exppro has not finished and that
c         it is requesting another matrix vector product before
c         continuing. The user must compute Ax where A is the matrix
c         and x is the vector provided by exppro, and return the 
c         result in y. Then exppro must be called again without
c         changing any other argument. typically this must be
c         implemented in a loop with exppro being called as long
c         indic is returned with a value .ne. 1.
c
c ierr  = error indicator. 
c         ierr = 1 means phipro was called with indic=1 (not allowed)
c         ierr = -1 means that the input is zero the solution has been 
c         unchanged.
c 
c NOTES:  m should not exceed 60 in this version  (see mmax below)
c-----------------------------------------------------------------------
c written by Y. Saad -- version feb, 1991.
c----------------------------------------------------------------------- 
c For reference see following papers : 
c (1) E. Gallopoulos and Y. Saad: Efficient solution of parabolic 
c     equations by Krylov approximation methods. RIACS technical
c     report 90-14.
c (2) Y.Saad: Analysis of some Krylov subspace approximations to the
c     matrix exponential operator. RIACS Tech report. 90-14 
c-----------------------------------------------------------------------
c local variables 
c 
      integer mmax
      parameter (mmax=60) 
      real*8 errst, tcur, told, dtl, beta, red, dabs, dble
      real*8 hh(mmax+2,mmax+1), z(mmax+1)
      complex*16   wkc(mmax+1) 
      integer ih, job
      logical verboz
      data verboz/.true./
      save
c-----------------------------------------------------------------------
c indic = 3  means  passing through only with result of y= Ax to exphes
c indic = 2  means exphes has finished its job
c indic = 1  means exppro has finished its job (real end)/
c-----------------------------------------------------------------------
      ierr = 0 
      if (indic .eq. 3) goto 101
      if (indic .eq. 1) then 
         ierr = 1
         return
      endif
c----- 
      ih = mmax 
      m  = min0(m,mmax) 
      tcur = 0.0d0
      dtl = tn-tcur
      job = -1
c-------------------- outer loop ----------------------------- 
 100  continue
      if(verboz) print *,'In EXPPRO, current time = ', tcur ,'---------'
c------------------------------------------------------------- 
c ---- call exponential propagator --------------------------- 
c------------------------------------------------------------- 
      told = tcur 
 101  continue
c     if (told + dtl .gt. tn) dtl = tn-told
      call  exphes (n,m,dtl,eps,u,w,job,z,wkc,beta,errst,hh,ih,
     *              x,y,indic,ierr) 
c-----------------------------------------------------------------------
      if (ierr .ne. 0) return
      if (indic .ge. 3) return
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
      call project(n,m,u,z,w)
c never go beyond tcur
      job = 0
      dtl = dmin1(dtl, tn-tcur)
      if (dabs(tcur+dtl) .gt. dabs(tn)) dtl = tn-tcur 
      if (dabs(tcur) .lt. dabs(tn)) goto 100
      indic = 1
c
      return
      end
c----------end-of-expro-------------------------------------------------
c-----------------------------------------------------------------------
      subroutine exphes (n,m,dt,eps,u,w,job,z,wkc,beta,errst,hh,ih,
     *                   x, y, indic,ierr) 
c     implicit  real*8 (a-h,o-z)
      integer n, m, job, ih, indic, ierr
      real*8 hh(ih+2,m+1), u(n,m+1), w(n), z(m+1), x(n), y(n)
      complex*16 wkc(m+1) 
      real*8 dt, eps, beta, errst
c-----------------------------------------------------------------------
c this subroutine computes the Arnoldi basis and the corresponding 
c coeffcient vector in the approximation 
c 
c        	w  ::= beta  Vm  ym 
c               where ym = exp(- Hm *dt) * e1
c
c to the vector exp(-A dt) w where A is an arbitary matrix and 
c w is a given input vector. In case job = 0 the arnoldi basis 
c is recomputed. Otherwise the
c code assumes assumes that  u(*) contains an already computed 
c arnoldi basis and computes only the y-vector (which is stored in v(*))
c-----------------------------------------------------------------------
c on entry:
c---------- 
c n	= dimension of matrix
c
c m	= dimension of Krylov subspace (= degree of polynomial 
c         approximation to the exponential used. )
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
c w	= real array of length n = input vector to  which exp(-A) is
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
c w2	= resulting vector w2 = exp(-A *dt) * w
c beta  = real equal to the 2-norm of w. Needed if exppro will
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
      if (indic .ge. 3) goto 60
c
c------input fraction expansion of rational function ------------------ 
c chebyshev (14,14) 
      ldg= 7
      alp0 =  0.183216998528140087E-11
      alp(1)=( 0.557503973136501826E+02,-0.204295038779771857E+03)
      rd(1)=(-0.562314417475317895E+01, 0.119406921611247440E+01)
      alp(2)=(-0.938666838877006739E+02, 0.912874896775456363E+02)
      rd(2)=(-0.508934679728216110E+01, 0.358882439228376881E+01)
      alp(3)=( 0.469965415550370835E+02,-0.116167609985818103E+02)
      rd(3)=(-0.399337136365302569E+01, 0.600483209099604664E+01)
      alp(4)=(-0.961424200626061065E+01,-0.264195613880262669E+01)
      rd(4)=(-0.226978543095856366E+01, 0.846173881758693369E+01)
      alp(5)=( 0.752722063978321642E+00, 0.670367365566377770E+00)
      rd(5)=( 0.208756929753827868E+00, 0.109912615662209418E+02)
      alp(6)=(-0.188781253158648576E-01,-0.343696176445802414E-01)
      rd(6)=( 0.370327340957595652E+01, 0.136563731924991884E+02)
      alp(7)=( 0.143086431411801849E-03, 0.287221133228814096E-03)
      rd(7)=( 0.889777151877331107E+01, 0.166309842834712071E+02)
c-----------------------------------------------------------------------
c     
c     if job .gt. 0 skip arnoldi process:
c     
      if (job .gt. 0) goto 2
c------normalize vector w and put in first column of u -- 
      beta = dsqrt(ddot(n,w,1,w,1)) 
c-----------------------------------------------------------------------
      if(verboz) print *, ' In EXPHES, beta ', beta 
      if (beta .eq. 0.0d0) then
         ierr = -1 
         indic = 1
         return
      endif
c
      t = 1.0d0/beta 
      do 25 j=1, n
         u(j,1) = w(j)*t 
 25   continue
c------------------Arnoldi loop ------------------------------------- 
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
c     
c switch  for Lanczos version 
c     i0 = max0(1, i-1)
      call mgsr (n, i0, i1, u, hh(1,i))
      fnorm = fnorm + ddot(i1, hh(1,i),1, hh(1,i),1)
      if (hh(i1,i) .eq. 0.0) m = i
      if  (i .lt. m) goto 58
c--------------done with arnoldi loop ---------------------------------
      rm = dble(m) 
      fnorm = dsqrt( fnorm / rm )
c-------get  : beta*e1 into z 
      m1 = m+1 
      do 4 i=1,m1
         hh(i,m1) = 0.0
 4    continue
c
c     compute initial dt when  job .lt. 1 
c
      if (job .ge. 0) goto 2 
c
c     t = eps / beta 
c     
      t = eps 
      do 41 k=1, m-1
         t = t*(1.0d0 - dble(m-k)/rm ) 
 41   continue
c
      t = 2.0d0*rm* (t**(1.0d0/rm) )  / fnorm 
      if(verboz) print *, ' t, dt = ', t, dt
      t = dmin1(dabs(dt),t) 
      dt = dsign(t, dt) 
c     
 2    continue 
      z(1) = beta 
      do 3 k=2, m1
         z(k) = 0.0d0 
 3    continue
c-------get  : exp(H) * beta*e1
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
      subroutine project(n,m,u,v,w)
      integer n, m
      real*8 u(n,m), v(m), w(n)
c
c     computes the vector w = u * v
c
c local variables 
c
      integer j, k

      do 1 k=1,n
         w(k) = 0.d0
 1    continue

      do 100 j=1,m
         do 99 k=1,n
            w(k) = w(k) + v(j) * u(k,j) 
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
c computes exp ( H dt) * y    (1)
c where H = Hessenberg matrix (hh)
c y	  = arbitrary vector.
c ----------------------------
c ndg	= number of poles as determined by getrat
c m1    = dimension of hessenberg matrix
c hh	= hessenberg matrix (real)
c ih+2	= first dimension of hh
c dt	= scaling factor used for hh (see (1)) 
c y	= real vector. on return exp(H dt ) y is computed
c         and overwritten on y.
c root  = poles of the rational approximation to exp as
c         computed by getrat
c coef, 
c coef0 = coefficients of partial fraction expansion 
c         
c  exp(t) ~ coef0 +  sum     Real [   coef(i) / (t - root(i)  ]
c                  i=1,ndg  
c
c valid for real t.
c coef0 is real, coef(*) is a complex array.
c
c--------------------------------------------------------------------  
c local variables 
c
      integer m1max
      parameter (m1max=61) 
      complex*16 hloc(m1max+1,m1max), t, zpiv, dcmplx
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


