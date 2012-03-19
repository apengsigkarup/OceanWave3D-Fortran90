      subroutine dc(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      real*8 rhs(n), sol(n), fpar(16), w(*)
c-----------------------------------------------------------------------
c     This a version of Defect Correction Method (DC) implemented 
c     with reverse communication.
c
c     the space of the `w' is used as follows:
c     (1) input
c     (2) output
c
c     TOTAL SIZE REQUIRED == 2*n
c
c     Implemented by Allan P. Engsig-Karup, apek@imm.dtu.dk.
c
c     SPARSKIT2 routine for GMRES method was used as template for
c     this implementation.
c
c-----------------------------------------------------------------------
c     external functions used
c
      real*8 distdot
      external distdot
c
      real*8 one, zero
      parameter(one=1.0D0, zero=0.0D0)
c
c     local variables, ptr and p2 are temporary pointers,
c     hess points to the Hessenberg matrix,
c     vc, vs point to the cosines and sines of the Givens rotations
c     vrn points to the vectors of residual norms, more precisely
c     the right hand side of the least square problem solved.
c
      integer i,ii,idx,k,m,ptr,p2,hess,vc,vs,vrn
      real*8 alpha, c, s
      logical lp, rp
      save
c
c     check the status of the call
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (10, 20) ipar(10)
c
c     initialization
c
      i = 2*n
      call bisinit(ipar,fpar,i,1,lp,rp,w)
      if (ipar(1).lt.0) return
c
c     request for matrix vector multiplication A*x in the initialization
c
      ipar(1)  = 1
      ipar(8)  = n+1
      ipar(9)  = 1
      ipar(10) = 1
c
c     copy initial guess to workspace to be input vector
c
      do i = 1, n
         w(n+i) = sol(i)
      enddo
      return
c
c     compute residual r = b - Ax
c
 10   ipar(7) = ipar(7) + 1
c
c     put residual in workspace as input vector
c
      do i = 1, n
         w(n+i) = rhs(i) - w(i)
      enddo
      fpar(11) = fpar(11) + n
      ipar(1) = 3
      ipar(10) = 2
      return
c
c     left preconditioning
c
 20   alpha = sqrt(distdot(n,w,1,w,1))
c      print *,alpha
      fpar(11) = fpar(11) + 2*n
      if (ipar(7).eq.1 .and. ipar(3).ne.999) then
c
c     compute initial error norm. Only, do this after first matvec.
c
         fpar(3) = alpha
c
c     compute target error norm (scaled by user-defined tolerance levels)
c
         fpar(4) = fpar(1) * alpha + fpar(2)
c      print *,'fpar(4)=',fpar(4)

      endif
c
c     compute CURRENT error norm 
c
      fpar(5) = alpha
c      print *,'current error norm = ',alpha
c
c     check if we satisfy the stopping condition by comparing current 
c     error norm with target error norm
c
      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
c
c        if full-filled, then terminate solver
c
c      print *,'we should stop because condition is satisfied'

         ipar(1) = 0
c
c        current error norm is stored
c
         fpar(6) = alpha
         goto 300
      endif
c
c     defect correction step
c
      do i = 1, n
         sol(i) = sol(i) + w(i)
      enddo
c
c     copy solution to workspace
c
      do i = 1, n
         w(n+i) = sol(i)
      enddo
      fpar(11) = fpar(11) + n
c
c     repeat defect correction procedure
c
      if (ipar(7)<ipar(6)) then
          ipar(1)  = 1
          ipar(10) = 1          
      else
c
c         terminate because of user-defined limit to number of iterations
c
          ipar(1) = -1
      endif
      return
c
c     process the complete stopping criteria
c
c      if (ipar(3).eq.999) then
c
c     FIXME HERE
c
c         print*,'FIXME, dc.f!!!'
c         ipar(1) = 10
c         ipar(8) = -1
c         ipar(9) = idx + 1
c         ipar(10) = 7
c         return
c      else if (ipar(3).lt.0) then
c         if (ipar(7).le.m+1) then
c            fpar(3) = abs(w(vrn+1))
c            if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
c         endif
c         fpar(6) = abs(w(vrn+k))
c      else
c         fpar(6) = fpar(5)
c      endif
c
c     termination, set error code, compute convergence rate
c
c      if (ipar(1).gt.0) then
c         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
c            ipar(1) = 0
c         else if (ipar(3).ne.999 .and. fpar(6).le.fpar(4)) then
c            ipar(1) = 0
c         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
c            ipar(1) = -1
c         else
c            ipar(1) = -10
c         endif
c      endif
 300  if (fpar(3).ne.zero .and. fpar(6).ne.zero .and.
     +     ipar(7).gt.ipar(13)) then
c
c    estimate convergence rate
c
         fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
      else
         fpar(7) = zero
      endif
      return
      end
c-----end-of-dc
