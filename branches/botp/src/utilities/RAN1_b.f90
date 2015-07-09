 FUNCTION RAN1_b(IDUM)
   implicit none
   integer idum,ia,im,iq,ir,ntab,ndiv
   real ran1_b,am,eps,rnmx
   parameter(ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,      &
        ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)

! Minimal random number generator of Park and Miller with Bays-Durham
! shuffle and added safeguards.  Returns a uniform random deviate between
! 0.0 and 1.0 (exclusive of the endpoint values.)  Call with idum a negative
! integer to initialize; thereafter, do not alter idum between successive
! deviates in a sequence.  tnmx should approximate the largest floating
! value that is less than 1.  From Numerical Recipes.

   integer j,k,iv(ntab),iy
   save iv,iy
   data iv /ntab*0/, iy /0/
   if(idum.le.0.or.iy.eq.0)then
      idum=max(-idum,1)
      do j=ntab+8,1,-1
         k=idum/iq
         idum=ia*(idum-k*iq)-ir*k
         if(idum.lt.0) idum=idum+im
         if(j.le.ntab) iv(j)=idum
      end do
      iy=iv(1)
   end if
   k=idum/iq
   idum=ia*(idum-k*iq)-ir*k
   if(idum.lt.0) idum=idum+im
   j=1+iy/ndiv
   iy=iv(j)
   iv(j)=idum
   ran1_b=min(am*iy,rnmx)
   return
 END FUNCTION RAN1_B
