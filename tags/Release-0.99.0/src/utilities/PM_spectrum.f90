 FUNCTION PM_spectrum(FREQ, Hs, Tp)
!
! This function returns the value of the Pierson-Moskowitz spectrum
! for the frequency FREQ, significant wave height Hs, and peak period
! Tp.

   IMPLICIT none
   integer, parameter :: long=selected_real_kind(12,99)
   REAL(kind=long) :: pm_spectrum, Hs, Tp, A, B, freq, T02,   &
        four=4._long, pi, zero=0._long, one=1._long
   pi = acos (-one)

   T02=.71*Tp
   A = Hs*Hs/(four*pi*T02**4)
   B = one/(pi*T02**4)
   IF (freq.le.zero) then
      pm_spectrum = zero
   ELSE
      pm_spectrum = A/freq**5 * exp (-B/freq**4)
   ENDIF
   RETURN
 END FUNCTION PM_spectrum
