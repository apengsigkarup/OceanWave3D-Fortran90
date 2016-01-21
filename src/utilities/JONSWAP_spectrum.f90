 FUNCTION JONSWAP_spectrum(FREQ, Hs, Tp, gamma)
!
! This function returns the value of the JONSWAP spectrum
! for the frequency FREQ, significant wave height Hs, and peak period
! Tp.

   IMPLICIT none
   integer, parameter :: long=selected_real_kind(12,99)
   REAL(kind=long) :: jonswap_spectrum, Hs, Tp, A, B, a_exp, alpha, gamma, &
        sigma, freq, f_p, four=4._long, pi, zero=0._long, one=1._long,     &
        two=2._long, half=.5_long, five=5._long, grav=9.82_long
   pi = acos (-one)
   f_p=one/Tp
   alpha=3.32285_long * Hs**2 * f_p**4
   A=alpha*grav**2/(two*pi)**4
   B=five/four * f_p**4

   IF (freq.le.zero) then 
      jonswap_spectrum = zero
   ELSEIF(freq<=f_p)THEN
      sigma=.07_long
      a_exp=exp(-half*((freq-f_p)/(sigma*f_p))**2)
      jonswap_spectrum = A/freq**5 * exp(-B/freq**4) * gamma**a_exp
   ELSE
      sigma=.09_long
      a_exp=exp(-half*((freq-f_p)/(sigma*f_p))**2)
      jonswap_spectrum = A/freq**5 * exp(-B/freq**4) * gamma**a_exp
   ENDIF
   RETURN
 END FUNCTION JONSWAP_spectrum
