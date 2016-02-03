 SUBROUTINE build_coeff (ETA, N, Hs, Tp, DELT, SEED, seed2, i_spec,gamma)
!
! This subroutine returns Z(w), the Fourier transform of
! the wave elevation, zeta(t).  The coefficients are based on
! either a P-M or a JONSWAP spectrum, depending on i_spec.
!
   IMPLICIT none
   integer, parameter :: long=selected_real_kind(12,99)
   INTEGER SEED, seed2, i_spec, n, i
   REAL gasdev, ran1_b
   REAL(kind=long) :: Hs, Tp, pi, delt, df, cst, freq, spec,           &
        pm_spectrum, jonswap_spectrum, zero=0._long,          &
        one=1._long, two=2._long, four=4._long, half=.5_long,          &
        twopi, ran1, psi, mag, gamma
   REAL(kind=long) :: ETA(n) ! GD : modif so that types are compatible... n has to be a power of 2
   COMPLEX(kind=long) :: phase
   EXTERNAL gasdev, ran1_b, pm_spectrum, jonswap_spectrum

   pi = acos(-one)
   twopi = two * pi
   df = one / (n * delt)
   cst = sqrt(df)

! Take care of zero and Nyquist frequencies.

   eta (1) = zero
   eta (2) = zero
   open(unit=78,file='spectrum',status='unknown')
   write(78,79)
79 format('# Spectrum for the random wave wavemaker signal. f, T, S(f)')

! Build the Fourier coefficients using the chosen Spectrum and a random
! phase.

   DO  i = 2, n / 2
      freq = (i - 1) * df
      psi=ran1_b(seed2)
!
! Choose here whether the Fourier amplitudes should be exactly on the spectrum
! or Gaussianly distributed around the spectral value.
!
      mag=one
!      mag=gasdev(seed)  ! This will give a Guassian distribution of magnitude.
!
      phase=mag*cmplx(cos(twopi*psi),sin(twopi*psi))

! Get the wave spectrum value at this frequency.
      IF(i_spec==0)THEN
         spec= pm_spectrum(freq, Hs, Tp)
      ELSE
         spec= jonswap_spectrum(freq, Hs, Tp,gamma)
      ENDIF

      write(78,*)freq,1/freq,spec

      !eta (i)   = sqrt (two * spec) * cst * phase
      ! Real part
      eta (2*i-1) = REAL(sqrt (two * spec) * cst * phase)
      ! Imaginary part
      eta (2*i)   = AIMAG(sqrt (two * spec) * cst * phase)

   END DO

   close(78)

   RETURN
 END SUBROUTINE BUILD_COEFF
