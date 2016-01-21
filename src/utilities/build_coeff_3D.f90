 SUBROUTINE build_coeff_3D (ETA, beta, N, Hs, Tp, DELT, SEED, seed2, beta0,gamma)
!
! This subroutine returns Z(w) and kvec, the Fourier transform of
! the wave elevation, zeta(t) and a direction beta(w) for each frequency.  
! The coefficients are based on a JONSWAP spectrum with a normal spreading 
! around the heading angle beta0.
!
   USE Error_Function
   IMPLICIT none
   integer, parameter :: long=selected_real_kind(12,99)
   INTEGER SEED, seed2, n
   REAL(kind=long) :: Hs, Tp, delt, beta0,gamma
   REAL(kind=long) :: ETA(n), beta(n)
   ! Locals
   Integer i
   REAL gasdev, ran1_b
   REAL(kind=long) ::  pi,  df, cst, freq, spec,                       &
        jonswap_spectrum, zero=0._long,                                &
        one=1._long, two=2._long, four=4._long, half=.5_long,          &
        twopi, ran1, psi, mag, deg2rad, fp, sigma, erf1, erf2,         &
        rad2deg
   COMPLEX(kind=long) :: phase
   EXTERNAL gasdev, ran1_b, jonswap_spectrum
!
! Open the file to write the spectral coeeficients out to 
!
   open(unit=78,file='spectrum',status='unknown')
   write(78,79)
79 format('# Spectrum for the directional random wave wavemaker signal. f, T, S(f), beta(f)')

   pi = acos(-one)
   twopi = two * pi
   deg2rad=pi/180._long; rad2deg=1/deg2rad
   df = one / (n * delt)
   cst = sqrt(df)
   fp=1/Tp

! Initialize
   beta(:)=beta0
   
! Take care of zero and Nyquist frequencies.
   eta (1) = zero; 
   eta (2) = zero

! Build the Fourier coefficients using the chosen Spectrum and a random
! phase.  Each frequency has one heading angle which is also based on this 
! random number and the accumulated PDF of the spectrum.  

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

! Get the wave spectrum value and the direction at this frequency.

      spec= jonswap_spectrum(freq, Hs, Tp, gamma)

      sigma=0.2_long*(freq/fp)**2.5_long
      erf1=Erf(pi/(two*sqrt(two*sigma)))
      erf2=InvErf((two*psi-one)*erf1)

      beta(i)=beta0+rad2deg*sqrt(2*sigma)*erf2

      write(78,*)freq,1/freq,spec,beta(i)
      !
      ! Load the Fourier coefficients for this frequency
      !
      ! Real part
      eta (2*i-1) = REAL(sqrt (two * spec) * cst * phase)
      ! Imaginary part
      eta (2*i)   = AIMAG(sqrt (two * spec) * cst * phase)

   END DO
   close(78)

   RETURN
 END SUBROUTINE BUILD_COEFF_3D
