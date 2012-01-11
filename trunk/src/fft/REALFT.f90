  SUBROUTINE REALFT (DATA, N, ISIGN)
!
! Numerical recipes FFT of a real function routine.
!
! Calculates the Fourier transform of a set of n real-valued data points.
! Replaces this data by the positive frequency half of its complex Fourier
! transform.  The real valued first and last components of the complex transform are
! returned as data(1) and data(2).  n must be a power of 2.  For isign=-1 the
! inverse transform is returned and must be multiplied by 2/n.
!
      REAL(8) WR, WI, WPR, WPI, WTEMP, THETA
      REAL(4) DATA ( * )
      THETA = 3.141592653589793D0 / DBLE (N / 2)
      C1 = 0.5
      IF (ISIGN.EQ.1) THEN
         C2 = - 0.5
         CALL FOUR1 (DATA, N / 2, + 1)
      ELSE
         C2 = 0.5
         THETA = - THETA
      ENDIF
      WPR = - 2.0D0 * DSIN (0.5D0 * THETA) **2
      WPI = DSIN (THETA)
      WR = 1.0D0 + WPR
      WI = WPI
      N2P3 = N + 3
      DO 11 I = 2, N / 4
         I1 = 2 * I - 1
         I2 = I1 + 1
         I3 = N2P3 - I2
         I4 = I3 + 1
         WRS = SNGL (WR)
         WIS = SNGL (WI)
         H1R = C1 * (DATA (I1) + DATA (I3) )
         H1I = C1 * (DATA (I2) - DATA (I4) )
         H2R = - C2 * (DATA (I2) + DATA (I4) )
         H2I = C2 * (DATA (I1) - DATA (I3) )
         DATA (I1) = H1R + WRS * H2R - WIS * H2I
         DATA (I2) = H1I + WRS * H2I + WIS * H2R
         DATA (I3) = H1R - WRS * H2R + WIS * H2I
         DATA (I4) = - H1I + WRS * H2I + WIS * H2R
         WTEMP = WR
         WR = WR * WPR - WI * WPI + WR
         WI = WI * WPR + WTEMP * WPI + WI
11    END DO
      IF (ISIGN.EQ.1) THEN
        H1R = DATA (1)
        DATA (1) = H1R + DATA (2)
        DATA (2) = H1R - DATA (2)
      ELSE
        H1R = DATA (1)
        DATA (1) = C1 * (H1R + DATA (2) )
        DATA (2) = C1 * (H1R - DATA (2) )
        CALL FOUR1 (DATA, N / 2, - 1)
      ENDIF
      RETURN
      END SUBROUTINE REALFT
