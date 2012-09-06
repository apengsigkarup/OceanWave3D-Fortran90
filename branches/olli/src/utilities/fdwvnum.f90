SUBROUTINE fdwvnum (hbot, omega, kinf, ulen, grav, nu)
!
! This subroutine takes the dimensional (SI) values of omega, infinite
! depth wavenumber kinf, and depth hbot, and returns the non-dimensional
! finite depth wavenumber nu*ulen.
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: long=SELECTED_REAL_KIND(12,99)
  REAL(kind=long) hbot, ulen, omega, grav, nu, r0, hi, hh, his, his2, his4,       &
       his8, h, twoh, fourh, kinf, rh, r0h, fh, sh, bigh, smlh, bigexp,  &
       e2h
  DATA smlh / 1.d-05 /, bigh / 1.d+05 /, bigexp / 30.0d00 /

  h = hbot * kinf

  IF (h.lt.smlh) THEN
     nu = ulen * omega / SQRT (grav * hbot)
  ELSEIF (h.gt.bigh) THEN
     nu = ulen * kinf
  ELSE
     hi = 1. / h
     twoh = h + h
     fourh = twoh + twoh
     hh = h * h
     his = hi * hi
     his2 = his + his
     his4 = his2 + his2
     his8 = his4 + his4
!
! Evaluation of real root R0H of  R0H*TANH(R0H)=H
!
     e2h = EXP ( - MIN (twoh, bigexp) )
     IF (h.le.2.d00) THEN
        r0h = SQRT (h) * ( (.031d00 * h + .169d00) * h + 1.d00)
     ELSE
        r0h = h + h * e2h * (2.d00 - 12.d00 * e2h)
     ENDIF
     IF (h.le.4.8d00) THEN
        rh = hh - r0h * r0h
        sh = 1.d00 / (h - rh)
        fh = sh * (.5d00 * LOG ( (r0h + h) / (r0h - h) ) - r0h)
        r0h = r0h - fh * rh * (1.d00 + fh * r0h * h * sh)
     ENDIF
     nu = r0h * hi * kinf * ulen
  ENDIF
  RETURN
END SUBROUTINE fdwvnum
