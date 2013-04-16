FUNCTION Erf(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE REAL ERROR FUNCTION
! Based upon a Fortran 66 routine in the Naval Surface Warfare Center's
! Mathematics Library (1993 version).
! Adapted by Alan.Miller @ vic.cmis.csiro.au
!-----------------------------------------------------------------------
! IMPLICIT NONE
INTEGER, PARAMETER :: dp = selected_real_kind(12,99)

REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables

REAL (dp), PARAMETER :: c = .564189583547756_dp, one = 1.0_dp, half = 0.5_dp, &
                        zero = 0.0_dp
REAL (dp), PARAMETER ::  &
           a(5) = (/ .771058495001320D-04, -.133733772997339D-02, &
                     .323076579225834D-01,  .479137145607681D-01, &
                     .128379167095513D+00 /),  &
           b(3) = (/ .301048631703895D-02,  .538971687740286D-01,  &
                     .375795757275549D+00 /),  &
           p(8) = (/ -1.36864857382717D-07, 5.64195517478974D-01,  &
                      7.21175825088309D+00, 4.31622272220567D+01,  &
                      1.52989285046940D+02, 3.39320816734344D+02,  &
                      4.51918953711873D+02, 3.00459261020162D+02 /), &
           q(8) = (/  1.00000000000000D+00, 1.27827273196294D+01,  &
                      7.70001529352295D+01, 2.77585444743988D+02,  &
                      6.38980264465631D+02, 9.31354094850610D+02,  &
                      7.90950925327898D+02, 3.00459260956983D+02 /), &
           r(5) = (/  2.10144126479064D+00, 2.62370141675169D+01,  &
                      2.13688200555087D+01, 4.65807828718470D+00,  &
                      2.82094791773523D-01 /),  &
           s(4) = (/  9.41537750555460D+01, 1.87114811799590D+02,  &
                      9.90191814623914D+01, 1.80124575948747D+01 /)
REAL (dp) :: ax, bot, t, top, x2
!-------------------------
ax = ABS(x)

IF (ax <= half) THEN
  t = x*x
  top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + one
  bot = ((b(1)*t + b(2))*t + b(3))*t + one
  fn_val = x*(top/bot)
  RETURN
END IF

IF (ax <= 4.0_dp) THEN
  top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
        + p(6))*ax + p(7))*ax + p(8)
  bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
        + q(6))*ax + q(7))*ax + q(8)
  fn_val = half + (half - EXP(-x*x)*top/bot)
  IF (x < zero) fn_val = -fn_val
  RETURN
END IF

IF (ax < 5.8_dp) THEN
  x2 = x*x
  t = one / x2
  top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
  bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + one
  fn_val = (c - top/(x2*bot)) / ax
  fn_val = half + (half - EXP(-x2)*fn_val)
  IF (x < zero) fn_val = -fn_val
  RETURN
END IF

fn_val = SIGN(one, x)
RETURN
END FUNCTION Erf
