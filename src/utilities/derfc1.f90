FUNCTION derfc1(ind, x) RESULT(fn_val)
!-----------------------------------------------------------------------

!         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION

!          DERFC1(IND,X) = ERFC(X)           IF IND = 0
!          DERFC1(IND,X) = EXP(X*X)*ERFC(X)  OTHERWISE

!-----------------------------------------------------------------------
   integer, parameter :: dp=selected_real_kind(12,99)
INTEGER, INTENT(IN)   :: ind
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: ax, t, w
INTEGER   :: i, k
REAL (dp), PARAMETER :: a(21) = (/ .1283791670955125738961589031215D+00,  &
        -.3761263890318375246320529677070D+00,  &
        .1128379167095512573896158902931D+00,  &
        -.2686617064513125175943235372542D-01,  &
        .5223977625442187842111812447877D-02,  &
        -.8548327023450852832540164081187D-03,  &
        .1205533298178966425020717182498D-03,  &
        -.1492565035840625090430728526820D-04,  &
        .1646211436588924261080723578109D-05,  &
        -.1636584469123468757408968429674D-06,  &
        .1480719281587021715400818627811D-07,  &
        -.1229055530145120140800510155331D-08,  &
        .9422759058437197017313055084212D-10,  &
        -.6711366740969385085896257227159D-11,  &
        .4463222608295664017461758843550D-12,  &
        -.2783497395542995487275065856998D-13,  &
        .1634095572365337143933023780777D-14,  &
        -.9052845786901123985710019387938D-16,  &
        .4708274559689744439341671426731D-17,  &
        -.2187159356685015949749948252160D-18,  &
        .7043407712019701609635599701333D-20 /)
!-------------------------------

!                     ABS(X) <= 1

ax = ABS(x)
IF (ax <= 1._dp) THEN
  t = x * x
  w = a(21)
  DO i = 1, 20
    k = 21 - i
    w = t * w + a(k)
  END DO
  fn_val = 0.5_dp + (0.5_dp-x*(1._dp+w))
  IF (ind /= 0) fn_val = EXP(t) * fn_val
  RETURN
END IF

!                       X < -1

IF (x <= 0._dp) THEN
  IF (x < -8.3_dp) GO TO 20
  IF (ind /= 0) THEN
    fn_val = 2._dp * EXP(x*x) - derfc0(ax)
    RETURN
  END IF
  fn_val = 2._dp - EXP(-x*x) * derfc0(ax)
  RETURN
END IF

!                       X > 1

IF (ind /= 0) THEN
  fn_val = derfc0(x)
  RETURN
END IF
fn_val = 0._dp
IF (x > 100._dp) RETURN
t = x * x
IF (t > -dxparg(1)) RETURN
fn_val = EXP(-t) * derfc0(x)
RETURN

!             LIMIT VALUE FOR LARGE NEGATIVE X

20 fn_val = 2._dp
IF (ind /= 0) fn_val = 2._dp * EXP(x*x)
RETURN
END FUNCTION derfc1
