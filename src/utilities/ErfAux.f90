MODULE Erf_Auxilliary
! Contains the NSWC functions SPMPAR, DPMPAR, EPSLN, DEPSLN, EXPARG & DXPARG
!-----------------------------------------------------------------------
!     WRITTEN using F90 intrinsics by
!        Alan Miller
!        CSIRO Mathematical & Information Sciences
!        CLAYTON, VICTORIA, AUSTRALIA 3169
!     Latest revision - 1 February 1997
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)

CONTAINS

FUNCTION ipmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
!     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
!     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...

!  INTEGERS.

!     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM

!               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )

!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.

!     IPMPAR(1) = A, THE BASE (radix).

!     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS (digits).

!     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE (huge).

!  FLOATING-POINT NUMBERS.

!     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
!     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
!     NONZERO NUMBERS ARE REPRESENTED IN THE FORM

!               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)

!               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
!               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.

!     IPMPAR(4) = B, THE BASE.

!  SINGLE-PRECISION

!     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.

!     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.

!     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.

!  DOUBLE-PRECISION

!     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.

!     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.

!     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.

!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
INTEGER             :: fn_val

SELECT CASE(i)
  CASE( 1)
    fn_val = RADIX(i)
  CASE( 2)
    fn_val = DIGITS(i)
  CASE( 3)
    fn_val = HUGE(i)
  CASE( 4)
    fn_val = RADIX(1.0)
  CASE( 5)
    fn_val = DIGITS(1.0)
  CASE( 6)
    fn_val = MINEXPONENT(1.0)
  CASE( 7)
    fn_val = MAXEXPONENT(1.0)
  CASE( 8)
    fn_val = DIGITS(1.0D0)
  CASE( 9)
    fn_val = MINEXPONENT(1.0D0)
  CASE(10)
    fn_val = MAXEXPONENT(1.0D0)
  CASE DEFAULT
    RETURN
END SELECT

RETURN
END FUNCTION ipmpar



FUNCTION spmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
REAL                :: fn_val

! Local variable
REAL                :: one = 1.0

SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
END SELECT

RETURN
END FUNCTION spmpar



FUNCTION dpmpar (i) RESULT(fn_val)
!-----------------------------------------------------------------------

!     DPMPAR PROVIDES THE DOUBLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     DOUBLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

!        DPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

!        DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

!        DPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!-----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: i
REAL (dp)           :: fn_val

! Local variable
REAL (dp)    :: one = 1._dp

SELECT CASE (i)
  CASE (1)
    fn_val = EPSILON(one)
  CASE (2)
    fn_val = TINY(one)
  CASE (3)
    fn_val = HUGE(one)
END SELECT

RETURN
END FUNCTION dpmpar


FUNCTION epsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.0 + EPS .GT. 1.0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
IMPLICIT NONE
REAL                :: fn_val

! Local variable
REAL                :: one = 1.0

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION epsln


FUNCTION exparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     EXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
REAL                :: fn_val

! Local variable
REAL                :: one = 1.0

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION exparg


FUNCTION depsln () RESULT(fn_val)
!--------------------------------------------------------------------
!     THE EVALUATION OF LN(EPS) WHERE EPS IS THE SMALLEST NUMBER
!     SUCH THAT 1.D0 + EPS .GT. 1.D0 .  L IS A DUMMY ARGUMENT.
!--------------------------------------------------------------------
IMPLICIT NONE
REAL (dp)           :: fn_val

! Local variable
REAL (dp)    :: one = 1._dp

fn_val = LOG( EPSILON(one) )
RETURN
END FUNCTION depsln


FUNCTION dxparg (l) RESULT(fn_val)
!--------------------------------------------------------------------
!     IF L = 0 THEN  DXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
!     DEXP(W) CAN BE COMPUTED.
!
!     IF L IS NONZERO THEN  DXPARG(L) = THE LARGEST NEGATIVE W FOR
!     WHICH THE COMPUTED VALUE OF DEXP(W) IS NONZERO.
!
!     NOTE... ONLY AN APPROXIMATE VALUE FOR DXPARG(L) IS NEEDED.
!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: l
REAL (dp)           :: fn_val

! Local variable
REAL (dp)    :: one = 1._dp

IF (l == 0) THEN
  fn_val = LOG( HUGE(one) )
ELSE
  fn_val = LOG( TINY(one) )
END IF
RETURN
END FUNCTION dxparg


FUNCTION derfc1(ind, x) RESULT(fn_val)
!-----------------------------------------------------------------------

!         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION

!          DERFC1(IND,X) = ERFC(X)           IF IND = 0
!          DERFC1(IND,X) = EXP(X*X)*ERFC(X)  OTHERWISE

!-----------------------------------------------------------------------
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

FUNCTION derfc0(x) RESULT(fn_val)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

!           EVALUATION OF EXP(X**2)*ERFC(X) FOR X >= 1

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!        APRIL 1992
!-------------------------------
REAL (dp)            :: t, u, v, z
REAL (dp), PARAMETER :: rpinv = .56418958354775628694807945156077259_dp
REAL (dp), PARAMETER :: p0 = .16506148041280876191828601D-03,  &
                        p1 =  .15471455377139313353998665D-03,  &
                        p2 =  .44852548090298868465196794D-04,  &
                        p3 = -.49177280017226285450486205D-05,  &
                        p4 = -.69353602078656412367801676D-05,  &
                        p5 = -.20508667787746282746857743D-05,  &
                        p6 = -.28982842617824971177267380D-06,  &
                        p7 = -.17272433544836633301127174D-07,  &
                        q1 =  .16272656776533322859856317D+01,  &
                        q2 =  .12040996037066026106794322D+01,  &
                        q3 =  .52400246352158386907601472D+00,  &
                        q4 =  .14497345252798672362384241D+00,  &
                        q5 =  .25592517111042546492590736D-01,  &
                        q6 =  .26869088293991371028123158D-02,  &
                        q7 =  .13133767840925681614496481D-03
REAL (dp), PARAMETER :: r0 =  .145589721275038539045668824025D+00,  &
                        r1 = -.273421931495426482902320421863D+00,  &
                        r2 =  .226008066916621506788789064272D+00,  &
                        r3 = -.163571895523923805648814425592D+00,  &
                        r4 =  .102604312032193978662297299832D+00,  &
                        r5 = -.548023266949835519254211506880D-01,  &
                        r6 =  .241432239725390106956523668160D-01,  &
                        r7 = -.822062115403915116036874169600D-02,  &
                        r8 =  .180296241564687154310619200000D-02
REAL (dp), PARAMETER :: a0 = -.45894433406309678202825375D-03,   &
                        a1 = -.12281298722544724287816236D-01,  &
                        a2 = -.91144359512342900801764781D-01,  &
                        a3 = -.28412489223839285652511367D-01,  &
                        a4 =  .14083827189977123530129812D+01,  &
                        a5 =  .11532175281537044570477189D+01,  &
                        a6 = -.72170903389442152112483632D+01,  &
                        a7 = -.19685597805218214001309225D+01,  &
                        a8 =  .93846891504541841150916038D+01,  &
                        b1 =  .25136329960926527692263725D+02,  &
                        b2 =  .15349442087145759184067981D+03,  &
                        b3 = -.29971215958498680905476402D+03,  &
                        b4 = -.33876477506888115226730368D+04,  &
                        b5 =  .28301829314924804988873701D+04,  &
                        b6 =  .22979620942196507068034887D+05,  &
                        b7 = -.24280681522998071562462041D+05,  &
                        b8 = -.36680620673264731899504580D+05,  &
                        b9 =  .42278731622295627627042436D+05,  &
                        b10=  .28834257644413614344549790D+03,  &
                        b11=  .70226293775648358646587341D+03
REAL (dp), PARAMETER :: c0 = -.7040906288250128001000086D-04,   &
                        c1 = -.3858822461760510359506941D-02,  &
                        c2 = -.7708202127512212359395078D-01,  &
                        c3 = -.6713655014557429480440263D+00,  &
                        c4 = -.2081992124162995545731882D+01,  &
                        c5 =  .2898831421475282558867888D+01,  &
                        c6 =  .2199509380600429331650192D+02,  &
                        c7 =  .2907064664404115316722996D+01,  &
                        c8 = -.4766208741588182425380950D+02,  &
                        d1 =  .5238852785508439144747174D+02,  &
                        d2 =  .9646843357714742409535148D+03,  &
                        d3 =  .7007152775135939601804416D+04,  &
                        d4 =  .8515386792259821780601162D+04,  &
                        d5 = -.1002360095177164564992134D+06,  &
                        d6 = -.2065250031331232815791912D+06,  &
                        d7 =  .5695324805290370358175984D+06,  &
                        d8 =  .6589752493461331195697873D+06,  &
                        d9 = -.1192930193156561957631462D+07
REAL (dp), PARAMETER :: e0 = .540464821348814822409610122136D+00,  &
                        e1 = -.261515522487415653487049835220D-01, &
                        e2 = -.288573438386338758794591212600D-02, &
                        e3 = -.529353396945788057720258856000D-03
REAL (dp), PARAMETER :: s1 = .75000000000000000000D+00,   &
        s2  = -.18750000000000000000D+01, s3  = .65625000000000000000D+01,  &
        s4  = -.29531250000000000000D+02, s5  = .16242187500000000000D+03,  &
        s6  = -.10557421875000000000D+04, s7  = .79180664062500000000D+04,  &
        s8  = -.67303564453125000000D+05, s9  = .63938386230468750000D+06,  &
        s10 = -.67135305541992187500D+07, s11 = .77205601373291015625D+08
!-------------------------------
!     RPINV = 1/SQRT(PI)
!-------------------------------

!                     1 <= X <= 2

IF (x <= 2._dp) THEN
  u = ((((((p7*x + p6)*x + p5)*x + p4)*x + p3)*x + p2)*x + p1) * x + p0
  v = ((((((q7*x + q6)*x + q5)*x + q4)*x + q3)*x + q2)*x + q1) * x + 1._dp

  t = (x-3.75_dp) / (x+3.75_dp)
  fn_val = (((((((((u/v)*t + r8)*t + r7)*t + r6)*t + r5)*t + r4)*t + r3)*t + &
           r2)*t + r1) * t + r0
  RETURN
END IF

!                     2 < X <= 4

IF (x <= 4._dp) THEN
  z = 1._dp / (2.5_dp + x*x)
  u = (((((((a8*z + a7)*z + a6)*z + a5)*z + a4)*z + a3)*z + a2)*z + a1) * z + a0
  v = ((((((((((b11*z + b10)*z + b9)*z + b8)*z + b7)*z + b6)*z + b5)*z +  &
      b4)*z + b3)*z + b2)*z + b1) * z + 1._dp

  t = 13._dp * z - 1._dp
  fn_val = ((((u/v)*t + e2)*t + e1)*t + e0) / x
  RETURN
END IF

!                     4 < X < 50

IF (x < 50._dp) THEN
  z = 1._dp / (2.5_dp + x*x)
  u = (((((((c8*z + c7)*z + c6)*z + c5)*z + c4)*z + c3)*z + c2)*z + c1) * z + &
      c0
  v = ((((((((d9*z + d8)*z + d7)*z + d6)*z + d5)*z + d4)*z + d3)*z + d2)*z +  &
      d1)*z + 1._dp

  t = 13._dp * z - 1._dp
  fn_val = (((((u/v)*t + e3)*t + e2)*t + e1)*t + e0) / x
  RETURN
END IF

!                        X >= 50

t = (1._dp/x) ** 2
z = (((((((((((s11*t + s10)*t + s9)*t + s8)*t + s7)*t + s6)*t + s5)*t +  &
    s4)*t + s3)*t + s2)*t + s1)*t - 0.5_dp) * t + 1._dp
fn_val = rpinv * (z/x)
RETURN
END FUNCTION derfc0

END MODULE Erf_Auxilliary

