SUBROUTINE eta_nonlinear_3D_standing(H,angle,Lx,Ly,xc,yc,omega,time,eta,etax,etay,etaxx,etayy,etat)
!
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
IMPLICIT NONE
REAL(KIND=long), INTENT(IN)  :: H,Lx,Ly,xc,yc,omega,time
REAL(KIND=long), INTENT(OUT) :: eta,etax,etay
REAL(KIND=long), INTENT(OUT) :: etaxx,etayy,etat
! REAL(KIND=long), OPTIONAL :: y,etay,etayy !Boundary fitted ?
! Local variables
REAL(KIND=long) :: kx, ky, half_H, k_omegat, sink_omegat
!
REAL(KIND=long) :: HH,k,w,angle,tmp,x,y,etaxy
!
x = xc*COS(-angle)-yc*SIN(-angle)
y = xc*SIN(-angle)+yc*COS(-angle)
!
kx    = two*pi/Lx/two; ky = zero*two*pi/Ly
k = SQRT(kx**2+ky**2)
w = omega
HH = H
!
!for now we use the linear solution...
IF (time==0.d0) THEN
!IF (0.d0==0.d0) THEN
eta = half*(0.8867_long*10.0**(-1)*COS(pi*X)+&
            SQRT(SQRT(two))*0.5243*10.0**(-2)*COS(two*pi*X)+&
            SQRT(SQRT(three))*0.4978*10.0**(-3)*COS(three*pi*X)+&
            SQRT(SQRT(four))*0.6542*10.0**(-4)*COS(four*pi*X)+&
            SQRT(SQRT(five))*0.1007*10.0**(-4)*COS(five*pi*X)+&
            SQRT(SQRT(six))*0.1653*10.0**(-5)*COS(six*pi*X)+&
            SQRT(SQRT(seven))*0.2753*10.0**(-6)*COS(seven*pi*X)+&
            SQRT(SQRT(eight))*0.4522*10.0**(-7)*COS(eight*pi*X))*COS(w*time)
etax = -half*(0.8867_long*10.0**(-1)*pi*SIN(pi*X)+&
            SQRT(SQRT(two))*0.5243*10.0**(-2)*(two*pi)*SIN(two*pi*X)+&
            SQRT(SQRT(three))*0.4978*10.0**(-3)*(three*pi)*SIN(three*pi*X)+&
            SQRT(SQRT(four))*0.6542*10.0**(-4)*(four*pi)*SIN(four*pi*X)+&
            SQRT(SQRT(five))*0.1007*10.0**(-4)*(five*pi)*SIN(five*pi*X)+&
            SQRT(SQRT(six))*0.1653*10.0**(-5)*(six*pi)*SIN(six*pi*X)+&
            SQRT(SQRT(seven))*0.2753*10.0**(-6)*(seven*pi)*SIN(seven*pi*X)+&
            SQRT(SQRT(eight))*0.4522*10.0**(-7)*(eight*pi)*SIN(eight*pi*X))*COS(w*time)
etaxx = -half*(0.8867_long*10.0**(-1)*pi**2*COS(pi*X)+&
            SQRT(SQRT(two))*0.5243*10.0**(-2)*(two*pi)**2*COS(two*pi*X)+&
            SQRT(SQRT(three))*0.4978*10.0**(-3)*(three*pi)**2*COS(three*pi*X)+&
            SQRT(SQRT(four))*0.6542*10.0**(-4)*(four*pi)**2*COS(four*pi*X)+&
            SQRT(SQRT(five))*0.1007*10.0**(-4)*(five*pi)**2*COS(five*pi*X)+&
            SQRT(SQRT(six))*0.1653*10.0**(-5)*(six*pi)**2*COS(six*pi*X)+&
            SQRT(SQRT(seven))*0.2753*10.0**(-6)*(seven*pi)**2*COS(seven*pi*X)+&
            SQRT(SQRT(eight))*0.4522*10.0**(-7)*(eight*pi)**2*COS(eight*pi*X))*COS(w*time)
etay  = zero
etayy = zero
etaxy = zero
etat  = -half*w*(0.8867_long*10.0**(-1)*COS(pi*X)+&
            SQRT(SQRT(two))*0.5243*10.0**(-2)*COS(two*pi*X)+&
            SQRT(SQRT(three))*0.4978*10.0**(-3)*COS(three*pi*X)+&
            SQRT(SQRT(four))*0.6542*10.0**(-4)*COS(four*pi*X)+&
            SQRT(SQRT(five))*0.1007*10.0**(-4)*COS(five*pi*X)+&
            SQRT(SQRT(six))*0.1653*10.0**(-5)*COS(six*pi*X)+&
            SQRT(SQRT(seven))*0.2753*10.0**(-6)*COS(seven*pi*X)+&
            SQRT(SQRT(eight))*0.4522*10.0**(-7)*COS(eight*pi*X))*SIN(w*time)
ELSE
  eta   =  half*0.08934064760240*COS(w*time)*COS(pi*X)
  etax  = -pi*half*0.08934064760240*COS(w*time)*SIN(pi*X)
  etaxx = -pi**2*half*0.08934064760240*COS(w*time)*COS(pi*X)
  etat  = -half*w*0.08934064760240*SIN(w*time)*COS(pi*X)
  etay  = zero
  etayy = zero
  etaxy = zero
ENDIF
!
tmp = etax
etax = tmp*COS(angle)-etay*SIN(angle)
etay = tmp*SIN(angle)+etay*COS(angle)
!
tmp = etaxx
etaxx = tmp*COS(angle)**2+etayy*SIN(angle)**2-2.0_long*COS(angle)*SIN(angle)*etaxy
etayy = tmp*SIN(angle)**2+etayy*COS(angle)**2+2.0_long*COS(angle)*SIN(angle)*etaxy
!xc = x*COS(angle)-y*SIN(angle)
!yc = x*SIN(angle)+y*COS(angle)


END SUBROUTINE eta_nonlinear_3D_standing
