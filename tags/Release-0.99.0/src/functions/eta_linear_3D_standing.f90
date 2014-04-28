SUBROUTINE eta_linear_3D_standing(H,angle,Lx,Ly,xc,yc,omega,time,eta,etax,etay,etaxx,etayy,etat)
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
REAL(KIND=long) :: HH,k,w,angle,tmp,x,y
!
x = xc*COS(-angle)-yc*SIN(-angle)
y = xc*SIN(-angle)+yc*COS(-angle)
!
kx    = two*pi/Lx; ky = two*pi/Ly
k = SQRT(kx**2+ky**2)
w = omega
HH = H
!
eta   =  half*HH*COS(w*time)*COS(kx*X)*COS(ky*Y)
etax  = -kx*half*HH*COS(w*time)*SIN(kx*X)*COS(ky*Y)
etaxx = -kx**2*half*HH*COS(w*time)*COS(kx*X)*COS(ky*Y)
etay  = -ky*half*HH*COS(w*time)*SIN(kx*X)*COS(ky*Y)
etayy = -ky**2*half*HH*COS(w*time)*COS(kx*X)*COS(ky*Y)
etat  = -half*w*HH*SIN(w*time)*COS(kx*X)*COS(ky*Y)
!
tmp = etax
etax = tmp*COS(angle)-etay*SIN(angle)
etay = tmp*SIN(angle)+etay*COS(angle)
!
tmp = etaxx
etaxx = tmp*COS(angle)-etayy*SIN(angle)
etayy = tmp*SIN(angle)+etayy*COS(angle)
!
!xc = x*COS(angle)-y*SIN(angle)
!yc = x*SIN(angle)+y*COS(angle)


END SUBROUTINE eta_linear_3D_standing
