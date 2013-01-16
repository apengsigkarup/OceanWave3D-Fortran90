SUBROUTINE linearstandingwave2D(g,z,time,sol,Lx,Ly,h,X,Y,pp,uu,vv,ww,eta,etax,etaxx,etay,etayy,Nfs)
! By Allan P. Engsig-Karup.
! FIXME: This implementation has in fact not been validated for all output
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: Nfs
TYPE (SFparam) :: sol
REAL(KIND=long) :: Lx, Ly, kx, ky, h, w, T, time, k, z, HH, c, g
REAL(KIND=long), DIMENSION(Nfs) :: eta, etax, etaxx, pp, uu, vv, ww, X,Y, etay, etayy
kx    = two*pi/Lx; ky = two*pi/Ly 
k = SQRT(kx**2+ky**2)
c = SQRT(g/k*TANH(k*h))
T = SQRT(  (two*pi)**2/(g*k*TANH(k*h)) ); 
sol%T=T
sol%k=k
sol%c=c
w = two*pi/T
!print * ,''
!print * ,'  Computed wave period, T = ',T
HH = sol%HH
eta   =  half*HH*COS(w*time)*COS(kx*X)*COS(ky*Y)
etax  = -kx*half*HH*COS(w*time)*SIN(kx*X)*COS(ky*Y)
etaxx = -kx**2*half*HH*COS(w*time)*COS(kx*X)*COS(ky*Y)
etay  = -ky*half*HH*COS(w*time)*SIN(kx*X)*COS(ky*Y)
etayy = -ky**2*half*HH*COS(w*time)*COS(kx*X)*COS(ky*Y)
pp    = -half*HH*c*COSH(k*(z+h))/SINH(k*h)*SIN(w*time)*COS(kx*X)*COS(ky*Y)
uu    =  kx*half*HH*c*COSH(k*(z+h))/SINH(k*h)*SIN(w*time)*SIN(kx*X)*COS(ky*Y)
vv    =  ky*half*HH*c*COSH(k*(z+h))/SINH(k*h)*SIN(w*time)*COS(kx*X)*SIN(ky*Y)
ww    = -half*k*HH*c*SINH(k*(z+h))/SINH(k*h)*SIN(w*time)*COS(kx*X)*COS(ky*Y)
END SUBROUTINE linearstandingwave2D
