!> Reference: Agnon & Glozman (1996) [p. 145].
!!    "Periodic solutions for a complex Hamiltonian system: New standing
!!    water-waves"
!!
!! Based on the assumptions:
!!    L = 2 m
!!    T = 1.13409 s
!!    deep water wave -> we set kh = 2*pi -> h = 2 m
!!    H  = 0.08934064760240;       % wave height
!!
!! By Allan P. Engsig-Karup.
!<
SUBROUTINE nonlinearstandingwave1D(k,h,X,uu, ww, eta, etax, etaxx, Nfs)
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: Nfs
REAL(KIND=long) :: L,k,h
REAL(KIND=long), DIMENSION(Nfs) :: eta, etax, etaxx, uu, ww, X
L = two*pi/k
IF (L/h < half) THEN ! check if conditions given result in a deep water wave
	PRINT *,'Conditions for a deep water wave is not given. Check ICs.'
    STOP
END IF
eta = half*(0.8867_long*10.0**(-1)*COS(pi*X)+&
            SQRT(SQRT(two))*0.5243*10.0**(-2)*COS(two*pi*X)+&
            SQRT(SQRT(three))*0.4978*10.0**(-3)*COS(three*pi*X)+&
            SQRT(SQRT(four))*0.6542*10.0**(-4)*COS(four*pi*X)+&
            SQRT(SQRT(five))*0.1007*10.0**(-4)*COS(five*pi*X)+&
            SQRT(SQRT(six))*0.1653*10.0**(-5)*COS(six*pi*X)+&
            SQRT(SQRT(seven))*0.2753*10.0**(-6)*COS(seven*pi*X)+&
            SQRT(SQRT(eight))*0.4522*10.0**(-7)*COS(eight*pi*X))
etax = -half*(0.8867_long*10.0**(-1)*pi*SIN(pi*X)+&
            SQRT(SQRT(two))*0.5243*10.0**(-2)*(two*pi)*SIN(two*pi*X)+&
            SQRT(SQRT(three))*0.4978*10.0**(-3)*(three*pi)*SIN(three*pi*X)+&
            SQRT(SQRT(four))*0.6542*10.0**(-4)*(four*pi)*SIN(four*pi*X)+&
            SQRT(SQRT(five))*0.1007*10.0**(-4)*(five*pi)*SIN(five*pi*X)+&
            SQRT(SQRT(six))*0.1653*10.0**(-5)*(six*pi)*SIN(six*pi*X)+&
            SQRT(SQRT(seven))*0.2753*10.0**(-6)*(seven*pi)*SIN(seven*pi*X)+&
            SQRT(SQRT(eight))*0.4522*10.0**(-7)*(eight*pi)*SIN(eight*pi*X))
etaxx = -half*(0.8867_long*10.0**(-1)*pi**2*COS(pi*X)+&
            SQRT(SQRT(two))*0.5243*10.0**(-2)*(two*pi)**2*COS(two*pi*X)+&
            SQRT(SQRT(three))*0.4978*10.0**(-3)*(three*pi)**2*COS(three*pi*X)+&
            SQRT(SQRT(four))*0.6542*10.0**(-4)*(four*pi)**2*COS(four*pi*X)+&
            SQRT(SQRT(five))*0.1007*10.0**(-4)*(five*pi)**2*COS(five*pi*X)+&
            SQRT(SQRT(six))*0.1653*10.0**(-5)*(six*pi)**2*COS(six*pi*X)+&
            SQRT(SQRT(seven))*0.2753*10.0**(-6)*(seven*pi)**2*COS(seven*pi*X)+&
            SQRT(SQRT(eight))*0.4522*10.0**(-7)*(eight*pi)**2*COS(eight*pi*X))
uu = zero;
ww = zero;
END SUBROUTINE nonlinearstandingwave1D
