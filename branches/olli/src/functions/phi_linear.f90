SUBROUTINE phi_linear(H,k,x,z,d,omega,time,grav,phix,phiz,phi,phit)
!
USE PRECISION
IMPLICIT NONE
REAL(KIND=long), INTENT(IN)  :: H,k,x,z,d,omega,time,grav
REAL(KIND=long), INTENT(OUT) :: phix,phiz
REAL(KIND=long), INTENT(OUT), OPTIONAL :: phi,phit
   ! Local variables
REAL(KIND=long) :: kx, sinkx, temp1, temp2
!
kx=k*x-omega*time
sinkx = sin(kx)
!
temp1 = COSH(k*(z+d))/COSH(k*d)
temp2 = H*0.5d0*grav/omega
!
phiz = temp2*sinkx*temp1*k*tanh(k*(z+d))
phix = temp2*k*cos(kx)*temp1
!
IF(present(phit)) THEN
    phit = -phix*omega/k !we assume that k/=0
ENDIF
IF(present(phi)) THEN
    phi  =  temp2*sinkx*temp1
ENDIF
END SUBROUTINE phi_linear
