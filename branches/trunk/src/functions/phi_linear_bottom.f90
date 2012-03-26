SUBROUTINE phi_linear_bottom(H,k,x,z,d,hx,hxx,int_kx,omega,time,grav,phix,phiz,phi,phit)
!
USE PRECISION
IMPLICIT NONE
REAL(KIND=long), INTENT(IN)  :: H,k,x,z,d,omega,time,grav,int_kx,hx,hxx
REAL(KIND=long), INTENT(OUT) :: phix,phiz
REAL(KIND=long), INTENT(OUT), OPTIONAL :: phi,phit
   ! Local variables
REAL(KIND=long) :: kx, sinkx, coskx, temp1, temp2, temp3
!
kx=int_kx-omega*time
sinkx = sin(kx)
coskx = cos(kx)
!
temp1 = COSH(k*(z+d))/COSH(k*d)
temp3 = SINH(k*(z+d))/COSH(k*d)
temp2 = H*0.5d0*grav/omega
!
!phiz = temp2*(sinkx*temp1*k*tanh(k*(z+d)))
phiz = temp2*(-sinkx*temp3*k-hx*coskx*temp1*k)
!phix = temp2*k*coskx*temp1
phix = temp2*(-k*coskx*temp1+k*hx*sinkx*temp3-hxx*coskx*temp3)
!
IF(present(phit)) THEN
    !phit = -phix*omega/k !we assume that k/=0
    phit = temp2*omega*(-coskx*temp1+hx*sinkx*temp3)
ENDIF
IF(present(phi)) THEN
    phi  =  temp2*(-sinkx*temp1-hx*coskx*temp3)
ENDIF
END SUBROUTINE phi_linear_bottom
