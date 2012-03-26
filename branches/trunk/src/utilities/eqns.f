C
C-----------------------------------------------------------------------
      SUBROUTINE EQNS(N,NUM,PI,HOVERD,HEIGHT,VALUE,DEPTH,CASE,CURRNT,
     $           Z,COSA,SINA,COEFF,SOL,Y,RHS)
C-----------------------------------------------------------------------
C
C     Evaluation of equations
C
C-----------------------------------------------------------------------
C
      implicit none

      CHARACTER*10 DEPTH,CASE,CURRNT
      INTEGER n,num
      REAL*8 PI,HOVERD,HEIGHT,VALUE,Z(num),COSA(0:num),SINA(0:num), 
     $     COEFF(num),SOL(num,2),Y(num)
c      IMPLICIT REAL*8 (A-H,K,L,O-Z)
      real*8 :: RHS(num) !GD change: RHS(*)
C      real*8 :: RHS(*) 
      integer :: i, it, nm, m, j
      real*8 :: e, c, s, psi, u, v
C
      IF(DEPTH.EQ.'finite') THEN
         RHS(1)=Z(2)-Z(1)*HOVERD
      ELSE
         RHS(1)=Z(1)+1.D0
      ENDIF
      IF(CASE.EQ.'wavelength') THEN
         RHS(2)=Z(2)-2.D0*PI*HEIGHT
      ELSE
         RHS(2)=Z(2)-HEIGHT*Z(3)**2
      ENDIF
      RHS(3)=Z(4)*Z(3)-PI-PI
      RHS(4)=Z(5)+Z(7)-Z(4)
      RHS(5)=Z(6)+Z(7)-Z(4)
      IF(DEPTH.EQ.'finite') THEN
         RHS(5)=RHS(5)-Z(8)/Z(1)
         DO 2 I=1,N
    2       COEFF(I)=Z(N+I+10)/DCOSH(I*Z(1))        
      ENDIF
      IT=6
      IF(CURRNT.EQ.'Euler') IT=5
      RHS(6)=Z(IT)-VALUE*DSQRT(Z(2))
      RHS(7)=Z(10)+Z(N+10)
      DO 1 I=1,N-1
    1   RHS(7)=RHS(7)+Z(10+I)+Z(10+I)
      RHS(8)=Z(10)-Z(N+10)-Z(2)
      DO 3 M=0,N
         PSI=0.D0
         U=0.D0
         V=0.D0
         IF(DEPTH.EQ.'finite') THEN
            DO 4 J=1,N
               NM=MOD(M*J,N+N)
               E=DEXP(J*(Z(1)+Z(10+M)))
               S=0.5D0*(E-1.D0/E)
               C=0.5D0*(E+1.D0/E)
               PSI=PSI+COEFF(J)*S*COSA(NM)
               U=U+J*COEFF(J)*C*COSA(NM)
               V=V+J*COEFF(J)*S*SINA(NM)
    4       CONTINUE
         ELSE
            DO 5 J=1,N
               NM=MOD(M*J,N+N)
               E=DEXP(J*Z(10+M))
               PSI=PSI+Z(N+J+10)*E*COSA(NM)
               U=U+J*Z(N+J+10)*E*COSA(NM)
    5          V=V+J*Z(N+J+10)*E*SINA(NM)
         ENDIF
         RHS(M+9)=PSI-Z(8)-Z(7)*Z(M+10)
         RHS(N+M+10)=0.5D0*((-Z(7)+U)**2+V**2)+Z(M+10)-Z(9)
    3 CONTINUE
      RETURN
      END
