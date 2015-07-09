C
C-----------------------------------------------------------------------
      SUBROUTINE INIT(N,NUM,PI,HOVERD,HEIGHT,VALUE,DEPTH,CASE,CURRNT,
     $           Z,COSA,SINA,COEFF,SOL,Y)  
C-----------------------------------------------------------------------
C
C     Calculate initial solution from linear wave theory
C
C-----------------------------------------------------------------------
C
      implicit none
      CHARACTER*10 DEPTH,CASE,CURRNT
      INTEGER n,num
      REAL*8 PI,HOVERD,HEIGHT,VALUE,Z(num),COSA(0:num),SINA(0:num), 
     $     COEFF(num),SOL(num,2),Y(num)

      integer :: i
      real*8 :: t, b, a
C
      IF (DEPTH.EQ.'finite') THEN
         IF (CASE.EQ.'period') THEN
            A=4.D0*PI*PI*HEIGHT/HOVERD
            B=A/DSQRT(DTANH(A))
            T=DTANH(B)
            Z(1)=B+(A-B*T)/(T+B*(1.D0-T*T))
         ELSE
            Z(1)=2.D0*PI*HEIGHT/HOVERD
         ENDIF
         Z(2)=Z(1)*HOVERD
         Z(4)=DSQRT(DTANH(Z(1)))
      ELSE
         Z(1)=-1.D0
         Z(4)=1.D0
         IF (CASE.EQ.'period') THEN
            Z(2)=4.D0*PI*PI*HEIGHT
         ELSE
            Z(2)=2.D0*PI*HEIGHT
         ENDIF
      ENDIF
      Z(3)=2.D0*PI/Z(4)
      IF (CURRNT.EQ.'Euler') THEN
         Z(5)=VALUE*DSQRT(Z(2))
         Z(6)=0.D0
      ELSE
         Z(5)=0.D0
         Z(6)=VALUE*DSQRT(Z(2))
      ENDIF
      Z(7)=Z(4)
      Z(8)=0.D0
      Z(9)=0.5D0*Z(7)**2
      Z(10)=0.5D0*Z(2)
      COSA(0)=1.D0
      SINA(0)=0.D0
      DO 1 I=1,N
         COSA(I)=DCOS(I*PI/N)
         COSA(I+N)=DCOS((I+N)*PI/N)
         SINA(I)=DSIN(I*PI/N)
         SINA(I+N)=DSIN((I+N)*PI/N)
         Z(N+I+10)=0.D0
         Z(I+10)=0.5D0*Z(2)*COSA(I)
    1 CONTINUE
      Z(N+11)=0.5D0*Z(2)/Z(7)
C
    2 FORMAT(//,'   INITIAL LINEAR SOLUTION',10(/,6E13.6))
      DO 3 I=1,9
    3     SOL(I,1)=Z(I)
      DO 4 I=10,NUM
    4    SOL(I,1)=0.D0
      DO 5 I=1,NUM
    5    SOL(I,2)=0.D0
      RETURN
      END
