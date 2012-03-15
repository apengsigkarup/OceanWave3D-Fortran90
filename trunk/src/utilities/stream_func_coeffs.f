      SUBROUTINE stream_func_coeffs(depth,hoverd,case,height,currnt,
     $     value,n,nstep,num,cosa,sina,coeff,sol,rhs1,rhs2,a,b,ipvt,z,y)
C-----------------------------------------------------------------------
C
C     Calculation of steady waves
C
C-----------------------------------------------------------------------

      IMPLICIT none
C Arguements
      CHARACTER*10 DEPTH,CASE,CURRNT
      integer  n, nstep, num
      integer IPVT(num)
      real*8  HOVERD,HEIGHT,VALUE
      real*8  Z(num),COSA(0:num),SINA(0:num),COEFF(num),SOL(num,2),
     $     Y(num),RHS1(num),RHS2(num),A(num,num),B(num)


      integer  j, iter, i, ns, m, info, number
      real*8  PI, h, criter, crit, dho, dhe, sum

C.....Input data-----
C
C     "depth" IS EITHER 'deep' OR 'finite'
C     "hoverd" IS WAVE HEIGHT/DEPTH
C
C     "case" IS EITHER 'period' OR 'wavelength'
C     "height" IS HEIGHT/LENGTH IF "case" IS 'wavelength'
C     "height" IS HEIGHT/(g*T**2) IF "case" IS 'period'
C
C     "currnt" IS EITHER 'Euler' OR 'Stokes'
C     "value" is the magnitude of the mean Eulerian or Stokes velocities
C             non-dimensionalized with respect to wave height
C
C
C     "n" is the number of terms in the Fourier series and the number
C         of intervals in half a wavelength
C     "nstep" is the number of steps in wave height
C
C     "number" is the number of iterations for each height step
C
      NUMBER=256
C
C     "crit" is the criterion for convergence. If the sum of magnitudes
C            of corrections is smaller than crit, the iteration stops
C      
      CRIT=1.D-6 !Originally 10^-4 GD
!      WRITE(6,20) DEPTH,HOVERD
!      WRITE(6,21) HEIGHT,CASE
!      WRITE(6,22) CURRNT,VALUE
!      WRITE(7,20) DEPTH,HOVERD
!      WRITE(7,21) HEIGHT,CASE
!      WRITE(7,22) CURRNT,VALUE
      PI=4.D0*DATAN(1.D0)
      DHE=HEIGHT/NSTEP
      DHO=HOVERD/NSTEP
C
C.....Commence stepping through steps in wave height
C
      DO 1 NS=1,NSTEP
         HEIGHT=NS*DHE
         HOVERD=NS*DHO
C
C.....Calculate initial linear solution
C
         IF(NS.LE.1) THEN
            CALL INIT(N,NUM,PI,HOVERD,HEIGHT,VALUE,DEPTH,CASE,CURRNT,
     $           Z,COSA,SINA,COEFF,SOL,Y)  
         ELSE
C
C.....or extrapolate for next wave height, is neccessary
C
            DO 3 I=1,NUM
    3          Z(I)=2.D0*SOL(I,2)-SOL(I,1)
         ENDIF
C
C.....Commence iterative solution
C
         DO 4 ITER=1,NUMBER
C
C.....Calculate right sides of equations and differentiate numerically
C.....to obtain Jacobian matrix
C
            CALL EQNS(N,NUM,PI,HOVERD,HEIGHT,VALUE,DEPTH,CASE,CURRNT,
     $           Z,COSA,SINA,COEFF,SOL,Y,RHS1)
            DO 5 I=1,NUM
               H=0.01D0*Z(I)
               IF(DABS(Z(I)).LT.1.D-4) H=1.D-5
               Z(I)=Z(I)+H
               CALL EQNS(N,NUM,PI,HOVERD,HEIGHT,VALUE,DEPTH,CASE,CURRNT,
     $           Z,COSA,SINA,COEFF,SOL,Y,RHS2)
               Z(I)=Z(I)-H
               B(I)=-RHS1(I)
               DO 6 J=1,NUM
    6             A(J,I)=(RHS2(J)-RHS1(J))/H
    5       CONTINUE
C
C.....Solve matrix equation and correct variables, using LINPACK 
C
C            CALL DGEFA(A,NUM,NUM,IPVT,INFO)
C            IF (INFO.NE.0) THEN
C               WRITE(6,27)
C               WRITE(7,27)
C               STOP
C            ENDIF
C            CALL DGESL(A,NUM,NUM,IPVT,B,0)
C.....Employ LAPACK instead... APEK 17-01-2009.
             CALL DGESV(NUM,1,A,NUM,IPVT,B,NUM,INFO)
C
C.....The b(i) are now corrections to each variable
C
            SUM=0.D0
            DO 7 I=1,NUM
               SUM=SUM+DABS(B(I))
    7          Z(I)=Z(I)+B(I)
            CRITER=CRIT
            IF(NS.EQ.NSTEP) CRITER=0.01D0*CRIT
            IF(SUM.LT.CRITER) GOTO 8
    4    CONTINUE
         WRITE(6,26) NUMBER
!         WRITE(7,26) NUMBER
         STOP
    8    IF(NS.EQ.1) THEN
            DO 9 I=1,NUM
    9          SOL(I,2)=Z(I)
         ELSE
            DO 10 I=1,NUM
               SOL(I,1)=SOL(I,2)
   10          SOL(I,2)=Z(I)
         ENDIF
    1 CONTINUE
C
C.....Calculate Fourier coefficients of surface elevation
C
      DO  J=1,N
         SUM=0.5D0*(Z(10)+Z(N+10)*(-1.D0)**J)
         DO  M=1,N-1
            SUM=SUM+Z(10+M)*COSA(MOD(M*J,N+N))
         END DO
         Y(J)=2.D0*SUM/N
      END DO

   20 FORMAT(//,'   Depth: ',A6,', Height/Depth',F7.4)
   21 FORMAT(/,'   Wave height',F9.6,', dimensionless with respect to ',
     .         A10)
   22 FORMAT(/,'   Current criterion ',A6,', Magnitude ', F5.2)
   23 FORMAT(//,'   Height step ',I2,' of ',I2)
   24 FORMAT(//,'   Iteration ',I3)
   25 FORMAT(//,'   Solution vector',10(/,6E13.6))
   26 FORMAT(/,' did not converge sufficiently after ',I3,' iterations')
   27 FORMAT(/,'   MATRIX SINGULAR')
      RETURN
      END
