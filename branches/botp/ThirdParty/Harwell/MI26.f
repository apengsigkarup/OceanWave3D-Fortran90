C COPYRIGHT (c) 1995 Council for the Central Laboratory
*                    of the Research Councils
C######DATE 29 March 2001
C  March 2001: threadsafe version of MI06 and extend control arguments

      SUBROUTINE MI26ID(ICNTL,CNTL,ISAVE,RSAVE)
C
C  MI26 solves the linear system Ax = b using the
C  BiConjugate Gradient Stabilized iterative method with
C  optionally using preconditioning
C                    ^                   ^
C           P(L)AP(R)x = P(L)b,  x = P(R)x

C  P(L), P(R) are the preconditioners, which are not passed to the code,

C  but each time a matrix-vector product with P = P(L)P(R)
C  is required, control is passed back to the user.
C
C  Similarly, the matrix A is not passed to the code, but when
C  a matrix-vector with A is required, control is
C  passed back to the user.
C
C  MI26I/ID is the initialisation routine for MI26A/AD and should
C  be called once prior to calls to MI26A/AD.
C
C  Argument list.
C
C  ICNTL   (output) INTEGER control array, dimension 8.
C          ICNTL(1) is the stream number for error messages.
C          On exit, ICNTL(2) = 6.
C          ICNTL(2) is the stream number for warning messages.
C          On exit, ICNTL(2) = 6.
C          ICNTL(3) indicates whether the user wishes to use a
C          preconditioner. If ICNTL(3) is nonzero, preconditioning.
C          On exit, ICNTL(3) = 0
C          ICNTL(4) indicates whether
C          the convergence test offered by MI26A/AD is to be used.
C          If ICNTL(4) = 0, the computed solution x is accepted
C          if ||Ax - b||/||r sub 0|| < CNTL(1) (r sub 0 is
C          initial residual).
C          Otherwise, the user may perform his/her
C          own test for convergence when IACT = 1 is returned.
C          On exit, ICNTL(4) = 0
C          ICNTL(5) indicates whether the user wishes to supply an
C          initial guess for the solution vector x.
C          If ICNTL(5) = 0, the user does not wish to supply
C          an initial guess and x = (0,0,...,0) will be used
C          as the initial guess. Otherwise, the user
C          must supply an intial guess on the first call to
C          MI26A/AD. On exit, ICNTL(5) = 0.
C          ICNTL(6) determines the maximum number of iterations
C          allowed. It has default value -1 and, in this case,
C          the maximum number will be N. If the user does
C          not want the maximum number to be N, ICNTL(6) should
C          be set to the maximum  number the user wishes
C          to allow.
C          ICNTRL(7) and ICNTRL(8) are spare and set to zero.
C
C  CNTL    (output)  REAL (DOUBLE PRECISION) control array, dimension 5.

C          CNTL(1) and CNTL(2) are convergence tolerances.
C          On exit, set to square root of machine precision and 0,
C          respectively.
C          If ICNTL(4) is nonzero, CNTL(1) and CNTL(2) are
C          not accessed by MI26A/AD.
C          CNTL(3) is tolerance used to check whether the algorithm
C          has broken down. On exit, set to machine precision
C          CNTRL(4) and CNTRL(5) are spare and set to zero.
C
C  ISAVE   (output) INTEGER ARRAY, length 14, used to hold the
C          routine's persistent integer data.  The array contents
C          must not be altered by the user.
C
C  RSAVE   (output) DOUBLE PRECISION ARRAY, length 9, used to hold the
C          routine's persistent real data.  The array contents
C          must not be altered by the user.
C
C     .. Array Arguments ..
      DOUBLE PRECISION CNTL(5)
      INTEGER ICNTL(8)
      INTEGER ISAVE(14)
      DOUBLE PRECISION RSAVE(9)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION FD05AD
      EXTERNAL FD05AD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 0
      ICNTL(4) = 0
      ICNTL(5) = 0
      ICNTL(6) = -1

      ICNTL(7) = 0
      ICNTL(8) = 0

      CNTL(1) = SQRT(FD05AD(1))
      CNTL(2) = ZERO
      CNTL(3) = FD05AD(1)

      CNTL(4) = ZERO
      CNTL(5) = ZERO

C  Initialize persistent data to avoid undefined assignment
      DO 10 I = 1, 14
      ISAVE(I) = 0
   10 CONTINUE
      DO 20 I = 1, 9
      RSAVE(I) = 0.0
   20 CONTINUE

      RETURN
      END
C
      SUBROUTINE MI26AD(IACT,N,W,LDW,LOCY,LOCZ,RESID,ICNTL,CNTL,INFO,
     +                 ISAVE,RSAVE)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION RESID
      INTEGER IACT,LDW,LOCY,LOCZ,N
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CNTL(5),W(LDW,8)
      INTEGER ICNTL(8),INFO(4)
      INTEGER ISAVE(14)
      DOUBLE PRECISION RSAVE(9)
C     ..
C
C  Argument list.
C
C  IACT    (input) INTEGER.
C          IACT must be set to 0 prior to first call.
C          On each exit, IACT indicates the action required by
C          the user. Possible values of IACT and the action
C          required are as follows:
C  -1      Fatal error (see INFO(1)). Terminate computation.
C   1      If ICNTL(4) = 0 (the default), convergence has been
C          achieved and the user should terminate the computation.
C          If ICNTL(4) is nonzero, the user should test the norm of
C          the residual in RESID for convergence and recall MI26A/AD
C          if convergence has not been achieved.
C   2      The user must perform the matrix-vector product
C          y := Az ,
C          and recall MI26A/AD. The vectors y and z
C          are  held in the first N entries of
C          columns LOCY and LOCZ of array W, respectively.
C   3      The user must perform the preconditioning operations
C          y := Pz,
C          where P = P sub L P sub R  is the preconditioner
C          and recall MI26A/AD. The vectors y and z
C          are  held in the first N entries of
C          columns LOCY and LOCZ of array W, respectively.
C
C  N       (input) INTEGER.
C          On entry, the dimension of the matrix.
C          Unchanged on exit.
C
C  W       (right-hand side, solution and workspace, input/output)
C          REAL (DOUBLE PRECISION) array, dimension (LDW,8).
C          Prior to the first call, the first N entries of column
C          1 must be set to hold the right-hand side vector b and,
C          if ICNTL(5) is nonzero, the first N entries of column 2
C          must be set to the initial estimate of the solution vector
C          x.  On exit with IACT = 1, the first N entries of column 1
C          hold the current residual vector r = b - Ax, and the
C          current estimate of the solution x is held in the
C          first N entries of column 2.  On exit
C          with IACT > 1, the user is required to perform a
C          computation with columns LOCZ and LOCZ of W
C          (see argument IACT).  The remaining columns of
C          W must not be altered by the user between calls to
C          MI26A/AD.
C
C  LDW     (input) INTEGER
C          The leading dimension of the array W. LDW >= max(1,N).
C
C LOCY, LOCZ (output) INTEGER  variables
C          On exit with IACT > 1, LOCY and LOCZ define
C          the columns of the array W which hold y and z.
C          (see IACT).
C
C  RESID   (output) REAL (DOUBLE PRECISION)
C          On exit with IACT = 1,
C          RESID holds ||b - Ax||, where x is the
C          iterated solution.
C          If ICNTL(4) is nonzero, on exit with IACT = 1,
C          the user may carry out his/her
C          own test for convergnce at this point.
C
C  ICNTL   (input) INTEGER control array of dimension 8.
C          ICNTL may be initialised by calling MI26I/ID.
C          See MI26I/ID for details.
C
C  CNTL    (input) REAL (DOUBLE PRECISION) control array of dimension 5.

C          CNTL may be initialised by calling MI26I/ID.
C          See MI26I/ID for details.
C
C  INFO    (output) INTEGER ARRAY , length 4.
C          If INFO(1) = 0 on exit, no errors or warnings issued.
C          If INFO(1) > 0 on exit, warning to user.
C          INFO(1) = 1, value of CNTL(1) is out-of-range (u,1.0)
C          The default sqrt(u) is used, u=machine precision.
C          If INFO(1) < 0 on exit, illegal input parameter,
C               or breakdown occured during iteration.
C
C                Error input data:
C
C                   -1: matrix dimension N < 0
C                   -2: LDW < N
C
C                BREAKDOWN: If RHO becomes too small,
C                   the program will terminate.
C
C                   -3: RHO < CNTL(3)*||R|*||RTLD||
C          R and RTLD have become  orthogonal.
C          Error -3 also returned if OMEGA becomes too small
C          (S and T have become orthogonal relative to T'*T).
C
C          On each exit, INFO(2) holds the number of
C          iterations performed so far.
C
C  ISAVE   (output) INTEGER ARRAY, length 14, used to hold the
C          routine's persistent integer data.  The array contents
C          must not be altered by the user.
C
C  RSAVE   (output) DOUBLE PRECISION ARRAY, length 9, used to hold the
C          routine's persistent real data.  The array contents
C          must not be altered by the user.
C
C  BLAS CALLS: DAXPY, DCOPY, DDOT, DNRM2, DSCAL
C***************************************************
C
C     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALPHA,BETA,BNRM2,OMEGA,RHO,RHO1,RNRM2,RTNRM2,
     +                 RSTOP,SNORM2,TNORM2
      INTEGER B,I,IPOS,ITMAX,P,PHAT,R,RTLD,S,SHAT,T,V,X
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT,DNRM2,FD05AD
      EXTERNAL DDOT,DNRM2,FD05AD
C     ..
C     .. Intrinsic Functions ..

      INTRINSIC ABS,MAX,SQRT
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DCOPY,DSCAL
C
C  Restore persistent data
C
      IPOS   = ISAVE(1)
      ITMAX  = ISAVE(2)
      B      = ISAVE(3)
      R      = ISAVE(4)
      X      = ISAVE(5)
      P      = ISAVE(6)
      S      = ISAVE(7)
      T      = ISAVE(8)
      V      = ISAVE(9)
      PHAT   = ISAVE(10)
      RTLD   = ISAVE(11)
      SHAT   = ISAVE(12)

      BNRM2  = RSAVE(1)
      ALPHA  = RSAVE(2)
      BETA   = RSAVE(3)
      RHO    = RSAVE(4)
      RHO1   = RSAVE(5)
      RSTOP  = RSAVE(6)
      OMEGA  = RSAVE(7)

C Jump to appropriate place in code
      IF (IACT.EQ.0) GO TO 10
C Immediate return if error on a previous call
      IF (IACT.LT.0) GO TO 1000
C Immediate return if convergence already achieved
      IF (IACT.EQ.1 .AND. ICNTL(4).EQ.0) GO TO 1000
      IF (IACT.EQ.1 .AND. BNRM2.EQ.ZERO) GO TO 1000
C
      IF (IPOS.EQ.1) GO TO 40
      IF (IPOS.EQ.2) GO TO 70
      IF (IPOS.EQ.3) GO TO 80
      IF (IPOS.EQ.4) GO TO 90
      IF (IPOS.EQ.5) GO TO 100
      IF (IPOS.EQ.6) GO TO 110
      IF (IPOS.EQ.7) GO TO 120


   10 CONTINUE
C
C  Initial call.
      INFO(1) = 0
C  Test the input parameters.
      IF (N.LE.0) THEN
         INFO(1) = -1
      ELSE IF (LDW.LT.MAX(1,N)) THEN
         INFO(1) = -2
      END IF
      IF (INFO(1).LT.0) THEN
         IACT = -1
         IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9000) INFO(1)
         GO TO 1000
      END IF

C
C Alias workspace columns.
C
      B = 1
      X = 2
      R = 1
      RTLD = 3
      P = 4
      V = 5
      T = 6
      PHAT = 7
      SHAT = 8
      S = 1
C Set INFO(2) and ITMAX.

      INFO(2) = 0
      ITMAX = N
      IF (ICNTL(6).GT.0) ITMAX = ICNTL(6)

C Compute ||b||
C
      BNRM2 = DNRM2(N,W(1,B),1)
C Immediate return if ||b|| = 0.
      IF (BNRM2.EQ.ZERO) THEN
         IACT = 1
         DO 20 I = 1,N
            W(I,X) = ZERO
            W(I,B) = ZERO
   20    CONTINUE
         RESID = ZERO
         GO TO 1000
      END IF
C
      IF (ICNTL(4).EQ.0) THEN
C Check value of CNTL(1)
         IF (CNTL(1).LT.FD05AD(1) .OR. CNTL(1).GT.ONE) THEN
            INFO(1) = 1
            IF (ICNTL(2).GT.0) THEN
               WRITE (ICNTL(2),FMT=9010) INFO(1)
               WRITE (ICNTL(2),FMT=9020)
            END IF
            CNTL(1) = SQRT(FD05AD(1))
         END IF
      END IF
C
C Compute initial residual.
C
C If the user has not supplied an initial guess, set x = 0
C as the initial guess.
      IF (ICNTL(5).EQ.0) THEN
         DO 30 I = 1,N
            W(I,X) = ZERO
   30    CONTINUE
         GO TO 50
      ELSE
C Initial guess supplied by user
C If initial guess for solution is x = 0 no action is required.
C (r = b)
         IF (DNRM2(N,W(1,X),1).EQ.ZERO) GO TO 50
C
C Otherwise, return to user to compute Ax.
C Column P can be used temporarily to hold Ax.
         IPOS = 1
         IACT = 2
         LOCY = P
         LOCZ = X
         GO TO 1000
      END IF
C
C Compute r = b - Ax
   40 CALL DAXPY(N,-ONE,W(1,P),1,W(1,R),1)

   50 CONTINUE
C
C  Compute the norm of the initial residual
C
      RSTOP = DNRM2(N,W(1,R),1)
C
C Choose RTLD such that initially, (R,RTLD) = RHO is not equal to 0.
C Here we choose RTLD = R.
C
      CALL DCOPY(N,W(1,R),1,W(1,RTLD),1)

C Perform BiConjugate Gradient Stabilized iteration.

   60 CONTINUE

C Update iteration count

      INFO(2) = INFO(2) + 1
C
Check maximum number of iterations has not been exceeded.
      IF (INFO(2).GT.ITMAX) THEN
         INFO(1) = -4
         IACT = -1
         IF (ICNTL(1).GT.0) THEN
            WRITE (ICNTL(1),FMT=9000) INFO(1)
            WRITE (ICNTL(1),FMT=9030) ITMAX
         END IF
         GO TO 1000
      END IF


      RHO = DDOT(N,W(1,RTLD),1,W(1,R),1)
C Check for breakdown.
      IF (ABS(RHO).LT.CNTL(3)*N) THEN
C RHO getting small. Carry out more rigorous test
         RNRM2 = DNRM2(N,W(1,R),1)
         RTNRM2 = DNRM2(N,W(1,RTLD),1)
         IF (ABS(RHO).LT.CNTL(3)*RNRM2*RTNRM2) THEN
            INFO(1) = -3
            IACT = -1
            IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9000) INFO(1)
            GO TO 1000
         END IF
      END IF
C
C Compute vector P.
C
      IF (INFO(2).GT.1) THEN
         BETA = (RHO/RHO1)* (ALPHA/OMEGA)
         CALL DAXPY(N,-OMEGA,W(1,V),1,W(1,P),1)
         CALL DSCAL(N,BETA,W(1,P),1)
         CALL DAXPY(N,ONE,W(1,R),1,W(1,P),1)
      ELSE
         CALL DCOPY(N,W(1,R),1,W(1,P),1)
      END IF
C
C Compute direction adjusting vector PHAT and scalar ALPHA.
C
      IF (ICNTL(3).NE.0) THEN
C Return to user for preconditioning.
         IPOS = 2
         IACT = 3
         LOCY = PHAT
         LOCZ = P
         GO TO 1000
      ELSE
C No preconditioning (i.e identity preconditioning)
         CALL DCOPY(N,W(1,P),1,W(1,PHAT),1)
      END IF

   70 CONTINUE
C
C Return to user for matrix-vector product
      IPOS = 3
      IACT = 2
      LOCY = V
      LOCZ = PHAT
      GO TO 1000

   80 CONTINUE

      ALPHA = RHO/DDOT(N,W(1,RTLD),1,W(1,V),1)
C
C Early check for tolerance.
C
      CALL DAXPY(N,-ALPHA,W(1,V),1,W(1,R),1)

C Note: R=1 and S=1, so no need to copy R into S
C        CALL DCOPY( N, W(1,R), 1, W(1,S), 1 )

      CALL DAXPY(N,ALPHA,W(1,PHAT),1,W(1,X),1)

      RESID = DNRM2(N,W(1,S),1)
      IPOS = 4
      IF (ICNTL(4).NE.0) THEN
C Return the residual to the user for convergence testing.
         IACT = 1
         GO TO 1000
      ELSE
C Test the scaled residual for convergence.
         IF (RESID.LE.MAX(CNTL(2),RSTOP*CNTL(1))) THEN
C Convergence achieved
            IACT = 1
            GO TO 1000
         END IF
      END IF
C
   90 CONTINUE

C
C  Compute stabilizer vector SHAT and scalar OMEGA.
C
      IF (ICNTL(3).NE.0) THEN
C Return to user for preconditioning.
         IPOS = 5
         IACT = 3
         LOCY = SHAT
         LOCZ = S
         GO TO 1000
      ELSE
C No preconditioning (i.e identity preconditioning)
         CALL DCOPY(N,W(1,S),1,W(1,SHAT),1)
      END IF

  100 CONTINUE

C Return to user for matrix-vector product
      IPOS = 6
      IACT = 2
      LOCY = T
      LOCZ = SHAT
      GO TO 1000

  110 CONTINUE

      OMEGA = DDOT(N,W(1,T),1,W(1,S),1)/DDOT(N,W(1,T),1,W(1,T),1)

C Check OMEGA is not too small.
      IF (ABS(OMEGA).LT.CNTL(3)*N) THEN
C OMEGA getting small. Carry out more rigorous test
         SNORM2 = DNRM2(N,W(1,S),1)
         TNORM2 = DNRM2(N,W(1,T),1)
         IF (ABS(RHO).LT.CNTL(3)*SNORM2/TNORM2) THEN
            INFO(1) = -3
            IACT = -1
            IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9000) INFO(1)
            GO TO 1000
         END IF
      END IF
C
C Compute new solution approximation vector X.
C
      CALL DAXPY(N,OMEGA,W(1,SHAT),1,W(1,X),1)

C  Compute residual R, check for tolerance.
C
      CALL DAXPY(N,-OMEGA,W(1,T),1,W(1,R),1)

      RESID = DNRM2(N,W(1,R),1)
      IPOS = 7
      IF (ICNTL(4).NE.0) THEN
C Return the residual to the user for convergence testing.
         IACT = 1
         GO TO 1000
      ELSE
C Test the scaled residual for convergence.
         IF (RESID.LE.MAX(CNTL(2),RSTOP*CNTL(1))) THEN
C Convergence achieved
            IACT = 1
            GO TO 1000
         END IF
      END IF
C
  120 CONTINUE

      RHO1 = RHO

C Next iteration
      GO TO 60
C
C  Save persistent data and return to caller
C
 1000 CONTINUE
      ISAVE(1)  = IPOS
      ISAVE(2)  = ITMAX
      ISAVE(3)  = B
      ISAVE(4)  = R
      ISAVE(5)  = X
      ISAVE(6)  = P
      ISAVE(7)  = S
      ISAVE(8)  = T
      ISAVE(9)  = V
      ISAVE(10) = PHAT
      ISAVE(11) = RTLD
      ISAVE(12) = SHAT

      RSAVE(1)  = BNRM2
      RSAVE(2)  = ALPHA
      RSAVE(3)  = BETA
      RSAVE(4)  = RHO
      RSAVE(5)  = RHO1
      RSAVE(6)  = RSTOP
      RSAVE(7)  = OMEGA
      RETURN
C
C End of MI26A/AD
C
 9000 FORMAT (/' Error message from MI26A/AD. INFO(1) = ',I4)
 9010 FORMAT (/' Warning message from MI26A/AD. INFO(1) = ',I4)
 9020 FORMAT (' Convergence tolerance out of range.')
 9030 FORMAT (' Number of iterations required exceeds the maximum of ',
     +       I8,/' allowed by ICNTL(6)')
      END
