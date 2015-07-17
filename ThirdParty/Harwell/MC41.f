* COPYRIGHT (c) 1993 AEA Technology
*######DATE 6 Jan 1993
C       Toolpack tool decs employed.
C       ZERO and ONE made PARAMETER.
C
      SUBROUTINE MC41AD(N,KASE,X,EST,W,IW)
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION EST
      INTEGER KASE,N
      DOUBLE PRECISION W(*),X(*)
      INTEGER IW(*)
      DOUBLE PRECISION ALTSGN,TEMP
      INTEGER I,ITER,J,JLAST,JUMP
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      INTRINSIC DABS,DBLE,DSIGN,IDNINT
      SAVE ITER,J,JLAST,JUMP
      IF (N.LE.0) THEN
        KASE = -1
        RETURN
      END IF
      IF (KASE.EQ.0) THEN
        DO 10 I = 1,N
          X(I) = ONE/DBLE(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        RETURN
      END IF
      GO TO (100,200,300,400,500) JUMP
  100 CONTINUE
      IF (N.EQ.1) THEN
        W(1) = X(1)
        EST = DABS(W(1))
        GO TO 510
      END IF
      DO 110 I = 1,N
        X(I) = DSIGN(ONE,X(I))
        IW(I) = IDNINT(X(I))
  110 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
  200 CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
  220 CONTINUE
      DO 230 I = 1,N
        X(I) = ZERO
  230 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      RETURN
  300 CONTINUE
      DO 310 I = 1,N
        W(I) = X(I)
  310 CONTINUE
      DO 320 I = 1,N
        IF (IDNINT(DSIGN(ONE,X(I))).NE.IW(I)) GO TO 330
  320 CONTINUE
      GO TO 410
  330 CONTINUE
      DO 340 I = 1,N
        X(I) = DSIGN(ONE,X(I))
        IW(I) = IDNINT(X(I))
  340 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
  400 CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF ((DABS(X(JLAST)).NE.DABS(X(J))) .AND. (ITER.LT.ITMAX)) THEN
        ITER = ITER + 1
        GO TO 220
      END IF
  410 CONTINUE
      EST = ZERO
      DO 420 I = 1,N
        EST = EST + DABS(W(I))
  420 CONTINUE
      ALTSGN = ONE
      DO 430 I = 1,N
        X(I) = ALTSGN* (ONE+DBLE(I-1)/DBLE(N-1))
        ALTSGN = -ALTSGN
  430 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
  500 CONTINUE
      TEMP = ZERO
      DO 520 I = 1,N
        TEMP = TEMP + DABS(X(I))
  520 CONTINUE
      TEMP = 2.0*TEMP/DBLE(3*N)
      IF (TEMP.GT.EST) THEN
        DO 530 I = 1,N
          W(I) = X(I)
  530   CONTINUE
        EST = TEMP
      END IF
  510 KASE = 0
      RETURN
      END


