C COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
C######DATE 20 May 1993
C 24 September 1993. Some IVDEP comments added for speed on the Cray,
C     some minor bugs fixed, defaults changed to BLAS 3 with block size
C     32 and no amalgamation of blocks of the block triangular form.
C 18 October 1993. Minor bug in MA48BD when ICNTL(3)=3 corrected.
C 10 December 1993. Minor bugs in printing corrected.
C 14 June 1994. Minor bug rank calculation corrected.
C 14 June 1994. Minor bugs in printing corrected.
C 12/12/94 Calls of MC13D and MC21A changed to MC13DD and MC21AD
C 4/10/95. Redundant variables ONE removed.
C 1/11/95  IFLAG in MA48DD made of length 7.
C 1/11/95  Temporary variable NP introduced to avoid PFORT warnings.
C 14/2/96  Third argument in calls to MA48DD changed to LA-NE.
C 13/8/96  Bug corrected in print statements at end of MA48BD and
C          start of MA48CD.
C 07/09/96 Resetting of IPTRD and IRN after call to MA50BD moved to
C          immediately after the call so that they will be reset
C          correctly in the event of an error return from MA50BD.
C 12/02/01 Corrections made to the copying of the permutations found
C          by MA50A/AD in rectangular case.
C          Commas added after 1P edit descriptors.
C          Length of ICNTL5 corrected in MA48D/DD

      SUBROUTINE MA48AD(M,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,IW,INFO,
     +                  RINFO)
      INTEGER M,N,NE,JOB,LA
      DOUBLE PRECISION A(LA)
      INTEGER IRN(LA),JCN(LA),KEEP(*)
      DOUBLE PRECISION CNTL(5)
      INTEGER ICNTL(9),IW(6*M+3*N),INFO(12)
      DOUBLE PRECISION RINFO
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.D0)
      DOUBLE PRECISION CNTL5(4)
      INTEGER EYE,HEADC,I,IB,ICNTL5(7),IDUMMY,INFO5(7),IP,IPTRA,IPTRD,
     +        IPTRO,IPTRP,IQ,ISW,IW13,IW21,IW50,J,JAY,JB,JFIRST,J1,J2,
     +        J3,K,KB,KBLOCK,KD,KK,KL,KO,L,LASTR,LASTC,LBLOCK
      LOGICAL LDUP
      INTEGER LENC,LENP,LENR,LP,MBLOCK,MINBLK,MP,NB,NBLOCK,NC,NDIAG,
     +        NEXTC,NEXTR,NP,NR,NXTEYE,NZA,NZB,NZD,PTRD,PTRO
      DOUBLE PRECISION RINFO5,TOL
      EXTERNAL MC13DD,MC21AD,MA50AD
      INTRINSIC ABS,MAX
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (M.LE.0 .OR. N.LE.0) THEN
         INFO(1) = -1
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9100) M,N
         GO TO 530
      END IF
      IF (NE.LE.0) THEN
         INFO(1) = -2
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9110) NE
         GO TO 530
      END IF
      IF (LA.LT.2*NE) THEN
         INFO(1) = -3
         INFO(3) = 2*NE
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9120) LA,INFO(3)
         GO TO 530
      END IF
      IF (JOB.LT.1 .OR. JOB.GT.3) THEN
         INFO(1) = -6
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9130) JOB
         GO TO 530
      END IF
      IF (JOB.EQ.2) THEN
         DO 10 I = 1,MAX(M,N)
            IW(N+I) = 0
   10    CONTINUE
         DO 20 I = 1,M
            J = KEEP(I)
            IF (J.LT.1 .OR. J.GT.M) GO TO 40
            IF (IW(N+J).EQ.1) GO TO 40
            IW(N+J) = 1
   20    CONTINUE
         DO 30 I = 1,N
            J = KEEP(M+I)
            IF (J.LT.1 .OR. J.GT.N) GO TO 40
            IF (IW(N+J).EQ.2) GO TO 40
            IW(N+J) = 2
   30    CONTINUE
         GO TO 50
   40    INFO(1) = -5
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9140)
         GO TO 530
      END IF
   50 INFO(1) = 0
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,
     +     '(/A/A,I7,A,I6,A,I7,A,I2,A,I7/A,1P,4D12.4/A,4I8/A,3I8)')
     +     ' Entering MA48AD with',' M =',M,'     N =',N,'     NE =',NE,
     +     '     JOB =',JOB,'     LA =',LA,' CNTL (1:4) =',
     +     (CNTL(I),I=1,4),' ICNTL(1:4) = ', (ICNTL(I),I=1,4),
     +     ' ICNTL(6:8) = ', (ICNTL(I),I=6,8)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9000) (A(K),IRN(K),JCN(K),K=1,NE)
 9000       FORMAT (' Entries:'/3 (1PD12.4,2I6))
         ELSE
            WRITE (MP,9000) (A(K),IRN(K),JCN(K),K=1,MIN(9,NE))
         END IF
         IF (JOB.EQ.2) THEN
            WRITE (MP,'(A)') ' Permutations input (JOB=2)'
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9010) (KEEP(I),I=1,M)
 9010          FORMAT (' Positions of original rows in the permuted ma',
     +                'trix'/ (10I6))
               WRITE (MP,9020) (KEEP(M+I),I=1,N)
 9020          FORMAT (' Positions of columns of permuted matrix ','in',
     +                ' or','iginal matrix '/ (10I6))
            ELSE
               WRITE (MP,9010) (KEEP(I),I=1,MIN(10,M))
               WRITE (MP,9020) (KEEP(M+I),I=1,MIN(10,N))
            END IF
         END IF
         IF (ICNTL(8).NE.0) THEN
            WRITE (MP,'(A,I6)')
     +        ' Value of IW entries on call with ICNTL(8) =',ICNTL(8)
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9030) (IW(I),I=1,N)
 9030          FORMAT (10I6)
            ELSE
               WRITE (MP,9030) (IW(I),I=1,MIN(10,N))
            END IF
         END IF
      END IF
      INFO(2) = 0
      INFO(3) = NE*2
      INFO(4) = NE
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
      INFO(8) = 0
      INFO(9) = 0
      INFO(10) = MIN(M,N)
      INFO(11) = 0
      INFO(12) = 0
      RINFO = ZERO
      DO 60 I = 1,4
         CNTL5(I) = CNTL(I)
         ICNTL5(I) = ICNTL(I)
   60 CONTINUE
      ICNTL5(3) = 0
      ICNTL5(5) = ICNTL(5)
      ICNTL5(6) = 0
      ICNTL5(7) = 0
      IF (JOB.EQ.2) THEN
         ICNTL5(7) = 2
         ICNTL5(4) = 1
      END IF
      IF (JOB.EQ.3) ICNTL5(7) = 1
      TOL = MAX(ZERO,CNTL(4))
      MINBLK = MAX(1,ICNTL(6))
      LDUP = .FALSE.
      IPTRD = M + 3*N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1
      HEADC = N + 1
      LASTC = HEADC + N
      DO 70 J = 1,N
         IW(HEADC+J) = 0
   70 CONTINUE
      DO 80 K = 1,NE
         I = IRN(K)
         J = JCN(K)
         IF (I.LT.1 .OR. I.GT.M .OR. J.LT.1 .OR. J.GT.N) THEN
            INFO(12) = INFO(12) + 1
            IF (MP.GT.0 .AND. INFO(12).LE.10 .AND.
     +          ICNTL(3).GE.2) WRITE (MP,'(A,I7,A,2I6)')
     +          ' Message from MA48A/AD .. indices for entry ',K,' are',
     +          I,J
            JCN(K) = 0
         ELSE
            JCN(K) = IW(HEADC+J)
            IW(HEADC+J) = K
         END IF
   80 CONTINUE
      IF (MINBLK.GE.N .OR. M.NE.N .OR. JOB.GT.1) GO TO 190
      IW21 = 2*N
      IW13 = IW21
      IPTRA = LASTC + N
      LENC = IPTRA + N
      IB = LENC
      IP = LENC + N
      IPTRP = IP + N
      LENP = IPTRP + N
      DO 90 I = 1,N
         IW(LASTC+I) = 0
   90 CONTINUE
      LDUP = .TRUE.
      K = 1
      DO 120 J = 1,N
         EYE = IW(HEADC+J)
         IW(IPTRA+J) = K
         DO 100 IDUMMY = 1,NE
            IF (EYE.EQ.0) GO TO 110
            I = IRN(EYE)
            IF (IW(LASTC+I).NE.J) THEN
               IW(LASTC+I) = J
               IRN(NE+K) = I
               K = K + 1
            ELSE
               INFO(11) = INFO(11) + 1
               IF (MP.GT.0 .AND. INFO(11).LE.10 .AND.
     +             ICNTL(3).GE.2) WRITE (MP,'(A,I7,A,2I6)')
     +             ' Message from MA48A/AD .. duplicate in position ',K,
     +             ' with indices',I,J
            END IF
            EYE = JCN(EYE)
  100    CONTINUE
  110    IW(LENC+J) = K - IW(IPTRA+J)
  120 CONTINUE
      NZA = K - 1
      CALL MC21AD(N,IRN(NE+1),NZA,IW(IPTRA+1),IW(LENC+1),IW(IP+1),NDIAG,
     +           KEEP(IW21+1))
      INFO(10) = NDIAG
      IF (NDIAG.LT.N) THEN
         IF (ICNTL(7).NE.0) THEN
            INFO(1) = -4
            IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,
     +          '(A,A/A,I7,A,I7)')
     +        ' Error return from MA48A/AD because matrix structurally '
     +          ,' singular',' order is ',N,' and structural rank',NDIAG
            GO TO 530
         END IF
         GO TO 190
      END IF
CDIR$ IVDEP
      DO 130 J = 1,N
         JAY = IW(IP+J)
         IW(IPTRP+J) = IW(IPTRA+JAY)
         IW(LENP+J) = IW(LENC+JAY)
  130 CONTINUE
      CALL MC13DD(N,IRN(NE+1),NZA,IW(IPTRP+1),IW(LENP+1),KEEP(M+1),
     +           IW(IB+1),NB,KEEP(IW13+1))
      DO 140 JB = 2,NB
         IW(IB+JB-1) = IW(IB+JB) - IW(IB+JB-1)
  140 CONTINUE
      IW(IB+NB) = N + 1 - IW(IB+NB)
      IF (IW(IB+1).EQ.1) IW(IB+1) = -1
      KB = 1
      DO 150 JB = 2,NB
         L = IW(IB+JB)
         IF (L.EQ.1 .AND. IW(IB+KB).LE.0) THEN
            IW(IB+KB) = IW(IB+KB) - 1
         ELSE
            KB = KB + 1
            IF (L.EQ.1) THEN
               IW(IB+KB) = -1
            ELSE
               IW(IB+KB) = L
            END IF
         END IF
  150 CONTINUE
      NB = KB
      KB = 1
      DO 160 JB = 2,NB
         IF (ABS(IW(IB+KB)).LT.MINBLK) THEN
            IW(IB+KB) = ABS(IW(IB+KB)) + ABS(IW(IB+JB))
         ELSE
            KB = KB + 1
            IW(IB+KB) = IW(IB+JB)
         END IF
  160 CONTINUE
      NB = KB
      DO 170 JB = 1,NB
         KEEP(NBLOCK+3*JB) = ABS(IW(IB+JB))
         KEEP(MBLOCK+3*JB) = IW(IB+JB)
  170 CONTINUE
      DO 180 J = 1,N
         KEEP(KEEP(M+J)) = J
         KEEP(M+J) = IW(IP+KEEP(M+J))
  180 CONTINUE
      GO TO 220
  190 NB = 1
      IF (JOB.EQ.1 .OR. JOB.EQ.3) THEN
         DO 200 I = 1,M
            KEEP(I) = I
  200    CONTINUE
         DO 210 I = 1,N
            KEEP(M+I) = I
  210    CONTINUE
      END IF
      KEEP(NBLOCK+3) = N
      KEEP(MBLOCK+3) = 0
  220 IF (ICNTL(8).NE.0) THEN
         LBLOCK = KBLOCK + 3*NB
         IF (JOB.EQ.2) THEN
            DO 230 I = 1,N
               IF (IW(I).EQ.0) GO TO 240
  230       CONTINUE
  240       KEEP(LBLOCK+1) = N - I + 1
         ELSE
            J = 1
            DO 270 JB = 1,NB
               KEEP(LBLOCK+JB) = 0
               J2 = J + KEEP(NBLOCK+3*JB) - 1
               J1 = J2
               IF (KEEP(MBLOCK+3*JB).LT.0) GO TO 260
  250          IF (J.EQ.J2) GO TO 260
               IF (IW(KEEP(M+J)).EQ.0) THEN
                  KEEP(LBLOCK+JB) = KEEP(LBLOCK+JB) + 1
                  ISW = KEEP(M+J2)
                  KEEP(M+J2) = KEEP(M+J)
                  KEEP(M+J) = ISW
                  J2 = J2 - 1
               ELSE
                  J = J + 1
               END IF
               GO TO 250
  260          J = J1 + 1
  270       CONTINUE
         END IF
      END IF
      LASTR = LASTC + M
      DO 280 I = 1,M
         IW(LASTC+I) = 0
  280 CONTINUE
      KEEP(KBLOCK+3) = NB
      K = 1
      KK = NE
      J2 = 0
      DO 310 JB = 1,NB
         J1 = J2 + 1
         J2 = J1 + KEEP(NBLOCK+3*JB) - 1
         DO 300 JAY = J1,J2
            J = KEEP(M+JAY)
            EYE = IW(HEADC+J)
            KEEP(IPTRD+JAY) = K
            KEEP(IPTRO+JAY) = KK
            IF (KEEP(MBLOCK+3*JB).LT.0) THEN
               IW(LASTC+JAY) = JAY
               A(NE+K) = ZERO
               IRN(NE+K) = JAY
               IW(LASTR+JAY) = K
               K = K + 1
            END IF
            DO 290 IDUMMY = 1,NE
               IF (EYE.EQ.0) GO TO 300
               NXTEYE = JCN(EYE)
               I = KEEP(IRN(EYE))
               IF (IW(LASTC+I).NE.JAY) THEN
                  IW(LASTC+I) = JAY
                  IF ((I.GE.J1.AND.I.LE.J2) .OR. (M.NE.N)) THEN
                     A(NE+K) = A(EYE)
                     IRN(NE+K) = I
                     IW(LASTR+I) = K
                     JCN(EYE) = K
                     K = K + 1
                  ELSE
                     A(NE+KK) = A(EYE)
                     IRN(NE+KK) = I
                     IW(LASTR+I) = KK
                     JCN(EYE) = KK
                     KK = KK - 1
                  END IF
               ELSE
                  IF (.NOT.LDUP) THEN
                     INFO(11) = INFO(11) + 1
                     IF (MP.GT.0 .AND. INFO(11).LE.10 .AND.
     +                   ICNTL(3).GE.2) WRITE (MP,'(A,I7,A,2I6)')
     +                ' Message from MA48A/AD .. duplicate in position '
     +                   ,EYE,' with indices',IRN(EYE),J
                  END IF
                  KL = IW(LASTR+I)
                  JCN(EYE) = KL
                  A(NE+KL) = A(NE+KL) + A(EYE)
               END IF
               EYE = NXTEYE
  290       CONTINUE
  300    CONTINUE
  310 CONTINUE
      KEEP(IPTRD+N+1) = K
      KEEP(IPTRO+N+1) = KK
      NZD = K - 1
      DO 320 I = 1,K-1
         IRN(I) = IRN(NE+I)
  320 CONTINUE
      DO 325 I = KK+1,NE
         IRN(I) = IRN(NE+I)
  325 CONTINUE
      DO 326 I = K,KK
         IRN(I) = 0
  326 CONTINUE
      IP = 0
      IQ = M
      PTRD = M + N
      JFIRST = M + 2*N
      LENR = 2*M + 2*N
      LASTR = 3*M + 2*N
      NEXTR = 4*M + 2*N
      IW50 = 5*M + 2*N
      NEXTC = 6*M + 2*N
      PTRO = NEXTC
      LASTC = M + N
      LENC = LASTC + N
      J1 = N + 1
      KB = 0
      DO 390 JB = NB,1,-1
         NC = KEEP(NBLOCK+3*JB)
         J2 = J1 - 1
         J1 = J2 + 1 - NC
         IF (KEEP(MBLOCK+3*JB).LT.0) THEN
            DO 330 J = J1,J2
               IF (ABS(A(NE+KEEP(IPTRD+J))).GT.TOL)
     +             INFO(5) = INFO(5) + 1
               IW(IP+J) = J
               IW(IQ+J) = J
  330       CONTINUE
         ELSE
            NZB = KEEP(IPTRD+J2+1) - KEEP(IPTRD+J1)
            DO 340 K = KEEP(IPTRD+J1),KEEP(IPTRD+J2+1) - 1
               IRN(NE+K) = IRN(K) - J1 + 1
  340       CONTINUE
            K = KEEP(IPTRD+J1) - 1
            DO 350 J = J1,J2
               IW(IQ+J) = KEEP(IPTRD+J) - K
  350       CONTINUE
            NR = NC
            IF (NB.EQ.1) NR = M
            IF (JOB.EQ.2) THEN
               DO 360 J = J1,J1 + NR - 1
                  IW(IP+J) = J - J1 + 1
  360          CONTINUE
               DO 370 J = J1,J2
                  IW(PTRD+J) = J - J1 + 1
  370          CONTINUE
            END IF
            INFO(7) = MAX(INFO(7),NR)
            INFO(8) = INFO(8) + NC
            INFO(9) = INFO(9) + NZB
            IF (ICNTL(8).NE.0) ICNTL5(6) = KEEP(LBLOCK+JB)
            CALL MA50AD(NR,NC,NZB,LA-NE-K,A(NE+K+1),IRN(NE+K+1),
     +                  JCN(NE+1),IW(IQ+J1),CNTL5,ICNTL5,IW(IP+J1),
     +                  NP,IW(JFIRST+1),IW(LENR+1),
     +                  IW(LASTR+1),IW(NEXTR+1),IW(IW50+1),IW(PTRD+J1),
     +                  KEEP(LENC+KB+1),KEEP(LASTC+KB+1),IW(NEXTC+1),
     +                  INFO5,RINFO5)
            KEEP(MBLOCK+3*JB) = NP
            DO 380 J = J1,J1+NR-1
               IW(IP+J) = IW(IP+J) + J1 - 1
  380       CONTINUE
            DO 385 J = J1,J2
               IW(IQ+J) = IW(IQ+J) + J1 - 1
  385       CONTINUE
            IF (INFO5(1).EQ.1) THEN
               IF (INFO(1).EQ.0 .OR. INFO(1).EQ.4) INFO(1) = INFO(1) + 2
            END IF
            IF (INFO5(1).EQ.2) THEN
               IF (INFO(1).EQ.0 .OR. INFO(1).EQ.2) INFO(1) = INFO(1) + 4
            END IF
            IF (INFO5(1).EQ.3 .AND. INFO(1).GE.0) INFO(1) = 6
            IF (INFO5(1).EQ.-3) INFO(1) = -3
            INFO(2) = INFO(2) + INFO5(2)
            INFO(3) = MAX(INFO(3),NE+K+INFO5(3))
            INFO(5) = INFO(5) + INFO5(5)
            INFO(6) = INFO(6) + INFO5(6)
            RINFO = RINFO + RINFO5
            KB = KB + 1
            KEEP(LENC+KB) = INFO5(4) - 2*NR
            KEEP(LASTC+KB) = INFO5(4)
         END IF
  390 CONTINUE
      INFO(4) = NE
      K = NE
      DO 400 JB = KB,1,-1
         INFO(4) = MAX(INFO(4),K+KEEP(LASTC+JB))
         K = K + KEEP(LENC+JB)
  400 CONTINUE
      IF (INFO(1).EQ.-3) THEN
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9120) LA,INFO(3)
         GO TO 530
      END IF
      DO 410 K = 1,NE
         IRN(NE+K) = IRN(K)
  410 CONTINUE
      DO 420 J = 1,N
         IW(PTRD+J) = KEEP(IPTRD+J)
         IW(PTRO+J) = KEEP(IPTRO+J+1) + 1
  420 CONTINUE
      KD = 1
      KO = NZD + 1
      DO 450 J = 1,N
         KEEP(IPTRD+J) = KD
         JAY = IW(IQ+J)
         KL = NZD
         IF (JAY.NE.N) KL = IW(PTRD+JAY+1) - 1
         DO 430 KK = IW(PTRD+JAY),KL
            IRN(KD) = IW(IP+IRN(NE+KK))
            IRN(NE+KK) = KD
            KD = KD + 1
  430    CONTINUE
         KEEP(IPTRO+J) = KO
         KL = NE
         IF (JAY.NE.1) KL = IW(PTRO+JAY-1) - 1
         DO 440 KK = IW(PTRO+JAY),KL
            IRN(KO) = IW(IP+IRN(NE+KK))
            IRN(NE+KK) = KO
            KO = KO + 1
  440    CONTINUE
  450 CONTINUE
      KEEP(IPTRO+N+1) = KO
      DO 460 I = 1,M
         KEEP(I) = IW(IP+KEEP(I))
  460 CONTINUE
      DO 465 I = 1,N
         IW(IQ+I) = KEEP(M+IW(IQ+I))
  465 CONTINUE
      DO 470 I = 1,N
         KEEP(M+I) = IW(IQ+I)
  470 CONTINUE
      IW(1) = IRN(NE)
      IRN(NE) = NE
      DO 480 K = 1,NE
         JCN(K) = IRN(NE+JCN(K))
  480 CONTINUE
      IRN(NE) = IW(1)
      IF (INFO(11).GT.0 .OR. INFO(12).GT.0) THEN
         INFO(1) = INFO(1) + 1
         IF (MP.GT.0 .AND. ICNTL(3).GE.2) THEN
            IF (INFO(11).GT.0) WRITE (MP,9150) INFO(11)
            IF (INFO(12).GT.0) WRITE (MP,9160) INFO(12)
         END IF
         IF (INFO(11).GT.0) JCN(1) = -JCN(1)
      END IF
      IF (INFO(10).LT.INFO(5)) THEN
         INFO(5) = INFO(10)
         IF (INFO(1).NE.2 .AND. INFO(1).NE.3 .AND.
     +       INFO(1).LT.6) INFO(1) = INFO(1) + 2
      END IF
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/A/A,I5,11I6/A,F12.1)') ' Leaving MA48AD with',
     +     ' INFO  =',INFO,' RINFO =',RINFO
         WRITE (MP,'(A)') ' Permuted matrix by blocks (IRN)'
         KB = NB
         IF (ICNTL(3).EQ.3) THEN
            WRITE (MP,'(A)')
     +        ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         END IF
         WRITE (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         DO 500 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (J1.LE.J2) WRITE (MP,'(A,I6)') ' Block',JB
            J3 = J2
            IF (ICNTL(3).EQ.3) J3 = J1
            DO 490 J = J1,J3
               WRITE (MP,'(A,I5,(T13,10I6))') ' Column',J,
     +           (IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
  490       CONTINUE
            J1 = J2 + 1
  500    CONTINUE
         IF (KEEP(IPTRO+N+1).GT.KEEP(IPTRD+N+1)) THEN
            WRITE (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            DO 520 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               IF (ICNTL(3).EQ.3) J3 = J1
               DO 510 J = J1,J3
                  IF (KEEP(IPTRO+J+1).GT.KEEP(IPTRO+J)) WRITE (MP,
     +                '(A,I5,(T13,10I6))') ' Column',J,
     +                (IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
  510          CONTINUE
               J1 = J2 + 1
  520       CONTINUE
         END IF
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9040) (JCN(K),K=1,NE)
 9040       FORMAT (' JCN (MAP) ='/ (6X,10I6))
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9010) (KEEP(I),I=1,M)
            WRITE (MP,9020) (KEEP(M+I),I=1,N)
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9050) (KEEP(IPTRD+J),J=1,N+1)
 9050       FORMAT (' IPTRD ='/ (8X,10I6))
            WRITE (MP,9060) (KEEP(IPTRO+J),J=1,N+1)
 9060       FORMAT (' IPTRO ='/ (8X,10I6))
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9070) (KEEP(NBLOCK+3*JB),JB=1,NB)
            WRITE (MP,9080) (KEEP(MBLOCK+3*JB),JB=1,NB)
 9070       FORMAT (' NBLOCK (order blocks) ='/ (8X,10I6))
 9080       FORMAT (' MBLOCK (triangular flag and number packed rows) ='
     +             / (8X,10I6))
 9090       FORMAT (' LBLOCK (number of changed columns) ='/ (8X,10I6))
            IF (ICNTL(8).NE.0) WRITE (MP,9090) (KEEP(LBLOCK+JB),JB=1,NB)
         ELSE
            WRITE (MP,9040) (JCN(K),K=1,MIN(10,NE))
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9010) (KEEP(I),I=1,MIN(10,M))
            WRITE (MP,9020) (KEEP(M+I),I=1,MIN(10,N))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9050) (KEEP(IPTRD+J),J=1,MIN(10,N+1))
            WRITE (MP,9060) (KEEP(IPTRO+J),J=1,MIN(10,N+1))
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9070) (KEEP(NBLOCK+3*JB),JB=1,MIN(10,NB))
            WRITE (MP,9080) (KEEP(MBLOCK+3*JB),JB=1,MIN(10,NB))
            IF (ICNTL(8).NE.0) WRITE (MP,9090) (KEEP(LBLOCK+JB),JB=1,
     +          MIN(10,NB))
         END IF
      END IF
  530 RETURN
 9100 FORMAT (' Error return from MA48A/AD because M =',I10,' and N =',
     +       I10)
 9110 FORMAT (' Error return from MA48A/AD because NE =',I10)
 9120 FORMAT (' Error return from MA48A/AD because LA is',I10/' and ',
     +       'must be at least',I10)
 9130 FORMAT (' Error return from MA48A/AD because ','JOB = ',I10)
 9140 FORMAT (' Error return from MA48A/AD because ','faulty permutati',
     +       'ons input when JOB = 2')
 9150 FORMAT (' Message from MA48A/AD ..',I8,' duplicates found')
 9160 FORMAT (' Message from MA48A/AD ..',I8,' out-of-range indices fo',
     +       'und')
      END
      SUBROUTINE MA48BD(M,N,NE,JOB,LA,A,IRN,JCN,KEEP,CNTL,ICNTL,W,IW,
     +                  INFO,RINFO)
      INTEGER M,N,NE,JOB,LA
      DOUBLE PRECISION A(LA)
      INTEGER IRN(LA),JCN(NE),KEEP(*)
      DOUBLE PRECISION CNTL(5)
      INTEGER ICNTL(9)
      DOUBLE PRECISION W(M)
      INTEGER IW(2*M+2*N),INFO(12)
      DOUBLE PRECISION RINFO
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.D0)
      DOUBLE PRECISION CNTL5(4)
      INTEGER I,ICNTL5(7),INFO5(7),IPTRD,IPTRL,IPTRO,IPTRU,IQB(1),J,JB,
     +        J1,J2,J3,K,KB,KBLOCK,KK,LBLOCK,LP,MBLOCK,MP,NB,NBLOCK,
     +        NEWNE,NC,NP,NR,NRF,NZB
      DOUBLE PRECISION RINFO5,TOL
      LOGICAL TRISNG
      EXTERNAL MA50BD
      INTRINSIC MAX
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (M.LE.0 .OR. N.LE.0) THEN
         INFO(1) = -1
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9160) M,N
         GO TO 240
      END IF
      IF (NE.LE.0) THEN
         INFO(1) = -2
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9170) NE
         GO TO 240
      END IF
      IF (LA.LT.2*NE) THEN
         INFO(1) = -3
         INFO(4) = 2*NE
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9180) LA,INFO(4)
         GO TO 240
      END IF
      IF (JOB.LT.1 .OR. JOB.GT.3) THEN
         INFO(1) = -6
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9190) JOB
         GO TO 240
      END IF
      INFO(1) = 0
      DO 10 I = 1,4
         CNTL5(I) = CNTL(I)
         ICNTL5(I) = ICNTL(I)
   10 CONTINUE
      ICNTL5(3) = 0
      ICNTL5(5) = ICNTL(5)
      ICNTL5(6) = 0
      ICNTL5(7) = 0
      IPTRL = M + N
      IPTRU = IPTRL + N
      IPTRD = IPTRU + N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1
      NB = KEEP(KBLOCK+3)
      LBLOCK = KBLOCK + 3*NB
      NEWNE = KEEP(IPTRO+N+1) - 1
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/A/3(A,I8),A,I2,A,I8/A,1P,3D12.4/A,3I8/A,I8/A,I8)')
     +     ' Entering MA48BD with',' M =',M,'     N =',N,'     NE =',NE,
     +     '     JOB =',JOB,'     LA =',LA,' CNTL (2:4) =',
     +     (CNTL(I),I=2,4),' ICNTL(1:3) =', (ICNTL(I),I=1,3),
     +     ' ICNTL(5)   =',ICNTL(5),' ICNTL(8)   =',ICNTL(8)
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9000) (A(K),K=1,NE)
         ELSE
            WRITE (MP,9000) (A(K),K=1,MIN(10,NE))
         END IF
 9000    FORMAT (' A ='/ (4X,1P,5D12.4))
         WRITE (MP,'(A)') ' Indices for permuted matrix by blocks'
         KB = NB
         IF (ICNTL(3).EQ.3) THEN
            WRITE (MP,'(A)')
     +        ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         END IF
         WRITE (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         DO 30 JB = 1,KB
            WRITE (MP,'(A,I6)') ' Block',JB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            J3 = J2
            IF(ICNTL(3).EQ.3) J3 = J1
            DO 20 J = J1,J3
               WRITE (MP,'(A,I5,(T13,10I6))') ' Column',J,
     +           (IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
   20       CONTINUE
            J1 = J2 + 1
   30    CONTINUE
         IF (KEEP(IPTRO+N+1).GT.KEEP(IPTRD+N+1)) THEN
            WRITE (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            DO 50 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               IF(ICNTL(3).EQ.3) J3 = J1
               DO 40 J = J1,J3
                  IF (KEEP(IPTRO+J+1).GT.KEEP(IPTRO+J)) WRITE (MP,
     +                '(A,I5,(T13,10I6))') ' Column',J,
     +                (IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
   40          CONTINUE
               J1 = J2 + 1
   50       CONTINUE
         END IF
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9010) (JCN(K),K=1,NE)
 9010       FORMAT (' JCN (MAP) ='/ (6X,10I6))
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9020) (KEEP(I),I=1,M)
            WRITE (MP,9030) (KEEP(M+I),I=1,N)
 9020       FORMAT (' Positions of original rows in the permuted matrix'
     +             / (10I6))
 9030       FORMAT (' Positions of columns of permuted matrix ','in ',
     +             'original matrix '/ (10I6))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9040) (KEEP(IPTRD+J),J=1,N+1)
 9040       FORMAT (' IPTRD ='/ (8X,10I6))
            WRITE (MP,9050) (KEEP(IPTRO+J),J=1,N+1)
 9050       FORMAT (' IPTRO ='/ (8X,10I6))
            IF (JOB.GT.1) THEN
               WRITE (MP,9060) (KEEP(IPTRL+J),J=1,N)
 9060          FORMAT (' IPTRL ='/ (8X,10I6))
               WRITE (MP,9070) (KEEP(IPTRU+J),J=1,N)
 9070          FORMAT (' IPTRU ='/ (8X,10I6))
            END IF
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,NB)
            WRITE (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,NB)
 9080       FORMAT (' NBLOCK (order blocks) ='/ (8X,10I6))
 9090       FORMAT (' MBLOCK (triangular flag and number packed rows) ='
     +             / (8X,10I6))
 9100       FORMAT (' KBLOCK (position of beginning of block) ='/
     +             (8X,10I6))
 9110       FORMAT (' LBLOCK (number of changed columns) ='/ (8X,10I6))
            IF (JOB.GT.1) THEN
               WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,NB)
               IF (ICNTL(8).NE.0) WRITE (MP,9110) (KEEP(LBLOCK+JB),JB=1,
     +             NB)
            END IF
         ELSE
            WRITE (MP,9010) (JCN(K),K=1,MIN(10,NE))
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9020) (KEEP(I),I=1,MIN(10,M))
            WRITE (MP,9030) (KEEP(M+I),I=1,MIN(10,N))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9040) (KEEP(IPTRD+J),J=1,MIN(10,N+1))
            WRITE (MP,9050) (KEEP(IPTRO+J),J=1,MIN(10,N+1))
            IF (JOB.GT.1) THEN
               WRITE (MP,9060) (KEEP(IPTRL+J),J=1,MIN(10,N))
               WRITE (MP,9070) (KEEP(IPTRU+J),J=1,MIN(10,N))
            END IF
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,MIN(10,NB))
            WRITE (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,MIN(10,NB))
            IF (JOB.GT.1) THEN
               WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,MIN(10,NB))
               IF (ICNTL(8).NE.0) WRITE (MP,9110) (KEEP(LBLOCK+JB),JB=1,
     +             MIN(10,NB))
            END IF
         END IF
      END IF
      INFO(4) = NE
      INFO(5) = 0
      INFO(6) = 0
      RINFO = ZERO
      TOL = MAX(ZERO,CNTL(4))
      IF (JCN(1).GT.0) THEN
         DO 60 K = 1,NE
            A(NE+K) = A(K)
   60    CONTINUE
CDIR$ IVDEP
         DO 70 K = 1,NE
            A(JCN(K)) = A(NE+K)
   70    CONTINUE
      ELSE
         DO 80 K = 1,NE
            A(NE+K) = A(K)
            A(K) = ZERO
   80    CONTINUE
         A(-JCN(1)) = A(NE+1)
         DO 90 K = 2,NE
            KK = JCN(K)
            A(KK) = A(KK) + A(NE+K)
   90    CONTINUE
      END IF
      IQB(1) = 0
      KK = 0
      J2 = 0
      DO 150 JB = 1,NB
         NC = KEEP(NBLOCK+3*JB)
         J1 = J2 + 1
         J2 = J1 + NC - 1
         KEEP(KBLOCK+3*JB) = 0
         IF (KEEP(MBLOCK+3*JB).LT.0) THEN
            TRISNG = .FALSE.
            DO 100 J = J1,J2
               IF (ABS(A(KEEP(IPTRD+J))).LE.TOL) TRISNG = .TRUE.
               KEEP(IPTRL+J) = 0
               KEEP(IPTRU+J) = 0
  100       CONTINUE
            IF (.NOT.TRISNG)  THEN
               INFO(5) = INFO(5) + NC
               GO TO 150
            ENDIF
            IF (JOB.EQ.2) THEN
               INFO(1) = -7
               IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,'(A)')
     +              ' Error return from MA48B/BD with JOB=2 because ',
     +              ' the matrix is incompatible with expectations'
               GO TO 240
            ELSE
               KEEP(MBLOCK+3*JB) = NC
            ENDIF
         END IF
         NR = NC
         IF (NB.EQ.1) NR = M
         NZB = KEEP(IPTRD+J2+1) - KEEP(IPTRD+J1)
         IF (ICNTL(8).NE.0) ICNTL5(6) = KEEP(LBLOCK+JB)
         DO 110 K = KEEP(IPTRD+J1),KEEP(IPTRD+J2+1) - 1
            IRN(K) = IRN(K) - J1 + 1
  110    CONTINUE
         K = KEEP(IPTRD+J1) - 1
         DO 120 J = J1,J1+NR-1
            KEEP(IPTRD+J) = KEEP(IPTRD+J) - K
            IW(J) = J - J1 + 1
  120    CONTINUE
         NP = KEEP(MBLOCK+3*JB)
         CALL MA50BD(NR,NC,NZB,JOB,A(K+1),IRN(K+1),KEEP(IPTRD+J1),CNTL5,
     +               ICNTL5,IW(J1),IQB,NP,LA-NEWNE-KK,
     +               A(NEWNE+KK+1),IRN(NEWNE+KK+1),KEEP(IPTRL+J1),
     +               KEEP(IPTRU+J1),W,IW(M+1),INFO5,RINFO5)
         DO 130 J = J1,J2
            KEEP(IPTRD+J) = KEEP(IPTRD+J) + K
  130    CONTINUE
         DO 140 K = KEEP(IPTRD+J1),KEEP(IPTRD+J2+1) - 1
            IRN(K) = IRN(K) + J1 - 1
  140    CONTINUE
         IF (INFO5(1).EQ.-6) THEN
            INFO(1) = -6
            IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,'(A)')
     +          ' Error return from MA48B/BD with JOB greater than 1',
     +          ' and entries dropped during previous factorization'
            GO TO 240
         END IF
         IF (INFO5(1).LT.-7) THEN
            INFO(1) = -7
            IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,'(A)')
     +          ' Error return from MA48B/BD with JOB=2 because ',
     +          ' the matrix is incompatible with expectations'
            GO TO 240
         END IF
         IF (INFO5(1).EQ.-3) INFO(1) = -3
         IF (INFO(1).EQ.-3) THEN
            INFO(4) = INFO(4) + INFO5(4)
            KK = 0
         ELSE
            INFO(4) = MAX(INFO(4),KK+NEWNE+INFO5(4))
            NRF = IRN(NEWNE+KK+2)
            KEEP(KBLOCK+3*JB) = KK + 1
            KK = KK + KEEP(IPTRL+J2) + MAX((NC-KEEP(MBLOCK+3*JB))*
     +           (NRF), (NC-KEEP(MBLOCK+3*JB))+(NRF))
         END IF
         IF (INFO5(1).EQ.1) THEN
            IF (INFO(1).NE.-3) INFO(1) = 2
         END IF
         RINFO = RINFO + RINFO5
         INFO(5) = INFO(5) + INFO5(5)
         INFO(6) = INFO(6) + INFO5(6)
  150 CONTINUE
      INFO(4) = MAX(NE*2,INFO(4))
      KEEP(KBLOCK+3) = NB
      IF (INFO(1).EQ.-3) THEN
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9180) LA,INFO(4)
         GO TO 240
      END IF
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/A/A,I6/A,3I6/A,F12.1)') ' Leaving MA48BD with',
     +     ' INFO(1)   = ',INFO(1),' INFO(4:6) = ', (INFO(I),I=4,6),
     +     ' RINFO     =',RINFO
         WRITE (MP,'(A)') ' Permuted matrix by blocks'
         KB = NB
         IF (ICNTL(3).EQ.3) THEN
            WRITE (MP,'(A)')
     +        ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         END IF
         WRITE (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         DO 170 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (J1.LE.J2) WRITE (MP,'(A,I6)') ' Block',JB
            J3 = J2
            IF (ICNTL(3).EQ.3) J3 = MIN(J1,J2)
            DO 160 J = J1,J3
               WRITE (MP,'(A,I5,(T13,3(1PD12.4,I5)))') ' Column',J,
     +           (A(I),IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
  160       CONTINUE
            J1 = J2 + 1
  170    CONTINUE
         IF (KEEP(IPTRO+N+1).GT.KEEP(IPTRD+N+1)) THEN
            WRITE (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            DO 190 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               IF (ICNTL(3).EQ.3) J3 = MIN(J1,J2)
               DO 180 J = J1,J3
                  IF (KEEP(IPTRO+J+1).GT.KEEP(IPTRO+J)) WRITE (MP,
     +                '(A,I5,(T13,3(1P,D12.4,I5)))') ' Column',J,
     +                (A(I),IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
  180          CONTINUE
               J1 = J2 + 1
  190       CONTINUE
         END IF
         WRITE (MP,'(A)') ' Factorized matrix by blocks'
         J1 = 1
         DO 230 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (KEEP(MBLOCK+3*JB).LT.0) GO TO 220
            NC = J2 - J1 + 1
            NR = NC
            IF (KB.EQ.1) NR = M
            WRITE (MP,'(A,I6)') ' Block',JB
            K = NEWNE
            IF (JB.GT.1) K = KEEP(KBLOCK+3*JB) + NEWNE - 1
            IF (KEEP(IPTRL+J1).GT.KEEP(IPTRU+J1)) WRITE (MP,
     +          '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J1,' of L',
     +          (A(K+I),IRN(K+I),I=KEEP(IPTRU+J1)+1,KEEP(IPTRL+J1))
            IF (ICNTL(3).EQ.3) GO TO 210
            DO 200 J = J1 + 1,J2
               IF (KEEP(IPTRU+J).GT.KEEP(IPTRL+J-1)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J,' of U',
     +              (A(K+I),IRN(K+I),I=KEEP(IPTRL+J-1)+1,KEEP(IPTRU+J))
               IF (KEEP(IPTRU+J).LT.KEEP(IPTRL+J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J,' of L',
     +              (A(K+I),IRN(K+I),I=KEEP(IPTRU+J)+1,KEEP(IPTRL+J))
  200       CONTINUE
  210       WRITE (MP,'(A)') ' Full block'
            WRITE (MP,'(A)') ' Row indices'
            NRF = IRN(K+2)
            K = K + KEEP(IPTRL+J2)
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9120) (IRN(K+I),I=1,NRF)
               WRITE (MP,'(A)') ' Column pivoting information'
               WRITE (MP,9120) (IRN(K+I),I=NRF+1,
     +           NRF+NC-KEEP(MBLOCK+3*JB))
               WRITE (MP,'(A)') ' Reals by columns'
               WRITE (MP,9130) (A(K+I),I=1,
     +           (NRF)* (NC-KEEP(MBLOCK+3*JB)))
 9120          FORMAT (10I6)
 9130          FORMAT (1P,5D12.4)
            ELSE
               WRITE (MP,9120) (IRN(K+I),I=1,MIN(10,NRF))
               WRITE (MP,'(A)') ' Column pivoting information'
               WRITE (MP,9120) (IRN(K+I),I=NRF+1,
     +           NRF+MIN(10,NC-KEEP(MBLOCK+3*JB)))
               WRITE (MP,'(A)') ' Reals by columns'
               WRITE (MP,9130) (A(K+I),I=1,
     +           MIN(10, (NRF)* (NC-KEEP(MBLOCK+3*JB))))
            END IF
  220       J1 = J2 + 1
  230    CONTINUE
         IF (JOB.EQ.1 .OR. JOB.EQ.3) THEN
            WRITE (MP,'(A)') ' Contents of KEEP array'
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9020) (KEEP(I),I=1,M)
               WRITE (MP,9030) (KEEP(M+I),I=1,N)
               WRITE (MP,'(A)') ' Pointer information from KEEP'
               WRITE (MP,9140) (KEEP(IPTRL+J),J=1,N)
 9140          FORMAT (' IPTRL ='/ (8X,10I6))
               WRITE (MP,9150) (KEEP(IPTRU+J),J=1,N)
 9150          FORMAT (' IPTRU ='/ (8X,10I6))
            ELSE
               WRITE (MP,9020) (KEEP(I),I=1,MIN(10,M))
               WRITE (MP,9030) (KEEP(M+I),I=1,MIN(10,N))
               WRITE (MP,'(A)') ' Pointer information from KEEP'
               WRITE (MP,9140) (KEEP(IPTRL+J),J=1,MIN(10,N))
               WRITE (MP,9150) (KEEP(IPTRU+J),J=1,MIN(10,N))
            END IF
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,KB)
         END IF
      END IF
  240 RETURN
 9160 FORMAT (' Error return from MA48B/BD because M =',I10,' and N =',
     +       I10)
 9170 FORMAT (' Error return from MA48B/BD because NE =',I10)
 9180 FORMAT (' Error return from MA48B/BD because LA is',I10/' and mu',
     +       'st be at least',I10)
 9190 FORMAT (' Error return from MA48B/BD because ','JOB = ',I10)
      END
      SUBROUTINE MA48CD(M,N,TRANS,JOB,LA,A,IRN,KEEP,CNTL,ICNTL,RHS,X,
     +                  ERROR,W,IW,INFO)
      INTEGER M,N
      LOGICAL TRANS
      INTEGER JOB,LA
      DOUBLE PRECISION A(LA)
      INTEGER IRN(LA),KEEP(*)
      DOUBLE PRECISION CNTL(5)
      INTEGER ICNTL(9)
      DOUBLE PRECISION RHS(*),X(*),ERROR(3),W(*)
      INTEGER IW(*),INFO(12)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.D0,ONE=1.0D0)
      DOUBLE PRECISION COND(2),CTAU,DXMAX
      INTEGER I,ICNTL5(7),IPTRD,IPTRL,IPTRO,IPTRU,IQB(1),J,JB,JJ,J1,J2,
     +        J3,K,KASE,KB,KBLOCK,KK
      LOGICAL LCOND(2)
      INTEGER LP,MBLOCK,MP,NB,NBLOCK,NC,NE,NEQ,NR,NRF,NVAR
      DOUBLE PRECISION OLDOMG(2),OMEGA(2),OM1,OM2,TAU
      EXTERNAL FD05AD,MA48DD,MA50CD,MC41AD
      DOUBLE PRECISION FD05AD
      INTRINSIC ABS,MAX
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (N.LE.0 .OR. M.LE.0) THEN
         INFO(1) = -1
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9140) M,N
         GO TO 380
      END IF
      IF (JOB.GT.4 .OR. JOB.LT.1) THEN
         INFO(1) = -6
         IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9150) JOB
         GO TO 380
      END IF
      INFO(1) = 0
      IPTRL = M + N
      IPTRU = IPTRL + N
      IPTRD = IPTRU + N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1
      NB = KEEP(KBLOCK+3)
      NE = KEEP(IPTRO+N+1) - 1
      OMEGA(1) = ZERO
      OMEGA(2) = ZERO
      ERROR(3) = ZERO
      CTAU = 1000.*FD05AD(1)
      IQB(1) = 0
      DO 10 I = 1,7
         ICNTL5(I) = 0
   10 CONTINUE
      ICNTL5(5) = ICNTL(5)
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,
     +     '(/A/3(A,I8),A,I2/A,L2,A,I7/A,1P,D12.4/A,3I6,2(/A,I6))')
     +     ' Entering MA48CD with',' M =',M,'     N =',N,'     LA =',LA,
     +     '      JOB =',JOB,'   TRANS =',TRANS,'      No. of blocks =',
     +     NB,'   CNTL(5)    = ',CNTL(5),'   ICNTL(1:3) = ',
     +     (ICNTL(I),I=1,3),'   ICNTL(5)   = ',ICNTL(5),
     +     '   ICNTL(9)   = ',ICNTL(9)
         WRITE (MP,'(A)') ' Permuted matrix by blocks'
         KB = NB
         IF (ICNTL(3).EQ.3) THEN
            WRITE (MP,'(A)')
     +        ' Only first column of up to 10 blocks printed'
            KB = MIN(10,NB)
         END IF
         WRITE (MP,'(A)') ' Diagonal blocks'
         J1 = 1
         DO 30 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (J1.LE.J2) WRITE (MP,'(A,I6)') ' Block',JB
            J3 = J2
            IF (ICNTL(3).EQ.3) J3 = MIN(J1,J2)
            DO 20 J = J1,J3
               WRITE (MP,'(A,I5,(T13,3(1P,D12.4,I5)))') ' Column',J,
     +           (A(I),IRN(I),I=KEEP(IPTRD+J),KEEP(IPTRD+J+1)-1)
   20       CONTINUE
            J1 = J2 + 1
   30    CONTINUE
         IF (KEEP(IPTRO+N+1).GT.KEEP(IPTRD+N+1)) THEN
            WRITE (MP,'(A)') ' Off-diagonal entries'
            J1 = 1
            DO 50 JB = 1,KB
               J2 = J1 + KEEP(NBLOCK+3*JB) - 1
               J3 = J2
               IF (ICNTL(3).EQ.3) J3 = MIN(J1,J2)
               DO 40 J = J1,J3
                  IF (KEEP(IPTRO+J+1).GT.KEEP(IPTRO+J)) WRITE (MP,
     +                '(A,I5,(T13,3(1P,D12.4,I5)))') ' Column',J,
     +                (A(I),IRN(I),I=KEEP(IPTRO+J),KEEP(IPTRO+J+1)-1)
   40          CONTINUE
               J1 = J2 + 1
   50       CONTINUE
         END IF
         WRITE (MP,'(A)') ' Factorized matrix by blocks'
         J1 = 1
         DO 90 JB = 1,KB
            J2 = J1 + KEEP(NBLOCK+3*JB) - 1
            IF (KEEP(MBLOCK+3*JB).LT.0) GO TO 80
            NC = J2 - J1 + 1
            NR = NC
            IF (KB.EQ.1) NR = M
            WRITE (MP,'(A,I6)') ' Block',JB
            K = NE
            IF (JB.GT.1) K = KEEP(KBLOCK+3*JB) + NE - 1
            IF (KEEP(IPTRL+J1).GT.KEEP(IPTRU+J1)) WRITE (MP,
     +          '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J1,' of L',
     +          (A(K+I),IRN(K+I),I=KEEP(IPTRU+J1)+1,KEEP(IPTRL+J1))
            IF (ICNTL(3).EQ.3) GO TO 70
            DO 60 J = J1 + 1,J2
               IF (KEEP(IPTRU+J).GT.KEEP(IPTRL+J-1)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J,' of U',
     +              (A(K+I),IRN(K+I),I=KEEP(IPTRL+J-1)+1,KEEP(IPTRU+J))
               IF (KEEP(IPTRU+J).LT.KEEP(IPTRL+J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,D12.4,I5)))') ' Column',J,' of L',
     +              (A(K+I),IRN(K+I),I=KEEP(IPTRU+J)+1,KEEP(IPTRL+J))
   60       CONTINUE
   70       WRITE (MP,'(A)') ' Full block'
            WRITE (MP,'(A)') ' Row indices'
            NRF = IRN(K+2)
            K = K + KEEP(IPTRL+J2)
            IF (ICNTL(3).GT.3) THEN
               WRITE (MP,9000) (IRN(K+I),I=1,NRF)
               WRITE (MP,'(A)') ' Column pivoting information'
               WRITE (MP,9000) (IRN(K+I),I=NRF+1,
     +           NRF+NC-KEEP(MBLOCK+3*JB))
               WRITE (MP,'(A)') ' Reals by columns'
               WRITE (MP,9010) (A(K+I),I=1,
     +           (NRF)* (NC-KEEP(MBLOCK+3*JB)))
 9000          FORMAT (10I6)
 9010          FORMAT (1P,5D12.4)
            ELSE
               WRITE (MP,9000) (IRN(K+I),I=1,MIN(10,NRF))
               WRITE (MP,'(A)') ' Column pivoting information'
               WRITE (MP,9000) (IRN(K+I),I=NRF+1,
     +           NRF+MIN(10,NC-KEEP(MBLOCK+3*JB)))
               WRITE (MP,'(A)') ' Reals by columns'
               WRITE (MP,9010) (A(K+I),I=1,
     +           MIN(10, (NRF)* (NC-KEEP(MBLOCK+3*JB))))
            END IF
   80       J1 = J2 + 1
   90    CONTINUE
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9020) (KEEP(I),I=1,M)
            WRITE (MP,9030) (KEEP(M+I),I=1,N)
 9020       FORMAT (' Positions of original rows in the permuted matrix'
     +             / (10I6))
 9030       FORMAT (' Positions of columns of permuted matrix ','in or',
     +             'ig','inal matrix '/ (10I6))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9040) (KEEP(IPTRD+J),J=1,N+1)
 9040       FORMAT (' IPTRD ='/ (8X,10I6))
            WRITE (MP,9050) (KEEP(IPTRO+J),J=1,N+1)
 9050       FORMAT (' IPTRO ='/ (8X,10I6))
            WRITE (MP,9060) (KEEP(IPTRL+J),J=1,N)
 9060       FORMAT (' IPTRL ='/ (8X,10I6))
            WRITE (MP,9070) (KEEP(IPTRU+J),J=1,N)
 9070       FORMAT (' IPTRU ='/ (8X,10I6))
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,NB)
            WRITE (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,NB)
            WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,NB)
 9080       FORMAT (' NBLOCK (order blocks) ='/ (8X,10I6))
 9090       FORMAT (' MBLOCK (triangular flag and number packed rows) ='
     +             / (8X,10I6))
 9100       FORMAT (' KBLOCK (position of beginning of block) ='/
     +             (8X,10I6))
            IF (TRANS) THEN
               WRITE (MP,9110) (RHS(I),I=1,N)
 9110          FORMAT (' RHS =  ',1P,5D12.4/ (8X,5D12.4))
            ELSE
               WRITE (MP,9110) (RHS(I),I=1,M)
            END IF
         ELSE
            WRITE (MP,'(A)') ' Contents of KEEP array'
            WRITE (MP,9020) (KEEP(I),I=1,MIN(10,M))
            WRITE (MP,9030) (KEEP(M+I),I=1,MIN(10,N))
            WRITE (MP,'(A)') ' Pointer information from KEEP'
            WRITE (MP,9040) (KEEP(IPTRD+K),K=1,MIN(10,N+1))
            WRITE (MP,9050) (KEEP(IPTRO+K),K=1,MIN(10,N+1))
            WRITE (MP,9060) (KEEP(IPTRL+J),J=1,MIN(10,N))
            WRITE (MP,9070) (KEEP(IPTRU+J),J=1,MIN(10,N))
            WRITE (MP,'(A)') ' Block structure information from KEEP'
            WRITE (MP,9080) (KEEP(NBLOCK+3*JB),JB=1,KB)
            WRITE (MP,9090) (KEEP(MBLOCK+3*JB),JB=1,KB)
            WRITE (MP,9100) (KEEP(KBLOCK+3*JB),JB=1,KB)
            IF (TRANS) THEN
               WRITE (MP,9110) (RHS(I),I=1,MIN(10,N))
            ELSE
               WRITE (MP,9110) (RHS(I),I=1,MIN(10,M))
            END IF
         END IF
      END IF
      IF (TRANS) THEN
         NEQ = N
         NVAR = M
         DO 100 I = 1,NEQ
            W(I) = RHS(KEEP(M+I))
  100    CONTINUE
      ELSE
         NEQ = M
         NVAR = N
         DO 110 I = 1,NEQ
            W(KEEP(I)) = RHS(I)
  110    CONTINUE
      END IF
      IF (JOB.EQ.1) THEN
         IF (NB.EQ.1 .AND. KEEP(MBLOCK+3).GE.0) THEN
            CALL MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),TRANS,LA-NE,
     +                  A(NE+1),IRN(NE+1),KEEP(IPTRL+1),KEEP(IPTRU+1),W,
     +                  X,W(NEQ+1),INFO)
         ELSE
            CALL MA48DD
     +           (N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,KEEP(IPTRD+1),
     +            KEEP(IPTRO+1),NB,KEEP(NBLOCK+3),KEEP(IPTRL+1),
     +            KEEP(IPTRU+1),W,X,TRANS,ICNTL5,W(NEQ+1))
         END IF
         GO TO 340
      END IF
      DO 120 I = 1,NVAR
         X(I) = ZERO
  120 CONTINUE
      DO 130 I = 1,NEQ
         RHS(I) = W(I)
  130 CONTINUE
      OM1 = ZERO
      DO 260 K = 1,ICNTL(9)
         IF (NB.EQ.1 .AND. KEEP(MBLOCK+3).GE.0) THEN
            CALL MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),TRANS,LA-NE,
     +                A(NE+1),IRN(NE+1),KEEP(IPTRL+1),KEEP(IPTRU+1),
     +                W,W(NEQ+1),W(M+N+1),INFO)
         ELSE
            CALL MA48DD
     +           (N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,KEEP(IPTRD+1),
     +            KEEP(IPTRO+1),NB,KEEP(NBLOCK+3),KEEP(IPTRL+1),
     +            KEEP(IPTRU+1),W,W(NEQ+1),TRANS,ICNTL5,W(M+N+1))
         END IF
         DO 140 I = 1,NVAR
            X(I) = X(I) + W(NEQ+I)
  140    CONTINUE
         DO 150 I = 1,NEQ
            W(I) = RHS(I)
            W(NEQ+I) = ZERO
            W(2*NEQ+I) = ZERO
  150    CONTINUE
         IF (TRANS) THEN
            DO 180 J = 1,N
CDIR$ IVDEP
               DO 160 JJ = KEEP(IPTRD+J),KEEP(IPTRD+J+1) - 1
                  I = IRN(JJ)
                  W(J) = W(J) - A(JJ)*X(I)
                  W(NEQ+J) = W(NEQ+J) + ABS(A(JJ)*X(I))
                  W(2*NEQ+J) = W(2*NEQ+J) + ABS(A(JJ))
  160          CONTINUE
CDIR$ IVDEP
               DO 170 JJ = KEEP(IPTRO+J),KEEP(IPTRO+J+1) - 1
                  I = IRN(JJ)
                  W(J) = W(J) - A(JJ)*X(I)
                  W(NEQ+J) = W(NEQ+J) + ABS(A(JJ)*X(I))
                  W(2*NEQ+J) = W(2*NEQ+J) + ABS(A(JJ))
  170          CONTINUE
  180       CONTINUE
         ELSE
            DO 210 J = 1,N
CDIR$ IVDEP
               DO 190 JJ = KEEP(IPTRD+J),KEEP(IPTRD+J+1) - 1
                  I = IRN(JJ)
                  W(I) = W(I) - A(JJ)*X(J)
                  W(NEQ+I) = W(NEQ+I) + ABS(A(JJ)*X(J))
                  W(2*NEQ+I) = W(2*NEQ+I) + ABS(A(JJ))
  190          CONTINUE
CDIR$ IVDEP
               DO 200 JJ = KEEP(IPTRO+J),KEEP(IPTRO+J+1) - 1
                  I = IRN(JJ)
                  W(I) = W(I) - A(JJ)*X(J)
                  W(NEQ+I) = W(NEQ+I) + ABS(A(JJ)*X(J))
                  W(2*NEQ+I) = W(2*NEQ+I) + ABS(A(JJ))
  200          CONTINUE
  210       CONTINUE
         END IF
         DXMAX = ZERO
         DO 220 I = 1,NVAR
            DXMAX = MAX(DXMAX,ABS(X(I)))
  220    CONTINUE
         OMEGA(1) = ZERO
         OMEGA(2) = ZERO
         DO 230 I = 1,NEQ
            TAU = (W(2*NEQ+I)*DXMAX+ABS(RHS(I)))*NVAR*CTAU
            IF ((W(NEQ+I)+ABS(RHS(I))).GT.TAU) THEN
               OMEGA(1) = MAX(OMEGA(1),ABS(W(I))/
     +                    (W(NEQ+I)+ABS(RHS(I))))
               IW(I) = 1
            ELSE
               IF (TAU.GT.ZERO) THEN
                  OMEGA(2) = MAX(OMEGA(2),ABS(W(I))/
     +                       (W(NEQ+I)+W(2*NEQ+I)*DXMAX))
               END IF
               IW(I) = 2
            END IF
  230    CONTINUE
         IF (JOB.EQ.2) GO TO 340
         OM2 = OMEGA(1) + OMEGA(2)
         IF ((OM2+ONE).LE.ONE) GO TO 270
         IF (K.GT.1 .AND. OM2.GT.OM1*CNTL(5)) THEN
            IF (OM2.GT.OM1) THEN
               OMEGA(1) = OLDOMG(1)
               OMEGA(2) = OLDOMG(2)
               DO 240 I = 1,NVAR
                  X(I) = W(3*NEQ+I)
  240          CONTINUE
            END IF
            GO TO 270
         END IF
         DO 250 I = 1,NVAR
            W(3*NEQ+I) = X(I)
  250    CONTINUE
         OLDOMG(1) = OMEGA(1)
         OLDOMG(2) = OMEGA(2)
         OM1 = OM2
  260 CONTINUE
      INFO(1) = -8
      IF (LP.GT.0 .AND. ICNTL(3).GE.1) WRITE (LP,9170) INFO(1),ICNTL(9)
      GO TO 340
  270 IF (JOB.LE.3) GO TO 340
      IF (M.NE.N) GO TO 340
      LCOND(1) = .FALSE.
      LCOND(2) = .FALSE.
      DO 280 I = 1,NEQ
         IF (IW(I).EQ.1) THEN
            W(I) = W(NEQ+I) + ABS(RHS(I))
            W(NEQ+I) = ZERO
            LCOND(1) = .TRUE.
         ELSE
            W(NEQ+I) = W(NEQ+I) + W(2*NEQ+I)*DXMAX
            W(I) = ZERO
            LCOND(2) = .TRUE.
         END IF
  280 CONTINUE
      KASE = 0
      DO 330 K = 1,2
         IF (LCOND(K)) THEN
            DO 310 KK = 1,40
               CALL MC41AD(N,KASE,W(3*NEQ+1),COND(K),RHS,IW)
               IF (KASE.EQ.0) GO TO 320
               IF (KASE.EQ.1) THEN
                  IF (NB.EQ.1 .AND. KEEP(MBLOCK+3).GE.0) THEN
                     CALL MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),
     +                           .NOT.TRANS,LA-NE,A(NE+1),IRN(NE+1),
     +                           KEEP(IPTRL+1),KEEP(IPTRU+1),W(3*NEQ+1),
     +                           W(2*NEQ+1),RHS,INFO)
                  ELSE
                     CALL MA48DD(N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,
     +                           KEEP(IPTRD+1),KEEP(IPTRO+1),NB,
     +                           KEEP(NBLOCK+3),KEEP(IPTRL+1),
     +                           KEEP(IPTRU+1),W(3*NEQ+1),W(2*NEQ+1),
     +                           .NOT.TRANS,ICNTL5,RHS)
                  END IF
                  DO 290 I = 1,M
                     W(3*NEQ+I) = W((K-1)*NEQ+I)*W(2*NEQ+I)
  290             CONTINUE
               END IF
               IF (KASE.EQ.2) THEN
                  DO 300 I = 1,N
                     W(2*NEQ+I) = W((K-1)*NEQ+I)*W(3*NEQ+I)
  300             CONTINUE
                  IF (NB.EQ.1 .AND. KEEP(MBLOCK+3).GE.0) THEN
                     CALL MA50CD(M,N,ICNTL5,IQB,KEEP(MBLOCK+3),TRANS,
     +                           LA-NE,A(NE+1),IRN(NE+1),KEEP(IPTRL+1),
     +                           KEEP(IPTRU+1),W(2*NEQ+1),W(3*NEQ+1),
     +                           RHS,INFO)
                  ELSE
                     CALL MA48DD(N,NE,LA-NE,A(NE+1),A,IRN(NE+1),IRN,
     +                           KEEP(IPTRD+1),KEEP(IPTRO+1),NB,
     +                           KEEP(NBLOCK+3),KEEP(IPTRL+1),
     +                           KEEP(IPTRU+1),W(2*NEQ+1),W(3*NEQ+1),
     +                           TRANS,ICNTL5,RHS)
                  END IF
               END IF
  310       CONTINUE
            INFO(1) = -9
            IF (LP.NE.0 .AND. ICNTL(3).GE.1) WRITE (LP,9160)
            GO TO 340
  320       IF (DXMAX.GT.ZERO) COND(K) = COND(K)/DXMAX
            ERROR(3) = ERROR(3) + OMEGA(K)*COND(K)
         END IF
  330 CONTINUE
  340 DO 350 I = 1,NVAR
         W(I) = X(I)
  350 CONTINUE
      IF (.NOT.TRANS) THEN
         DO 360 I = 1,NVAR
            X(KEEP(M+I)) = W(I)
  360    CONTINUE
      ELSE
         DO 370 I = 1,NVAR
            X(I) = W(KEEP(I))
  370    CONTINUE
      END IF
      IF (JOB.GE.2) THEN
         ERROR(1) = OMEGA(1)
         ERROR(2) = OMEGA(2)
      END IF
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/A,I6)') ' Leaving MA48CD with INFO(1) =',INFO(1)
         IF (JOB.GT.1) THEN
            K = 2
            IF (JOB.EQ.4 .AND. INFO(1).NE.-9) K = 3
            WRITE (MP,9120) (ERROR(I),I=1,K)
 9120       FORMAT (' ERROR =',1P,3D12.4)
         END IF
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,9130) (X(I),I=1,NVAR)
 9130       FORMAT (' X =    ',1P,5D12.4:/ (8X,5D12.4))
         ELSE
            WRITE (MP,9130) (X(I),I=1,MIN(10,NVAR))
         END IF
      END IF
  380 RETURN
 9140 FORMAT (' Error return from MA48C/CD because M =',I10,' and N =',
     +       I10)
 9150 FORMAT (' Error return from MA48C/CD because ','JOB = ',I10)
 9160 FORMAT (' Error return from MA48C/CD because of ','error in MC41',
     +       'A/AD'/' ERROR(3) not calculated')
 9170 FORMAT (' Error return from MA48C/CD because of ','nonconvergenc',
     +       'e of iterative refinement'/' Error INFO(1) = ',I2,'  wit',
     +       'h ICNTL','(9) = ',I10)
      END
      SUBROUTINE MA48DD(N,NE,LA,A,AA,IRN,IRNA,IPTRD,IPTRO,NB,IBLOCK,
     +                  IPTRL,IPTRU,RHS,X,TRANS,ICNTL5,W)
      INTEGER N,NE,LA
      DOUBLE PRECISION A(LA),AA(NE)
      INTEGER IRN(LA),IRNA(NE),IPTRD(N+1),IPTRO(N+1),NB
      INTEGER IBLOCK(3,NB),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION RHS(N),X(N)
      LOGICAL TRANS
      INTEGER ICNTL5(7)
      DOUBLE PRECISION W(N)
      INTEGER I,IFLAG(7),IQB(1),J,JB,JJ,J1,K1,K2,NC,NUMB
      EXTERNAL MA50CD
      IQB(1) = 0
      NUMB = IBLOCK(3,1)
      IF (.NOT.TRANS) THEN
         K1 = N + 1
         DO 50 JB = NUMB,1,-1
            NC = IBLOCK(1,JB)
            K2 = K1 - 1
            K1 = K1 - NC
            IF (IBLOCK(2,JB).LT.0) THEN
               DO 20 J = K2,K1,-1
                  X(J) = RHS(J)/AA(IPTRD(J))
CDIR$ IVDEP
                  DO 10 JJ = IPTRD(J) + 1,IPTRD(J+1) - 1
                     I = IRNA(JJ)
                     RHS(I) = RHS(I) - AA(JJ)*X(J)
   10             CONTINUE
   20          CONTINUE
            ELSE
               J1 = 1
               IF (JB.GT.1) J1 = IBLOCK(3,JB)
               CALL MA50CD(NC,NC,ICNTL5,IQB,IBLOCK(2,JB),TRANS,LA+1-J1,
     +                     A(J1),IRN(J1),IPTRL(K1),IPTRU(K1),RHS(K1),
     +                     X(K1),W,IFLAG)
            END IF
            IF (JB.EQ.1) GO TO 50
            DO 40 J = K1,K2
CDIR$ IVDEP
               DO 30 JJ = IPTRO(J),IPTRO(J+1) - 1
                  I = IRNA(JJ)
                  RHS(I) = RHS(I) - AA(JJ)*X(J)
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      ELSE
         K2 = 0
         DO 100 JB = 1,NUMB
            NC = IBLOCK(1,JB)
            K1 = K2 + 1
            K2 = K2 + NC
            IF (JB.GT.1) THEN
               DO 70 J = K1,K2
                  DO 60 JJ = IPTRO(J),IPTRO(J+1) - 1
                     I = IRNA(JJ)
                     RHS(J) = RHS(J) - AA(JJ)*X(I)
   60             CONTINUE
   70          CONTINUE
            END IF
            IF (IBLOCK(2,JB).LT.0) THEN
               DO 90 J = K1,K2
                  DO 80 JJ = IPTRD(J) + 1,IPTRD(J+1) - 1
                     I = IRNA(JJ)
                     RHS(J) = RHS(J) - AA(JJ)*X(I)
   80             CONTINUE
                  X(J) = RHS(J)/AA(IPTRD(J))
   90          CONTINUE
            ELSE
               J1 = 1
               IF (JB.GT.1) J1 = IBLOCK(3,JB)
               CALL MA50CD(NC,NC,ICNTL5,IQB,IBLOCK(2,JB),TRANS,LA+1-J1,
     +                     A(J1),IRN(J1),IPTRL(K1),IPTRU(K1),RHS(K1),
     +                     X(K1),W,IFLAG)
            END IF
  100    CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE MA48ID(CNTL,ICNTL)
      DOUBLE PRECISION CNTL(5)
      INTEGER ICNTL(9)
      CNTL(1) = 0.5D0
      CNTL(2) = 0.1D0
      CNTL(3) = 0.0D0
      CNTL(4) = 0.0D0
      CNTL(5) = 0.5D0
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 2
      ICNTL(4) = 3
      ICNTL(5) = 32
      ICNTL(6) = 1
      ICNTL(7) = 1
      ICNTL(8) = 0
      ICNTL(9) = 10
      END


