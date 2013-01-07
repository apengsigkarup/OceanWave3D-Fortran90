C COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
C######DATE 20 May 1993
C 24 September 1993 Some IVDEP comments added for speed on the Cray,
C     some minor bugs fixed, default changed to BLAS 3 with block size
C     32.
C 6 December 1993. Minor bug fixed re threshold test for pivots.
C 4/10/95. IQ in MA50BD made assumed size
C 1/11/95. IQ in MA50CD made assumed size
C 1/11/95. DTRSV not called for zero-sized array.
C 14/2/96. NP initialized to 0.
C 13/11/97 INFO(4) and INFO(6) in MA50AD made to reflect the situation
C          at the point of failure in the case of insufficient storage.
C 17/3/98  In MA50AD, copy a row forward if there is space at its front,
C          rather than put new entry at front. Makes the result
C          repeatable as LA is altered.
      SUBROUTINE MA50AD(M,N,NE,LA,A,IRN,JCN,IQ,CNTL,ICNTL,IP,NP,JFIRST,
     +                  LENR,LASTR,NEXTR,IW,IFIRST,LENC,LASTC,NEXTC,
     +                  INFO,RINFO)
      INTEGER M,N,NE,LA
      DOUBLE PRECISION A(LA)
      INTEGER IRN(LA),JCN(LA),IQ(N)
      DOUBLE PRECISION CNTL(4)
      INTEGER ICNTL(7),IP(M),NP,JFIRST(M),LENR(M),LASTR(M),NEXTR(M),
     +        IW(M),IFIRST(N),LENC(N),LASTC(N),NEXTC(N),INFO(7)
      DOUBLE PRECISION RINFO
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      EXTERNAL MA50DD
      INTRINSIC ABS,MAX,MIN
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0D0,ONE=1.0D0)
      DOUBLE PRECISION ALEN,AMULT,ANEW,ASW,AU,COST,CPIV
      INTEGER DISPC,DISPR,EYE,I,IDROP,IDUMMY,IEND,IFILL,IFIR,II,IJ,
     +        IJPOS,IOP,IPIV,IPOS,ISRCH,I1,I2,J,JBEG,JEND,JJ,JLAST,
     +        JMORE,JNEW,JPIV,JPOS,J1,J2,L,LC,LEN,LENPIV,LP,LR
      DOUBLE PRECISION MAXENT
      INTEGER MINC,MORD,MP,MSRCH,NC,NDROP,NEFACT,NEPR,NERED,NE1,NORD,
     +        NORD1,NR,NULLC,NULLI,NULLJ,NULLR,PIVBEG,PIVCOL,PIVEND,
     +        PIVOT
      DOUBLE PRECISION PIVR,PIVRAT,U
      LP = ICNTL(1)
      IF (ICNTL(3).LE.0) LP = 0
      MP = ICNTL(2)
      IF (ICNTL(3).LE.1) MP = 0
      INFO(1) = 0
      INFO(2) = 0
      INFO(3) = NE
      INFO(4) = NE
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
      RINFO = ZERO
      NP = 0
      IF (M.LT.1 .OR. N.LT.1) GO TO 690
      IF (NE.LT.1) GO TO 700
      IF (LA.LT.NE) THEN
         INFO(3) = NE
         GO TO 710
      END IF
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(/2(A,I6),A,I8,A,I8/A,1P,4E10.2/A,7I4)')
     +     ' Entering MA50AD with M =',M,' N =',N,' NE =',NE,' LA =',LA,
     +     ' CNTL =',CNTL,' ICNTL =',ICNTL
         IF (N.EQ.1 .OR. ICNTL(3).GT.3) THEN
            DO 10 J = 1,N - 1
               IF (IQ(J).LT.IQ(J+1)) WRITE (MP,
     +             '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J,
     +             (A(II),IRN(II),II=IQ(J),IQ(J+1)-1)
   10       CONTINUE
            IF (IQ(N).LE.NE) WRITE (MP,'(A,I5,(T13,3(1P,E12.4,I5)))')
     +          ' Column',N, (A(II),IRN(II),II=IQ(N),NE)
         ELSE
            IF (IQ(1).LT.IQ(2)) WRITE (MP,
     +          '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',1,
     +          (A(II),IRN(II),II=IQ(1),IQ(2)-1)
         END IF
         IF (ICNTL(7).EQ.2) THEN
            WRITE (MP,'(A,(T10,10(I7)))') ' IP = ',IP
            WRITE (MP,'(A,(T10,10(I7)))') ' IFIRST = ',IFIRST
         END IF
      END IF
      MINC = 1
      NERED = NE
      U = MIN(CNTL(2),ONE)
      U = MAX(U,ZERO)
      MSRCH = ICNTL(4)
      IF (MSRCH.EQ.0) MSRCH = N
      JLAST = N - ICNTL(6)
      IF (JLAST.LT.1 .OR. JLAST.GT.N) JLAST = N
      NULLI = 0
      NULLJ = 0
      MORD = 0
      NORD = 0
      NDROP = 0
      NEFACT = 0
      DO 20 I = 1,N - 1
         LENC(I) = IQ(I+1) - IQ(I)
   20 CONTINUE
      LENC(N) = NE + 1 - IQ(N)
      IF (CNTL(3).GT.ZERO) THEN
         NERED = 0
         DO 40 J = 1,N
            I = IQ(J)
            IQ(J) = NERED + 1
            DO 30 II = I,I + LENC(J) - 1
               IF (ABS(A(II)).GE.CNTL(3)) THEN
                  NERED = NERED + 1
                  A(NERED) = A(II)
                  IRN(NERED) = IRN(II)
               ELSE
                  INFO(6) = INFO(6) + 1
               END IF
   30       CONTINUE
            LENC(J) = NERED + 1 - IQ(J)
   40    CONTINUE
      END IF
      IF (ICNTL(7).EQ.2) THEN
         DO 50 I = 1,M
            NEXTR(I) = IP(I)
   50    CONTINUE
         IF (ICNTL(4).NE.1) GO TO 740
      END IF
      DISPR = NERED + 1
      DISPC = NERED + 1
      DO 60 I = 1,M
         IW(I) = 0
         LENR(I) = 0
         JFIRST(I) = 0
   60 CONTINUE
      DO 70 II = 1,NERED
         I = IRN(II)
         LENR(I) = LENR(I) + 1
   70 CONTINUE
      IP(1) = LENR(1) + 1
      DO 80 I = 2,M
         IP(I) = IP(I-1) + LENR(I)
   80 CONTINUE
      DO 100 J = 1,N
         I = IQ(J)
         DO 90 II = I,I + LENC(J) - 1
            I = IRN(II)
            IF (IW(I).EQ.J) GO TO 720
            IW(I) = J
            IPOS = IP(I) - 1
            JCN(IPOS) = J
            IP(I) = IPOS
   90    CONTINUE
  100 CONTINUE
      DO 110 I = 1,M
         IW(I) = 0
  110 CONTINUE
      IF (ICNTL(4).LE.0) THEN
         DO 120 I = 1,N
            IFIRST(I) = 0
  120    CONTINUE
         DO 130 I = M,1,-1
            NE1 = LENR(I)
            IF (NE1.GT.0) THEN
               IFIR = IFIRST(NE1)
               IFIRST(NE1) = I
               LASTR(I) = 0
               NEXTR(I) = IFIR
               IF (IFIR.GT.0) LASTR(IFIR) = I
            ELSE
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            END IF
  130    CONTINUE
      ELSE
         DO 140 I = M,1,-1
            NE1 = LENR(I)
            IF (NE1.EQ.0) THEN
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            END IF
  140    CONTINUE
      END IF
      DO 150 J = N,1,-1
         NE1 = LENC(J)
         IF (NE1.EQ.0) THEN
            IF (ICNTL(7).NE.2) THEN
               IF (J.LE.JLAST) THEN
                  NORD = NORD + 1
                  IQ(J) = -NORD
                  IF (NORD.EQ.JLAST) THEN
                     NORD = NORD + NULLJ
                     JLAST = N
                     NULLJ = 0
                  END IF
               ELSE
                  NULLJ = NULLJ + 1
                  IQ(J) = - (JLAST+NULLJ)
               END IF
               LASTC(J) = 0
               NEXTC(J) = 0
            END IF
         ELSE
            IFIR = JFIRST(NE1)
            JFIRST(NE1) = J
            NEXTC(J) = IFIR
            LASTC(J) = 0
            IF (IFIR.GT.0) LASTC(IFIR) = J
         END IF
  150 CONTINUE
      IF (INFO(6).EQ.0) THEN
         NULLC = NORD + NULLJ
         NULLR = NULLI
      END IF
      DO 630 PIVOT = 1,N
         IF (NERED.GE. (MIN(CNTL(1),ONE)*(N-NORD))*
     +       (M-MORD)) GO TO 640
         IF (ICNTL(7).EQ.2) THEN
            IPIV = 0
            J = IFIRST(PIVOT)
            IF (J.LT.1 .OR. J.GT.N) GO TO 730
            IF (IQ(J).LT.0) GO TO 730
            LEN = LENC(J)
            IF (LEN.LE.0) GO TO 320
            ALEN = LEN - 1
            I1 = IQ(J)
            I2 = I1 + LEN - 1
            II = IDAMAX(LEN,A(I1),1)
            MAXENT = ABS(A(I1+II-1))
            IF (MAXENT.LE.CNTL(4)) GO TO 320
            AU = MAX(MAXENT*U,CNTL(4))
            DO 160 II = I1,I2
               IF (ABS(A(II)).LT.AU) GO TO 160
               I = IRN(II)
               IF (IPIV.NE.0) THEN
                  IF (NEXTR(I).GE.NEXTR(IPIV)) GO TO 160
               END IF
               CPIV = ALEN*(LENR(I)-1)
               IJPOS = II
               IPIV = I
               JPIV = J
  160       CONTINUE
            GO TO 330
         END IF
         LEN = MINC
         DO 170 MINC = LEN,M - MORD
            IF (JFIRST(MINC).NE.0) GO TO 180
            IF (ICNTL(4).LE.0) THEN
               IF (IFIRST(MINC).NE.0) GO TO 180
            END IF
  170    CONTINUE
  180    CPIV = M
         CPIV = CPIV*N
         PIVRAT = ZERO
         ISRCH = 0
         DO 300 LEN = MINC,M - MORD
            ALEN = LEN - 1
            IF (CPIV.LE.ALEN**2 .AND. ICNTL(4).LE.0) GO TO 310
            IJ = JFIRST(LEN)
            DO 220 IDUMMY = 1,N
               IF (IJ.LE.0) GO TO 230
               J = IJ
               IJ = NEXTC(J)
               IF (J.GT.JLAST) GO TO 220
               MAXENT = ZERO
               I1 = IQ(J)
               I2 = I1 + LEN - 1
               II = IDAMAX(LEN,A(I1),1)
               MAXENT = ABS(A(I1+II-1))
               IF (MAXENT.LE.CNTL(4)) GO TO 320
               AU = MAX(MAXENT*U,CNTL(4))
               IF (ICNTL(7).EQ.1) THEN
                  DO 190 II = I1,I2
                     IF (IRN(II).EQ.J) GO TO 200
  190             CONTINUE
                  GO TO 220
  200             I1 = II
                  I2 = II
               END IF
               DO 210 II = I1,I2
                  IF (ABS(A(II)).LT.AU) GO TO 210
                  I = IRN(II)
                  COST = ALEN*(LENR(I)-1)
                  IF (COST.GT.CPIV) GO TO 210
                  PIVR = ABS(A(II))/MAXENT
                  IF (COST.EQ.CPIV) THEN
                     IF (PIVR.LE.PIVRAT) GO TO 210
                  END IF
                  CPIV = COST
                  IJPOS = II
                  IPIV = I
                  JPIV = J
                  IF (CPIV.LE.ALEN**2 .AND. ICNTL(4).LE.0) GO TO 330
                  PIVRAT = PIVR
  210          CONTINUE
               ISRCH = ISRCH + 1
               IF (ISRCH.GE.MSRCH) THEN
                  IF (PIVRAT.GT.ZERO) GO TO 330
               END IF
  220       CONTINUE
  230       IF (ICNTL(4).GT.0) GO TO 300
            IF (CPIV.LE.ALEN*(ALEN+1)) GO TO 310
            IF (LEN.GT.N-NORD) GO TO 300
            IJ = IFIRST(LEN)
            DO 290 IDUMMY = 1,M
               IF (IJ.EQ.0) GO TO 300
               I = IJ
               IJ = NEXTR(IJ)
               J1 = IP(I)
               J2 = J1 + LEN - 1
               IF (ICNTL(7).EQ.1) THEN
                  DO 240 JJ = J1,J2
                     IF (JCN(JJ).EQ.I) GO TO 250
  240             CONTINUE
                  GO TO 290
  250             J1 = JJ
                  J2 = JJ
               END IF
               DO 280 JJ = J1,J2
                  J = JCN(JJ)
                  IF (J.GT.JLAST) GO TO 280
                  COST = ALEN*(LENC(J)-1)
                  IF (COST.GE.CPIV) GO TO 280
                  I1 = IQ(J)
                  I2 = I1 + LENC(J) - 1
                  II = IDAMAX(LENC(J),A(I1),1)
                  MAXENT = ABS(A(I1+II-1))
                  DO 260 II = I1,I2 - 1
                     IF (IRN(II).EQ.I) GO TO 270
  260             CONTINUE
  270             JPOS = II
                  IF (MAXENT.LE.CNTL(4)) GO TO 320
                  IF (ABS(A(JPOS)).LT.MAXENT*U) GO TO 280
                  CPIV = COST
                  IPIV = I
                  JPIV = J
                  IJPOS = JPOS
                  PIVRAT = ABS(A(JPOS))/MAXENT
                  IF (CPIV.LE.ALEN*(ALEN+1)) GO TO 330
  280          CONTINUE
  290       CONTINUE
  300    CONTINUE
  310    IF (PIVRAT.GT.ZERO) GO TO 330
         INFO(1) = INFO(1) + 2
         IF (MP.GT.0) WRITE (MP,'(A/A)')
     +       ' Warning message from MA50AD: no suitable diagonal pivot',
     +       ' found, so switched to full matrix processing.'
         GO TO 640
  320    IPIV = 0
         JPIV = J
  330    NEFACT = NEFACT + LENC(JPIV)
         PIVBEG = IQ(JPIV)
         PIVEND = PIVBEG + LENC(JPIV) - 1
         NORD = NORD + 1
         NORD1 = NORD
         IF (NORD.EQ.JLAST) THEN
            NORD = NORD + NULLJ
            JLAST = N
            NULLJ = 0
         END IF
         IF (ICNTL(4).LE.0) THEN
            DO 340 II = PIVBEG,PIVEND
               I = IRN(II)
               LR = LASTR(I)
               NR = NEXTR(I)
               IF (NR.NE.0) LASTR(NR) = LR
               IF (LR.EQ.0) THEN
                  NE1 = LENR(I)
                  IFIRST(NE1) = NR
               ELSE
                  NEXTR(LR) = NR
               END IF
  340       CONTINUE
         END IF
         IF (IPIV.GT.0) THEN
            NEPR = LENR(IPIV) - 1
            NEFACT = NEFACT + NEPR
            RINFO = RINFO + CPIV*2 + LENR(IPIV)
            J1 = IP(IPIV)
            DO 350 JJ = J1,J1 + NEPR
               J = JCN(JJ)
               LC = LASTC(J)
               NC = NEXTC(J)
               IF (NC.NE.0) LASTC(NC) = LC
               IF (LC.EQ.0) THEN
                  NE1 = LENC(J)
                  JFIRST(NE1) = NC
               ELSE
                  NEXTC(LC) = NC
               END IF
  350       CONTINUE
            IF (PIVBEG.NE.IJPOS) THEN
               ASW = A(PIVBEG)
               A(PIVBEG) = A(IJPOS)
               A(IJPOS) = ASW
               IRN(IJPOS) = IRN(PIVBEG)
               IRN(PIVBEG) = IPIV
            END IF
         ELSE
            NEPR = 0
            NE1 = LENC(JPIV)
            IF (CNTL(3).GT.ZERO) NDROP = NDROP + NE1
            IF (NE1.GT.0) THEN
               LC = LASTC(JPIV)
               NC = NEXTC(JPIV)
               IF (NC.NE.0) LASTC(NC) = LC
               IF (LC.EQ.0) THEN
                  JFIRST(NE1) = NC
               ELSE
                  NEXTC(LC) = NC
               END IF
            END IF
         END IF
         DO 360 II = PIVBEG + 1,PIVEND
            I = IRN(II)
            IW(I) = II - PIVBEG
  360    CONTINUE
         LENPIV = PIVEND - PIVBEG
         DO 390 II = PIVBEG,PIVEND
            I = IRN(II)
            LENR(I) = LENR(I) - 1
            J1 = IP(I)
            J2 = J1 + LENR(I)
            DO 370 JJ = J1,J2 - 1
               IF (JCN(JJ).EQ.JPIV) GO TO 380
  370       CONTINUE
  380       JCN(JJ) = JCN(J2)
            JCN(J2) = 0
  390    CONTINUE
         DO 600 EYE = 1,NEPR
            J = JCN(IP(IPIV)+EYE-1)
            IDROP = 0
            JBEG = IQ(J)
            JEND = JBEG + LENC(J) - 1
            DO 400 II = JBEG,JEND - 1
               IF (IRN(II).EQ.IPIV) GO TO 410
  400       CONTINUE
  410       AMULT = -A(II)/A(IQ(JPIV))
            A(II) = A(JEND)
            IRN(II) = IRN(JEND)
            LENC(J) = LENC(J) - 1
            IRN(JEND) = 0
            JEND = JEND - 1
            IF (LENPIV.EQ.0) GO TO 600
            IOP = 0
CDIR$ IVDEP
            DO 420 II = JBEG,JEND
               I = IRN(II)
               IF (IW(I).GT.0) THEN
                  IOP = IOP + 1
                  PIVCOL = IQ(JPIV) + IW(I)
                  IW(I) = -IW(I)
                  A(II) = A(II) + AMULT*A(PIVCOL)
               END IF
  420       CONTINUE
            IF (CNTL(3).GT.ZERO) THEN
               JNEW = JBEG
               DO 450 II = JBEG,JEND
                  IF (ABS(A(II)).GE.CNTL(3)) THEN
                     A(JNEW) = A(II)
                     IRN(JNEW) = IRN(II)
                     JNEW = JNEW + 1
                  ELSE
                     I = IRN(II)
                     J1 = IP(I)
                     J2 = J1 + LENR(I) - 1
                     DO 430 JJ = J1,J2 - 1
                        IF (JCN(JJ).EQ.J) GO TO 440
  430                CONTINUE
  440                JCN(JJ) = JCN(J2)
                     JCN(J2) = 0
                     LENR(I) = LENR(I) - 1
                  END IF
  450          CONTINUE
               DO 460 II = JNEW,JEND
                  IRN(II) = 0
  460          CONTINUE
               IDROP = JEND + 1 - JNEW
               JEND = JNEW - 1
               LENC(J) = LENC(J) - IDROP
               NERED = NERED - IDROP
               INFO(6) = INFO(6) + IDROP
            END IF
            IFILL = LENPIV - IOP
            NERED = NERED + IFILL
            INFO(3) = MAX(INFO(3),NERED+LENC(J))
            IF (IFILL.EQ.0) THEN
CDIR$ IVDEP
               DO 470 II = PIVBEG + 1,PIVEND
                  I = IRN(II)
                  IW(I) = -IW(I)
  470          CONTINUE
               GO TO 600
            END IF
            DO 480 IPOS = JEND + 1,MIN(JEND+IFILL,DISPC-1)
               IF (IRN(IPOS).NE.0) GO TO 490
  480       CONTINUE
            IF (IPOS.EQ.JEND+IFILL+1) GO TO 540
            IF (JEND+IFILL+1.LE.LA+1) THEN
               DISPC = JEND + IFILL + 1
               GO TO 540
            END IF
            IPOS = LA
            DISPC = LA + 1
  490       JMORE = JEND + IFILL - IPOS + 1
            DO 500 IPOS = JBEG - 1,MAX(JBEG-JMORE,1),-1
               IF (IRN(IPOS).NE.0) GO TO 510
  500       CONTINUE
            IPOS = IPOS + 1
            IF (IPOS.EQ.JBEG-JMORE) GO TO 520
  510       IF (DISPC+LENC(J)+IFILL.GT.LA+1) THEN
               INFO(2) = INFO(2) + 1
               CALL MA50DD(LA,A,IRN,IQ,N,DISPC,.TRUE.)
               JBEG = IQ(J)
               JEND = JBEG + LENC(J) - 1
               PIVBEG = IQ(JPIV)
               PIVEND = PIVBEG + LENC(JPIV) - 1
               IF (DISPC+LENC(J)+IFILL.GT.LA+1) GO TO 705
            END IF
            IPOS = DISPC
            DISPC = DISPC + LENC(J) + IFILL
  520       IQ(J) = IPOS
            DO 530 II = JBEG,JEND
               A(IPOS) = A(II)
               IRN(IPOS) = IRN(II)
               IPOS = IPOS + 1
               IRN(II) = 0
  530       CONTINUE
            JBEG = IQ(J)
            JEND = IPOS - 1
  540       IDROP = 0
            DO 580 II = PIVBEG + 1,PIVEND
               I = IRN(II)
               INFO(3) = MAX(INFO(3),NERED+LENR(I)+1)
               IF (IW(I).LT.0) THEN
                  IW(I) = -IW(I)
                  GO TO 580
               END IF
               ANEW = AMULT*A(II)
               IF (ABS(ANEW).LT.CNTL(3)) THEN
                  IDROP = IDROP + 1
               ELSE
                  JEND = JEND + 1
                  A(JEND) = ANEW
                  IRN(JEND) = I
                  IEND = IP(I) + LENR(I)
                  IF (IEND.LT.DISPR) THEN
                     IF (JCN(IEND).EQ.0) GO TO 560
                  ELSE
                     IF (DISPR.LE.LA) THEN
                        DISPR = DISPR + 1
                        GO TO 560
                     END IF
                  END IF
                  IF (IP(I).GT.1) THEN
                     IF (JCN(IP(I)-1).EQ.0) THEN
                        IEND = IEND - 1
                        DO 545 JJ = IP(I),IEND
                           JCN(JJ-1) = JCN(JJ)
  545                   CONTINUE
                        IP(I) = IP(I) - 1
                        GO TO 560
                     END IF
                  END IF
                  IF (DISPR+LENR(I).GT.LA) THEN
                     INFO(2) = INFO(2) + 1
                     CALL MA50DD(LA,A,JCN,IP,M,DISPR,.FALSE.)
                     IF (DISPR+LENR(I).GT.LA) GO TO 705
                  END IF
                  J1 = IP(I)
                  J2 = IP(I) + LENR(I) - 1
                  IP(I) = DISPR
                  DO 550 JJ = J1,J2
                     JCN(DISPR) = JCN(JJ)
                     JCN(JJ) = 0
                     DISPR = DISPR + 1
  550             CONTINUE
                  IEND = DISPR
                  DISPR = IEND + 1
  560             JCN(IEND) = J
                  LENR(I) = LENR(I) + 1
               END IF
  580       CONTINUE
            INFO(6) = INFO(6) + IDROP
            NERED = NERED - IDROP
            DO 590 II = 1,IDROP
               IRN(JEND+II) = 0
  590       CONTINUE
            LENC(J) = LENC(J) + IFILL - IDROP
  600    CONTINUE
         DO 610 EYE = 1,NEPR
            JJ = IP(IPIV) + EYE - 1
            J = JCN(JJ)
            JCN(JJ) = 0
            NE1 = LENC(J)
            LASTC(J) = 0
            IF (NE1.GT.0) THEN
               IFIR = JFIRST(NE1)
               JFIRST(NE1) = J
               NEXTC(J) = IFIR
               IF (IFIR.NE.0) LASTC(IFIR) = J
               MINC = MIN(MINC,NE1)
            ELSE IF (ICNTL(7).NE.2) THEN
               IF (INFO(6).EQ.0) NULLC = NULLC + 1
               IF (J.LE.JLAST) THEN
                  NORD = NORD + 1
                  IQ(J) = -NORD
                  IF (NORD.EQ.JLAST) THEN
                     NORD = NORD + NULLJ
                     JLAST = N
                     NULLJ = 0
                  END IF
               ELSE
                  NULLJ = NULLJ + 1
                  IQ(J) = - (JLAST+NULLJ)
               END IF
            END IF
  610    CONTINUE
         NERED = NERED - NEPR
         IF (IPIV.NE.0) THEN
            LENR(IPIV) = 0
            IW(IPIV) = 0
            IRN(PIVBEG) = 0
            MORD = MORD + 1
            PIVBEG = PIVBEG + 1
            IP(IPIV) = -MORD
         END IF
         NERED = NERED - LENPIV - 1
         DO 620 II = PIVBEG,PIVEND
            I = IRN(II)
            IW(I) = 0
            IRN(II) = 0
            NE1 = LENR(I)
            IF (NE1.EQ.0) THEN
               IF (INFO(6).EQ.0) NULLR = NULLR + 1
               IP(I) = -M + NULLI
               NULLI = NULLI + 1
            ELSE IF (ICNTL(4).LE.0) THEN
               IFIR = IFIRST(NE1)
               LASTR(I) = 0
               NEXTR(I) = IFIR
               IFIRST(NE1) = I
               IF (IFIR.NE.0) LASTR(IFIR) = I
               MINC = MIN(MINC,NE1)
            END IF
  620    CONTINUE
         IQ(JPIV) = -NORD1
  630 CONTINUE
  640 INFO(5) = MORD + MIN(M-MORD-NULLI,N-NORD-NULLJ)
      DO 650 L = 1,MIN(M-MORD,N-NORD)
         RINFO = RINFO + M - MORD - L + 1 + REAL(M-MORD-L)*(N-NORD-L)*2
  650 CONTINUE
      NP = NORD
      INFO(4) = 2 + NEFACT + M*2 + MAX(N-NORD+M-MORD,
     +          (N-NORD)*(M-MORD))
      INFO(6) = INFO(6) + NDROP
      INFO(7) = M - MORD
      DO 660 L = 1,M
         IF (IP(L).LT.0) THEN
            IP(L) = -IP(L)
         ELSE
            MORD = MORD + 1
            IP(L) = MORD
         END IF
  660 CONTINUE
      DO 670 L = 1,N
         IF (IQ(L).LT.0) THEN
            LASTC(L) = -IQ(L)
         ELSE
            IF (NORD.EQ.JLAST) NORD = NORD + NULLJ
            NORD = NORD + 1
            LASTC(L) = NORD
         END IF
  670 CONTINUE
      DO 680 L = 1,N
         IQ(LASTC(L)) = L
  680 CONTINUE
      IF (INFO(5).LT.MIN(M,N)) INFO(1) = INFO(1) + 1
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(A,I6,A,F12.1/A,7I8)') ' Leaving MA50AD with NP =',
     +     NP,' RINFO =',RINFO,' INFO =',INFO
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            WRITE (MP,'(A,(T6,10(I7)))') ' IQ = ',IQ
         END IF
      END IF
      GO TO 750
  690 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,'(/A/(2(A,I8)))')
     +    ' **** Error return from MA50AD ****',' M =',M,' N =',N
      GO TO 750
  700 INFO(1) = -2
      IF (LP.GT.0) WRITE (LP,'(/A/(A,I10))')
     +    ' **** Error return from MA50AD ****',' NE =',NE
      GO TO 750
  705 INFO(4) =  NEFACT + NERED
      INFO(6) = INFO(6) + NDROP
  710 INFO(1) = -3
      IF (LP.GT.0) WRITE (LP,'(/A/A,I9,A,I9)')
     +    ' **** Error return from MA50AD ****',
     +    ' LA  must be increased from',LA,' to at least',INFO(3)
      GO TO 750
  720 INFO(1) = -4
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50AD ****',' Entry in row',I,
     +    ' and column',J,' duplicated'
      GO TO 750
  730 INFO(1) = -5
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50AD ****',' Fault in component ',
     +    PIVOT,' of column permutation given in IFIRST'
      GO TO 750
  740 INFO(1) = -6
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50AD ****',' ICNTL(4) = ',ICNTL(4),
     +    ' when ICNTL(6) = 2'
  750 END
      SUBROUTINE MA50BD(M,N,NE,JOB,AA,IRNA,IPTRA,CNTL,ICNTL,IP,IQ,NP,
     +                  LFACT,FACT,IRNF,IPTRL,IPTRU,W,IW,INFO,RINFO)
      INTEGER M,N,NE,JOB
      DOUBLE PRECISION AA(NE)
      INTEGER IRNA(NE),IPTRA(N)
      DOUBLE PRECISION CNTL(4)
      INTEGER ICNTL(7),IP(M),IQ(*),NP,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION W(M)
      INTEGER IW(M+2*N),INFO(7)
      DOUBLE PRECISION RINFO
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0D0,ONE=1.0D0)
      DOUBLE PRECISION AMULT,ASW
      INTEGER BEGCOL
      LOGICAL DROP
      INTEGER ENDCOL,EYE,EYE1,I,IA1,IA2,IF1,IF2,II,IL1,IL2,IPIV,IQPIV,
     +        IU1,IU2,ISW,J,JDUMMY,JJ,JLAST,K,LP
      DOUBLE PRECISION MAXENT
      INTEGER MF,MORD,MP,NEU,NF,NULLC
      DOUBLE PRECISION PIVLIM
      INTEGER RANK
      DOUBLE PRECISION U
      EXTERNAL MA50ED,MA50FD,MA50GD
      INTRINSIC ABS,MAX,MIN
      INFO(1) = 0
      INFO(4) = 0
      INFO(5) = 0
      INFO(6) = 0
      INFO(7) = 0
      RINFO = ZERO
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (ICNTL(3).LE.0) LP = 0
      IF (ICNTL(3).LE.1) MP = 0
      IF (M.LT.1 .OR. N.LT.1) THEN
         INFO(1) = -1
         IF (LP.GT.0) WRITE (LP,'(/A/A,I8,A,I8)')
     +       ' **** Error return from MA50BD ****',' M =',M,' N =',N
         GO TO 550
      END IF
      IF (NE.LE.0) THEN
         INFO(1) = -2
         IF (LP.GT.0) WRITE (LP,'(/A/A,I6)')
     +       ' **** Error return from MA50BD ****',' NE =',NE
         GO TO 550
      END IF
      IF (NP.LT.0 .OR. NP.GT.N) THEN
         INFO(1) = -7
         IF (LP.GT.0) WRITE (LP,'(/A/A,I8,A,I8)')
     +       ' **** Error return from MA50BD ****',' NP =',NP,' N =',N
         GO TO 550
      END IF
      IF (LFACT.LT.MAX(M,NE+2)) THEN
         INFO(4) = MAX(M,NE+2)
         GO TO 520
      END IF
      IF (JOB.EQ.1) THEN
      ELSE IF (JOB.EQ.2 .OR. JOB.EQ.3) THEN
         IF (IRNF(1).NE.0) THEN
            INFO(1) = -6
            IF (LP.GT.0) WRITE (LP,'(/A/A,I1,A)')
     +          ' **** Error return from MA50BD ***',' Call with JOB=',
     +          JOB,' follows JOB=1 call in which entries were dropped'
            GO TO 550
         END IF
      ELSE
         INFO(1) = -5
         IF (LP.GT.0) WRITE (LP,'(/A/A,I2)')
     +       ' **** Error return from MA50BD ****',' JOB =',JOB
         GO TO 550
      END IF
      IF (MP.GT.0) THEN
         IF (ICNTL(3).GT.2) WRITE (MP,
     +       '(/2(A,I6),A,I8,A,I3/A,I8,A,I7/A,1P,4E10.2/A,7I8)')
     +       ' Entering MA50BD with M =',M,' N =',N,' NE =',NE,' JOB =',
     +       JOB,' LFACT =',LFACT,' NP =',NP,' CNTL =',CNTL,' ICNTL =',
     +       ICNTL
         IF (ICNTL(3).GT.3) THEN
            WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            IF (IQ(1).GT.0) THEN
               WRITE (MP,'(A,(T6,10(I7)))') ' IQ = ',(IQ(J),J=1,N)
            ELSE
               WRITE (MP,'(A,(T6,I7))') ' IQ = ',IQ(1)
            END IF
            DO 10 J = 1,N - 1
               IF (IPTRA(J).LT.IPTRA(J+1)) WRITE (MP,
     +             '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',J,
     +             (AA(II),IRNA(II),II=IPTRA(J),IPTRA(J+1)-1)
   10       CONTINUE
            IF (IPTRA(N).LE.NE) WRITE (MP,
     +          '(A,I5,(T13,3(1P,E12.4,I5)))') ' Column',N,
     +          (AA(II),IRNA(II),II=IPTRA(N),NE)
         END IF
      END IF
      JLAST = 0
      NULLC = 0
      IF (JOB.GT.1 .AND. ICNTL(6).GT.0 .AND.
     +    ICNTL(6).LT.N) JLAST = MIN(NP,N-ICNTL(6))
      U = MIN(CNTL(2),ONE)
      U = MAX(U,ZERO)
      DO 20 I = 1,M
         IW(I+N) = 0
         W(I) = ZERO
   20 CONTINUE
      MORD = 0
      IF1 = LFACT + 1
      IF2 = 0
      NF = N - NP
      MF = 0
      IL2 = 2
      IF (JLAST.GT.0) IL2 = IPTRL(JLAST)
      NEU = 0
      IF (JOB.EQ.2) GO TO 370
      IF (JOB.EQ.3) THEN
         DO 30 J = 1,NP
            IA1 = IPTRU(J) + 1
            IF (IA1.GT.IPTRL(J)) GO TO 30
            IF (J.LE.JLAST) THEN
               MORD = MORD + 1
               IP(IRNF(IA1)) = -J
            ELSE
               IP(IRNF(IA1)) = J
            END IF
   30    CONTINUE
         MF = IRNF(2)
         IA1 = IPTRL(N)
         DO 40 J = 1,MF
            IP(IRNF(IA1+J)) = NP + J
   40    CONTINUE
      END IF
      DO 50 K = 1,JLAST
         IW(M+N+K) = IPTRL(K)
   50 CONTINUE
      DO 310 K = JLAST + 1,N
         DROP = .FALSE.
         IF (K.EQ.NP+1) THEN
            MF = M - MORD
            IF1 = LFACT + 1 - MF
            II = 0
            DO 60 I = 1,M
               IF (IP(I).GT.0) THEN
                  IW(I+N) = N
                  IRNF(IF1+II) = I
                  II = II + 1
                  IP(I) = NP + II
               END IF
   60       CONTINUE
            IF1 = LFACT + 1 - MAX(MF*NF,MF+NF)
            IF2 = IF1 - 1 + MF*MAX(0,JLAST-NP)
         END IF
         J = K
         IF (IQ(1).GT.0) J = IQ(K)
         IA1 = IPTRA(J)
         IA2 = NE
         IF (J.NE.N) IA2 = IPTRA(J+1) - 1
         IU1 = IL2 + 1
         IU2 = IU1 - 1
         IL1 = IF1 - 1 + IA1 - IA2
         IL2 = IL1 - 1
         INFO(4) = MAX(INFO(4),NEU+LFACT-IL1+IU2+M+1)
         IF (IL1-IU2.LE.M) THEN
            IF (INFO(1).NE.-3) THEN
               INFO(1) = -3
               NEU = IL2 + LFACT + 1 - MF - IF1
               IF1 = LFACT + 1 - MF
               IF2 = IF1 - 1
               IL2 = 0
               EYE = 0
               DO 80 J = 1,MIN(K-1,NP)
                  IU2 = IPTRU(J)
                  IPTRU(J) = EYE
                  IL2 = IPTRL(J)
                  NEU = NEU + IU2 - IL2
                  DO 70 II = IU2 + 1,IL2
                     EYE = EYE + 1
                     IRNF(EYE) = IRNF(II)
                     FACT(EYE) = FACT(II)
   70             CONTINUE
                  IPTRL(J) = EYE
                  IW(M+N+J) = EYE
   80          CONTINUE
               IU1 = EYE + 1
               IU2 = EYE
               IL1 = IF1 - 1 + IA1 - IA2
               IL2 = IL1 - 1
            END IF
            IF (IL1-IU2.LE.M) GO TO 480
         END IF
         EYE = IL1
         DO 90 II = IA1,IA2
            I = IRNA(II)
            IF (IW(I+N).EQ.-1) GO TO 540
            IW(I+N) = -1
            W(I) = AA(II)
            IRNF(EYE) = I
            EYE = EYE + 1
   90    CONTINUE
         IPTRL(K) = EYE - 1
         IW(M+N+K) = EYE - 1
         IW(K) = IL1
         J = K
         DO 120 JDUMMY = 1,2*K
            DO 100 II = IW(J),ABS(IW(M+N+J))
               I = IRNF(II)
               IF (IW(I+N).GE.K) GO TO 100
               IF (IP(I).LE.0) GO TO 110
               IW(I+N) = K
               IL1 = IL1 - 1
               IRNF(IL1) = I
  100       CONTINUE
            IF (J.EQ.K) GO TO 130
            IU2 = IU2 + 1
            I = IRNF(IPTRU(J)+1)
            IRNF(IU2) = I
            J = -IW(I+N)
            IW(I+N) = K
            GO TO 120
  110       IW(I+N) = -J
            IW(J) = II + 1
            J = -IP(I)
            IW(J) = IPTRU(J) + 2
  120    CONTINUE
  130    DO 150 II = IU2,IU1,-1
            I = IRNF(II)
            J = -IP(I)
            EYE1 = IPTRU(J) + 1
            IF (ABS(W(I)).LT.CNTL(3)) GO TO 150
            AMULT = -W(I)*FACT(EYE1)
            W(I) = AMULT
            DO 140 EYE = EYE1 + 1,IPTRL(J)
               I = IRNF(EYE)
               W(I) = W(I) + AMULT*FACT(EYE)
  140       CONTINUE
            RINFO = RINFO + ONE + 2*(IPTRL(J)-EYE1)
  150    CONTINUE
         IF (CNTL(3).GT.ZERO) THEN
            EYE = IU1
            DO 160 II = IU1,IU2
               I = IRNF(II)
               IF (ABS(W(I)).LT.CNTL(3)) THEN
                  INFO(6) = INFO(6) + 1
               ELSE
                  IRNF(EYE) = -IP(I)
                  FACT(EYE) = W(I)
                  EYE = EYE + 1
               END IF
               W(I) = ZERO
  160       CONTINUE
            IU2 = EYE - 1
         ELSE
            DO 170 II = IU1,IU2
               I = IRNF(II)
               IRNF(II) = -IP(I)
               FACT(II) = W(I)
               W(I) = ZERO
  170       CONTINUE
         END IF
         IF (INFO(1).EQ.-3) THEN
            NEU = NEU + IU2 - IU1 + 1
            IU2 = IU1 - 1
         END IF
         IPTRU(K) = IU2
         IF (K.LE.NP) THEN
            MAXENT = ZERO
            IF (CNTL(3).GT.ZERO) THEN
               EYE = IL1
               DO 180 II = IL1,IL2
                  I = IRNF(II)
                  IF (ABS(W(I)).LT.CNTL(3)) THEN
                     INFO(6) = INFO(6) + 1
                     W(I) = ZERO
                     DROP = .TRUE.
                  ELSE
                     IRNF(EYE) = I
                     EYE = EYE + 1
                     MAXENT = MAX(ABS(W(I)),MAXENT)
                  END IF
  180          CONTINUE
               IL2 = EYE - 1
            ELSE
               DO 190 II = IL1,IL2
                  MAXENT = MAX(ABS(W(IRNF(II))),MAXENT)
  190          CONTINUE
            END IF
            PIVLIM = U*MAXENT
            EYE = IU2
            IQPIV = M + N
            IF (IL1.GT.IL2) NULLC = NULLC + 1
            DO 200 II = IL1,IL2
               I = IRNF(II)
               EYE = EYE + 1
               IRNF(EYE) = I
               FACT(EYE) = W(I)
               W(I) = ZERO
               IF (ABS(FACT(EYE)).GE.PIVLIM) THEN
                  IF (ABS(FACT(EYE)).GT.CNTL(4)) THEN
                     IF (IP(I).LT.IQPIV) THEN
                        IQPIV = IP(I)
                        IPIV = EYE
                     END IF
                  END IF
               END IF
  200       CONTINUE
            IL1 = IU2 + 1
            IL2 = EYE
            IF (IL1.LE.IL2) THEN
               IF (IQPIV.EQ.M+N) THEN
                  IF (CNTL(3).GT.ZERO) INFO(6) = INFO(6) + EYE - IU2
                  IL2 = IU2
               ELSE
                  IF (IL1.NE.IPIV) THEN
                     ASW = FACT(IPIV)
                     FACT(IPIV) = FACT(IL1)
                     FACT(IL1) = ASW
                     ISW = IRNF(IL1)
                     IRNF(IL1) = IRNF(IPIV)
                     IRNF(IPIV) = ISW
                  END IF
                  INFO(5) = INFO(5) + 1
                  FACT(IL1) = ONE/FACT(IL1)
                  RINFO = RINFO + ONE
                  MORD = MORD + 1
                  IP(IRNF(IL1)) = -K
               END IF
            END IF
         ELSE
            IL2 = IPTRU(K)
CDIR$ IVDEP
            DO 210 II = LFACT - MF + 1,LFACT
               I = IRNF(II)
               IF2 = IF2 + 1
               FACT(IF2) = W(I)
               W(I) = ZERO
  210       CONTINUE
            IF (INFO(1).EQ.-3) IF2 = IF2 - MF
         END IF
         IW(M+N+K) = IL2
         IPTRL(K) = IL2
         IF (DROP) GO TO 310
         DO 300 II = IU1,IU2
            I = IRNF(II)
            IF (IW(M+N+I).LT.0) GO TO 300
            BEGCOL = IPTRU(I) + 2
            ENDCOL = IPTRL(I)
            IF (K.LE.NP) THEN
               DO 220 JJ = BEGCOL,ENDCOL
                  IF (IP(IRNF(JJ)).EQ.-K) GO TO 230
  220          CONTINUE
               GO TO 300
            END IF
  230       DO 280 JDUMMY = BEGCOL,ENDCOL
               JJ = BEGCOL
               DO 240 BEGCOL = JJ,ENDCOL
                  IF (IP(IRNF(BEGCOL)).GT.0) GO TO 250
  240          CONTINUE
               GO TO 290
  250          JJ = ENDCOL
               DO 260 ENDCOL = JJ,BEGCOL,-1
                  IF (IP(IRNF(ENDCOL)).LT.0) GO TO 270
  260          CONTINUE
               GO TO 290
  270          ASW = FACT(BEGCOL)
               FACT(BEGCOL) = FACT(ENDCOL)
               FACT(ENDCOL) = ASW
               J = IRNF(BEGCOL)
               IRNF(BEGCOL) = IRNF(ENDCOL)
               IRNF(ENDCOL) = J
               BEGCOL = BEGCOL + 1
               ENDCOL = ENDCOL - 1
  280       CONTINUE
  290       IW(M+N+I) = -ENDCOL
  300    CONTINUE
  310 CONTINUE
      IF (N.EQ.NP) THEN
         MF = M - MORD
         IF1 = LFACT + 1 - MF
         II = 0
         DO 320 I = 1,M
            IF (IP(I).GT.0) THEN
               IW(I+N) = N
               IRNF(IF1+II) = I
               II = II + 1
               IP(I) = NP + II
            END IF
  320    CONTINUE
         IF1 = LFACT + 1 - MAX(MF*NF,MF+NF)
         IF2 = IF1 - 1 + MF*MAX(0,JLAST-NP)
      END IF
      IF (INFO(5).EQ.MIN(M,N)) THEN
         DO 330 I = 1,M
            IP(I) = ABS(IP(I))
  330    CONTINUE
      ELSE
         MORD = NP
         DO 340 I = 1,M
            IF (IP(I).LT.0) THEN
               IP(I) = -IP(I)
            ELSE
               MORD = MORD + 1
               IP(I) = MORD
            END IF
  340    CONTINUE
      END IF
      IRNF(1) = INFO(6)
      IRNF(2) = MF
      INFO(7) = MF
      FACT(1) = CNTL(3)
      FACT(2) = CNTL(4)
      IF (INFO(1).EQ.-3) GO TO 520
      IF2 = IF2 - MF*NF
      DO 350 II = 1,MF*NF
         FACT(IL2+II) = FACT(IF1-1+II)
  350 CONTINUE
      DO 360 II = 1,MF
         IRNF(IL2+II) = IRNF(LFACT-MF+II)
  360 CONTINUE
      IF1 = IL2 + 1
      GO TO 440
  370 MF = IRNF(2)
      IF1 = IPTRL(N) + 1
      IF2 = IF1 - 1
      DO 430 K = JLAST + 1,N
         J = K
         IF (IQ(1).GT.0) J = IQ(K)
         IA1 = IPTRA(J)
         IA2 = NE
         IF (J.NE.N) IA2 = IPTRA(J+1) - 1
         IU1 = IL2 + 1
         IU2 = IPTRU(K)
         IL1 = IU2 + 1
         IL2 = IPTRL(K)
         DO 380 II = IA1,IA2
            W(IRNA(II)) = AA(II)
  380    CONTINUE
         DO 400 II = IU2,IU1,-1
            J = IRNF(II)
            I = IRNF(IPTRU(J)+1)
            EYE1 = IPTRU(J) + 1
            AMULT = -W(I)*FACT(EYE1)
            FACT(II) = AMULT
            W(I) = ZERO
            DO 390 EYE = EYE1 + 1,IPTRL(J)
               I = IRNF(EYE)
               W(I) = W(I) + AMULT*FACT(EYE)
  390       CONTINUE
            RINFO = RINFO + ONE + 2*(IPTRL(J)-EYE1)
  400    CONTINUE
         IF (K.LE.NP) THEN
            IF (IL1.LE.IL2) THEN
CDIR$ IVDEP
               DO 410 II = IL1,IL2
                  I = IRNF(II)
                  FACT(II) = W(I)
                  W(I) = ZERO
  410          CONTINUE
               IF (ABS(FACT(IL1)).LE.CNTL(4)) THEN
                  GO TO 530
               ELSE
                  INFO(5) = INFO(5) + 1
                  FACT(IL1) = ONE/FACT(IL1)
                  RINFO = RINFO + ONE
               END IF
            END IF
         ELSE
            DO 420 II = IF1,IF1 + MF - 1
               I = IRNF(II)
               IF2 = IF2 + 1
               FACT(IF2) = W(I)
               W(I) = ZERO
  420       CONTINUE
         END IF
  430 CONTINUE
      INFO(4) = MAX(IF1+MF+NF-1,IF2)
  440 IF (MF.GT.0 .AND. NF.GT.0) THEN
         IF (ICNTL(5).GT.1) CALL MA50GD(MF,NF,FACT(IF1),MF,ICNTL(5),
     +                                  CNTL(4),IRNF(IF1+MF),RANK)
         IF (ICNTL(5).EQ.1) CALL MA50FD(MF,NF,FACT(IF1),MF,CNTL(4),
     +                                  IRNF(IF1+MF),RANK)
         IF (ICNTL(5).LE.0) CALL MA50ED(MF,NF,FACT(IF1),MF,CNTL(4),
     +                                  IRNF(IF1+MF),RANK)
         INFO(5) = INFO(5) + RANK
         DO 450 I = 1,MIN(MF,NF)
            RINFO = RINFO + MF - I + 1 + REAL(MF-I)*(NF-I)*2
  450    CONTINUE
      END IF
      IF (INFO(5).LT.MIN(M,N)) INFO(1) = 1
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) THEN
         WRITE (MP,'(A,I6,A,F12.1/A,I3,A,4I8)')
     +     ' Leaving MA50BD with IRNF(2) =',IRNF(2),' RINFO =',RINFO,
     +     ' INFO(1) =',INFO(1),' INFO(4:7) =', (INFO(J),J=4,7)
         IF (ICNTL(3).GT.3) THEN
            IF (JOB.NE.2) WRITE (MP,'(A,(T6,10(I7)))') ' IP = ',IP
            DO 460 J = 1,N
               IF (J.GT.1) THEN
                  IF (IPTRL(J-1).LT.IPTRU(J)) WRITE (MP,
     +                '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,
     +                ' of U', (FACT(II),IRNF(II),II=IPTRL(J-1)+1,
     +                IPTRU(J))
               END IF
               IF (IPTRU(J).LT.IPTRL(J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L',
     +              (FACT(II),IRNF(II),II=IPTRU(J)+1,IPTRL(J))
  460       CONTINUE
            WRITE (MP,'(A)') ' Full part'
            WRITE (MP,'((6I12))') (IRNF(IF1+MF+J),J=0,NF-1)
            DO 470 I = 0,MF - 1
               WRITE (MP,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') IRNF(IF1+I),
     +            (FACT(IF1+I+J*MF),J=0,NF-1)
  470       CONTINUE
         END IF
      END IF
      GO TO 550
  480 DO 490 I = 1,M
         IW(I) = 0
  490 CONTINUE
      DO 500 I = 1,M
         IF (IP(I).GT.0) THEN
            IW(IP(I)) = I
         ELSE
            IP(I) = -IP(I)
         END IF
  500 CONTINUE
      DO 510 I = 1,M
         IF (IW(I).GT.0) THEN
            IP(IW(I)) = K
            K = K + 1
         END IF
  510 CONTINUE
  520 INFO(1) = -3
      IF (LP.GT.0) THEN
         WRITE (LP,'(/A/A,I7,A,I7)')
     +     ' **** Error return from MA50BD **** ',
     +     ' LFACT must be increased from',LFACT,' to at least',INFO(4)
      END IF
      GO TO 550
  530 INFO(1) = - (7+K)
      IF (LP.GT.0) WRITE (LP,'(/A/A,I6,A)')
     +    ' **** Error return from MA50BD **** ',
     +    ' Small pivot found in column',K,' of the permuted matrix.'
      GO TO 550
  540 INFO(1) = -4
      IF (LP.GT.0) WRITE (LP,'(/A/(3(A,I9)))')
     +    ' **** Error return from MA50BD ****',' Entry in row',I,
     +    ' and column',J,' duplicated'
  550 END
      SUBROUTINE MA50CD(M,N,ICNTL,IQ,NP,TRANS,LFACT,FACT,IRNF,IPTRL,
     +                  IPTRU,B,X,W,INFO)
      INTEGER M,N,ICNTL(7),IQ(*),NP
      LOGICAL TRANS
      INTEGER LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER IRNF(LFACT),IPTRL(N),IPTRU(N)
      DOUBLE PRECISION B(*),X(*),W(*)
      INTEGER INFO(7)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0D0)
      INTEGER I,II,IA1,IF1,J,LP,MF,MP,NF
      DOUBLE PRECISION PROD
      EXTERNAL MA50HD
      LP = ICNTL(1)
      MP = ICNTL(2)
      IF (ICNTL(3).LE.0) LP = 0
      IF (ICNTL(3).LE.1) MP = 0
      IF (M.LT.1 .OR. N.LT.1) GO TO 250
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) WRITE (MP,
     +    '(/2(A,I6),A,I4,A,L2)') ' Entering MA50CD with M=',M,' N =',N,
     +    ' NP =',NP,' TRANS =',TRANS
      IF1 = IPTRL(N) + 1
      MF = IRNF(2)
      NF = N - NP
      IF (MP.GT.0 .AND. ICNTL(3).GT.2) WRITE (MP,
     +    '(A,I5,A,I5)') ' Size of full submatrix',MF,' by',NF
      IF (MP.GT.0 .AND. ICNTL(3).GT.3) THEN
         DO 10 J = 1,N
            IF (J.GT.1) THEN
               IF (IPTRL(J-1).LT.IPTRU(J)) WRITE (MP,
     +             '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of U',
     +              (FACT(II),IRNF(II),II=IPTRL(J-1)+1,IPTRU(J))
            END IF
            IF (IPTRU(J).LT.IPTRL(J)) WRITE (MP,
     +          '(A,I5,A,(T18,3(1P,E12.4,I5)))') ' Column',J,' of L',
     +          (FACT(II),IRNF(II),II=IPTRU(J)+1,IPTRL(J))
   10    CONTINUE
         WRITE (MP,'(A)') ' Full part'
         WRITE (MP,'((6I12))') (IRNF(IF1+MF+J),J=0,NF-1)
         DO 20 I = 0,MF - 1
            WRITE (MP,'(I4,1P,6E12.4:/(4X,1P,6E12.4))') IRNF(IF1+I),
     +        (FACT(IF1+I+J*MF),J=0,NF-1)
   20    CONTINUE
      END IF
      IF (TRANS) THEN
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A4,5F10.4:/(4X,5F10.4))') ' B =', (B(I),I=1,N)
         IF (IQ(1).GT.0) THEN
            DO 30 I = 1,N
               W(I) = B(IQ(I))
   30       CONTINUE
         ELSE
            DO 40 I = 1,N
               W(I) = B(I)
   40       CONTINUE
         END IF
         DO 50 I = 1,M
            X(I) = ZERO
   50    CONTINUE
         DO 70 I = 2,N
            PROD = ZERO
            DO 60 II = IPTRL(I-1) + 1,IPTRU(I)
               PROD = PROD + FACT(II)*W(IRNF(II))
   60       CONTINUE
            W(I) = W(I) + PROD
   70    CONTINUE
         DO 80 I = 1,NF
            X(I) = W(NP+I)
   80    CONTINUE
         IF (MF.GT.0 .AND. NF.GT.0) THEN
            CALL MA50HD(TRANS,MF,NF,FACT(IF1),MF,IRNF(IF1+MF),X,
     +                  ICNTL(5))
         ELSE
            DO 90 I = 1,MF
               X(I) = ZERO
   90       CONTINUE
         END IF
         DO 100 I = MF,1,-1
            J = IRNF(IF1+I-1)
            IF (J.NE.I) X(J) = X(I)
  100    CONTINUE
         DO 120 I = NP,1,-1
            IA1 = IPTRU(I) + 1
            IF (IA1.GT.IPTRL(I)) GO TO 120
            PROD = ZERO
            DO 110 II = IA1 + 1,IPTRL(I)
               PROD = PROD + FACT(II)*X(IRNF(II))
  110       CONTINUE
            X(IRNF(IA1)) = (W(I)-PROD)*FACT(IA1)
  120    CONTINUE
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A/(4X,5F10.4))') ' Leaving MA50CD with X =', (X(I),I=1,M)
      ELSE
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A4,5F10.4:/(4X,5F10.4))') ' B =', (B(I),I=1,M)
         DO 130 I = 1,M
            W(I) = B(I)
  130    CONTINUE
         DO 150 I = 1,NP
            IA1 = IPTRU(I) + 1
            IF (IA1.LE.IPTRL(I)) THEN
               X(I) = W(IRNF(IA1))*FACT(IA1)
               IF (X(I).NE.ZERO) THEN
CDIR$ IVDEP
                  DO 140 II = IA1 + 1,IPTRL(I)
                     W(IRNF(II)) = W(IRNF(II)) - FACT(II)*X(I)
  140             CONTINUE
               END IF
            END IF
  150    CONTINUE
         IF (MF.GT.0 .AND. NF.GT.0) THEN
            DO 160 I = 1,MF
               W(I) = W(IRNF(IF1+I-1))
  160       CONTINUE
            CALL MA50HD(TRANS,MF,NF,FACT(IF1),MF,IRNF(IF1+MF),W,
     +                  ICNTL(5))
            DO 170 I = 1,NF
               X(NP+I) = W(I)
  170       CONTINUE
         ELSE
            DO 180 I = 1,NF
               X(NP+I) = ZERO
  180       CONTINUE
         END IF
         DO 200 J = N,MAX(2,NP+1),-1
            PROD = X(J)
CDIR$ IVDEP
            DO 190 II = IPTRL(J-1) + 1,IPTRU(J)
               X(IRNF(II)) = X(IRNF(II)) + FACT(II)*PROD
  190       CONTINUE
  200    CONTINUE
         DO 220 J = NP,2,-1
            IA1 = IPTRU(J)
            IF (IA1.GE.IPTRL(J)) THEN
               X(J) = ZERO
            ELSE
               PROD = X(J)
CDIR$ IVDEP
               DO 210 II = IPTRL(J-1) + 1,IA1
                  X(IRNF(II)) = X(IRNF(II)) + FACT(II)*PROD
  210          CONTINUE
            END IF
  220    CONTINUE
         IF (NP.GE.1 .AND. IPTRU(1).GE.IPTRL(1)) X(1) = ZERO
         IF (IQ(1).GT.0) THEN
            DO 230 I = 1,N
               W(I) = X(I)
  230       CONTINUE
            DO 240 I = 1,N
               X(IQ(I)) = W(I)
  240       CONTINUE
         END IF
         IF (MP.GT.0 .AND. ICNTL(3).GT.3) WRITE (MP,
     +       '(A/(4X,5F10.4))') ' Leaving MA50CD with X =', (X(I),I=1,N)
      END IF
      RETURN
  250 INFO(1) = -1
      IF (LP.GT.0) WRITE (LP,'(/A/2(A,I8))')
     +    ' **** Error return from MA50CD ****',' M =',M,' N =',N
      END
      SUBROUTINE MA50DD(LA,A,IND,IPTR,N,DISP,REALS)
      INTEGER LA,N,DISP
      DOUBLE PRECISION A(LA)
      INTEGER IPTR(N)
      LOGICAL REALS
      INTEGER IND(LA)
      INTEGER J,K,KN
      DO 10 J = 1,N
         K = IPTR(J)
         IF (K.GT.0) THEN
            IPTR(J) = IND(K)
            IND(K) = -J
         END IF
   10 CONTINUE
      KN = 0
      DO 20 K = 1,DISP - 1
         IF (IND(K).EQ.0) GO TO 20
         KN = KN + 1
         IF (REALS) A(KN) = A(K)
         IF (IND(K).LE.0) THEN
            J = -IND(K)
            IND(K) = IPTR(J)
            IPTR(J) = KN
         END IF
         IND(KN) = IND(K)
   20 CONTINUE
      DISP = KN + 1
      END
      SUBROUTINE MA50ED(M,N,A,LDA,PIVTOL,IPIV,RANK)
**
      INTEGER LDA,M,N,RANK
      DOUBLE PRECISION PIVTOL
      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N)
**
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER I,J,JP,K
      LOGICAL PIVOT
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      EXTERNAL DAXPY,DSCAL,DSWAP
      INTRINSIC ABS
      J = 1
      DO 30 K = 1,N
         DO 10 I = 1,J - 1
            IF (M.GT.I) CALL DAXPY(M-I,-A(I,J),A(I+1,I),1,A(I+1,J),1)
   10    CONTINUE
         IF (J.LE.M) THEN
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .GT. PIVTOL
         ELSE
            PIVOT = .FALSE.
         END IF
         IF (PIVOT) THEN
            IF (JP.NE.J) CALL DSWAP(N+J-K,A(J,1),LDA,A(JP,1),LDA)
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)
            J = J + 1
         ELSE
            DO 20 I = J,M
               A(I,J) = ZERO
   20       CONTINUE
            IF (K.LT.N) CALL DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IPIV(N-K+J) = -J
         END IF
   30 CONTINUE
      RANK = J - 1
      END
      SUBROUTINE MA50FD(M,N,A,LDA,PIVTOL,IPIV,RANK)
      INTEGER LDA,M,N,RANK
      DOUBLE PRECISION PIVTOL
      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N)
**
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER I,J,JP,K
      LOGICAL PIVOT
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      EXTERNAL DGEMV,DSCAL,DSWAP
      INTRINSIC ABS
      J = 1
      DO 20 K = 1,N
         IF (J.LE.M) THEN
            CALL DGEMV('No transpose',M-J+1,J-1,-ONE,A(J,1),LDA,A(1,J),
     +                 1,ONE,A(J,J),1)
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .GT. PIVTOL
         ELSE
            PIVOT = .FALSE.
         END IF
         IF (PIVOT) THEN
            IF (JP.NE.J) CALL DSWAP(N+J-K,A(J,1),LDA,A(JP,1),LDA)
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)
            IF (J.LT.N) THEN
               CALL DGEMV('Transpose',J-1,N-J,-ONE,A(1,J+1),LDA,A(J,1),
     +                    LDA,ONE,A(J,J+1),LDA)
            END IF
            J = J + 1
         ELSE
            DO 10 I = J,M
               A(I,J) = ZERO
   10       CONTINUE
            IF (K.LT.N) CALL DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IPIV(N-K+J) = -J
         END IF
   20 CONTINUE
      RANK = J - 1
      END
      SUBROUTINE MA50GD(M,N,A,LDA,NB,PIVTOL,IPIV,RANK)
      INTEGER LDA,M,N,NB,RANK
      DOUBLE PRECISION PIVTOL
      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N)
**
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER I,J,JJ,JP,J1,J2,K
      LOGICAL PIVOT
      DOUBLE PRECISION TEMP
      EXTERNAL DGEMM,DGEMV,DSWAP,DSCAL,DTRSM,DTRSV
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      INTRINSIC ABS,MIN
      J = 1
      J1 = 1
      J2 = MIN(N,NB)
      DO 70 K = 1,N
         IF (J.LE.M) THEN
            CALL DGEMV('No transpose',M-J+1,J-J1,-ONE,A(J,J1),LDA,
     +                 A(J1,J),1,ONE,A(J,J),1)
            JP = J - 1 + IDAMAX(M-J+1,A(J,J),1)
            IPIV(J) = JP
            PIVOT = ABS(A(JP,J)) .GT. PIVTOL
         ELSE
            PIVOT = .FALSE.
         END IF
         IF (PIVOT) THEN
            IF (JP.NE.J) CALL DSWAP(J2-J1+1,A(J,J1),LDA,A(JP,J1),LDA)
            IF (J.LT.M) CALL DSCAL(M-J,ONE/A(J,J),A(J+1,J),1)
            IF (J+1.LE.J2) THEN
               CALL DGEMV('Transpose',J-J1,J2-J,-ONE,A(J1,J+1),LDA,
     +                    A(J,J1),LDA,ONE,A(J,J+1),LDA)
            END IF
            J = J + 1
         ELSE
            DO 10 I = J,M
               A(I,J) = ZERO
   10       CONTINUE
            IPIV(N-K+J) = -J
            IF (K.NE.N) CALL DSWAP(M,A(1,J),1,A(1,N-K+J),1)
            IF (N-K+J.GT.J2) THEN
               DO 20 I = J1,J - 1
                  JP = IPIV(I)
                  TEMP = A(I,J)
                  A(I,J) = A(JP,J)
                  A(JP,J) = TEMP
   20          CONTINUE
               IF(J.GT.J1) CALL DTRSV('Lower','No transpose','Unit',
     +                                J-J1,A(J1,J1),LDA,A(J1,J),1)
            ELSE
               J2 = J2 - 1
            END IF
         END IF
         IF (J.GT.J2) THEN
            DO 40 JJ = 1,J1 - 1
               DO 30 I = J1,J2
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   30          CONTINUE
   40       CONTINUE
            DO 60 JJ = J2 + 1,N - K + J - 1
               DO 50 I = J1,J2
                  JP = IPIV(I)
                  TEMP = A(I,JJ)
                  A(I,JJ) = A(JP,JJ)
                  A(JP,JJ) = TEMP
   50          CONTINUE
   60       CONTINUE
            IF (K.NE.N) THEN
               CALL DTRSM('Left','Lower','No transpose','Unit',J2-J1+1,
     +                    N-K,ONE,A(J1,J1),LDA,A(J1,J2+1),LDA)
               IF (M.GT.J2) CALL DGEMM('No transpose','No transpose',
     +                                 M-J2,N-K,J2-J1+1,-ONE,A(J2+1,J1),
     +                                 LDA,A(J1,J2+1),LDA,ONE,
     +                                 A(J2+1,J2+1),LDA)
            END IF
            J1 = J2 + 1
            J2 = MIN(J2+NB,N-K+J-1)
         END IF
   70 CONTINUE
      RANK = J - 1
      END
      SUBROUTINE MA50HD(TRANS,M,N,A,LDA,IPIV,B,ICNTL5)
      LOGICAL TRANS
      INTEGER LDA,M,N,ICNTL5
      INTEGER IPIV(N)
      DOUBLE PRECISION A(LDA,N),B(*)
      INTEGER I,K,RANK
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0D0)
      DOUBLE PRECISION TEMP
      INTRINSIC MIN
      EXTERNAL DAXPY,DDOT,DTRSV
      DOUBLE PRECISION DDOT
      RANK = 0
      DO 10 RANK = MIN(M,N),1,-1
         IF (IPIV(RANK).GT.0) GO TO 20
   10 CONTINUE
   20 IF (.NOT.TRANS) THEN
         DO 30 I = 1,RANK
            K = IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   30    CONTINUE
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('L','NoTrans','Unit',RANK,A,LDA,B,
     +                                1)
         ELSE
            DO 40 K = 1,RANK - 1
               IF (B(K).NE.ZERO) CALL DAXPY(RANK-K,-B(K),A(K+1,K),1,
     +                                B(K+1),1)
   40       CONTINUE
         END IF
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('U','NoTrans','NonUnit',RANK,A,
     +                                LDA,B,1)
         ELSE
            DO 50 K = RANK,2,-1
               IF (B(K).NE.ZERO) THEN
                  B(K) = B(K)/A(K,K)
                  CALL DAXPY(K-1,-B(K),A(1,K),1,B(1),1)
               END IF
   50       CONTINUE
            IF (RANK.GT.0) B(1) = B(1)/A(1,1)
         END IF
         DO 60 K = RANK + 1,N
            B(K) = ZERO
   60    CONTINUE
         DO 70 I = RANK + 1,N
            K = -IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   70    CONTINUE
      ELSE
         DO 80 I = N,RANK + 1,-1
            K = -IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
   80    CONTINUE
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('U','Trans','NonUnit',RANK,A,LDA,
     +                                B,1)
         ELSE
            IF (RANK.GT.0) B(1) = B(1)/A(1,1)
            DO 90 I = 2,RANK
               TEMP = B(I) - DDOT(I-1,A(1,I),1,B(1),1)
               B(I) = TEMP/A(I,I)
   90       CONTINUE
         END IF
         IF (ICNTL5.GT.0) THEN
            IF (RANK.GT.0) CALL DTRSV('L','Trans','Unit',RANK,A,LDA,B,1)
         ELSE
            DO 100 I = RANK - 1,1,-1
               B(I) = B(I) - DDOT(RANK-I,A(I+1,I),1,B(I+1),1)
  100       CONTINUE
         END IF
         DO 110 I = RANK + 1,M
            B(I) = ZERO
  110    CONTINUE
         DO 120 I = RANK,1,-1
            K = IPIV(I)
            TEMP = B(I)
            B(I) = B(K)
            B(K) = TEMP
  120    CONTINUE
      END IF
      END
      SUBROUTINE MA50ID(CNTL,ICNTL)
      DOUBLE PRECISION CNTL(4)
      INTEGER ICNTL(7)
      CNTL(1) = 0.5D0
      CNTL(2) = 0.1D0
      CNTL(3) = 0.0D0
      CNTL(4) = 0.0D0
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 1
      ICNTL(4) = 3
      ICNTL(5) = 32
      ICNTL(6) = 0
      ICNTL(7) = 0
      END


