C COPYRIGHT (c) 1995 Timothy A. Davis, Patrick Amestoy and Council for
C                    the Central Laboratory of the Research Councils
C######DATE 30 November 1995
      SUBROUTINE MC47AD(N, NE, PE, IW, IWLEN, MP, INFO)
      INTEGER N, NE, PE(N+1), IWLEN, IW(IWLEN), MP, INFO(7)
      INTEGER DEGREE
      DOUBLE PRECISION DUMMY(1)
      INTEGER ELEN,HEAD,I,IFLAG,II,I1,I2,J,LAST,LEN,LENIW,NCMPA,
     *        NEXT,NV,PFREE,W
      EXTERNAL MC49AD,MC34AD,MC47BD
      INTEGER LP,MP49,IOUT,JOUT,IDUP,NZOUT
      COMMON /MC49ED/ LP, MP49, IOUT, JOUT, IDUP, NZOUT
      SAVE /MC49ED/
      INFO(1) = 0
      IF (N.LT.1) THEN
        INFO(1) = -1
        GO TO 1000
      ENDIF
      IF (PE(1).LT.1) THEN
        IF (2*NE+N.GT.IWLEN) THEN
          INFO(1) = -2
          GO TO 1000
        ENDIF
      ELSE
        IF (NE+N.GT.IWLEN) THEN
          INFO(1) = -2
          GO TO 1000
        ENDIF
      ENDIF
      IF (MP.GT.0) THEN
        WRITE(MP,'(/A)') 'Entry to MC47A/AD'
        WRITE(MP,'(A,I10,A,I10,A)') 'Matrix of order',N,' with',NE,
     *                            ' entries'
        IF (PE(1).LT.0)  THEN
          WRITE(MP,'(A)') 'Matrix input in coordinate form'
          WRITE(MP,'(A/(4(I8,I8,4X)))') 'Row and column indices',
     *          (IW(I),IW(NE+I),I=1,NE)
        ELSE
          WRITE(MP,'(A)') 'Matrix input by columns'
          DO 10 J=1,N
            WRITE(MP,'(A,I4/(10I8))') 'Column',J,
     *                                (IW(I),I=PE(J),PE(J+1)-1)
   10     CONTINUE
        ENDIF
      ENDIF
      LAST   = IWLEN  - N + 1
      ELEN   = LAST   - N
      NV     = ELEN   - N
      W      = NV     - N
      DEGREE = W      - N
      HEAD   = DEGREE - N
      NEXT   = HEAD   - N
      LEN    = NEXT   - N
      LENIW = LEN-1
      INFO(6) = 0
      INFO(7) = 0
      IF (PE(1).LT.0) THEN
        DO 20 I=1,NE
          IF (IW(I).LE.IW(NE+I)) THEN
            IF (IW(I).EQ.IW(NE+I) .AND. IW(I).NE.0) THEN
              INFO(7) = INFO(7) + 1
            ELSE
              IF (IW(I).GT.0) INFO(6) = INFO(6) + 1
            ENDIF
            IW(I)=0
          ENDIF
   20   CONTINUE
        MP49 = 0
        CALL MC49AD(2,N,N,NE,IW,IW(NE+1),.FALSE.,1,DUMMY,N+1,PE,N+1,
     *              IW(2*NE+1),IFLAG)
      ELSE
        IDUP = 0
        IOUT = 0
        JOUT = 0
        DO 30 I = 1,N
          IW(NE+I) = 0
   30   CONTINUE
        DO 50 J=1,N
          I1 = PE(J)
          PE(J) = I1-(IOUT+IDUP)
          I2 = PE(J+1)-1
          IF (I2.LT.I1-1) THEN
            INFO(1) = -3
            GO TO 1000
          ENDIF
          DO 40 II = I1,I2
            I = IW(II)
            IF (I.LE.J .OR. I.GT.N) THEN
              IF (I.EQ.J) INFO(7) = INFO(7) + 1
              IF (I.GT.0 .AND. I.LT.J) INFO(6) = INFO(6) + 1
              IOUT = IOUT + 1
            ELSE
              IF (IW(NE+I).EQ.J) THEN
                IDUP = IDUP + 1
              ELSE
                IW(NE+I)=J
                IW(II-(IOUT+IDUP)) = I
              ENDIF
            ENDIF
   40     CONTINUE
   50   CONTINUE
        PE(N+1) = NE - (IOUT+IDUP) + 1
      ENDIF
      IF (IDUP.GT.0) THEN
        INFO(1) = 1
        INFO(4) = IDUP
      ELSE
        INFO(4) = 0
      ENDIF
      IF (IOUT.GT.0 .OR. JOUT.GT.0) THEN
        INFO(1) = 1
        INFO(5) = IOUT + JOUT - INFO(7)
      ELSE
        INFO(5) = 0
      ENDIF
      IF (INFO(6).GT.0 .OR. INFO(7).GT.0) INFO(1) = 1
      IF (NE-(IOUT+IDUP).EQ.0) THEN
        INFO(1) = -4
        GO TO 1000
      ENDIF
      IF (LENIW.LT.2*(PE(N+1)-1)) THEN
        INFO(1) = -2
        GO TO 1000
      ENDIF
      CALL MC34AD(N,IW,PE,.FALSE.,DUMMY,IW(W))
      PFREE = PE(N+1)
      DO 60 I=1,N
        IW(LEN+I-1) = PE(I+1) - PE(I)
   60 CONTINUE
      CALL MC47BD(N,LENIW,PE,PFREE,IW(LEN),IW,IW(NV),IW(ELEN),
     *            IW(LAST),NCMPA,IW(DEGREE),IW(HEAD),IW(NEXT),IW(W))
      INFO(2) = NCMPA
      INFO(3) = PFREE+8*N
      IF (MP.GT.0) THEN
        WRITE(MP,'(/A)') 'Exit from MC47A/AD'
        WRITE(MP,'(A/7I10)') 'INFO(1-7):',(INFO(I),I=1,7)
        WRITE(MP,'(A/(8I10))') 'Parent array',(PE(I),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Permutation',(IW(ELEN+I-1),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Inverse permutation',
     *                         (IW(LAST+I-1),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Degree array',(IW(NV+I-1),I=1,N)
      ENDIF
 1000 RETURN
      END
      SUBROUTINE MC47BD (N, IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     $                   LAST, NCMPA, DEGREE, HEAD, NEXT, W)
      INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N),
     $        ELEN(N), LAST(N), NCMPA, DEGREE(N), HEAD(N), NEXT(N),
     $        W(N)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, HASH, HMOD, I,
     $        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     $        LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     $        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      INTRINSIC MAX, MIN, MOD
C=======================================================================
C=======================================================================
      WFLG = 2
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = MAX (1, N-1)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      DO 10 I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
        DEGREE (I) = LEN (I)
   10 CONTINUE
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG .GT. 0) THEN
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          HEAD (DEG) = I
        ELSE
          NEL = NEL + 1
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
C=======================================================================
C=======================================================================
   30 IF (NEL .LT. N) THEN
C=======================================================================
C=======================================================================
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        INEXT = NEXT (ME)
        IF (INEXT .NE. 0) LAST (INEXT) = 0
        HEAD (DEG) = INEXT
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
C=======================================================================
C=======================================================================
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (DEGREE (I)) = INEXT
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          DO 120 KNT1 = 1, ELENME + 1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            DO 110 KNT2 = 1, LN
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = IW (PN)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = PE (J)
                      PE (J) = PDST
                      PDST = PDST + 1
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (DEGREE (I)) = INEXT
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = MAX (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = PME2 - PME1 + 1
        IF (WFLG+N .LE. WFLG) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C=======================================================================
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
C=======================================================================
C=======================================================================
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
            ELSE IF (DEXT .EQ. 0) THEN
              PE (E) = -ME
              W (E) = 0
            ENDIF
  160     CONTINUE
          ELEN (I) = PN - P1 + 1
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              DEG = DEG + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + J
            ENDIF
  170     CONTINUE
          IF (DEG .EQ. 0) THEN
            PE (I) = -ME
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE (I) = MIN (DEGREE (I), DEG)
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = PN - P1 + 1
            HASH = MOD (HASH, HMOD) + 1
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = HASH
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = MAX (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG+N .LE. WFLG) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C=======================================================================
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF (NV (I) .LT. 0) THEN
            HASH = LAST (I)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                PE (J) = -I
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
C=======================================================================
C=======================================================================
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            DEG = MIN (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
            MINDEG = MIN (MINDEG, DEG)
            DEGREE (I) = DEG
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
C=======================================================================
C=======================================================================
        NV (ME) = NVPIV + DEGME
        LEN (ME) = P - PME1
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
C=======================================================================
      GO TO 30
      ENDIF
C=======================================================================
C=======================================================================
C=======================================================================
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = -PE (I)
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = -PE (J)
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = -PE (J)
              PE (J) = -E
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = ABS (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
C=======================================================================
C=======================================================================
      PFREE = MAXMEM
      RETURN
      END



