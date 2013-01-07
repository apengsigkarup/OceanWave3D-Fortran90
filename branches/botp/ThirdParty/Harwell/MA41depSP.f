C *******************************************************************
C COPYRIGHT (c) 1995 Timothy A. Davis, Patrick Amestoy and Council for
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
C Licence, see http://hsl.rl.ac.uk/acuk/cou.html
C
C Please note that for a UK ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C                    the Central Laboratory of the Research Councils
C######DATE 30 November 1995
C  April 2001: call to MC49 changed to MC59 to make routine threadsafe
C 20/2/02 Cosmetic changes applied to reduce single/double differences

      SUBROUTINE MC47A(N, NE, PE, IW, IWLEN, MP, INFO)
      INTEGER N, NE, PE(N+1), IWLEN, IW(IWLEN), MP, INFO(7)
      INTEGER DEGREE
      REAL DUMMY(1)
      INTEGER ELEN,HEAD,I,IFLAG,II,I1,I2,J,LAST,LEN,LENIW,NCMPA,
     *        NEXT,NV,PFREE,W
      INTEGER ICT59(10),INFO59(10),IOUT,JOUT,IDUP,NZOUT
      EXTERNAL MC59A,MC34A,MC47B
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
        ICT59(1) = 0
        ICT59(2) = 1
        ICT59(3) = 1
        ICT59(4) = MP
        ICT59(5) = -1
        ICT59(6) = 0
        CALL MC59A(ICT59,N,N,NE,IW,NE,IW(NE+1),1,DUMMY,
     *             N+1,PE,N+1,IW(2*NE+1),INFO59)
        IFLAG = INFO59(1)
        IDUP  = INFO59(3)
        IOUT  = INFO59(4)
        JOUT  = INFO59(5)
        NZOUT = INFO59(6)
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
      CALL MC34A(N,IW,PE,.FALSE.,DUMMY,IW(W))
      PFREE = PE(N+1)
      DO 60 I=1,N
        IW(LEN+I-1) = PE(I+1) - PE(I)
   60 CONTINUE
      CALL MC47B(N,LENIW,PE,PFREE,IW(LEN),IW,IW(NV),IW(ELEN),
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
      SUBROUTINE MC47B (N, IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
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
* *******************************************************************
* COPYRIGHT (c) 1977 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
* Licence, see http://hsl.rl.ac.uk/acuk/cou.html
*
* Please note that for a UK ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 8 Oct 1992
C######8/10/92 Toolpack tool decs employed.
C######8/10/92 D version created by name change only.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
      SUBROUTINE MC21A(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
      INTEGER LICN,N,NUMNZ
      INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
      EXTERNAL MC21B
      CALL MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),IW(1,3),
     +           IW(1,4))
      RETURN
      END
      SUBROUTINE MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
      INTEGER LICN,N,NUMNZ
      INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
      DO 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      NUMNZ = 0
      DO 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = ICN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = -1
   30     CONTINUE
          OUT(J) = LENR(J) - 1
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = ICN(II)
              IF (CV(I).EQ.JORD) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 100
   60     CONTINUE
   70   CONTINUE
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   CONTINUE
  100 CONTINUE
      IF (NUMNZ.EQ.N) RETURN
      DO 110 I = 1,N
        ARP(I) = 0
  110 CONTINUE
      K = 0
      DO 130 I = 1,N
        IF (IPERM(I).NE.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130
  120   CONTINUE
        J = IPERM(I)
        ARP(J) = I
  130 CONTINUE
      K = 0
      DO 140 I = 1,N
        IF (ARP(I).NE.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 CONTINUE
      RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1993 Hyprotech UK and
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
* Licence, see http://hsl.rl.ac.uk/acuk/cou.html
*
* Please note that for a UK ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
* Council for the Central Laboratory of the Research Councils
      SUBROUTINE MC29A(M,N,NE,A,IRN,ICN,R,C,W,LP,IFAIL)
      INTEGER M,N,NE
      REAL A(NE)
      INTEGER IRN(NE),ICN(NE)
      REAL R(M),C(N),W(M*2+N*3)
      INTEGER LP,IFAIL
      INTRINSIC LOG,ABS,MIN
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      REAL ONE,SMIN,ZERO
      PARAMETER (ONE=1.0,SMIN=0.1,ZERO=0.0)
      INTEGER I,I1,I2,I3,I4,I5,ITER,J,K
      REAL E,E1,EM,Q,Q1,QM,S,S1,SM,U,V
      IFAIL = 0
      IF (M.LT.1 .OR. N.LT.1) THEN
         IFAIL = -1
         GO TO 220
      ELSE IF (NE.LE.0) THEN
         IFAIL = -2
         GO TO 220
      END IF
      I1 = 0
      I2 = M
      I3 = M + N
      I4 = M + N*2
      I5 = M + N*3
      DO 10 I = 1,M
         R(I) = ZERO
         W(I1+I) = ZERO
   10 CONTINUE
      DO 20 J = 1,N
         C(J) = ZERO
         W(I2+J) = ZERO
         W(I3+J) = ZERO
         W(I4+J) = ZERO
   20 CONTINUE
      DO 30 K = 1,NE
         U = ABS(A(K))
         IF (U.EQ.ZERO) GO TO 30
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 30
         U = LOG(U)
         W(I1+I) = W(I1+I) + ONE
         W(I2+J) = W(I2+J) + ONE
         R(I) = R(I) + U
         W(I3+J) = W(I3+J) + U
   30 CONTINUE
      DO 40 I = 1,M
         IF (W(I1+I).EQ.ZERO) W(I1+I) = ONE
         R(I) = R(I)/W(I1+I)
         W(I5+I) = R(I)
   40 CONTINUE
      DO 50 J = 1,N
         IF (W(I2+J).EQ.ZERO) W(I2+J) = ONE
         W(I3+J) = W(I3+J)/W(I2+J)
   50 CONTINUE
      SM = SMIN*NE
      DO 60 K = 1,NE
         IF (A(K).EQ.ZERO) GO TO 60
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 60
         R(I) = R(I) - W(I3+J)/W(I1+I)
   60 CONTINUE
      E = ZERO
      Q = ONE
      S = ZERO
      DO 70 I = 1,M
         S = S + W(I1+I)*R(I)**2
   70 CONTINUE
      IF (S.LE.SM) GO TO 160
      DO 150 ITER = 1,MAXIT
         DO 80 K = 1,NE
            IF (A(K).EQ.ZERO) GO TO 80
            J = ICN(K)
            I = IRN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 80
            C(J) = C(J) + R(I)
   80    CONTINUE
         S1 = S
         S = ZERO
         DO 90 J = 1,N
            V = -C(J)/Q
            C(J) = V/W(I2+J)
            S = S + V*C(J)
   90    CONTINUE
         E1 = E
         E = Q*S/S1
         Q = ONE - E
         IF (S.LE.SM) E = ZERO
         DO 100 I = 1,M
            R(I) = R(I)*E*W(I1+I)
  100    CONTINUE
         IF (S.LE.SM) GO TO 180
         EM = E*E1
         DO 110 K = 1,NE
            IF (A(K).EQ.ZERO) GO TO 110
            I = IRN(K)
            J = ICN(K)
            IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 110
            R(I) = R(I) + C(J)
  110    CONTINUE
         S1 = S
         S = ZERO
         DO 120 I = 1,M
            V = -R(I)/Q
            R(I) = V/W(I1+I)
            S = S + V*R(I)
  120    CONTINUE
         E1 = E
         E = Q*S/S1
         Q1 = Q
         Q = ONE - E
         IF (S.LE.SM) Q = ONE
         QM = Q*Q1
         DO 130 J = 1,N
            W(I4+J) = (EM*W(I4+J)+C(J))/QM
            W(I3+J) = W(I3+J) + W(I4+J)
  130    CONTINUE
         IF (S.LE.SM) GO TO 160
         DO 140 J = 1,N
            C(J) = C(J)*E*W(I2+J)
  140    CONTINUE
  150 CONTINUE
  160 DO 170 I = 1,M
         R(I) = R(I)*W(I1+I)
  170 CONTINUE
  180 DO 190 K = 1,NE
         IF (A(K).EQ.ZERO) GO TO 190
         I = IRN(K)
         J = ICN(K)
         IF (MIN(I,J).LT.1 .OR. I.GT.M .OR. J.GT.N) GO TO 190
         R(I) = R(I) + W(I3+J)
  190 CONTINUE
      DO 200 I = 1,M
         R(I) = R(I)/W(I1+I) - W(I5+I)
  200 CONTINUE
      DO 210 J = 1,N
         C(J) = -W(I3+J)
  210 CONTINUE
      RETURN
  220 IF (LP.GT.0) WRITE (LP,'(/A/A,I3)')
     +    ' **** Error return from MC29A ****',' IFAIL =',IFAIL
      END
* *******************************************************************
* COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
* Licence, see http://hsl.rl.ac.uk/acuk/cou.html
*
* Please note that for a UK ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
C 29 January 2001. Modified from MC49 to be threadsafe.

      SUBROUTINE MC59A(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP,
     &                  LIW,IW,INFO)
      INTEGER LA,LIP,LIW,LJCN,NC,NE,NR
      REAL A(LA)
      INTEGER ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
      INTEGER I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
      INTEGER IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
      LOGICAL LCHECK
      EXTERNAL MC59B,MC59C,MC59D,MC59E,MC59F
      INTRINSIC MAX
      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE
      ICNTL1 = ICNTL(1)
      ICNTL2 = ICNTL(2)
      ICNTL3 = ICNTL(3)
      ICNTL6 = ICNTL(6)
      LCHECK = (ICNTL1.EQ.0)
      LP = ICNTL(4)
      MP = ICNTL(5)
      IF (ICNTL2.GT.2 .OR. ICNTL2.LT.0) THEN
         INFO(1) = -1
         INFO(2) = ICNTL2
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9010) ICNTL2
         END IF
         GO TO 70
      END IF
      IF (ICNTL6.GT.2 .OR. ICNTL6.LT.-2) THEN
         INFO(1) = -11
         INFO(2) = ICNTL6
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9150) ICNTL6
         END IF
         GO TO 70
      END IF
      IF (NC.LT.1) THEN
        INFO(1) = -2
        INFO(2) = NC
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9020) NC
        END IF
        GO TO 70
      END IF
      IF (NR.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9030) NR
        END IF
        GO TO 70
      END IF
      IF (ICNTL6.NE.0 .AND. NR.NE.NC) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9035) NC,NR
        END IF
        GO TO 70
      END IF
      IF (NE.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NE
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9040) NE
        END IF
        GO TO 70
      END IF
      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.1) THEN
        IF (LJCN.LT.NE) THEN
          INFO(1) = -5
          INFO(2) = NE
        END IF
      ELSE
        IF (LJCN.LT.1) THEN
          INFO(1) = -5
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-5) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9050) LJCN,INFO(2)
         END IF
         GO TO 70
      END IF
      IF (ICNTL3.EQ.0) THEN
        IF (LA.LT.NE) THEN
          INFO(1) = -6
          INFO(2) = NE
        END IF
      ELSE
        IF (LA.LT.1) THEN
          INFO(1) = -6
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-6) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9060) LA,INFO(2)
         END IF
         GO TO 70
      END IF
      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.2) THEN
        IF (LIP.LT.NC+1) THEN
          INFO(1) = -7
          INFO(2) = NC+1
        END IF
      ELSE IF (LIP.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -7
        INFO(2) = MAX(NR,NC)+1
      END IF
      IF (INFO(1).EQ.-7) THEN
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9065) LIP,INFO(2)
        END IF
        GO TO 70
      END IF
      IF (LIW.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -8
        INFO(2) = MAX(NR,NC)+1
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9070) LIW,INFO(2)
        END IF
        GO TO 70
      END IF
      LAA = NE
      IF (ICNTL3.NE.0) LAA = 1
      IOUT = 0
      JOUT = 0
      IDUP = 0
      IUP = 0
      PART = 0
      IF (ICNTL6.NE.0) PART = 1
      IF (ICNTL2.EQ.0) THEN
        CALL MC59B(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW,
     +              IOUT,JOUT,KNE)
        IF (KNE.EQ.0) GO TO 50
        IF (LCHECK) CALL MC59E(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP,
     &                          KNE,ICNTL6)
      ELSE IF (ICNTL2.EQ.1) THEN
        IF (ICNTL6.NE.0) PART = -1
        CALL MC59B(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP,
     +              JOUT,IOUT,KNE)
        IF (KNE.EQ.0) GO TO 50
        IF (LCHECK) CALL MC59E(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP,
     &                          IDUP,KNE,ICNTL6)
        CALL MC59C(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)
      ELSE IF (ICNTL2.EQ.2) THEN
        IF (LCHECK) THEN
          CALL MC59F(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP,
     +                IOUT,IUP,KNE,ICNTL6,INFO)
          IF (INFO(1).EQ.-9) GO TO 40
          IF (KNE.EQ.0) GO TO 50
        ELSE
           KNE = NE
        END IF
        CALL MC59D(NC,KNE,IRN,IP,LAA,A)
      END IF
      INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(6) = KNE
      INFO(7) = IUP
      IF (IDUP.GT.0) INFO(1) = INFO(1) + 1
      IF (IOUT.GT.0) INFO(1) = INFO(1) + 2
      IF (JOUT.GT.0) INFO(1) = INFO(1) + 4
      IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
        WRITE (MP,FMT=9080) INFO(1)
        IF (IOUT.GT.0) WRITE (MP,FMT=9090) IOUT
        IF (JOUT.GT.0) WRITE (MP,FMT=9110) JOUT
        IF (IDUP.GT.0) WRITE (MP,FMT=9100) IDUP
        IF (IUP.GT.0)  WRITE (MP,FMT=9130) IUP
      END IF
      GO TO 70
   40 INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(7) = IUP
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9140)
      END IF
      GO TO 70
   50 INFO(1) = -10
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(2) = IOUT + JOUT
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9120)
      END IF
   70 RETURN
 9000 FORMAT (/,' *** Error return from MC59A *** INFO(1) = ',I3)
 9010 FORMAT (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 FORMAT (1X,'NC = ',I6,' is out of range')
 9030 FORMAT (1X,'NR = ',I6,' is out of range')
 9035 FORMAT (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 FORMAT (1X,'NE = ',I10,' is out of range')
 9050 FORMAT (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 FORMAT (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 FORMAT (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 FORMAT (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 FORMAT (/,' *** Warning message from MC59A *** INFO(1) = ',I3)
 9090 FORMAT (1X,I8,' entries in IRN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9100 FORMAT (1X,I8,' duplicate entries were supplied by the user')
 9110 FORMAT (1X,I8,' entries in JCN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9120 FORMAT (1X,'All entries out of range')
 9130 FORMAT (1X,I8,' of these entries were in the upper triangular ',
     +       /,'       part of matrix')
 9140 FORMAT (1X,'Entries in IP are not monotonic increasing')
 9150 FORMAT (1X,'ICNTL(6) = ',I2,' is out of range')
      END
C***********************************************************************
      SUBROUTINE MC59B(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT,
     +                  JOUT,KNE)
      INTEGER LA,NC,NE,NR,IOUT,JOUT,KNE,PART
      LOGICAL LCHECK
      REAL A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
      REAL ACE,ACEP
      INTEGER I,ICE,ICEP,J,JCE,JCEP,K,L,LOC
      DO 10 J = 1,NC + 1
        IW(J) = 0
   10 CONTINUE
      KNE = 0
      IOUT = 0
      JOUT = 0
      IF (LCHECK) THEN
        IF (LA.GT.1) THEN
          IF (PART.EQ.0) THEN
            DO 20 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                A(KNE) = A(K)
                IW(J) = IW(J) + 1
              END IF
   20       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 21 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   21       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 22 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   22       CONTINUE
          END IF
        ELSE
          IF (PART.EQ.0) THEN
            DO 25 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                IW(J) = IW(J) + 1
              END IF
   25       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 26 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   26       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 27 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   27       CONTINUE
          END IF
        END IF
        IF (KNE.EQ.0) GO TO 130
      ELSE
        KNE = NE
        IF (PART.EQ.0) THEN
          DO 30 K = 1,NE
            J = JCN(K)
            IW(J) = IW(J) + 1
   30     CONTINUE
        ELSE IF (PART.EQ.1) THEN
          DO 35 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            IF (I.LT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   35     CONTINUE
        ELSE IF (PART.EQ.-1) THEN
          DO 36 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            IF (I.GT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   36     CONTINUE
        END IF
      END IF
      IP(1) = 1
      DO 37 J = 2,NC + 1
        IP(J) = IW(J-1) + IP(J-1)
        IW(J-1) = IP(J-1)
   37 CONTINUE
      IF (LA.EQ.1) THEN
        DO 70 L = 1,NC
          DO 60 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            DO 40 J = 1,NE
              IF (JCE.EQ.L) GO TO 50
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
   40       CONTINUE
   50       JCN(K) = JCE
            IRN(K) = ICE
   60     CONTINUE
   70   CONTINUE
      ELSE
        DO 120 L = 1,NC
          DO 110 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            ACE = A(K)
            DO 90 J = 1,NE
              IF (JCE.EQ.L) GO TO 100
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
              ACEP = A(LOC)
              A(LOC) = ACE
              ACE = ACEP
   90       CONTINUE
  100       JCN(K) = JCE
            IRN(K) = ICE
            A(K) = ACE
  110     CONTINUE
  120   CONTINUE
      END IF
  130 CONTINUE
      RETURN
      END
C**********************************************************
      SUBROUTINE MC59C(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
      INTEGER LA,NC,NE,NR
      REAL A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
      REAL ACE,ACEP
      INTEGER I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP
      DO 10 J = 1,NC
        IP(J) = 0
   10 CONTINUE
      IF (LA.GT.1) THEN
        DO 20 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
          IRN(K) = JCN(K)
   20   CONTINUE
        IP(NC+1) = NE + 1
        IP(1) = IP(1) + 1
        DO 30 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
   30   CONTINUE
        DO 50 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 40 J = J1,J2
            K = IRN(J)
            L = IP(K) - 1
            JCN(J) = L
            IRN(J) = I
            IP(K) = L
   40     CONTINUE
   50   CONTINUE
        IP(NC+1) = NE + 1
        DO 70 J = 1,NE
          LOC = JCN(J)
          IF (LOC.EQ.0) GO TO 70
          ICE = IRN(J)
          ACE = A(J)
          JCN(J) = 0
          DO 60 K = 1,NE
            LOCP = JCN(LOC)
            ICEP = IRN(LOC)
            ACEP = A(LOC)
            JCN(LOC) = 0
            IRN(LOC) = ICE
            A(LOC) = ACE
            IF (LOCP.EQ.0) GO TO 70
            ICE = ICEP
            ACE = ACEP
            LOC = LOCP
   60     CONTINUE
   70   CONTINUE
      ELSE
        DO 90 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
   90   CONTINUE
        IP(NC+1) = NE + 1
        IP(1) = IP(1) + 1
        DO 100 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
  100   CONTINUE
        DO 120 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 110 J = J1,J2
            K = JCN(J)
            L = IP(K) - 1
            IRN(L) = I
            IP(K) = L
  110     CONTINUE
  120   CONTINUE
      END IF
      RETURN
      END
C**********************************************************
      SUBROUTINE MC59D(NC,NE,IRN,IP,LA,A)
      INTEGER LA,NC,NE
      REAL A(LA)
      INTEGER IRN(NE),IP(NC)
      REAL ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
      INTRINSIC ABS
      IF (LA.GT.1) THEN
        KMAX = NE
        DO 50 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 40
          KOR = KMAX
          DO 30 KDUMMY = KLO,KMAX
            ACE = A(KOR-1)
            ICE = IRN(KOR-1)
            DO 10 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 20
              IRN(K-1) = IK
              A(K-1) = A(K)
   10       CONTINUE
            K = KMAX + 1
   20       IRN(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
   30     CONTINUE
   40     KMAX = KLO - 2
   50   CONTINUE
      ELSE
        KMAX = NE
        DO 150 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 140
          KOR = KMAX
          DO 130 KDUMMY = KLO,KMAX
            ICE = IRN(KOR-1)
            DO 110 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 120
              IRN(K-1) = IK
  110       CONTINUE
            K = KMAX + 1
  120       IRN(K-1) = ICE
            KOR = KOR - 1
  130     CONTINUE
  140     KMAX = KLO - 2
  150   CONTINUE
      END IF
      END
C***********************************************************************
      SUBROUTINE MC59E(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)
      INTEGER ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
      REAL A(LA)
      INTEGER IRN(NE),IP(LIP),IW(NR)
      INTEGER I,J,K,KSTART,KSTOP,NZJ
      IDUP = 0
      KNE = 0
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE
      KSTART = IP(1)
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
              IDUP = IDUP + 1
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE
      ELSE
        DO 50 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE MC59F(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT,
     +                  IUP,KNE,ICNTL6,INFO)
      INTEGER ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
      REAL A(LA)
      INTEGER IRN(NE),IP(LIP),IW(LIW),INFO(2)
      INTEGER I,J,K,KSTART,KSTOP,NZJ,LOWER
      IDUP = 0
      IOUT = 0
      IUP = 0
      KNE = 0
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE
      KSTART = IP(1)
      LOWER = 1
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.LT.J) IUP = IUP + 1
            ELSE IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
              IDUP = IDUP + 1
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE
      ELSE
        DO 50 J = 1,NC
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO  40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.GT.1) IUP = IUP + 1
            ELSE IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF
      RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1987 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of a UK ACADEMIC
* Licence, see http://hsl.rl.ac.uk/acuk/cou.html
*
* Please note that for a UK ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 10 Feb 1993
C       Toolpack tool decs employed.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

      SUBROUTINE MC34A(N,IRN,JCOLST,YESA,A,IW)
      INTEGER N
      LOGICAL YESA
      REAL A(*)
      INTEGER IRN(*),IW(*),JCOLST(*)
      INTEGER CKP1,I,I1,I2,II,IPKP1,IPOS,J,JSTART,LENK,NDIAG,NEWTAU,
     +        OLDTAU
      OLDTAU = JCOLST(N+1) - 1
      DO 5 I = 1,N
        IW(I) = 0
    5 CONTINUE
      NDIAG = 0
      DO 20 J = 1,N
        I1 = JCOLST(J)
        I2 = JCOLST(J+1) - 1
        IW(J) = IW(J) + I2 - I1 + 1
        DO 10 II = I1,I2
          I = IRN(II)
          IF (I.NE.J) THEN
            IW(I) = IW(I) + 1
          ELSE
            NDIAG = NDIAG + 1
          END IF
   10   CONTINUE
   20 CONTINUE
      NEWTAU = 2*OLDTAU - NDIAG
      IPKP1 = OLDTAU + 1
      CKP1 = NEWTAU + 1
      DO 40 J = N,1,-1
        I1 = JCOLST(J)
        I2 = IPKP1
        LENK = I2 - I1
        JSTART = CKP1
        IPKP1 = I1
        I2 = I2 - 1
        DO 30 II = I2,I1,-1
          JSTART = JSTART - 1
          IF (YESA) A(JSTART) = A(II)
          IRN(JSTART) = IRN(II)
   30   CONTINUE
        JCOLST(J) = JSTART
        CKP1 = CKP1 - IW(J)
        IW(J) = LENK
   40 CONTINUE
      DO 80 J = N,1,-1
        I1 = JCOLST(J)
        I2 = JCOLST(J) + IW(J) - 1
        DO 60 II = I1,I2
          I = IRN(II)
          IF (I.EQ.J) GO TO 60
          JCOLST(I) = JCOLST(I) - 1
          IPOS = JCOLST(I)
          IF (YESA) A(IPOS) = A(II)
          IRN(IPOS) = J
   60   CONTINUE
   80 CONTINUE
      JCOLST(N+1) = NEWTAU + 1
      RETURN
      END
