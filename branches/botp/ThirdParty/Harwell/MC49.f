
* COPYRIGHT (c) 1993 AEA Technology
*######DATE 10 Feb 1993
C       Toolpack tool decs employed.
C
      SUBROUTINE MC34AD(N,IRN,JCOLST,YESA,A,IW)
      INTEGER N
      LOGICAL YESA
      DOUBLE PRECISION A(*)
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



* COPYRIGHT (c) 1992 Council for the Central Laboratory
*                    of the Research Councils
*######DATE 21 Dec 1992
C       Toolpack tool decs employed.
C
C  EAT 21/6/93 EXTERNAL statement put in for block data on VAXs.
C  JAS 29/4/98 Bug corrected. If (abs(IND).eq.1) IW must be length
C              MAX(NC,NR)+1
C
C
      SUBROUTINE MC49AD(IND,NC,NR,NNZ,IRN,JCN,YESA,LA,A,LIP,IP,LIW,IW,
     +                  IFLAG)
      INTEGER IFLAG,IND,LA,LIP,LIW,NC,NNZ,NR
      LOGICAL YESA
      DOUBLE PRECISION A(LA)
      INTEGER IP(LIP),IRN(NNZ),IW(LIW),JCN(NNZ)
      INTEGER I,J,K,KSTART,KSTOP,NZJ
      EXTERNAL MC49BD,MC49CD
      EXTERNAL MC49DD
      INTRINSIC ABS,MAX
      COMMON /MC49ED/LP,MP,IOUT,JOUT,IDUP,NZOUT
      INTEGER IDUP,IOUT,JOUT,LP,MP,NZOUT
      SAVE /MC49ED/
      IFLAG = 0
      NZOUT = 0
      IOUT = 0
      JOUT = 0
      IDUP = 0
      IF (IND.GT.2 .OR. IND.LT.-2 .OR. IND.EQ.0) THEN
        IFLAG = -1
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) IFLAG
          WRITE (LP,FMT=9010) IND
        END IF
        GO TO 70
      END IF
      IF (NC.LT.1 .OR. NR.LT.1 .OR. NNZ.LT.1) THEN
        IFLAG = -2
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) IFLAG
          WRITE (LP,FMT=9020)
          WRITE (LP,FMT=9030) NC,NR,NNZ
        END IF
        GO TO 70
      END IF
      IF (YESA) THEN
        IF (LA.LT.NNZ) THEN
          IFLAG = -3
          IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) IFLAG
            WRITE (LP,FMT=9040) LA,NNZ
          END IF
          GO TO 70
        END IF
      ELSE
        IF (LA.LT.1) THEN
          IFLAG = -3
          IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) IFLAG
            WRITE (LP,FMT=9050) LA
          END IF
          GO TO 70
        END IF
      END IF
      IF (ABS(IND).EQ.1 .AND. LIW.LT.MAX(NR,NC)+1) THEN
        IFLAG = -4
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) IFLAG
          WRITE (LP,FMT=9060) LIW,MAX(NR,NC) + 1
        END IF
        GO TO 70
      ELSE IF (ABS(IND).EQ.2 .AND. LIW.LT.NR+1) THEN
        IFLAG = -4
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) IFLAG
          WRITE (LP,FMT=9060) LIW,NR + 1
        END IF
        GO TO 70
      END IF
      IF (ABS(IND).EQ.1 .AND. LIP.LT.NC+1) THEN
        IFLAG = -5
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) IFLAG
          WRITE (LP,FMT=9070) LIP,NC + 1
        END IF
        GO TO 70
      ELSE IF (ABS(IND).EQ.2 .AND. LIP.LT.MAX(NR,NC)+1) THEN
        IFLAG = -5
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) IFLAG
          WRITE (LP,FMT=9070) LIP,MAX(NR,NC) + 1
        END IF
        GO TO 70
      END IF
      IF (IND.LT.0) THEN
        NZOUT = NNZ
        GO TO 20
      END IF
      DO 10 K = 1,NNZ
        I = IRN(K)
        J = JCN(K)
        IF (I.GT.NR .OR. I.LT.1) THEN
          IOUT = IOUT + 1
        ELSE IF (J.GT.NC .OR. J.LT.1) THEN
          JOUT = JOUT + 1
        ELSE
          NZOUT = NZOUT + 1
          IRN(NZOUT) = I
          JCN(NZOUT) = J
          IF (YESA) A(NZOUT) = A(K)
        END IF
   10 CONTINUE
      IF (IOUT.GT.0) THEN
        IFLAG = 2
        IF (MP.GT.0) THEN
          WRITE (MP,FMT=9080) IFLAG
          WRITE (MP,FMT=9090) IOUT
        END IF
      END IF
      IF (JOUT.GT.0) THEN
        IFLAG = 3
        IF (MP.GT.0) THEN
          WRITE (MP,FMT=9080) IFLAG
          WRITE (MP,FMT=9110) JOUT
        END IF
      END IF
      IF (IOUT+JOUT.EQ.NNZ) THEN
        NZOUT = 0
        GO TO 70
      END IF
   20 CONTINUE
      IF (ABS(IND).EQ.1) THEN
        CALL MC49BD(NC,NZOUT,IRN,JCN,YESA,LA,A,IP,IW)
      ELSE
        CALL MC49BD(NR,NZOUT,JCN,IRN,YESA,LA,A,IW,IP)
        CALL MC49CD(NC,NR,NZOUT,IRN,JCN,YESA,LA,A,IP,IW)
      END IF
      IF (IND.GT.0) THEN
        NZOUT = 0
        KSTART = 1
        NZJ = 0
        DO 30 I = 1,NR
          IW(I) = 0
   30   CONTINUE
        DO 50 J = 1,NC
          KSTOP = IP(J+1) - 1
          IP(J+1) = IP(J)
          DO 40 K = KSTART,KSTOP
            I = IRN(K)
            IF (IW(I).LE.NZJ) THEN
              NZOUT = NZOUT + 1
              IRN(NZOUT) = I
              IF (YESA) A(NZOUT) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = NZOUT
            ELSE
              IDUP = IDUP + 1
              IF (YESA) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   40     CONTINUE
          KSTART = KSTOP + 1
          NZJ = NZOUT
   50   CONTINUE
        IF (IDUP.GT.0) THEN
          IFLAG = 1
          IF (MP.GT.0) THEN
            WRITE (MP,FMT=9080) IFLAG
            WRITE (MP,FMT=9100) IDUP
          END IF
        END IF
      END IF
   70 RETURN
 9000 FORMAT (/,' *** ERROR RETURN FROM MC49A/AD *** IFLAG = ',I2)
 9010 FORMAT (1X,'IND=',I2,' IS OUT OF RANGE')
 9020 FORMAT (1X,'NC, NR, OR, NNZ IS OUT OF RANGE')
 9030 FORMAT (1X,'NC=',I6,' NR=',I6,' NNZ=',I10)
 9040 FORMAT (1X,'INCREASE LA FROM',I8,' TO AT LEAST ',I8)
 9050 FORMAT (1X,'INCREASE LA FROM',I8,' TO AT LEAST 1')
 9060 FORMAT (1X,'INCREASE LIW FROM',I8,' TO AT LEAST ',I8)
 9070 FORMAT (1X,'INCREASE LIP FROM',I8,' TO AT LEAST ',I8)
 9080 FORMAT (/,' *** WARNING MESSAGE FROM MC49A/AD *** IFLAG = ',I2)
 9090 FORMAT (1X,I6,' ENTRIES IN IRN SUPPLIED BY THE USER WERE OUT OF ',
     +       /,'       RANGE AND WERE IGNORED BY THE ROUTINE')
 9100 FORMAT (1X,I6,' DUPLICATE ENTRIES WERE SUPPLIED BY THE USER')
 9110 FORMAT (1X,I6,' ENTRIES IN JCN SUPPLIED BY THE USER WERE OUT OF ',
     +       /,'       RANGE AND WERE IGNORED BY THE ROUTINE')
      END
C***********************************************************************
      SUBROUTINE MC49BD(NC,NNZ,IRN,JCN,YESA,LA,A,IP,IW)
      INTEGER LA,NC,NNZ
      LOGICAL YESA
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NNZ),IW(NC+1),JCN(NNZ)
      DOUBLE PRECISION ACE,ACEP
      INTEGER ICE,ICEP,J,JCE,JCEP,K,L,LOC
      DO 10 J = 1,NC + 1
        IW(J) = 0
   10 CONTINUE
C**      COUNT THE NUMBER OF ELEMENTS IN EACH COLUMN AND STORE IN IW.
      DO 20 K = 1,NNZ
        J = JCN(K)
        IW(J) = IW(J) + 1
   20 CONTINUE
C**      PUT INTO IP AND IW THE POSITIONS WHERE EACH COLUMN
      IP(1) = 1
      DO 30 J = 2,NC + 1
        IP(J) = IW(J-1) + IP(J-1)
        IW(J-1) = IP(J-1)
   30 CONTINUE
C******  REORDER THE ELEMENTS INTO COLUMN ORDER.
      IF (YESA) GO TO 80
      DO 70 L = 1,NC
        DO 60 K = IW(L),IP(L+1) - 1
          ICE = IRN(K)
          JCE = JCN(K)
          DO 40 J = 1,NNZ
            IF (JCE.EQ.L) GO TO 50
            LOC = IW(JCE)
            JCEP = JCN(LOC)
            ICEP = IRN(LOC)
            IW(JCE) = LOC + 1
            JCN(LOC) = JCE
            IRN(LOC) = ICE
            JCE = JCEP
            ICE = ICEP
   40     CONTINUE
   50     JCN(K) = JCE
          IRN(K) = ICE
   60   CONTINUE
   70 CONTINUE
      GO TO 130
   80 CONTINUE
      DO 120 L = 1,NC
        DO 110 K = IW(L),IP(L+1) - 1
          ICE = IRN(K)
          JCE = JCN(K)
          ACE = A(K)
          DO 90 J = 1,NNZ
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
   90     CONTINUE
  100     JCN(K) = JCE
          IRN(K) = ICE
          A(K) = ACE
  110   CONTINUE
  120 CONTINUE
  130 CONTINUE
      RETURN
      END
C**********************************************************
      SUBROUTINE MC49CD(NC,NR,NNZ,IRN,JCN,YESA,LA,A,IP,IW)
      INTEGER LA,NC,NNZ,NR
      LOGICAL YESA
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NNZ),IW(NR+1),JCN(NNZ)
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,K,L,LOC,LOCP
      DO 10 J = 1,NC
        IP(J) = 0
   10 CONTINUE
      IF (.NOT.YESA) GO TO 80
      DO 20 K = 1,NNZ
        I = JCN(K)
        IP(I) = IP(I) + 1
        IRN(K) = JCN(K)
   20 CONTINUE
      IP(NC+1) = NNZ + 1
      IP(1) = IP(1) + 1
      DO 30 J = 2,NC
        IP(J) = IP(J) + IP(J-1)
   30 CONTINUE
      DO 50 I = NR,1,-1
        DO 40 J = IW(I),IW(I+1) - 1
          K = IRN(J)
          L = IP(K) - 1
          JCN(J) = L
          IP(K) = L
          IRN(J) = I
   40   CONTINUE
   50 CONTINUE
      IP(NC+1) = NNZ + 1
      DO 70 J = 1,NNZ
        LOC = JCN(J)
        IF (LOC.EQ.0) GO TO 70
        ICE = IRN(J)
        ACE = A(J)
        JCN(J) = 0
        DO 60 K = 1,NNZ
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
   60   CONTINUE
   70 CONTINUE
      GO TO 130
   80 CONTINUE
      DO 90 K = 1,NNZ
        I = JCN(K)
        IP(I) = IP(I) + 1
   90 CONTINUE
      IP(NC+1) = NNZ + 1
      IP(1) = IP(1) + 1
      DO 100 J = 2,NC
        IP(J) = IP(J) + IP(J-1)
  100 CONTINUE
      DO 120 I = NR,1,-1
        DO 110 J = IW(I),IW(I+1) - 1
          K = JCN(J)
          L = IP(K) - 1
          IP(K) = L
          IRN(L) = I
  110   CONTINUE
  120 CONTINUE
  130 RETURN
      END
C**********************************************************************
      BLOCK DATA MC49DD
      COMMON /MC49ED/LP,MP,IOUT,JOUT,IDUP,NZOUT
      INTEGER IDUP,IOUT,JOUT,LP,MP,NZOUT
      SAVE /MC49ED/
      DATA LP/6/,MP/6/
      END



