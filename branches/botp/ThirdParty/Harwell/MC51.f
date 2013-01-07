C COPYRIGHT (c) 1995 Council for the Central Laboratory
*                    of the Research Councils
C######DATE 30 November 1995
CCC
C--------------------------------------------------------------------
C            HARWELL SUBROUTINE LIBRARY Release 12 (1995)
C        --
C-             Copyright Rutherford Appleton Laboratory
C        --
C--------------------------------------------------------------------
C********************************************
C SCALING : driver plus routines doing various
C           scalings
C********************************************
       SUBROUTINE MC51AD (N, NZ, NSCA, ASPK, IRN, ICN, COLSCA, ROWSCA,
     *                    S, MAXS, ICNTL, INFO)
CCCC*************************************************************
C**********************************
      INTEGER N, NZ, NSCA, MAXS
      INTEGER IRN(NZ), ICN(NZ)
      INTEGER ICNTL(20), INFO(20)
      DOUBLE PRECISION    ASPK(NZ), COLSCA(*), ROWSCA(*)
      DOUBLE PRECISION    S(MAXS)
C***** working variables
      INTEGER MPRINT,LP, MP
      INTEGER ISPW1, IWNOR
      INTEGER I, K, ITOT
      LOGICAL PROK
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      EXTERNAL MC51VD, MC51WD, MC51XD, MC51YD
C***************************************************************
      LP     = ICNTL(1)
      MP     = ICNTL(2)
      MPRINT = ICNTL(3)
      PROK   = (MPRINT.GE.0)
      IF (PROK) WRITE(MPRINT,101)
 101    FORMAT(/'****** Scaling of original matrix '/)
CC
        IF (NSCA.EQ.1) THEN
         IF (PROK)
     *    WRITE (MPRINT,'(A)') 'Diagonal scaling '
        ELSEIF (NSCA.EQ.2) THEN
         IF (PROK)
     *   WRITE (MPRINT,'(A)') 'Scaling based on MC29'
        ELSEIF (NSCA.EQ.3) THEN
         IF (PROK)
     *   WRITE (MPRINT,'(A)') 'Column scaling'
        ELSEIF (NSCA.EQ.4) THEN
         IF (PROK)
     *   WRITE (MPRINT,'(A)') 'Row and column scaling'
        ELSEIF (NSCA.EQ.5) THEN
         IF (PROK)
     *   WRITE (MPRINT,'(A)') 'MC29 followed by row and column scaling'
        ELSEIF (NSCA.EQ.6) THEN
         IF (PROK)
     *   WRITE (MPRINT,'(A)') 'MC29 followed by column scaling'
        ENDIF
        DO 10 I=1,N
            COLSCA(I) = ONE
            ROWSCA(I) = ONE
 10     CONTINUE
C********************************
        IF ((NSCA.EQ.5).OR.(NSCA.EQ.6)) THEN
          ITOT = 5*N + NZ
          IF (ITOT.GT.MAXS) GOTO 400
          ISPW1 = MAXS - NZ
          DO 15 K=1,NZ
           S(ISPW1+K-1) = ASPK(K)
  15      CONTINUE
        ELSE
          ISPW1 = MAXS
          ITOT  = 5*N
          IF (ITOT.GT.MAXS) GOTO 400
        ENDIF
        IWNOR = ISPW1 - 5*N
          IF (NSCA.EQ.1) THEN
            CALL MC51VD(N,NZ,ASPK,IRN,ICN,S(IWNOR),
     *        COLSCA,ROWSCA,MPRINT)
          ELSEIF (NSCA.EQ.2) THEN
            CALL MC51WD(N,NZ,ASPK,IRN,ICN,
     *      ROWSCA,COLSCA,S(IWNOR),MPRINT,MP,NSCA)
          ELSEIF (NSCA.EQ.3) THEN
            CALL MC51YD(N,NZ,ASPK,IRN,ICN,S(IWNOR),COLSCA,
     *      MPRINT)
          ELSEIF (NSCA.EQ.4) THEN
            CALL MC51XD (N,NZ,IRN,ICN,ASPK,
     *      S(IWNOR),S(IWNOR+N),COLSCA,ROWSCA,MPRINT)
          ELSEIF (NSCA.EQ.5) THEN
            CALL MC51WD(N,NZ,S(ISPW1),IRN,ICN,
     *      ROWSCA,COLSCA,S(IWNOR),MPRINT,MP,NSCA)
            CALL MC51YD(N,NZ,S(ISPW1),IRN,ICN,S(IWNOR),COLSCA,
     *          MPRINT)
          ELSEIF (NSCA.EQ.6) THEN
            CALL MC51WD(N,NZ,S(ISPW1),IRN,ICN,
     *      ROWSCA,COLSCA,S(IWNOR),MPRINT,MP,NSCA)
            CALL MC51XD (N,NZ,IRN,ICN,S(ISPW1),
     *      S(IWNOR),S(IWNOR+N),COLSCA,ROWSCA,MPRINT)
          ENDIF
      GOTO 500
C** ERROR return message
 400  INFO(1) = -5
      INFO(2) = ITOT
      IF ((LP.GE.0).AND.(ICNTL(4).GE.1))
     * WRITE(LP,'(/A)') '*** Error: Not enough space to scale matrix'
 500  RETURN
      END
CCC
      SUBROUTINE MC51BD (N, KASE, X, EST, W, IW, ITER, J, JLAST,
     *                   JUMP)
      INTEGER N, IW(N), KASE
      DOUBLE PRECISION W(N), X(N), EST
      INTRINSIC DBLE, ABS, NINT, SIGN
      INTEGER IDAMAX
      INTEGER ITMAX
      PARAMETER (ITMAX = 5)
      INTEGER I, ITER, J, JLAST, JUMP
      DOUBLE PRECISION ALTSGN, TEMP, ZERO, ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      EXTERNAL IDAMAX
      IF (KASE .EQ. 0) THEN
         DO 10 I = 1,N
            X(I) = ONE/DBLE(N)
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      ENDIF
      GOTO (100, 200, 300, 400, 500) JUMP
  100 CONTINUE
      IF (N .EQ. 1) THEN
         W(1) = X(1)
         EST = ABS(W(1))
         GOTO 510
      ENDIF
      DO 110 I = 1,N
         X(I) = SIGN(ONE,X(I))
         IW(I) = NINT(X(I))
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
         IF ( NINT( SIGN(ONE,X(I)) ) .NE. IW(I) ) GOTO 330
  320 CONTINUE
      GOTO 410
  330 CONTINUE
      DO 340 I = 1,N
         X(I) = SIGN(ONE,X(I))
         IW(I) = NINT(X(I))
  340 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
  400 CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF (   (  ABS(X(JLAST)) .NE. ABS(X(J)) ) .AND.
     +      (ITER .LT. ITMAX)   ) THEN
         ITER = ITER + 1
         GOTO 220
      ENDIF
  410 CONTINUE
      EST = ZERO
      DO 420 I = 1,N
         EST = EST + ABS(W(I))
  420 CONTINUE
      ALTSGN = ONE
      DO 430 I = 1,N
         X(I) = ALTSGN * (ONE+DBLE(I-1)/DBLE(N-1))
         ALTSGN = -ALTSGN
  430 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
  500 CONTINUE
      TEMP = ZERO
      DO 520 I = 1,N
         TEMP = TEMP + ABS(X(I))
  520 CONTINUE
      TEMP = 2.0*TEMP/DBLE(3*N)
      IF (TEMP. GT. EST) THEN
         DO 530 I = 1,N
            W(I) = X(I)
  530    CONTINUE
         EST = TEMP
      ENDIF
  510 KASE = 0
      RETURN
      END
CCC
CC
CC
      SUBROUTINE MC51CD
      RETURN
      END
      SUBROUTINE MC51DD
      RETURN
      END
CCC
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MC51ED(NTASKS,NPROCW,N,IW,LIW,A,LA,
     *                NSTK,PERM,IFLAG,ND,FILS,FRERE,MAXFRT,
     *                NTOTPV,PTRIST,PTRAST,PTRARW,PTRAIW,
     *                IPTA,TLKJ,ALOCLU,AGAIN,ITLOC,
     *                IERROR,IPOOL, LPOOL, IPOOLB, LPOOLB,
     *                NIRBDU,IPTBEG,LENBUD,
     *                LENA, IPTEND,RINFO,LCKS,
     *                ISPON1,ISPON2,ISPON3,ISPON4,NBNODE,NBACTI,
     *                ENDTAG, POSFAC ,IWPOS, LRLU, IPTRLU, IPTINI,
     *                LRLUS, NBLEC, III, IIIB, LEAF, LEAFB, NBROOT,
     *                UU, ICNTL, PTLUST, NSTEPS, ISTEP, INFO, KEEP)
C**********************************************
      INTEGER NBUD
      PARAMETER (NBUD=29)
      INTEGER NTASKS,NPROCW,N,IFLAG,NTOTPV,MAXFRT,LA,LIW,
     *        IERROR,NIRBDU,IPTBEG,LENA, IPTEND,
     *        NSTEPS, ISTEP, INFO(20)
      INTEGER ISPON1,ISPON2,ISPON3,ISPON4,NBNODE,NBACTI
      DOUBLE PRECISION A(LA)
      INTEGER TLKJ(N),ITLOC(N), KEEP(50)
      INTEGER IW(LIW), PERM(N), NSTK(N), IPTA(NBUD)
      INTEGER PTRARW(N), PTRAIW(N), ND(N)
      INTEGER FILS(N),FRERE(N),PTRIST(N),PTRAST(N)
      INTEGER RC,POSELT,LENBUD,ICNTL(20),PTLUST(NSTEPS)
      INTEGER LPOOL, LPOOLB, IPOOL(LPOOL), IPOOLB(LPOOLB)
      INTEGER LCKS(20)
      DOUBLE PRECISION RINFO(20)
      LOGICAL ENDTAG
      INTEGER POSFAC,IWPOS,LRLU,
     *    IPTRLU, IPTINI, LRLUS,
     *    NBLEC, III, IIIB,
     *    LEAF, LEAFB, NBROOT
      DOUBLE PRECISION UU
      LOGICAL IOK,AGAIN,ALOCLU
      INTRINSIC MOD
C***************************************
C***
C**************************************
C*********************
C*********************
      DOUBLE PRECISION ZERO
      INTEGER I,NFRONT,NPIV,INODE,INOPV,IOLDPS
      INTEGER LKJIB,IFINB,NASS,NEL1,IROW,NPIVB
      INTEGER IN,IF,NPIVE
      EXTERNAL MC51FD, MC51GD, MC51HD, MC51ID, MC51MD, MC51ND,
     *         MC51OD, MC51PD, MC51QD, MC51RD, MC51SD, MC51TD
      INTEGER MAXFRW, NPVW, NOFFW, NELVAW, NCMPBW
      DOUBLE PRECISION    OPASSW, OPELIW
      LOGICAL LASTSON
      PARAMETER (ZERO=0.0D0)
CFPP$ NOAUTOEXPAND R
      DO 1234 I=1,N
       ITLOC(I)=0
 1234 CONTINUE
      MAXFRW = 0
      NPVW   = 0
      NOFFW  = 0
      NELVAW = 0
      NCMPBW  = 0
      OPASSW = ZERO
      OPELIW = ZERO
 123  CONTINUE
      IF (IFLAG.LT.0) GOTO 500
C***************************************************
CC
CC Suppress the following 5 lines  in a uniproc version of the code
       IF (NBNODE.LE.-1) THEN
          GOTO 500
       ENDIF
      IF (IIIB.NE.LEAFB) THEN
        I = IIIB
        RC = 1
        IIIB = IIIB + 1
        IF (IIIB.GT.LPOOLB) IIIB = 1
      ELSE
        IF (III.NE.LEAF) THEN
          I = III
          NFRONT = ND(IPOOL(I))
          IF (NBNODE.GE.NPROCW) THEN
            RC = -1
          ELSE
            RC = 0
            III = III + 1
            IF (III.GT.LPOOL) III = 1
             NBNODE = NBNODE+1
          ENDIF
         ELSE
          RC = -1
         ENDIF
      ENDIF
C********************************************************
      IF (RC.GE.0) GO TO 124
      IF (ENDTAG) GO TO 500
      GO TO 123
 124  IF (RC.EQ.1) THEN
        INODE=IPOOLB(I)
      ELSE
        INODE=IPOOL(I)
      ENDIF
      IF (INODE.GT.N) GO TO 55
      IF (INODE.LT.-N) GO TO 555
C*********************************
C*********************************
 234  CALL MC51FD(NTASKS,N,INODE,IW,LIW,A,LA,
     *        IFLAG,IERROR,ND,
     *        FILS,FRERE,MAXFRW,OPASSW,
     *        IPTA,PTRIST,PTRAST,PTRARW,PTRAIW,ITLOC,
     *        TLKJ(INODE),NIRBDU,AGAIN,ALOCLU,LENBUD,
     *        LENA, IPTEND,IPTBEG,NBNODE,
     *        NBACTI, LCKS, LRLU, IPTRLU, IPTINI,
     *        IWPOS, POSFAC, LRLUS, NBLEC,
     *        ICNTL, KEEP)
      IF (IFLAG.LT.0) GOTO 500
      IF (INODE.LT.0) THEN
        IPOOL(LEAF) = -INODE
        LEAF = LEAF + 1
        IF (LEAF.GT.LPOOL) LEAF = 1
       GOTO 123
      ELSE
         ISTEP = ISTEP + 1
         PTLUST(ISTEP) = PTRIST(INODE)
      ENDIF
 20   CALL MC51ID(N,INODE,IW,LIW,A,LA,INOPV,NOFFW,
     *                IFLAG,PTRIST,PTRAST,UU)
      IF (IFLAG.LT.0) GOTO 500
      IF (INOPV.EQ.1) GO TO 50
      IF (INOPV.EQ.2) THEN
CC NO MORE PIVOTS CAN BE SELECTED BUT UPDATING OPERATIONS HAVE
CC STILL TO BE PERFORMED. WE MEET AN UNSTABLE LINE FOR PIVOTING
CC WHICH IS NOT IN THE LAST BLOCK OF THE KJI FACTORIZATION.
CC SO WE HAVE TO UPDATE THE REMAINING BLOCKS WITH RESPECT TO THE
CC SELECTED VARIABLES OF THE CURRENT BLOCK.
CC
         CALL MC51SD(N,INODE,IW,LIW,A,LA,
     *            PTRIST,PTRAST,TLKJ(INODE),KEEP(4))
         GOTO 20
      ENDIF
C***
      NPVW = NPVW + 1
C**************************************
      IOLDPS = PTRIST(INODE)
      IF (IW(IOLDPS+2).LE.1) THEN
       CALL MC51OD(N,INODE,IW,LIW,A,LA,
     *                 PTRIST,PTRAST)
       IW(IOLDPS+1) = IW(IOLDPS+1) + 1
       GO TO 61
      ENDIF
C*******
C*******
         LKJIB = IABS(TLKJ(INODE))
       CALL MC51MD(N,INODE,IW,LIW,A,LA,
     *             PTRIST,PTRAST,IFINB,LKJIB,KEEP(4))
       IW(IOLDPS+1) = IW(IOLDPS+1) + 1
C*******
C*******
       IF (IFINB.EQ.0) GOTO 20
       IF (IFINB.EQ.(-1)) GOTO 50
C*******
C******
       POSELT = PTRAST(INODE)
       NFRONT = IW(IOLDPS)
       NPIV   = IW(IOLDPS+1)
       NASS   = IW(IOLDPS+2)
       NEL1   = NASS - NPIV
       LKJIB  = TLKJ(INODE)
       IF ((LKJIB.LT.0).OR.(NEL1.LT.ISPON1)) THEN
CC
CC we don't parallelize this phase
CC we update the remaining block of rows
CC so that we can eliminate another block of rows
CC
         LKJIB = IABS(LKJIB)
         CALL MC51QD(A,LA,
     *           NFRONT,NPIV,NASS,POSELT,LKJIB)
         GO TO 20
       ENDIF
CC
CC
      PERM(INODE)  = NEL1/ISPON2
      IF (MOD(NEL1,ISPON2).NE.0) PERM(INODE)=PERM(INODE)+1
C*****************************************
C**************************************
      DO 556 I = 1,NEL1,ISPON2
        IPOOLB(LEAFB) = -(N+1)*I-INODE
        LEAFB        = LEAFB +1
        IF (LEAFB.GT.LPOOLB) LEAFB = 1
 556  CONTINUE
C***************************************
C**************************************
      GO TO 123
 555  INODE  = -INODE
      IROW   = INODE/(N+1)
      INODE  = INODE - (N+1)*IROW
         LKJIB = IABS(TLKJ(INODE))
      CALL MC51RD(N,INODE,IW,LIW,A,LA,IROW,
     *       PTRIST,PTRAST,LKJIB,ISPON2)
      IOK = .FALSE.
C**************************************
C***************************************
      PERM(INODE) = PERM(INODE) - 1
      IF (PERM(INODE).EQ.0) IOK = .TRUE.
C***************************************
C**************************************
      IF (.NOT.IOK) GO TO 123
      GOTO 20
 50   IOLDPS = PTRIST(INODE)
      POSELT = PTRAST(INODE)
      NFRONT = IW(IOLDPS)
      NPIV   = IW(IOLDPS+1)
      NASS   = IW(IOLDPS+2)
      IF (NPIV.LE.0) GO TO 60
      NEL1   = NFRONT - NASS
C----
C----
      IF (NEL1.LE.0) GO TO 60
C***
      IF ((NEL1.LT.ISPON3).OR.(NPIV.LT.KEEP(6))) THEN
        IROW = 1
        CALL MC51PD(A,LA,IROW,NFRONT,
     *      NPIV,NASS,POSELT,ISPON3,ISPON4,KEEP(6))
       GO TO 60
      ENDIF
      PERM(INODE) = NEL1/ISPON4
      IF (MOD(NEL1,ISPON4).NE.0) PERM(INODE)=PERM(INODE)+1
C*****************************************
C**************************************
      DO 56 I = 1,NEL1,ISPON4
        IPOOLB(LEAFB) = (N+1)*I+INODE
        LEAFB        = LEAFB +1
        IF (LEAFB.GT.LPOOLB) LEAFB = 1
 56   CONTINUE
C***************************************
C**************************************
      GO TO 123
 55   IROW   = INODE/(N+1)
      INODE  = INODE - (N+1)*IROW
      IOLDPS = PTRIST(INODE)
      POSELT = PTRAST(INODE)
      NFRONT = IW(IOLDPS)
      NPIV   = IW(IOLDPS+1)
      NASS   = IW(IOLDPS+2)
      CALL MC51PD(A,LA,IROW,NFRONT,
     *     NPIV,NASS,POSELT,ISPON3,ISPON4,KEEP(6))
      IOK = .FALSE.
C**************************************
C***************************************
      PERM(INODE) = PERM(INODE) - 1
      IF (PERM(INODE).EQ.0) IOK = .TRUE.
C***************************************
C**************************************
      IF (.NOT.IOK) GO TO 123
CCCCCC
CCCCCC
 60   IOLDPS = PTRIST(INODE)
      NPIV   = IW(IOLDPS+1)
      NASS   = IW(IOLDPS+2)
CCC
CCC
      IW(IOLDPS+4) = NPIV
      IF (NASS.EQ.NPIV) GOTO 61
 62   CALL MC51HD(N,INODE,IW,LIW,A,LA,INOPV,NOFFW,
     *                PTRIST,PTRAST,UU)
      IF (INOPV.NE.1) THEN
       NPVW = NPVW + 1
C***************************************
       CALL MC51ND(N,INODE,IW,LIW,A,LA,
     *                 PTRIST,PTRAST,IFINB)
       IW(IOLDPS+1) = IW(IOLDPS+1) + 1
       IF (IFINB.EQ.0) GOTO 62
      ENDIF
      IOLDPS = PTRIST(INODE)
      POSELT = PTRAST(INODE)
      NFRONT = IW(IOLDPS)
      NPIV   = IW(IOLDPS+1)
      NASS   = IW(IOLDPS+2)
      NPIVB  = IW(IOLDPS+4)
      NPIVE  = NPIV - NPIVB
      NEL1   = NFRONT - NASS
      IF ((NPIVE.LE.0).OR.(NEL1.EQ.0)) GO TO 61
        CALL MC51TD(A,LA,NPIVB,
     *                NFRONT,NPIV,NASS,POSELT)
C=====================
C=====================
 61   CALL MC51GD(NTASKS,N,INODE,A,LA,IW,LIW,A,
     *      IFLAG,IERROR,OPELIW,NELVAW,
     *       PTRIST,PTRAST,IPTA,ALOCLU,AGAIN,
     *       LENA,IPTEND,IPTBEG,LCKS,NBACTI,
     *       POSFAC,LRLU,IPTRLU,IPTINI,LRLUS,NBLEC,
     *       NCMPBW,KEEP(17),KEEP(22),ICNTL)
      IF (IFLAG.LT.0) GOTO 500
C============================
C============================
      IN = INODE
 30   IN = FRERE(IN)
      IF (IN.GT.0) GO TO 30
      IF (IN.EQ.0) THEN
        NBROOT = NBROOT - 1
        IF (NBROOT.EQ.0) ENDTAG = .TRUE.
         NBACTI = NBACTI-1
         NBNODE = NBNODE-1
       GOTO 123
      ENDIF
      IF = -IN
        NSTK(IF)=NSTK(IF)-1
        LASTSON = (NSTK(IF).EQ.0)
      IF (LASTSON) THEN
       IF ((III.EQ.LEAF).AND.(NBNODE.EQ.1).AND.
     *   (NPROCW.GT.1))  THEN
          NBACTI = NBACTI-1
          NBNODE = -1
          IPOOL(LEAF) = IF
          LEAF = LEAF +1
          IF (LEAF.GT.LPOOL) LEAF = 1
         GOTO 500
       ELSE
         INODE = IF
           NBACTI = NBACTI-1
         GOTO 234
       ENDIF
      ENDIF
C==================================
C==================================
        NBACTI = NBACTI-1
        NBNODE = NBNODE-1
      GO TO 123
500   CONTINUE
C======================================
C======================================
       MAXFRT       = MAX0(MAXFRT,MAXFRW)
       NTOTPV       = NTOTPV + NPVW
       INFO(12)   = INFO(12) + NOFFW
       RINFO(2)     = RINFO(2)  + OPASSW
       RINFO(3)     = RINFO(3)  + OPELIW
       INFO(13)   = INFO(13) + NELVAW
       INFO(14)     = INFO(14) + NCMPBW
      RETURN
C====================
C====================
      END
CCC
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MC51FD(NTASKS,N,INODE,IW,LIW,
     *       A,LA,IFLAG,IERROR,ND, FILS,FRERE,
     *       MAXFRW,OPASSW, IPTA,PTRIST,PTRAST,PTRARW,
     *       PTRAIW,ITLOC, LKJIB, NIRBDU,AGAIN,ALOCLU,
     *       LENBUD, LENA, IPTEND, IPTBEG,
     *       NBNODE, NBACTI, LCKS, LRLU, IPTRLU,
     *       IPTINI, IWPOS, POSFAC, LRLUS, NBLEC,
     *       ICNTL, KEEP)
C==================================
      INTEGER NBUD,ICNTL(20)
      PARAMETER (NBUD=29)
      INTEGER N,LIW,LA
      INTEGER KEEP(50)
      INTEGER NTASKS,IFLAG,IERROR,LKJIB,INODE,MAXFRW,
     *        LENBUD,LENA,IPTEND,IPTBEG,NBNODE,NBACTI,
     *        LRLU, IPTRLU, IPTINI, IWPOS, LRLUS,
     *        NBLEC, POSFAC, SIZFR
      INTEGER IW(LIW), ITLOC(N),NIRBDU,
     *        IPTA(NBUD), PTRARW(N), PTRAIW(N), ND(N),
     *        FILS(N), FRERE(N), PTRIST(N), PTRAST(N)
      INTEGER LCKS(20)
      DOUBLE PRECISION A(LA), OPASSW
      LOGICAL AGAIN,ALOCLU
      INTEGER LP
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER NFRONT,ISTCHK,LSTK,LREQ,LAELL,IPT,LSIZ,IRES,K
      INTEGER LAPOS2,NEWEL1,INEW1,J1,JT1,INEW,NTOTFS,J2
      INTEGER NFS,NELIM,JJ,JJ1,JJ2,J3,J,
     *        NEWEL,IBROT,IORG,I,IP1
      INTEGER IP2,K1,K2,IASTK,IACHK,NFSON1,NFSON,JPOS,ICT11
      INTEGER JK,IJROW,NBCOL,ICT13,NUMORG,IOLDPS,IOLDP2,J4
      INTEGER APOS, APOS2, AINPUT, POSELT, POSEL1, ICT12
      LOGICAL ISONLU,IOK
      INTRINSIC DBLE
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      EXTERNAL MC51KD, MC51UD, MC51ZD
C***************************************
        IN = INODE
        NUMORG = 0
 3      NUMORG = NUMORG + 1
        IN = FILS(IN)
        IF (IN.GT.0) GO TO 3
        NUMSTK =  0
        NASS =  0
        IFSON =  -IN
        ISON = IFSON
        IF (ISON.EQ.0) GO TO 12
 8      NUMSTK = NUMSTK + 1
        NASS = NASS + IW(PTRIST(ISON)+1)
        ISON=FRERE(ISON)
        IF (ISON.GT.0) GO TO 8
 12     NFRONT = ND(INODE) + NASS
        MAXFRW = MAX0(MAXFRW,NFRONT)
        NASS1 = NASS + NUMORG
        LAELL = NFRONT*NFRONT
        IF (NTASKS.EQ.1) THEN
C-------------------------------
C------------------------------
         IF (LRLU.LT.LAELL) GOTO 630
         LRLU     = LRLU - LAELL
         POSELT   = POSFAC
         POSFAC   = POSFAC + LAELL
        ELSE
C---------------------------
C--------------------------
         LREQ = LAELL + 3
CC**********
CC Allocate space of length LREQ
         CALL MC51ZD(A, LA, IPTA, LREQ, IPT, LSIZ, IRES,
     *    LENA,NBACTI,LCKS)
CC*********
         IF (IRES.LT.0) THEN
            NBNODE = NBNODE -1
         ENDIF
         IF (IRES.EQ.-2) GOTO 625
         IF (IRES.EQ.-1) GOTO 620
         POSELT = IPT+3
C----------------------------
C----------------------------
        ENDIF
C**************************************
             NBACTI = NBACTI+1
         POSEL1 = POSELT - NFRONT
         LAPOS2 = POSELT + LAELL - 1
         DO 230 K=POSELT,LAPOS2
           A(K) = ZERO
  230    CONTINUE
CC
C=================================================
         LREQ = 2*NFRONT+ 5
         IOK = .FALSE.
        IF ((IWPOS +LREQ).LE.NIRBDU) IOK = .TRUE.
        IOLDPS = IWPOS
        IWPOS  = IWPOS + LREQ
        IF (.NOT.IOK) GO TO 610
C=================================================
        IOLDP2 = IOLDPS + 4
        NEWEL = IOLDP2 + NASS1
        NEWEL1= NASS1
        IW(IOLDPS) = NFRONT
C==========================================
C==========================================
        IN = INODE
        INEW = IOLDPS + 5
        INEW1= 1
 73     J1 = PTRAIW(IN)+2
        JT1        = IW(J1)
        IW(J1)     = INEW1
        ITLOC(JT1) = INEW1
        IW(INEW)   = JT1
        INEW = INEW + 1
        INEW1= INEW1+ 1
        IN = FILS(IN)
        IF (IN.GT.0) GO TO 73
        IF (NUMSTK.NE.0) THEN
         NTOTFS = NUMORG
         ISON =  IFSON
         ICT11  = IOLDP2 + NFRONT
         DO 100 IELL=1,NUMSTK
          J2    = PTRIST(ISON)
          ISON  = FRERE(ISON)
          LSTK  = IW(J2)
          NFS   = IW(J2+1)
          NELIM = IW(J2+3)
          NFSON = NELIM + LSTK
          J1 = J2 + NFSON + 5 + NELIM
          J2 = J1 + LSTK - 1
          J3 = J1 + NFS -1
          IF (NFS.EQ.0) GO TO 75
          DO 74 JJ=J1,J3
            NTOTFS            = NTOTFS + 1
            JT1               = IW(JJ)
            IW(ICT11+NTOTFS)  = JT1
            ITLOC(JT1)        = NTOTFS
            IW(JJ)            = NTOTFS
            IW(IOLDP2+NTOTFS) = IW(JJ-NFSON)
  74      CONTINUE
  75      J1 = J3 + 1
          IF (NASS1.NE.NFRONT) THEN
           DO 90 JJ=J1,J2
            J = IW(JJ)
            IF (ITLOC(J).EQ.0) THEN
             NEWEL      = NEWEL + 1
             NEWEL1     = NEWEL1+ 1
             IW(NEWEL)  = J
             IW(JJ)     = NEWEL1
             ITLOC(J)   = NEWEL1
            ELSE
             IW(JJ)     = ITLOC(J)
            ENDIF
   90      CONTINUE
          ELSE
           DO 92 JJ=J1,J2
            IW(JJ) = ITLOC(IW(JJ))
   92      CONTINUE
          ENDIF
  100    CONTINUE
        ENDIF
        IBROT = INODE
        DO 180 IORG=1,NUMORG
          J1 = PTRAIW(IBROT)+2
          IBROT = FILS(IBROT)
CCC
CC with new MA41H/HD iw(j1-1) is always the number of nonzeros
CC in the row
CCC
          J2 = J1 + IW(J1-2) - IW(J1-1)
          J1 = J1 +1
          IF (J1.LE.J2) THEN
           DO 170 JJ=J1,J2
            J   = IW(JJ)
            IF (ITLOC(J).EQ.0) THEN
             NEWEL      = NEWEL + 1
             NEWEL1     = NEWEL1+ 1
             IW(NEWEL)  = J
             IW(JJ)     = NEWEL1
             ITLOC(J)   = NEWEL1
            ELSE
             IW(JJ)     = ITLOC(J)
            ENDIF
  170      CONTINUE
          ENDIF
  180   CONTINUE
        IP1   = IOLDPS + NASS1+ 5
        IP2   = IOLDPS + 5 + NFRONT - 1
        DO 183 I=IP1,IP2
          IW(I+NFRONT) = IW(I)
  183   CONTINUE
        IP1   = IOLDPS + 5
        IP2   = IOLDPS + 5 + NUMORG - 1
        DO 184 I=IP1,IP2
          IW(I+NFRONT) = IW(I)
  184   CONTINUE
C*********************************************
C** SET THE MODIFIED VERSION OF ITLOC TO ZERO
C*********************************************
        K1 = IOLDPS + 5 + NUMORG
        K2 = K1 + NFRONT - 1 + NASS
        DO 554 K=K1,K2
         I        = IW(K)
         ITLOC(I) = 0
  554   CONTINUE
C===================================
C===================================
        IF (NUMSTK.EQ.0) GO TO 290
        ISON   = IFSON
        DO 280 IELL=1,NUMSTK
          ISTCHK = PTRIST(ISON)
          LSTK   = IW(ISTCHK)
          SIZFR  = LSTK*LSTK
          OPASSW = OPASSW + DBLE(SIZFR)
          NFS    = IW(ISTCHK+1)
          NELIM  = IW(ISTCHK+3)
          NFSON  = NELIM + LSTK
          J1     = ISTCHK + NFSON + 5 + NELIM
          J2     = J1 + LSTK - 1
          IF (NTASKS.EQ.1) THEN
C-------------------
C-------------------
             IACHK  = PTRAST(ISON)
             NFSON1 = LSTK
          ELSE
C-------------------
C-------------------
            ISONLU = (PTRAST(ISON).LE.IPTBEG)
            IF (ISONLU) THEN
                NBLEC = NBLEC +1
            ENDIF
            IASTK  = PTRAST(ISON)
            IF (ISONLU) THEN
             IACHK =  IASTK
             NFSON1 = LSTK
            ELSE
             IACHK  = IASTK + NFSON*NELIM + NELIM
             NFSON1 = NFSON
            ENDIF
          ENDIF
C----------------------------------------------------
C----------------------------------------------------
          DO 270 JJ=J1,J2
            APOS = POSEL1 + IW(JJ)*NFRONT
            DO 268 JJ1=1,LSTK
             JJ2 = APOS+IW(J1+JJ1-1)-1
             A(JJ2) = A(JJ2) + A(IACHK+JJ1-1)
  268       CONTINUE
            IACHK = IACHK + NFSON1
  270     CONTINUE
C****************************************
          J3 = J1 + NFS
          DO 275 JJ=J3,J2
           IW(JJ) = IW(JJ-NFSON)
  275     CONTINUE
          IF (NFS.NE.0) THEN
            J3 = J3 -1
            DO 278 JJ=J1,J3
             JPOS = IW(JJ) + ICT11
             IW(JJ) = IW(JPOS)
  278       CONTINUE
          ENDIF
C==========================================================
C==========================================================
      IF (NTASKS.EQ.1) THEN
C-------------------
C-------------------
         IPTRLU = IPTRLU + SIZFR
         LRLU   = LRLU + SIZFR
      ELSE
C-------------------
C-------------------
        IF (ISONLU) THEN
          CALL MC51KD(N, A, ISON, PTRAST, IPTBEG,
     *       LRLU, LRLUS, IPTRLU, IPTINI,LCKS)
              NBLEC = NBLEC -1
        ELSE
          IPT = IASTK - 3
          CALL MC51UD(A,LA,IPTA, IPT,
     *       IPTEND, IPTBEG,LCKS)
        ENDIF
       ENDIF
          ISON   = FRERE(ISON)
  280   CONTINUE
C====================================================
C====================================================
C====================================================
C====================================================
  290   IBROT = INODE
        DO 320 IORG=1,NUMORG
          JK = PTRAIW(IBROT)
          AINPUT = PTRARW(IBROT)
          IBROT = FILS(IBROT)
          JJ = JK + 1
          J1 = JJ + 1
          J2 = J1 + IW(JK)
          J3 = J2 + 1
          J4 = J2 - IW(JJ)
          IJROW = IW(J1)
          ICT12 = POSELT - NFRONT+ IJROW - 1
C*****************************
C*****************************
          DO 300 JJ=J1,J2
            APOS2 = ICT12 + IW(JJ)*NFRONT
            A(APOS2) = A(APOS2) + A(AINPUT)
            AINPUT = AINPUT + 1
  300     CONTINUE
          IF (J3.GT.J4) GO TO 320
          ICT13 = POSELT + (IJROW-1)*NFRONT
          NBCOL = J4 - J3 + 1
C*****************************
C*****************************
          DO 310 JJ=1,NBCOL
           JJ1 = ICT13+IW(J3+JJ-1)-1
           A(JJ1) = A(JJ1) + A(AINPUT+JJ-1)
  310     CONTINUE
  320   CONTINUE
C====================================================
C====================================================
        NASS = NASS1
        PTRAST(INODE) = POSELT
        PTRIST(INODE) = IOLDPS
        IW(IOLDPS+2) = NASS1
        IW(IOLDPS+1) = 0
        IW(IOLDPS+3) = -NASS1
        IW(IOLDPS+4) = 1
        IF (NASS.GT.KEEP(3)) THEN
          LKJIB = KEEP(6)
        ELSE
          LKJIB = -KEEP(5)
        ENDIF
      GOTO 640
C*************
C*************
  610 CONTINUE
         NBACTI = NBACTI-1
         NBNODE = NBNODE-1
       IFLAG = -8
       IERROR = LREQ
       IF ((ICNTL(1).GE.0).AND.(ICNTL(4).GE.1)) THEN
       LP = ICNTL(1)
       WRITE(LP,'(/A)') 'Failure in integer allocation during assembly'
       ENDIF
       GO TO 640
  620 CONTINUE
CC
      IF ((LREQ.GT.LENA).AND.(LREQ.GT.LENBUD)) THEN
        IFLAG = -9
        IERROR = LREQ
C***************************************
        IF ((ICNTL(1).GE.0).AND.(ICNTL(4).GE.1)) THEN
        LP = ICNTL(1)
C**************************************
        WRITE(LP,'(/A)') 'Failure in real allocation during assembly'
        WRITE(LP,'(A,I8)') 'Space required in active area was : ',LREQ
C***************************************
        ENDIF
        GOTO 640
      ENDIF
C**************************************
      AGAIN  = .TRUE.
      ALOCLU = .TRUE.
      INODE = -INODE
      GOTO 640
  625 CONTINUE
       IFLAG = -9
       IERROR = LREQ
       GOTO 640
  630  CONTINUE
       IFLAG = -9
       IERROR = LAELL-LRLU
 640  RETURN
      END
CCC
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MC51GD(NTASKS, N, INODE,
     *         AFACT, LA, IW, LIW, A,
     *         IFLAG,IERROR,OPELIW,NELVAW,
     *         PTRIST,PTRAST,IPTA,ALOCLU,AGAIN,
     *         LENA,IPTEND,IPTBEG,LCKS,NBACTI,
     *         POSFAC,LRLU,IPTRLU,IPTINI,LRLUS,NBLEC,
     *         NCMPB,PTLUAC,LINAC,ICNTL)
C**********************************************
C******************************************
      INTEGER NBUD
      PARAMETER (NBUD=29)
      INTEGER NTASKS, N, LA, LIW, INODE,IFLAG,IERROR
      INTEGER LENA,IPTEND,IPTBEG,NCMPB, NBACTI
      INTEGER POSFAC, LRLU, IPTRLU, IPTINI, LRLUS, NBLEC
      INTEGER IW(LIW), IPTA(NBUD), PTLUAC, LINAC
      INTEGER LCKS(20)
      INTEGER ICNTL(20)
      INTEGER PTRIST(N), PTRAST(N), NELVAW
      DOUBLE PRECISION    A(LA), AFACT(LA), OPELIW
      LOGICAL ALOCLU,AGAIN
      INTEGER APOS, POSELT, OPSFAC, LSIZ, LP
      INTEGER IOLDPS,NFRONT,NPIV,NASS,LREQCB,LCONT
      INTEGER IPT,I,ISP,IRES,NPOS,J,J1,J2,LREQLU
      DOUBLE PRECISION    FLOP1
      LOGICAL IOK, FRESP, IOKLU
C***************************************
      INTRINSIC DBLE
      EXTERNAL MC51JD, MC51LD, MC51UD, MC51ZD
      IOLDPS = PTRIST(INODE)
      NFRONT = IW(IOLDPS)
      NPIV = IW(IOLDPS+1)
      NASS = IW(IOLDPS+2)
      LCONT   = NFRONT-NPIV
      POSELT = PTRAST(INODE)
      NELVAW = NELVAW + NASS - NPIV
C=================================
C=================================
      FLOP1  = DBLE(2*NFRONT*NPIV)*DBLE(NFRONT-NPIV-1)+
     *       DBLE(NPIV*(NPIV+1))*DBLE(2*NPIV+1)/DBLE(3)
      FLOP1  = FLOP1 + DBLE(((2*NFRONT-NPIV-1)*NPIV)/2)
      OPELIW = OPELIW + FLOP1
C------------------------------------------------------
      ISP = 2*NPIV*NFRONT - NPIV*NPIV
C------------------------------------
C----------------------------------
C----------------------------------
      IF (NTASKS.EQ.1) THEN
       IW(IOLDPS+4) = POSELT
       IF (LCONT.EQ.0) GOTO 610
       LREQCB = LCONT*LCONT
       IF (LREQCB.GT.LRLU) GOTO 630
       IPTRLU = IPTRLU - LREQCB
       POSFAC = POSFAC - LREQCB
       NPOS   = IPTRLU+1
       PTRAST(INODE) = NPOS
       OPSFAC = POSELT + NPIV*NFRONT + NPIV
       APOS = OPSFAC
       DO 12 I=1, LCONT
        J1= APOS
        DO 11 J=1,LCONT
            A(NPOS) = A(J1)
            NPOS    = NPOS + 1
            J1      = J1 + 1
 11     CONTINUE
        APOS = APOS + NFRONT
 12    CONTINUE
       IF ((NPIV.EQ.0).OR.(LCONT.EQ.1)) GOTO 610
       APOS = POSELT + (NPIV+1)*NFRONT
          DO 14 I=1,LCONT-1
          DO 13 J=1,NPIV
            A(OPSFAC) = A(APOS)
            OPSFAC=OPSFAC+1
            APOS  = APOS+1
 13       CONTINUE
          APOS = APOS + LCONT
 14      CONTINUE
       GOTO 610
      ENDIF
C----------------------------------
C----------------------------------
      IF (NPIV.EQ.0) GOTO 525
      IPT    = POSELT -3
      IOKLU = .FALSE.
         IF (LRLUS.GE.ISP) THEN
            LRLUS = LRLUS -ISP
            IOKLU = .TRUE.
         ENDIF
      IF (IOKLU) THEN
        IOK = .FALSE.
          IF (LRLU.GE.ISP) THEN
             LRLU  = LRLU - ISP
             IOK = .TRUE.
          ENDIF
        IF (.NOT.IOK) THEN
         CALL MC51JD(N,A,PTRAST,ISP,IRES,IPTBEG,
     *       LRLU,IPTRLU,IPTINI,NBLEC,LCKS)
         NCMPB = NCMPB + 1
         IF (IRES.EQ.-1) THEN
           IERROR = ISP
           GOTO 640
         ENDIF
        ENDIF
        IOK = .TRUE.
          OPSFAC = POSFAC
          POSFAC = POSFAC + ISP
          IF ((POSFAC-1).GT.IPTRLU) IOK = .FALSE.
        IF (.NOT.IOK) THEN
          IERROR = ISP
          GO TO 640
        ENDIF
      ELSE
        IF (NPIV.EQ.NFRONT) THEN
           LINAC = LINAC + ISP
           IW (IOLDPS+4) = POSELT
         GOTO 610
        ENDIF
        LREQLU = ISP + 3 + 3
        CALL MC51ZD(A, LA, IPTA, LREQLU, OPSFAC, LSIZ, IRES,
     *    LENA,NBACTI,LCKS)
        IF (IRES.LT.0) GOTO 635
         LINAC = LINAC + ISP
         A(OPSFAC+3) = DBLE(PTLUAC)
         A(OPSFAC+4) = DBLE(ISP)
         A(OPSFAC+5) = DBLE(INODE)
         PTLUAC      = OPSFAC
        OPSFAC = OPSFAC+6
      ENDIF
      IW(IOLDPS+4) = OPSFAC
C------------------------------------------------
C------------------------------------------------
      APOS = POSELT
      J1 = APOS
      J2 = APOS + NPIV*NFRONT - 1
      DO 500 J=J1,J2
            AFACT(OPSFAC) = A(J)
            OPSFAC        = OPSFAC + 1
  500 CONTINUE
      IF (LCONT.GT.0) THEN
          J1 = J2 + 1
          DO 520 I=1,LCONT
            DO 510 J=1,NPIV
             AFACT(OPSFAC) = A(J1)
             J1            = J1+1
             OPSFAC        = OPSFAC + 1
  510       CONTINUE
            J1 = J1 + LCONT
  520     CONTINUE
      ENDIF
C==========================================
C==========================================
  525   LREQCB   = LCONT*LCONT + 4
        IPT    = POSELT -3
        IF (LCONT.NE.0) THEN
C********************************************
          IF (.NOT.ALOCLU) GOTO 610
          IOK = .FALSE.
C*********************************************
CC*********************************************
CCCCCCCCCCCCCCCCCCCCCCCCCC
            IF ((LRLU.GE.LREQCB).AND.(LRLUS.GE.LREQCB)) THEN
               LRLUS = LRLUS -LREQCB
               LRLU  = LRLU - LREQCB
               IOK = .TRUE.
            ENDIF
CC
CCC IF (.NOT.AGAIN).AND.(.NOT.IOK)) contribution block
CC  is not allocated in LU area
CC
           IF ((.NOT.AGAIN).AND.(.NOT.IOK)) GOTO 610
           IF ((AGAIN).AND.(.NOT.IOK)) THEN
CC
            FRESP = .FALSE.
            IF (LRLUS.GE.LREQCB) THEN
             LRLUS = LRLUS -LREQCB
             FRESP = .TRUE.
            ENDIF
            IF (.NOT.FRESP) GOTO 610
            CALL MC51JD(N,A,PTRAST,LREQCB,IRES,IPTBEG,
     *       LRLU,IPTRLU,IPTINI,NBLEC,LCKS)
            NCMPB = NCMPB + 1
            IF (IRES.EQ.-1) THEN
             IERROR = LREQCB
             GOTO 640
            ENDIF
           ENDIF
           CALL MC51LD(N,A,LREQCB,INODE,PTRAST,IPTBEG,
     *         IPTRLU,IPTINI,LCKS)
C**********************************
           NBLEC = NBLEC +1
C**********************************
           APOS = POSELT + NPIV*NFRONT + NPIV
           NPOS = PTRAST(INODE)
           DO 710 I=1,LCONT
              J1 = APOS
              DO 700 J=1,LCONT
                A(NPOS) = A(J1)
                NPOS    = NPOS + 1
                J1    = J1 + 1
  700         CONTINUE
              APOS = APOS + NFRONT
  710      CONTINUE
C**********************************
            NBLEC = NBLEC -1
        ENDIF
C**********************************
CC
CC we free the working space in buddy system
CC
            CALL MC51UD(A,LA,IPTA,IPT, IPTEND, IPTBEG,LCKS)
CC
  610 IW(IOLDPS+3) = NPIV
      IW(IOLDPS+2) = NFRONT
      IW(IOLDPS+1) = NASS - NPIV
      IW(IOLDPS)   = LCONT
      GO TO 650
C**************************************
  630 CONTINUE
C***************************************
      IFLAG = -9
      IERROR = LREQCB-LRLU
      LP     = ICNTL(1)
      IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
      WRITE(LP,'(/A)') 'Error: During stack'
      WRITE(LP,'(A)') 'Failure in reserving real space during stack'
      ENDIF
      GOTO 650
C**************************************
  635 CONTINUE
C***************************************
      IFLAG = -9
      LP     = ICNTL(1)
      IERROR = LREQLU
      IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
      WRITE(LP,'(/A)') 'Error: during stack'
      WRITE(LP,'(A)') 'Failure in reserving real space for stacking '
      WRITE(LP,'(A)') 'LU factors in active area '
      ENDIF
      GOTO 650
C***************************************
  640 CONTINUE
C***************************************
CC
CC
      IFLAG = -9
      LP     = ICNTL(1)
      IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
      WRITE(LP,'(/A)') 'Error: during stack-memory scheme'
      WRITE(LP,'(A)') 'Failure in reserving real space for stacking '
      WRITE(LP,'(A)') 'contribution block after a compress'
      ENDIF
C***************************************
  650 RETURN
      END
CCC
      SUBROUTINE MC51HD(N,INODE,IW,LIW,A,LA,
     *   INOPV,NOFFW,PTRIST,PTRAST,UU)
      INTEGER N,LIW,LA,INODE,INOPV
      DOUBLE PRECISION UU, A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION ZERO, RMAX, AMROW
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION  SWOP
      INTEGER APOS, POSELT, PTRIST(N), PTRAST(N)
      INTEGER NFRONT,NOFFW,IOLDPS,NPIV,NASS,IPIV
      INTEGER NPIVP1,JMAX,J1,J3,JJ,J2,IDIAG,ISW,ISWPS1
      INTEGER ISWPS2,KSW
      INTEGER IDAMAX
      INTRINSIC MAX
      EXTERNAL IDAMAX
        INOPV   = 0
        POSELT  = PTRAST(INODE)
        IOLDPS  = PTRIST(INODE)
        NFRONT  = IW(IOLDPS)
        NPIV    = IW(IOLDPS+1)
        NPIVP1  = NPIV + 1
        NASS    = IW(IOLDPS+2)
C***************
          DO 460 IPIV=NPIVP1,NASS
            APOS = POSELT + NFRONT*NPIV + (IPIV-1)
            JMAX = 1
            AMROW = ZERO
            J1 = APOS
            J3    = NASS -NPIV
            JMAX  = IDAMAX(J3,A(J1),NFRONT)
            JJ    = J1 + (JMAX-1)*NFRONT
            AMROW = ABS(A(JJ))
            RMAX = AMROW
            J1 = APOS +  (NASS-NPIV) * NFRONT
            J3 = NFRONT - NASS
            IF (J3.EQ.0) GOTO 370
            DO 360 JJ=1,J3
              RMAX = MAX(ABS(A(J1)),RMAX)
              J1 = J1 + NFRONT
  360       CONTINUE
  370       IF (RMAX.EQ.ZERO) GO TO 460
            IDIAG = APOS + (IPIV - NPIVP1)*NFRONT
            IF (ABS(A(IDIAG)).GE.UU*RMAX) JMAX = IPIV - NPIV
            IF (ABS(A(IDIAG)).GE.UU*RMAX) GO TO 380
            IF (AMROW.LT.UU*RMAX) GO TO 460
C***************************************
            NOFFW = NOFFW + 1
C**************************************
  380       IF (IPIV.EQ.NPIVP1) GO TO 400
            J1 = POSELT + NPIV
            J3 = POSELT + (IPIV-1)
            DO 390 JJ= 1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J3)
              A(J3) = SWOP
              J1 = J1 + NFRONT
              J3 = J3 + NFRONT
  390       CONTINUE
            ISWPS1 = IOLDPS + 4 + NPIVP1 + NFRONT
            ISWPS2 = IOLDPS + 4 + IPIV + NFRONT
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            J1 = POSELT + NPIV*NFRONT
            J2 = POSELT + (NPIV + JMAX - 1)*NFRONT
            DO 410 KSW=1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + 1
              J2 = J2 + 1
  410       CONTINUE
            ISWPS1 = IOLDPS + 4 + NPIV + 1
            ISWPS2 = IOLDPS + 4 + NPIV + JMAX
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460     CONTINUE
C*************************************************
       INOPV = 1
  420 RETURN
      END
CCC
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
C********************************************************
C***********PIVOT********
      SUBROUTINE MC51ID(N,INODE,IW,LIW,A,LA,
     *    INOPV,NOFFW,IFLAG,PTRIST,PTRAST,UU)
      INTEGER N,LA,LIW,INODE,IFLAG,INOPV,NOFFW
      DOUBLE PRECISION A(LA)
      DOUBLE PRECISION UU
      INTEGER IW(LIW), PTRIST(N), PTRAST(N)
      DOUBLE PRECISION SWOP
      INTEGER APOS, POSELT
      DOUBLE PRECISION ZERO, RMAX, AMROW
      INTEGER NFRONT,IOLDPS,NPIV,NASS,NASSW,IPIV
      INTEGER NPIVP1,JMAX,J1,J3,JJ,J2,IDIAG,ISW,ISWPS1
      INTEGER ISWPS2,KSW
      INTEGER IDAMAX
      INTRINSIC MAX
      PARAMETER (ZERO=0.0D0)
      EXTERNAL IDAMAX
        INOPV   = 0
        POSELT  = PTRAST(INODE)
        IOLDPS  = PTRIST(INODE)
        NFRONT  = IW(IOLDPS)
        NPIV    = IW(IOLDPS+1)
        NPIVP1  = NPIV + 1
        NASS    = IW(IOLDPS+2)
C*****
C*****
        NASSW   = IABS(IW(IOLDPS+3))
C***************
          DO 460 IPIV=NPIVP1,NASSW
            APOS = POSELT + NFRONT*(IPIV-1) + NPIV
            JMAX = 1
            IF (UU.GT.ZERO) GO TO 340
            IF (A(APOS).EQ.ZERO) GO TO 630
            GO TO 380
  340       AMROW = ZERO
            J1 = APOS
            J2 = APOS - NPIV + NASS - 1
CCC
             J3    = NASS -NPIV
             JMAX  = IDAMAX(J3,A(J1),1)
             JJ    = JMAX + J1 - 1
             AMROW = ABS(A(JJ))
CCC
            RMAX = AMROW
            J1 = J2 + 1
            J2 = APOS - NPIV + NFRONT - 1
            IF (J2.LT.J1) GO TO 370
            DO 360 JJ=J1,J2
              RMAX = MAX(ABS(A(JJ)),RMAX)
  360       CONTINUE
  370       IF (RMAX.EQ.ZERO) GO TO 460
            IDIAG = APOS + IPIV - NPIVP1
            IF (ABS(A(IDIAG)).GE.UU*RMAX) JMAX = IPIV - NPIV
            IF (ABS(A(IDIAG)).GE.UU*RMAX) GO TO 380
            IF (AMROW.LT.UU*RMAX) GO TO 460
C**************************************
            NOFFW = NOFFW + 1
C***************************************
  380       IF (IPIV.EQ.NPIVP1) GO TO 400
            J1 = POSELT + NPIV*NFRONT
            J2 = J1 + NFRONT - 1
            J3 = POSELT + (IPIV-1)*NFRONT
            DO 390 JJ=J1,J2
              SWOP = A(JJ)
              A(JJ) = A(J3)
              A(J3) = SWOP
              J3 = J3 + 1
  390       CONTINUE
            ISWPS1 = IOLDPS + 4 + NPIVP1
            ISWPS2 = IOLDPS + 4 + IPIV
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            J1 = POSELT + NPIV
            J2 = POSELT + NPIV + JMAX - 1
            DO 410 KSW=1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + NFRONT
              J2 = J2 + NFRONT
  410       CONTINUE
            ISWPS1 = IOLDPS + 4 + NFRONT + NPIV + 1
            ISWPS2 = IOLDPS + 4 + NFRONT + NPIV + JMAX
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460     CONTINUE
C*************************************************
      IF (NASSW.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
C***************************************
      IFLAG = -10
C**************************************
  420 CONTINUE
      RETURN
      END
CCC
      SUBROUTINE MC51JD(N,A,PTRAST,ISPA,
     *   IRES,IPTBEG, LRLU,IPTRLU,IPTINI,NBLEC,LCKS)
      INTEGER N,LRLU,IPTRLU,IPTINI,NBLEC
      INTEGER IPTBEG,IPREV,INODE,LREQ,IPTNEW
      INTEGER I,II,ILAST,IPTCU,IOLD,LBLOCK,ICB
      DOUBLE PRECISION A(IPTBEG)
      INTEGER  IRES,PTRAST(N),ISPA
      LOGICAL IOK
      INTEGER LCKS(20)
      INTRINSIC DBLE
C***
C***
      LCKS(1) = 0
CCC
CCC
      IRES = 1
      IOK  = .FALSE.
C*************************************************
 220   CONTINUE
      IF (NBLEC.GT.0) THEN
CC
CC
         GOTO 220
      ENDIF
C*****
           IF (LRLU.GE.ISPA) THEN
            IOK = .TRUE.
            LRLU = LRLU -ISPA
           ENDIF
      IF (IOK) GOTO 410
C*************************
CCCCC
CCCCC
CCC
      IF (IPTINI.EQ.IPTBEG) THEN
       IRES =-1
       GOTO 410
      ENDIF
C***************************
      LBLOCK = IPTBEG - IPTINI + 1
      IPREV  =  INT (A(IPTINI +2) + 0.5)
      INODE  =  INT (A(IPTINI +3) + 0.5)
      LREQ   = INT (A(IPTINI) + 0.5)
      IPTNEW = IPTBEG - LREQ + 1
      IF (LREQ.LT.LBLOCK) THEN
       ILAST  = IPTINI + LREQ - 1
       II     = IPTBEG
       IPTINI = IPTNEW
       IF (2*LREQ.LT.LBLOCK) THEN
        DO 10 I=1,LREQ
         A(II) = A(ILAST)
         ILAST = ILAST -1
         II    = II - 1
 10     CONTINUE
       ELSE
        DO 20 I=1,LREQ
         A(II) = A(ILAST)
         ILAST = ILAST -1
         II    = II - 1
 20     CONTINUE
       ENDIF
       PTRAST(INODE) = IPTNEW+4
      ENDIF
      DO 50 ICB=1,N
       IF (IPREV.EQ.0) GOTO 400
       A(IPREV+1)    = DBLE(IPTNEW)
       IPTCU = IPREV
       IOLD   = IPTNEW
       IPTNEW = IPTCU
       LBLOCK = IOLD - IPTCU
       IPREV  =  INT (A(IPTCU+2) + 0.5)
       INODE  =  INT (A(IPTCU+3) + 0.5)
       LREQ   =  INT (A(IPTCU) + 0.5)
       IF (LREQ.GE.LBLOCK) GOTO 50
       IPTNEW = IOLD - LREQ
       ILAST  = IPTCU + LREQ - 1
       II     = IPTNEW + LREQ -1
       IF (2*LREQ.LT.LBLOCK) THEN
        DO 110 I=1,LREQ
          A(II) = A(ILAST)
          ILAST = ILAST -1
          II    = II - 1
 110    CONTINUE
       ELSE
        DO 120 I=1,LREQ
         A(II) = A(ILAST)
         ILAST = ILAST -1
         II    = II - 1
 120    CONTINUE
       ENDIF
       PTRAST(INODE) = IPTNEW+4
       A(IOLD+2) = DBLE(IPTNEW)
  50  CONTINUE
      IRES =-1
 400  CONTINUE
       LBLOCK = IPTNEW - IPTRLU - 1
       IPTRLU = IPTNEW - 1
        IF (LRLU+LBLOCK .GE.ISPA) THEN
           LRLU   = LRLU + LBLOCK - ISPA
         ELSE
           IRES=-1
         ENDIF
 410  CONTINUE
C***************************
C***************************
C***************************
C***************************
 500  RETURN
      END
CCC
      SUBROUTINE MC51KD(N,A,INODE,
     * PTRAST,IPTBEG,LRLU,LRLUS,IPTRLU,IPTINI,LCKS)
C******************************************************
C************************************************************
      INTEGER N
      INTEGER INODE,IPTBEG,IPREV,INEXT
      INTEGER LRLU, LRLUS, IPTRLU, IPTINI
      DOUBLE PRECISION A(IPTBEG)
      INTEGER  IPT, PTRAST(N),SIZFR
      DOUBLE PRECISION ZERO
      INTEGER LCKS(20)
      INTRINSIC DBLE, INT
      PARAMETER (ZERO=0.0D0)
C***
C***
      LCKS(1) = 0
CCC
CCC
        IPT = PTRAST(INODE) - 4
        SIZFR = INT(A(IPT) + 0.5)
           LRLUS = LRLUS + SIZFR
        IPREV = INT(A(IPT+2) + 0.5)
        INEXT = INT(A(IPT+1) + 0.5)
        IF (IPREV.EQ.0) THEN
CC
CC IPT points to the top of the stack
CC we compress the top of the stack
CC
         IF (INEXT.EQ.0) THEN
          IPTINI = IPTBEG
          IPTRLU = IPTBEG
          LRLU   = LRLU + IPTBEG - IPT + 1
         ELSE
          A(INEXT+2) = ZERO
          IPTRLU     = INEXT - 1
          LRLU       = LRLU + INEXT - IPT
         ENDIF
       ELSE
         IF (INEXT.EQ.0) THEN
          IPTINI = IPREV
          A(IPREV+1) = ZERO
         ELSE
          A(IPREV+1) = DBLE(INEXT)
          A(INEXT+2) = DBLE(IPREV)
         ENDIF
       ENDIF
      RETURN
      END
CCC
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
C**
      SUBROUTINE MC51LD(N,A,LREQ,INODE,
     *   PTRAST,IPTBEG, IPTRLU,IPTINI,LCKS)
C************************************************************
      INTEGER LREQ,IPTBEG,N
      DOUBLE PRECISION A(IPTBEG)
      INTEGER IPT,INODE,PTRAST(N)
      DOUBLE PRECISION ZERO
      INTEGER LCKS(20)
      INTEGER IPTRLU, IPTINI
      INTRINSIC DBLE
      PARAMETER (ZERO=0.0D0)
C***
C***
      LCKS(1) = 0
C*********************************
        IF (IPTRLU.EQ.IPTBEG) THEN
         IPTRLU = IPTBEG - LREQ
         IPT    = IPTRLU + 1
         IPTINI = IPT
         A(IPT) = DBLE(LREQ)
         A(IPT+1) = ZERO
         A(IPT+2) = ZERO
         A(IPT+3) = DBLE(INODE)
        ELSE
         IPT = IPTRLU -LREQ+1
         A(IPT)   = DBLE(LREQ)
         A(IPT+1) = DBLE(IPTRLU +1)
         A(IPT+2) = ZERO
         A(IPT+3) = DBLE(INODE)
         A(IPTRLU +3) = DBLE(IPT)
         IPTRLU = IPT - 1
        ENDIF
         PTRAST(INODE) = IPT+4
C**********************************
      RETURN
      END
CCC
      SUBROUTINE MC51MD(N,INODE,IW,LIW,A,LA,
     *     PTRIST, PTRAST,IFINB,LKJIB,LKJIT)
C****
      INTEGER N,LA,LIW,INODE,IFINB,LKJIB
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    ALPHA, VALPIV
      INTEGER APOS, POSELT, UUPOS, PTRIST(N), PTRAST(N)
      INTEGER LKJIT
      DOUBLE PRECISION ONE
      INTEGER IOLDPS,NFRONT,NPIV,NASS,JROW2
      INTEGER NEL2,NPIVP1,KROW,LPOS,NEL
      PARAMETER (ONE=1.0D0)
      EXTERNAL DGER
        POSELT = PTRAST(INODE)
        IOLDPS = PTRIST(INODE)
        NFRONT = IW(IOLDPS)
        NPIV   = IW(IOLDPS+1)
        NASS   = IW(IOLDPS+2)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        IFINB  = 0
C****
        IF (IW(IOLDPS+3).LE.0) THEN
          IF (NASS.LT.LKJIT) THEN
           IW(IOLDPS+3) = NASS
          ELSE
           IW(IOLDPS+3) = MIN0(NASS,LKJIB)
          ENDIF
        ENDIF
        JROW2 = IW(IOLDPS+3)
        NEL2   = JROW2 - NPIVP1
        IF (NEL2.EQ.0) THEN
         IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
         ELSE
          IFINB        = 1
          IW(IOLDPS+3) = MIN0(JROW2+LKJIB,NASS)
          IW(IOLDPS+4) = NPIVP1+1
         ENDIF
        ELSE
         APOS   = POSELT + NPIV*(NFRONT + 1)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + NFRONT
         DO 541 KROW = 1,NEL2
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT
 541     CONTINUE
         LPOS   = APOS + NFRONT
         UUPOS  = APOS+1
         ALPHA = -1.0D0
         CALL DGER(NEL,NEL2,ALPHA,A(UUPOS),1,A(LPOS),NFRONT,
     *              A(LPOS+1),NFRONT)
        ENDIF
        RETURN
        END
CCC
      SUBROUTINE MC51ND(N,INODE,IW,LIW,A,LA,
     *       PTRIST,PTRAST,IFINB)
      INTEGER N,LA,LIW,INODE,IFINB
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    ALPHA,VALPIV
      INTEGER APOS, POSELT,UUPOS,
     *        PTRIST(N), PTRAST(N)
      INTEGER IOLDPS,NFRONT,NPIV,NASS,KROW
      INTEGER NEL,LPOS,ICOL,NEL2,IRWPOS
      INTEGER NPIVP1
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      EXTERNAL DAXPY
        POSELT = PTRAST(INODE)
        IOLDPS = PTRIST(INODE)
        NFRONT = IW(IOLDPS)
        NPIV   = IW(IOLDPS+1)
        NASS   = IW(IOLDPS+2)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        NEL2   = NASS - NPIVP1
        IFINB  = 0
        IF (NPIVP1.EQ.NASS) IFINB = 1
        APOS   = POSELT + NPIV*(NFRONT + 1)
CPA WE DEFINE VALPIV= -1/ (THE VALUE OF THE PIVOT) = 1/A (APOS)
        VALPIV = ONE/A(APOS)
        LPOS   = APOS + NFRONT
        DO 541 KROW = 1,NEL
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT
 541    CONTINUE
        LPOS   = APOS + NFRONT
        UUPOS  = APOS+1
        DO 440 ICOL = 1,NEL
             IRWPOS  = LPOS + 1
             ALPHA   = -A(LPOS)
             CALL DAXPY(NEL2,ALPHA,A(UUPOS),1,A(IRWPOS),1)
             LPOS    = LPOS + NFRONT
  440   CONTINUE
        RETURN
        END
CCC************* MC51OD***************
      SUBROUTINE MC51OD(N,INODE,IW,LIW,A,LA,PTRIST,PTRAST)
      INTEGER N,INODE,LA,LIW
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    ALPHA,VALPIV
      INTEGER APOS, POSELT, UUPOS, PTRIST(N), PTRAST(N)
      INTEGER IOLDPS,NFRONT,NPIV,NEL
      INTEGER LPOS,JROW,IRWPOS
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      EXTERNAL DAXPY
        POSELT = PTRAST(INODE)
        IOLDPS = PTRIST(INODE)
        NFRONT = IW(IOLDPS)
        NPIV   = IW(IOLDPS+1)
        NEL    = NFRONT - NPIV - 1
        APOS   = POSELT + (NPIV)*NFRONT + NPIV
        IF (NEL.EQ.0) GO TO 650
CPA WE DEFINE VALPIV= -1/ (THE VALUE OF THE PIVOT) = -A (APOS)
        VALPIV = ONE/A(APOS)
        LPOS   = APOS + NFRONT
        DO 340 JROW = 1,NEL
            A(LPOS) = VALPIV*A(LPOS)
            LPOS    = LPOS + NFRONT
  340   CONTINUE
        LPOS   = APOS + NFRONT
        UUPOS  = APOS+1
        DO 440 JROW = 1,NEL
             IRWPOS  = LPOS + 1
             ALPHA   = -A(LPOS)
             CALL DAXPY(NEL,ALPHA,A(UUPOS),1,A(IRWPOS),1)
             LPOS    = LPOS + NFRONT
  440   CONTINUE
  650   RETURN
        END
CCC
C*******************************************
C*******************************************
      SUBROUTINE MC51PD(A,LA,JROW1,NFRONT,
     *       NPIV,NASS,POSELT,ISPON3,ISPON4,LKJPAR)
      INTEGER LA,POSELT,LKJPAR
      DOUBLE PRECISION    A(LA)
      INTEGER JROW1, NFRONT, NPIV, NASS, ISPON3, ISPON4
      DOUBLE PRECISION ONE
      INTEGER NEL1,NEL11,LPOS2,LPOS1,LPOS
      DOUBLE PRECISION ALPHA
      INTEGER JROW2,JROW3
      PARAMETER (ONE=1.0D0)
      EXTERNAL DTRSM, DGEMM
C************************************
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        JROW2 = MIN0(JROW1+ISPON4-1,NEL1)
        IF ((NEL1.LT.ISPON3).OR.(NPIV.LT.LKJPAR)) JROW2 = NEL1
CPA JROW3 = NB OF ROWS ON WHICH WE DO THE ELIMINATION
        JROW3  = JROW2-JROW1+1
        NEL1   = JROW3
CPA LPOS  = POSITION IN A OF THE FIRST ELEMENT OF HE BLOCK
        LPOS2  = POSELT + (NASS+JROW1-1)*NFRONT
C*** WE HAVE FIRST TO COMPUTE THE L FACTORS (BELOW DIAGONAL BLOCK)
        CALL DTRSM('L','L','N','N',NPIV,JROW3,ONE,A(POSELT),NFRONT,
     *              A(LPOS2),NFRONT)
C******
C******
        LPOS   = LPOS2 + NPIV
        LPOS1  = POSELT + NPIV
        ALPHA  = -1.0D0
        CALL DGEMM('N','N',NEL11,NEL1,NPIV,ALPHA,A(LPOS1),
     *          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
  500   RETURN
        END
CCC
      SUBROUTINE MC51QD(A,LA,NFRONT,NPIV,NASS,POSELT,LKJIB)
      INTEGER LA, NFRONT, NPIV, NASS, LKJIB
      DOUBLE PRECISION    A(LA)
      INTEGER POSELT
      DOUBLE PRECISION    ALPHA
      INTEGER NEL1, NEL11, NPBEG, LPOS, LPOS1, LPOS2
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      EXTERNAL DTRSM, DGEMM
CC in this case JROW1 = 1
CC and JROW2 = NEL1
        NEL1   = NASS - NPIV
        NPBEG  = NPIV - LKJIB + 1
        NEL11  = NFRONT - NPIV
        LPOS2  = POSELT + NPIV*NFRONT + NPBEG - 1
C*** WE HAVE FIRST TO COMPUTE THE L FACTORS (BELOW DIAGONAL BLOCK)
C***
        POSELT = POSELT + (NPBEG-1)*NFRONT + NPBEG - 1
        CALL DTRSM('L','L','N','N',LKJIB,NEL1,ONE,A(POSELT),
     *               NFRONT,A(LPOS2),NFRONT)
        LPOS   = LPOS2 + LKJIB
        LPOS1  = POSELT + LKJIB
        ALPHA  = -1.0D0
        CALL DGEMM('N','N',NEL11,NEL1,LKJIB,ALPHA,A(LPOS1),
     *       NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
  500   RETURN
        END
CCC
      SUBROUTINE MC51RD(N,INODE,IW,LIW,A,LA,
     *  JROW1, PTRIST,PTRAST,LKJIB,ISPON2)
      INTEGER N,LA,LIW
      INTEGER IOLDPS, NFRONT, NASS, NPIV, NPBEG,JROW2
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    ALPHA
      INTEGER INODE,POSELT,JROW1,
     *        PTRIST(N), PTRAST(N), LKJIB, ISPON2
      INTEGER NEL1, NEL11, JROW3, LPOS, LPOS1, LPOS2
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      EXTERNAL DTRSM, DGEMM
C***************************
C********************************
        POSELT = PTRAST(INODE)
        IOLDPS = PTRIST(INODE)
        NFRONT = IW(IOLDPS)
        NASS   = IW(IOLDPS+2)
        NPIV   = IW(IOLDPS+1)
        NEL1   = NASS - NPIV
        NPBEG  = NPIV - LKJIB + 1
        NEL11  = NFRONT - NPIV
        JROW2 = MIN0(JROW1+ISPON2-1,NEL1)
CPA JROW3 = NB OF ROWS ON WHICH WE DO THE ELIMINATION
        JROW3  = JROW2-JROW1+1
        NEL1   = JROW3
CPA LPOS  = POSITION IN A OF THE FIRST ELEMENT OF HE BLOCK
        LPOS2  = POSELT + (NPIV+JROW1-1)*NFRONT + NPBEG - 1
C*** WE HAVE FIRST TO COMPUTE THE L FACTORS (BELOW DIAGONAL BLOCK)
C***
        POSELT = POSELT + (NPBEG-1)*NFRONT + NPBEG - 1
        CALL DTRSM('L','L','N','N',LKJIB,JROW3,ONE,A(POSELT),
     *               NFRONT,A(LPOS2),NFRONT)
        LPOS   = LPOS2 + LKJIB
        LPOS1  = POSELT + LKJIB
        ALPHA  = -1.0D0
        CALL DGEMM('N','N',NEL11,NEL1,LKJIB,ALPHA,A(LPOS1),
     *          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
  500   RETURN
        END
CCC
C**** ROUTINE MC51SD **************
C*************************************
      SUBROUTINE MC51SD(N,INODE,IW,LIW,A,LA,
     *    PTRIST,PTRAST,LKJIB,LKJIT)
      INTEGER N,LA,LIW
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW), LKJIB, INODE
      DOUBLE PRECISION    ALPHA, ONE
      INTEGER POSELT, PTRIST(N), PTRAST(N)
      INTEGER IOLDPS, NFRONT, NPIV, NASS, JROW2, NPBEG
      INTEGER NONEL, LKABS, LKJIW, NEL1, NEL11
      INTEGER LBP, IPOS, IPOSEL, KPOS, LPOS2
      INTEGER LPOS1,LPOS,LBPT,I1,K1,II,ISWOP,LBP1
      INTEGER LKJIT
      PARAMETER (ONE=1.0D0)
      EXTERNAL DTRSM, DGEMM, DSWAP
        POSELT = PTRAST(INODE)
        IOLDPS = PTRIST(INODE)
        IPOSEL = POSELT
        NFRONT = IW(IOLDPS)
        NPIV   = IW(IOLDPS+1)
        NASS   = IW(IOLDPS+2)
        JROW2  = IABS(IW(IOLDPS+3))
        NPBEG  = IW(IOLDPS+4)
C*******
C******
        NONEL         = JROW2 - NPIV + 1
        IF (LKJIB.GT.0) NONEL = -NONEL
        LKABS         = IABS(LKJIB)
        IF ((NASS-NPIV).GE.LKJIT) THEN
          LKJIB       = LKJIB - NONEL
          LKABS       = IABS(LKJIB)
          IW(IOLDPS+3)= MIN0(NPIV+LKABS,NASS)
        ELSE
          IW(IOLDPS+3) = NASS
        ENDIF
        IW(IOLDPS+4) = NPIV + 1
        NEL1   = NASS - JROW2
        LKJIW  = NPIV - NPBEG + 1
        NEL11  = NFRONT - NPIV
        IF ((NEL1.EQ.0).OR.(LKJIW.EQ.0)) GO TO 500
        LPOS2  = POSELT + JROW2*NFRONT + NPBEG - 1
C*** WE HAVE FIRST TO COMPUTE THE L FACTORS (BELOW DIAGONAL BLOCK)
C***
         POSELT = POSELT + (NPBEG-1)*NFRONT + NPBEG - 1
         CALL DTRSM('L','L','N','N',LKJIW,NEL1,ONE,A(POSELT),NFRONT,
     *               A(LPOS2),NFRONT)
C******
C******
        LPOS   = LPOS2 + LKJIW
        LPOS1  = POSELT + LKJIW
        ALPHA  = -1.0D0
        CALL DGEMM('N','N',NEL11,NEL1,LKJIW,ALPHA,A(LPOS1),
     *          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
  500   LBP = JROW2 - NPIV
C***
C***
         LBPT  = LBP + LBP
        IF  ((NEL1.GE.LBPT).AND.(NEL1.GE.LKABS)) THEN
         I1 = IOLDPS + 5 + NPIV
         K1 = IOLDPS + 5 + NASS - LBP
         DO 10 II=1,LBP
          ISWOP  = IW(I1)
          IW(I1) = IW(K1)
          IW(K1) = ISWOP
          I1     = I1 +1
          K1     = K1 + 1
  10     CONTINUE
         IPOS = IPOSEL + NPIV*NFRONT
         KPOS = IPOSEL + (NASS-LBP)*NFRONT
         LBP1 = LBP * NFRONT
         CALL DSWAP(LBP1,A(IPOS),1,A(KPOS),1)
        ENDIF
  700   RETURN
        END
CCC
      SUBROUTINE MC51TD(A,LA,NPIVB,NFRONT,
     *                             NPIV,NASS,POSELT)
      INTEGER NPIVB,NASS,LA
      DOUBLE PRECISION    A(LA)
      DOUBLE PRECISION    ALPHA, ONE
      INTEGER APOS, POSELT
      INTEGER NFRONT, NPIV, NASSL
      INTEGER LPOS, LPOS1, LPOS2, NEL1, NEL11, NPIVE
      PARAMETER (ONE=1.0D0)
      EXTERNAL DTRSM, DGEMM
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        NPIVE  = NPIV - NPIVB
        NASSL  = NASS - NPIVB
        APOS   = POSELT + NPIVB*NFRONT + NPIVB
        LPOS2  = APOS + NASSL
C*** WE HAVE FIRST TO COMPUTE THE U FACTORS
        CALL DTRSM('R','U','N','U',NEL1,NPIVE,ONE,A(APOS),NFRONT,
     *              A(LPOS2),NFRONT)
C******
C******
        LPOS   = LPOS2 + NFRONT*NPIVE
        LPOS1  = APOS + NFRONT*NPIVE
        ALPHA  = -1.0D0
        CALL DGEMM('N','N',NEL1,NEL11,NPIVE,ALPHA,A(LPOS2),
     *          NFRONT,A(LPOS1),NFRONT,ONE,A(LPOS),NFRONT)
  500   RETURN
        END
CCC
      SUBROUTINE MC51UD(A, LA, IPTA, IPT,
     *       IPTEND, IPTBEG, LCKS)
C************************************************************
C************************************************************
      INTEGER NBUD
      PARAMETER (NBUD=29)
      INTEGER LA,IPTA(NBUD),IPTEND, IPTBEG
      INTEGER LCKS(20)
      DOUBLE PRECISION A(LA)
C****************
C****************
      INTEGER IPT,IPTOLD,LSIZ,IFW,LSIZP1,IPT1,IPTBUD
      INTEGER IBW,IHEAD
      DOUBLE PRECISION ZERO
      INTRINSIC DBLE, INT
      PARAMETER (ZERO=0.0D0)
C***
C***
      LCKS(1) = 0
C****************************
C****************************
C****************************
      IF (IPT.GE.IPTEND) THEN
       IPTOLD = IPTA(NBUD)
       IPTA(NBUD) = IPT
       A(IPT+1) = ZERO
       A(IPT+2) = DBLE(IPTOLD)
       IF (IPTOLD.GT.0) A(IPTOLD+1) = DBLE(IPT)
      ELSE
  45  LSIZ = INT(A(IPT)+0.5)
      IFW    = -1
      LSIZP1 = LSIZ + 1
      IPT1   = IPT - IPTBEG
      IF (((IPT1+2**LSIZP1-1)/2**LSIZP1)*2**LSIZP1.EQ.
     *      IPT1+2**LSIZP1-1) IFW = 1
      IPTBUD = IPT + IFW*2**LSIZ
      IF (IPTBUD.LT.IPTEND .AND. IPTBUD.GT.IPTBEG .AND.
     *    INT(A(IPTBUD+1)+0.5).GE.0 .AND.
     *    INT(A(IPTBUD)+0.5).EQ.LSIZ) THEN
        IFW = INT(A(IPTBUD+2)+0.5)
        IBW = INT(A(IPTBUD+1)+0.5)
        IF (IBW.EQ.0) IPTA(LSIZ) = IFW
        IF (IBW.GT.0) A(IBW+2)   = DBLE(IFW)
        IF (IFW.GT.0) A(IFW+1)   = DBLE(IBW)
        IPT    = MIN0(IPT,IPTBUD)
        A(IPT) = DBLE(LSIZ + 1)
        GO TO 45
      ELSE
        IHEAD = IPTA(LSIZ)
        IF (IHEAD.GT.0) A(IHEAD+1) = DBLE(IPT)
        A(IPT+1)   = ZERO
        A(IPT+2)   = DBLE(IHEAD)
        IPTA(LSIZ) = IPT
      ENDIF
      ENDIF
C****************************
C****************************
      RETURN
      END
CCC
      SUBROUTINE MC51VD(N,NZ,VAL,IRN,ICN,RNOR,
     *      COLSCA,ROWSCA,MPRINT)
C******************************
C******************************
      INTEGER   N, NZ
      DOUBLE PRECISION      VAL(NZ),RNOR(N),COLSCA(N)
      DOUBLE PRECISION      ROWSCA(N)
      INTEGER   IRN(NZ),ICN(NZ)
      DOUBLE PRECISION      VDIAG
      INTEGER   MPRINT,I,J,K
      INTRINSIC SQRT
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
CC
CC
      DO 10 I=1,N
       RNOR(I)   = 1.
  10  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          IF ((I.GT.N).OR.(I.LE.0)) GOTO 100
          J = ICN(K)
          IF (I.EQ.J) THEN
            VDIAG = ABS(VAL(K))
            IF (VDIAG.GT.ZERO) THEN
              RNOR(J) = 1./(SQRT(VDIAG))
            ENDIF
          ENDIF
 100   CONTINUE
       DO 110 I=1,N
        COLSCA(I) = RNOR(I)
        ROWSCA(I) = RNOR(I)
 110   CONTINUE
      IF (MPRINT.GE.0) WRITE(MPRINT,'(A)') 'End of diagonal scaling'
      RETURN
      END
CCC
      SUBROUTINE MC51WD(N, NZ, VAL, ROWIND, COLIND,
     *                  RNOR, CNOR, WNOR, MPRINT, MP, NSCA)
C*****************
C******************************
C*************
C***************
C***************
C******************************
      INTEGER N, NZ
      DOUBLE PRECISION    VAL(NZ),RNOR(N),CNOR(N),
     *        WNOR(5*N)
      INTEGER COLIND(NZ),ROWIND(NZ)
      INTEGER J,I,K
      INTEGER MPRINT,MP,NSCA
      INTEGER IFAIL9
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      EXTERNAL MC29AD
CC
CC
      DO 15 I=1,N
       RNOR(I)   = ZERO
       CNOR(I)   = ZERO
  15  CONTINUE
      CALL MC29AD(N,N,NZ,VAL,ROWIND,COLIND,
     *   RNOR,CNOR,WNOR, MP,IFAIL9)
CCVD$ NODEPCHK
CCVD$ VECTOR
CCVD$ CONCUR
      DO 30 I=1,N
       CNOR(I) = EXP(CNOR(I))
       RNOR(I) = EXP(RNOR(I))
  30  CONTINUE
C*********************************************
      IF ((NSCA.EQ.5).OR.(NSCA.EQ.6)) THEN
        DO 100 K=1,NZ
          I   = ROWIND(K)
          J   = COLIND(K)
          IF (MIN(I,J).LT.1 .OR. I.GT.N .OR. J.GT.N) GOTO 100
          VAL(K) = VAL(K) * CNOR(J) * RNOR(I)
 100    CONTINUE
      ENDIF
      IF (MPRINT.GE.0)
     *   WRITE(MPRINT,'(A)') 'End of scaling using MC29'
      RETURN
      END
CCC
      SUBROUTINE MC51XD(N,NZ,IRN,ICN,VAL,
     *    RNOR,CNOR,COLSCA,ROWSCA,MPRINT)
C*******************************************************
C*******************************************************
      INTEGER N, NZ
      DOUBLE PRECISION    VAL(NZ),RNOR(N),CNOR(N)
      DOUBLE PRECISION    COLSCA(N),ROWSCA(N)
      DOUBLE PRECISION    CMIN,CMAX,RMIN,ARNOR,ACNOR
      INTEGER IRN(NZ), ICN(NZ)
      DOUBLE PRECISION    VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
CC
CC
      DO 50 J=1,N
       CNOR(J)   = ZERO
       RNOR(J)   = ZERO
  50  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I.LE.0).OR.(I.GT.N).OR.
     *        (J.LE.0).OR.(J.GT.N)) GOTO 100
            VDIAG = ABS(VAL(K))
            IF (VDIAG.GT.CNOR(J)) THEN
              CNOR(J) =     VDIAG
            ENDIF
            IF (VDIAG.GT.RNOR(I)) THEN
              RNOR(I) =     VDIAG
            ENDIF
 100   CONTINUE
      IF (MPRINT.GE.0) THEN
       CMIN = CNOR(1)
       CMAX = CNOR(1)
       RMIN = RNOR(1)
C** STATISTICS OF MATRIX PRIOR ROW AND COL. SCALING
       DO 111 I=1,N
        ARNOR = ABS(RNOR(I))
        ACNOR = ABS(CNOR(I))
        IF (ACNOR.GT.CMAX) CMAX=ACNOR
        IF (ACNOR.LT.CMIN) CMIN=ACNOR
        IF (ARNOR.LT.RMIN) RMIN=ARNOR
 111   CONTINUE
       WRITE(MPRINT,'(/A,A)') '**** Statistics of matrix prior to ',
     *                 'row and column scaling'
       WRITE(MPRINT,'(A,1PD12.4)') 'Maximum max-norm of columns:',CMAX
       WRITE(MPRINT,'(A,1PD12.4)') 'Minimum max-norm of columns:',CMIN
       WRITE(MPRINT,'(A,1PD12.4)') 'Minimum max-norm of rows   :',RMIN
      ENDIF
C** END OF STATISTICS**
      DO 120 J=1,N
       IF (CNOR(J).LE.ZERO) THEN
         CNOR(J)   = 1.
       ELSE
         CNOR(J)   = 1./CNOR(J)
       ENDIF
 120  CONTINUE
      DO 130 J=1,N
       IF (RNOR(J).LE.ZERO) THEN
         RNOR(J)   = 1.
       ELSE
         RNOR(J)   = 1./RNOR(J)
       ENDIF
 130  CONTINUE
       DO 110 I=1,N
        ROWSCA(I) = ROWSCA(I)* RNOR(I)
        COLSCA(I) = COLSCA(I) * CNOR(I)
 110   CONTINUE
      IF (MPRINT.GE.0)
     *  WRITE(MPRINT,'(A)') 'End of scaling by max in row and column'
      RETURN
      END
CCC
      SUBROUTINE MC51YD(N,NZ,VAL,IRN,ICN,
     *       CNOR,COLSCA,MPRINT)
C*******************************************************
C*******************************************************
C******************************
      INTEGER N,NZ
      DOUBLE PRECISION VAL(NZ),CNOR(N),COLSCA(N)
      INTEGER IRN(NZ), ICN(NZ)
      DOUBLE PRECISION VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DO 10 J=1,N
       CNOR(J)   = ZERO
  10  CONTINUE
      DO 100 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.LE.0).OR.(I.GT.N).OR.
     *      (J.LE.0).OR.(J.GT.N)) GOTO 100
        VDIAG = ABS(VAL(K))
        IF (VDIAG.GT.CNOR(J)) THEN
           CNOR(J) =     VDIAG
        ENDIF
 100  CONTINUE
      DO 110 J=1,N
       IF (CNOR(J).LE.ZERO) THEN
         CNOR(J)   = 1.
       ELSE
         CNOR(J)   = 1./CNOR(J)
       ENDIF
 110  CONTINUE
       DO 215 I=1,N
        COLSCA(I) = COLSCA(I) * CNOR(I)
 215   CONTINUE
      IF (MPRINT.GE.0) WRITE(MPRINT,'(A)') 'End of column scaling'
      RETURN
      END
CCC
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
C***********************************************************
C**********************************************************
      SUBROUTINE MC51ZD(A, LA, IPTA, LREQ, IPT, LSIZ, IRES,
     *                  LENA,NBACTI,LCKS)
C*************************************************************
C************************************************************
      INTEGER NBUD
      PARAMETER (NBUD=29)
      INTEGER LA, IPTA(NBUD), IRES, NBACTI
      INTEGER LENA,LSIZ,LREQ,IPT
      INTEGER LCKS(20)
      INTEGER ISIZ,IFW,KSIZ,IPTBUD
      LOGICAL IOK
      DOUBLE PRECISION A(LA)
      DOUBLE PRECISION ZERO
      INTRINSIC DBLE, INT
      PARAMETER (ZERO=0.0D0)
C***
C***
      LCKS(1) = 0
      LSIZ = 0
      LSIZ = LOG(DBLE(LREQ))/LOG(2.0) + 1
C***************************************
C**************************************
      IOK = .FALSE.
      IF (LREQ.LE.LENA) THEN
        IPT = IPTA(NBUD)
        IF (IPT.NE.0) THEN
          IOK  = .TRUE.
          IRES = 1
          IPTA(NBUD) = INT(A(IPT+2)+0.5)
          IF (IPTA(NBUD).NE.0) A(IPTA(NBUD)+1) = 0
        ENDIF
      ENDIF
      IF (IOK) GOTO 400
      DO 300 ISIZ=LSIZ,NBUD-1
        IF (IPTA(ISIZ).NE.0) THEN
         IPT  = IPTA(ISIZ)
         IRES = 2
         IFW  = INT(A(IPT+2)+0.5)
         IPTA(ISIZ) = IFW
         IF (IFW.NE.0) A(IFW+1) = ZERO
         A(IPT+1) = DBLE(-2)
         A(IPT)   = DBLE(LSIZ)
         IF (ISIZ.EQ.LSIZ) GO TO 400
          DO 28 KSIZ = ISIZ-1,LSIZ,-1
            IPTBUD     = IPT + 2**KSIZ
            IPTA(KSIZ) = IPTBUD
            A(IPTBUD)  = DBLE(KSIZ)
            A(IPTBUD+1) = ZERO
            A(IPTBUD+2) = ZERO
 28       CONTINUE
          GOTO 400
        ENDIF
 300  CONTINUE
C=================================================================
C=================================================================
      IF (NBACTI.EQ.0) THEN
         IRES = -2
      ELSE
         IRES = -1
      ENDIF
C**************************************
 400    CONTINUE
C**************************************
      RETURN
      END


* COPYRIGHT (c) 1993 AEA Technology and
* Council for the Central Laboratory of the Research Councils
      SUBROUTINE MC29AD(M,N,NE,A,IRN,ICN,R,C,W,LP,IFAIL)
      INTEGER M,N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IRN(NE),ICN(NE)
      DOUBLE PRECISION R(M),C(N),W(M*2+N*3)
      INTEGER LP,IFAIL
      INTRINSIC LOG,ABS,MIN
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      DOUBLE PRECISION ONE,SMIN,ZERO
      PARAMETER (ONE=1D0,SMIN=0.1,ZERO=0D0)
      INTEGER I,I1,I2,I3,I4,I5,ITER,J,K
      DOUBLE PRECISION E,E1,EM,Q,Q1,QM,S,S1,SM,U,V
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
     +    ' **** Error return from MC29AD ****',' IFAIL =',IFAIL
      END


