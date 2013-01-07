C *******************************************************************
C COPYRIGHT (c) 1995  P.R. Amestoy and Council for the Central
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
C                     Laboratory of the Research Councils
C######DATE 30 November 1995
C  April 2001: threadsafe version of MA41
C  January 2002: features deleted in Fortran 95 removed.
      SUBROUTINE MA41I(CNTL, ICNTL, KEEP)
C****************************************************************
      REAL    CNTL(10)
      INTEGER ICNTL(20), KEEP(50)
C===========================================
C===========================================
C=========================
C========================
C-----
C-----
C----------------------------------------------------
C----------------------------------------------------
C===========================
C===========================
      INTEGER I
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
C===================================
C==================================
      CNTL(1)   = 0.01
      CNTL(2)   = 0.1E-4
      ICNTL(1)  = 6
      ICNTL(2)  = -1
      ICNTL(3)  = -1
      ICNTL(4)  = 2
      ICNTL(5)  = 1
      ICNTL(6)  = 0
      ICNTL(7)  = 0
      ICNTL(8)  = 0
      ICNTL(9)  = 1
      ICNTL(10)  = 0
      ICNTL(11)  = 0
      DO 110 I=12,20
        ICNTL(I) = 0
  110 CONTINUE
      DO 111 I=3,10
        CNTL(I) = ZERO
  111 CONTINUE
C===================================
C===================================
      KEEP(2) = -1
      KEEP(12)  = 10
      KEEP(24) = 0
      KEEP(1) = 8
      KEEP(3)  = 96
      KEEP(4)  = 32
      KEEP(5)  = 16
      KEEP(6) = 32
      KEEP(7) = 32
      KEEP(8) = 16
      KEEP(9) = 96
      KEEP(10) = 32
      RETURN
      END
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA41A (JOB, N, NE, IRN, JCN, ASPK,
     *   RHS, COLSCA, ROWSCA, KEEP, IS, MAXIS,
     *   S, MAXS, CNTL, ICNTL, INFO, RINFO)
C============================================
C============================================
C********
      INTEGER MAXS, JOB, N,NE, MAXIS
      INTEGER IS(MAXIS), IRN(NE), JCN(NE)
      INTEGER INFO(20), KEEP(50)
      REAL    ASPK(NE), RHS(N),
     *        COLSCA(*), ROWSCA(*)
      REAL    RINFO(20), S(MAXS)
      REAL    CNTL(10)
      INTEGER ICNTL(20)
      INTEGER MP
      LOGICAL LANAL, LFACTO, LSOLVE, PROK
      EXTERNAL MA41E
       INFO(1) = 0
       INFO(2) = 0
       MP      = ICNTL(2)
       PROK    = ((MP.GE.0).AND.(ICNTL(4).GE.3))
      IF ((N.LE.0).OR.((N+N+N)/3.NE.N)) GO TO 300
      IF (NE.LE.0)                      GO TO 350
      IF ((JOB.LT.1).OR.(JOB.GT.6))      GO TO 400
      LANAL  = .FALSE.
      LFACTO = .FALSE.
      LSOLVE = .FALSE.
      IF ((JOB.EQ.1).OR.(JOB.EQ.4).OR.
     *    (JOB.EQ.6))               LANAL  = .TRUE.
      IF ((JOB.EQ.2).OR.(JOB.EQ.4).OR.
     *    (JOB.EQ.5).OR.(JOB.EQ.6)) LFACTO = .TRUE.
      IF ((JOB.EQ.3).OR.(JOB.EQ.5).OR.
     *    (JOB.EQ.6))               LSOLVE = .TRUE.
C********
C********
        IF (PROK) WRITE(MP,99996) JOB,N,NE,MAXIS,MAXS
        CALL MA41E(N,NE,IRN,JCN,ASPK,IS,MAXIS,
     *   S,MAXS,RHS,COLSCA,ROWSCA,
     *   LANAL,LFACTO,LSOLVE,CNTL,ICNTL,
     *   INFO,RINFO,KEEP)
        IF (INFO(1).LT.0) GO TO 499
        GO TO 500
C=================
C=================
  300  INFO(1) = -1
       INFO(2) = N
       GO TO 499
  350  INFO(1) = -2
       INFO(2) = NE
       GO TO 499
  400  INFO(1) = -3
       INFO(2) = JOB
  499  PROK  = ((ICNTL(1).GE.0).AND.(ICNTL(4).GE.1))
       IF (PROK) WRITE (ICNTL(1),99995) INFO(1)
       IF (PROK) WRITE (ICNTL(1),99994) INFO(2)
500   RETURN
99996 FORMAT (/'Entering driver (MA41A) with ...........'/
     * '   JOB       N         NE       MAXIS          MAXS'/,
     * 4X, I2, I8, I11, I12, I14)
99995 FORMAT (/'** Error return ** from MA41A   INFO(1)=', I3)
99994 FORMAT ('INFO(2)=', I10)
      END
CCC
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA41B(N, NZ, NSTEPS, IRN, ICN, ASPK,
     * A, LA, IW, LIW, IKEEP,NFSIZ, FILS, FRERE, PTRAR, PTRIST,
     * PTLUST, IW1, NPROCS, IW2, IW3, LIW3, LSCAL,COLSCA, ROWSCA,
     * IPTA, CNTL, ICNTL, INFO, RINFO, KEEP)
      INTEGER NBUD
      PARAMETER (NBUD=29)
      INTEGER N,NZ,NSTEPS,LA,LIW,LIW3,NPROCS
      REAL A(LA),ASPK(NZ)
      REAL ROWSCA(*),COLSCA(*)
      REAL RINFO(20)
      LOGICAL LSCAL
      REAL CNTL(10)
      INTEGER ICNTL(20), IPTA(NBUD)
      INTEGER   INFO(20), KEEP(50)
      INTEGER   IRN(NZ), ICN(NZ), IW(LIW), IKEEP(N,3), FILS(N),
     *          FRERE(N), NFSIZ(N)
      INTEGER   PTRAR(N,4), PTRIST(N), PTLUST(NSTEPS)
      INTEGER   IW1(4*N), IW2(N,NPROCS), IW3(LIW3)
      REAL UULOC
      INTEGER MPRINT, LP, MP, LDIAG
      INTEGER NSTK,PTRAST,TLKJ, PERMW,
     *        KZ,K,I,KBLK,IPOS,IBLK,NCOLS,LPOOL
      INTEGER NPIV,IAPOS,JJ,J1,J2,LEN,IPIV,IROW,
     *        NIRBDU
      LOGICAL PROK
      REAL ZERO
      PARAMETER (ZERO = 0.0E0)
      EXTERNAL MA41H, MA41Z
      KEEP(17)= 0
      KEEP(22)= 0
      INFO(11)= 0
      INFO(12)= 0
      INFO(13)= 0
      INFO(14)  = 0
      RINFO(2) = ZERO
      RINFO(3) = ZERO
      MPRINT = ICNTL(3)
      LP     = ICNTL(1)
      MP     = ICNTL(2)
      LDIAG  = ICNTL(4)
      PROK   = (ICNTL(3).GE.0)
      IF ( (KEEP(26)+KEEP(14).GT.KEEP(15)) .OR.
     *   (KEEP(25).GT.KEEP(16)) .OR.
     *   (KEEP(16).LE.0) ) GO TO 50
      IF (KEEP(12).LT.0) KEEP(12)=10
      UULOC = CNTL(1)
      IF (UULOC.GT.1.0)   UULOC=1.0
      IF (UULOC.LT.ZERO)  UULOC=ZERO
      PTRAST = 1
      NSTK   = PTRAST + N
      TLKJ   = NSTK   + N
      PERMW  = TLKJ   + N
C*****************************
      IF (PROK) THEN
        WRITE (MPRINT,99999) N, NZ,
     *           KEEP(25), KEEP(26), LA, LIW, KEEP(16),
     *           KEEP(28), ICNTL(8), UULOC
         WRITE(MPRINT,99979) NPROCS
      ENDIF
      IF (LDIAG.GT.2 .AND. MP.GE.0) THEN
        IF (.NOT.PROK)
     *  WRITE (MP,99999) N, NZ,
     *           KEEP(25), KEEP(26), LA, LIW, KEEP(16),
     *           KEEP(28), ICNTL(8), UULOC
        KZ = MIN0(10,NZ)
        IF (LDIAG.EQ.4) KZ = NZ
        IF (NZ.GT.0) WRITE (MP,99998)
     *                 (ASPK(K),IRN(K),ICN(K),K=1,KZ)
        K = MIN0(10,N)
        IF (LDIAG.EQ.4) K = N
        IF (K.GT.0) WRITE (MP,99997) (IKEEP(I,1),I=1,K)
        IF (K.GT.0) WRITE (MP,99996) (IKEEP(I,2),I=1,K)
        IF (K.GT.0) WRITE (MP,99995) (IKEEP(I,3),I=1,K)
        IF (K.GT.0) WRITE (MP,99966) (PTRAR(I,1),I=1,K)
        IF (K.GT.0) WRITE (MP,99965) (PTRAR(I,2),I=1,K)
        IF (K.GT.0) WRITE (MP,99964) (PTRAR(I,3),I=1,K)
        IF (K.GT.0) WRITE (MP,99963) (PTRAR(I,4),I=1,K)
        IF (K.GT.0) WRITE (MP,99967) (NFSIZ(I),I=1,K)
      ENDIF
C**********************************
C**********************************
      CALL MA41H(N,NZ,KEEP(13),ASPK,A, LA,IRN, ICN, IW, LIW,
     *  IKEEP, IW1(TLKJ), PTRAR, PTRAR(1,4),
     *  PTRAR(1,3), LSCAL, COLSCA, ROWSCA)
C**********************************
C**********************************
      LPOOL   = KEEP(28)+1
      NIRBDU  = LIW - KEEP(14)
      CALL MA41Z(N, KEEP(13), NSTEPS, A, LA, IW,
     * LIW, IKEEP(1,2), IKEEP(1,3), IW1(NSTK),
     * INFO(11),INFO(1),INFO(2), NFSIZ,FILS,FRERE,
     * PTRAR(1,4),PTRAR(1,3),IW1(PTRAST),PTRIST,PTLUST,
     * IW1(TLKJ),IW1(PERMW),NPROCS,IW2,NIRBDU,
     * IPTA, IW3(1), LPOOL ,IW3(LPOOL+1),
     * UULOC,ICNTL, INFO, RINFO, KEEP)
      IF (INFO(1).LT.0) GO TO 100
      GO TO 100
C************************
C************************
   50 INFO(1) = -3
C****************
C****************
 100  IF (PROK) THEN
        WRITE (MPRINT,99980) INFO(1), INFO(2),
     *       KEEP(28), INFO(9), INFO(10), INFO(11), INFO(12),
     *       INFO(13), INFO(14), RINFO(2), RINFO(3)
      ENDIF
      IF (INFO(1).LT.0) GO TO 500
      IF (LDIAG.LE.0 .OR. MP.LT.0) GO TO 500
      IF (LDIAG.LE.2) GO TO 500
      IF (.NOT.PROK) THEN
        WRITE (MP,99980) INFO(1), INFO(2),
     *       KEEP(28), INFO(9), INFO(10), INFO(11), INFO(12),
     *       INFO(13), INFO(14),  RINFO(2), RINFO(3)
      ENDIF
      KBLK = KEEP(28)
      IF (LDIAG.EQ.3) KBLK = 1
      DO 120 IBLK=1,KBLK
        IPOS  = PTLUST(IBLK) + 2
        NCOLS = IW(IPOS)
        NPIV = IW(IPOS+1)
        IAPOS = IW(IPOS+2)
        J1 = IPOS + 3
        WRITE (MP,99985) IBLK, NPIV, NCOLS, IAPOS
        J2 = J1 + NCOLS - 1
        WRITE (MP,99984) (IW(JJ),JJ=J1,J2)
        J1 = J2 + 1
        J2 = J1 + NCOLS - 1
        WRITE (MP,99983) (IW(JJ),JJ=J1,J2)
        WRITE (MP,99982)
        LEN  = NCOLS
        IROW = 1
        IF (NPIV.EQ.0) GO TO 120
        DO 110 IPIV=1,NPIV
          J1 = IAPOS
          J2 = IAPOS + LEN - 1
          WRITE (MP,99981) (A(JJ),JJ=J1,J2)
          LEN = LEN - 1
          IAPOS = J2 + IROW
          IROW  = IROW + 1
          J1 = IAPOS
          J2 = IAPOS + (LEN - 1)*NCOLS
          IF (J2.GE.J1) WRITE(MP,99981)(A(JJ),JJ=J1,J2,NCOLS)
          IAPOS = IAPOS + 1
  110   CONTINUE
  120 CONTINUE
  500 RETURN
99999 FORMAT (/'Entering factorization step (MA41B) with ...'/
     1        'Order of input matrix (N)                     =',I12/
     2        'Entries in input matrix (NE)                  =',I12/
     3        'Entries   in factors (forecast)               =',I12/
     3        'Integers  in factors (forecast)               =',I12/
     4        'Size of real working space (LA)               =',I12/
     5        'Size of integer working space (LIW)           =',I12/
     6        'Estimated size  OF LU area                    =',I12/
     7        'Number nodes in elimination tree              =',I12/
     8        'Scaling parameter (ICNTL(8))                  =',I12/
     9        'Value of threshold parameter                  =',F12.5)
99998 FORMAT ('Matrix entries:   ASPK() IRN() ICN()'/
     * (1PE14.4, 2I8, 5X, 1PE14.4, 2I8))
99997 FORMAT ('IKEEP(.,1)=', 10I6/(12X, 10I6))
99996 FORMAT ('IKEEP(.,2)=', 10I6/(12X, 10I6))
99995 FORMAT ('IKEEP(.,3)=', 10I6/(12X, 10I6))
99967 FORMAT ('NFSIZ(.)  =', 10I6/(12X, 10I6))
99966 FORMAT ('PTRAR(.,1)=', 10I6/(12X, 10I6))
99965 FORMAT ('PTRAR(.,2)=', 10I6/(12X, 10I6))
99964 FORMAT ('PTRAR(.,3)=', 10I6/(12X, 10I6))
99963 FORMAT ('PTRAR(.,4)=', 10I6/(12X, 10I6))
99994 FORMAT (/'*** Warning message from routine MA41B **',
     *         '   INFO(1) =',I2/5X, 'Matrix is singular, rank=', I5)
99992 FORMAT ('Value of N out-of-range ... =', I10)
99991 FORMAT ('Value of NZ out-of-range .. =', I10)
99985 FORMAT (/'Block pivot =', I7, '   NPIV =', I7, '  NFRONT =', I7/
     *        'IAPOS =', I7)
99984 FORMAT ('Row indices    =', 9I6/(17X, 9I6))
99983 FORMAT ('Column indices =', 9I6/(17X, 9I6))
99982 FORMAT ('Real entries',
     *        ' .. each row and column starts on a new line')
99981 FORMAT (1P,5E14.6)
99980 FORMAT (/'Leaving factorization phase (MA41B) with ...'/
     1      'INFO (1)                                      =',I12/
     2      ' --- (2)                                      =',I12/
     5      '          Number of nodes in the tree         =',I12/
     6      'INFO (9)  Real space for factors              =',I12/
     7      ' --- (10) Integer space for factors           =',I12/
     8      ' --- (11) Maximum size of frontal matrices    =',I12/
     9      ' --- (12) Number of off diagonal pivots       =',I12/
     1      ' --- (13) Number of delayed pivots            =',I12/
     2      ' --- (14) Number of memory compresses         =',I12/
     3  'RINFO(2)  Operations during node assembly     =  ',1PE10.3/
     4  '-----(3)  Operations during node elimination  =  ',1PE10.3)
99979 FORMAT( 'Number of parallel tasks                      =',I12)
      END
CC
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA41C(N, A, LA, IW, LIW, W, MAXFRT,
     * RHS, PTLUST, NSTEPS, W2, MTYPE, ICNTL)
C************************
      INTEGER LA,MAXFRT,N,LIW,NSTEPS,MTYPE
      INTEGER PTLUST(NSTEPS), ICNTL(20)
      INTEGER IW(LIW)
      REAL    A(LA), W(MAXFRT), RHS(N),
     *        W2(N)
      INTEGER MP, LDIAG
      INTEGER IPOS,IBLK,NCOLS,NPIV,IAPOS,JJ,J1,J2
      INTEGER KBLK, LEN,K,I,KL,IPIV,IROW
      EXTERNAL MA41R, MA41S, MA41T, MA41U
      MP      = ICNTL(2)
      LDIAG   = ICNTL(4)
      IF (LDIAG.GT.2 .AND. MP.GE.0) THEN
C***************************************************
C***************************************************
        WRITE (MP,99999) N, LA, LIW, MAXFRT, NSTEPS
        KBLK = NSTEPS
        IF (LDIAG.EQ.3) KBLK = 1
        DO 20 IBLK=1,KBLK
          IPOS  = PTLUST(IBLK)+2
          NCOLS = IW(IPOS)
          NPIV = IW(IPOS+1)
          IAPOS = IW(IPOS+2)
          J1 = IPOS + 3
          WRITE (MP,99998) IBLK, NPIV, NCOLS, IAPOS
          J2 = J1 + NCOLS - 1
          WRITE (MP,99997) (IW(JJ),JJ=J1,J2)
          J1 = J2 + 1
          J2 = J1 + NCOLS - 1
          WRITE (MP,99996) (IW(JJ),JJ=J1,J2)
          WRITE (MP,99995)
          LEN = NCOLS
          IROW = 1
          IF (NPIV.EQ.0) GO TO 20
          DO 10 IPIV=1,NPIV
            J1 = IAPOS
            J2 = IAPOS + LEN - 1
            WRITE (MP,99994) (A(JJ),JJ=J1,J2)
            LEN = LEN - 1
            IAPOS = J2 + IROW
            IROW  = IROW + 1
            J1 = IAPOS
            J2 = IAPOS + (LEN - 1)*NCOLS
            IF (J2.GE.J1) WRITE(MP,99994)(A(JJ),JJ=J1,J2,NCOLS)
            IAPOS = IAPOS + 1
   10     CONTINUE
   20   CONTINUE
        K = MIN0(10,N)
        IF (LDIAG.EQ.4) K = N
        IF (N.GT.0) WRITE (MP,99993) (RHS(I),I=1,K)
      ENDIF
C*************
C*************
      IF (MTYPE.EQ.1) THEN
        CALL MA41R(N, A, LA, IW(1), LIW, W, MAXFRT, RHS,
     *             PTLUST, NSTEPS )
        CALL MA41S(N, A, LA, IW(1), LIW, W, MAXFRT, RHS,
     *              PTLUST, NSTEPS, W2)
      ELSE
        CALL MA41T(N, A, LA, IW(1), LIW, W, MAXFRT, RHS,
     *              PTLUST, NSTEPS )
        CALL MA41U(N, A, LA, IW(1), LIW, W, MAXFRT, RHS,
     *              PTLUST, NSTEPS, W2)
      ENDIF
      DO 100 KL=1,N
        RHS(KL) = W2(KL)
  100 CONTINUE
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
       K = MIN0(10,N)
       IF (LDIAG.EQ.4) K = N
       WRITE (MP,99992)
       IF (N.GT.0) WRITE (MP,99993) (RHS(I),I=1,K)
      ENDIF
      RETURN
99999 FORMAT (/'Entering solve phase (MA41C) with ......'/
     * '        N         LA         LIW     MAXFRT     NSTEPS'/,
     * 1X, I8, I11, I12, 2I11)
99998 FORMAT (/'Block pivot =', I7, '   NPIV =', I7, '  NFRONT =', I7/
     *        'IAPOS = ', I7)
99997 FORMAT ('Row indices    =', 9I6/(17X, 9I6))
99996 FORMAT ('Column indices =', 9I6/(17X, 9I6))
99995 FORMAT ('Real entries .. each row/column starts on a new line')
99994 FORMAT (1P,5E14.6)
99993 FORMAT ('RHS'/(1P,5E14.6))
99992 FORMAT (//'Leaving solve (MA41C) with')
      END
CCC
      SUBROUTINE MA41D(N, IPE, IW, LW, IWFR,NCMPA)
      INTEGER N,LW,IWFR,NCMPA
      INTEGER IPE(N)
      INTEGER   IW(LW)
      INTEGER I,K1,LWFR,IR,K,K2
      NCMPA = NCMPA + 1
      DO 10 I=1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE
      IWFR = 1
      LWFR = IWFR
      DO 60 IR=1,N
        IF (LWFR.GT.LW) GO TO 70
        DO 20 K=LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50
        DO 40 K=K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
      END
CCC
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA41E(N,NZ,IRN,ICN,ASPK,
     *   IS, MAXIS, S, MAXS, RHS, COLSCA,ROWSCA,
     *   LANAL, LFACTO, LSOLVE, CNTL, ICNTL,
     *   INFO, RINFO, KEEP)
C*******
      INTEGER NBUD
      PARAMETER (NBUD=29)
      INTEGER N,NZ,MAXIS,MAXS
      INTEGER IRN(NZ), ICN(NZ), INFO(20),
     *        IS(MAXIS), ICNTL(20)
      LOGICAL LANAL, LFACTO, LSOLVE
      REAL S(MAXS),ASPK(NZ),RHS(N),
     *     COLSCA(*),ROWSCA(*), RSOL(1), RINFO(20)
      REAL CNTL(10)
      INTEGER KEEP(50)
      INTEGER FILS, FRERE, PTRIST, PTRWB, PTRAR,
     *        PTLUST, IKEEP, ND, LA, LIW, IPTA
      INTEGER I, K, JJ, JPERM, LENA, NPROCS,
     *        ISTIW, IPOOL, LPOOL, ITLOC, IRW1,
     *        IRW2, ISTW1, ISTW3, NSTEPS, MAXFRT
      INTEGER IBEG1, IAV1, IPTBEG, IPT,
     *        LREQ, LSIZ, IPTEND, NBACTI, LCKS(20), IRES
      INTEGER MPRINT,LP,MP
      LOGICAL GIVSOL, ERANAL, LSCAL
      LOGICAL PROK, PB, PBSOLV, PBIR
      INTRINSIC INT
      EXTERNAL MA41B, MA41C, MA41F, MA41O, MA41Q, MA41V, 
     *        MC51A, MC51U, MC51Z
      MPRINT = ICNTL(3)
      LP     = ICNTL(1)
      MP     = ICNTL(2)
      PROK   = (MPRINT.GE.0)
      NPROCS = ICNTL(5)
CPA FINAL Cif defined(1)
CPA FINAL       NPROCS = 1
CPA FINAL Celse
      IF (ICNTL(5).LE.0) NPROCS=1
CPA FINAL Cendif
      IF (LANAL) THEN
       IF (PROK) WRITE(MPRINT,101)
 101    FORMAT(/'****** Analysis step ********'/)
       KEEP(23) = ICNTL(6)
       IF ((KEEP(23).EQ.1).AND.(ICNTL(7).EQ.1)) KEEP(23) = 0
       IF (KEEP(23).EQ.1) THEN
          LIW = MAXIS - N
          IF (LIW.LT.6*N+NZ) GO TO 460
          CALL MA41O(N, NZ, KEEP(23), IS(1),
     *       IRN, ICN, IS(N+1), LIW, ICNTL, INFO)
          IF (INFO(1).LT.0) GO TO 500
       ENDIF
       IF (KEEP(23).EQ.1) THEN
        IKEEP = N + 1
       ELSE
        IKEEP = 1
       ENDIF
       FILS  = 3*N+IKEEP
       FRERE = FILS+N
       PTRAR = FRERE + N
       ISTIW = PTRAR + 4*N
       LIW    = MAXIS - ISTIW + 1
       IF ((ICNTL(7).EQ.1).AND.(LIW.LT.NZ+3*N+1)) GO TO 475
       IF ((ICNTL(7).NE.1).AND.(LIW.LT.2*NZ+2*N+1)) GO TO 480
        CALL MA41F(N,NZ,IRN,ICN,IS(ISTIW),LIW,
     1   IS(IKEEP),IS(PTRAR),ICNTL(7), IS(ISTIW), IS(FILS),
     3   IS(FRERE), ICNTL, INFO, RINFO, KEEP)
        KEEP(25) = INFO(3)
        KEEP(26) = INFO(4)
        KEEP(27) = INFO(5)
        KEEP(28) = INFO(6)
        IF (INFO(1).LT.0) GO TO 500
        LPOOL      = INFO(6)+N+1
        INFO(7) = ISTIW-1 + INFO(6)+ 6*N+ N*NPROCS
     *               + LPOOL + KEEP(15) + NBUD
        KEEP(29) = INFO(7)
        KEEP(30) = INFO(8)
        IF (PROK) WRITE(MPRINT,99989) INFO(7), INFO(8)
      ENDIF
      IF (LFACTO) THEN
        IF (PROK) WRITE(MPRINT,102)
 102    FORMAT(/'****** Factorization step ********'/)
        IF (MAXS.LT.KEEP(30))  GO TO 485
        IF (MAXIS.LT.KEEP(29)) GO TO 490
        IF (PROK)
     *    WRITE(MPRINT,99990) MAXS, MAXIS
        IF (KEEP(23).EQ.1) THEN
          IKEEP = N+1
        ELSE
          IKEEP = 1
        ENDIF
        LSCAL = ((ICNTL(8).GT.0).AND.(ICNTL(8).LE.6))
        IF (LSCAL) THEN
         CALL MC51A(N, NZ, ICNTL(8), ASPK, IRN, ICN,
     *       COLSCA, ROWSCA, S, MAXS, ICNTL, INFO)
         IF (INFO(1).LT.0) GO TO 500
        ENDIF
        LSCAL = (LSCAL.OR.(ICNTL(8).EQ.-1))
        NSTEPS = KEEP(28)
        LPOOL  = NSTEPS+1+N
        FILS   = IKEEP+ 3*N
        FRERE  = FILS  + N
        PTRAR  = FRERE + N
        ND     = PTRAR + 4*N
        PTLUST = ND + N
        IPTA   = PTLUST + NSTEPS
        PTRIST = IPTA + NBUD
        PTRWB  = PTRIST  + N
        ITLOC  = PTRWB + 4*N
        IPOOL  = ITLOC + N*NPROCS
        ISTIW  = IPOOL + LPOOL
        KEEP(24) = ISTIW
        LIW    = MAXIS - ISTIW + 1
        CALL MA41B(N,NZ,NSTEPS,IRN,ICN,ASPK,
     1       S,MAXS,IS(ISTIW),LIW,IS(IKEEP),
     2       IS(ND),IS(FILS),IS(FRERE),IS(PTRAR),
     3       IS(PTRIST),IS(PTLUST),IS(PTRWB),
     4       NPROCS,IS(ITLOC),IS(IPOOL),LPOOL,
     5       LSCAL,COLSCA,ROWSCA,IS(IPTA),
     6       CNTL,ICNTL,INFO,RINFO,KEEP)
        IF (INFO(1).EQ.-9)  INFO(2)=INFO(2)+MAXS
        IF (INFO(1).EQ.-8)  INFO(2)=INFO(2)+MAXIS
        KEEP(31) = INFO(9)
        KEEP(32) = INFO(10)
        KEEP(33) = INFO(11)
        IF (INFO(1).LT.0) GO TO 500
      ENDIF
      IF (LSOLVE) THEN
        IF (PROK) WRITE(MPRINT,99998)
        IF (PROK) WRITE(MPRINT,99999) MAXS, ICNTL(9), ICNTL(10),
     *            ICNTL(11)
C****************************************
C=====================================================
       ERANAL = ((ICNTL(11).GT.0).OR.(ICNTL(10).GT.0))
       NSTEPS = KEEP(28)
       MAXFRT = KEEP(33)
       LA     = KEEP(31)-KEEP(22)
       LIW    = KEEP(32)
       IF (KEEP(23).EQ.1) THEN
         IKEEP = N+1
       ELSE
         IKEEP = 1
       ENDIF
       PBSOLV =  .FALSE.
       PBIR   =  .FALSE.
       IF (KEEP(22).EQ.0) THEN
        IBEG1   = LA + 1
       ELSE
        LA      = MAXS
        IBEG1   = KEEP(18) + KEEP(19) + KEEP(20) + 1
       ENDIF
       IAV1   = MAXS - IBEG1  + 1
       IRW1   = IBEG1
       IRW2   = IRW1 + N
       IF (N+MAXFRT.GT.IAV1) PBSOLV = .TRUE.
       IF (ERANAL) THEN
          IF (IAV1.GE.6*N) THEN
             ISTW3 = IBEG1
             IRW1  = ISTW3 + N
             IRW2  = IRW1 + N
          ELSE
             PBIR  = .TRUE.
          ENDIF
       ENDIF
       PB = (PBSOLV.OR.PBIR)
       IF ((KEEP(22).NE.0).AND.(PB)) THEN
C======================================================
C======================================================
         NBACTI = 0
         IPTA   = 10*N+IKEEP+NSTEPS
         LENA   = KEEP(21)
         IPTBEG = KEEP(18)
         IPTEND = IPTBEG+KEEP(19)+1
         IF (PBIR) THEN
          LREQ  = 5*N
          CALL MC51Z(S, MAXS, IS(IPTA), LREQ, IPT,
     *                  LSIZ, IRES,
     *                  LENA,NBACTI,LCKS)
          IF (IRES.GE.0) THEN
           ISTW3 = MAXS - KEEP(13) + 1
           IRW1 = IPT
           CALL MC51U(S, MAXS, IS(IPTA), IRW1,
     *       IPTEND, IPTBEG, LCKS)
           IRW2  = IRW1 + N
           PBIR  = .FALSE.
           PBSOLV = .FALSE.
          ENDIF
         ENDIF
         IF (PBSOLV) THEN
          LREQ = MAXFRT
          CALL MC51Z(S, MAXS, IS(IPTA), LREQ, IPT,
     *                  LSIZ, IRES,
     *                  LENA,NBACTI,LCKS)
          IF (IRES.GE.0) THEN
           IRW1 = MAXS - KEEP(13) + 1
           IRW2 = IPT
           CALL MC51U(S, MAXS, IS(IPTA), IRW2,
     *       IPTEND, IPTBEG, LCKS)
           PBSOLV = .FALSE.
          ENDIF
         ENDIF
       ENDIF
       IF (PBSOLV) GO TO 450
       LSCAL = ( ((ICNTL(8).GT.0).AND.(ICNTL(8).LE.6)).OR.
     *          (ICNTL(8).EQ.-1) )
       PTLUST = 10*N+IKEEP
       IPTA   = PTLUST+NSTEPS
       ISTIW  = KEEP(24)
       IF (KEEP(23).EQ.1) THEN
         IF (ICNTL(9).NE.1) THEN
           JJ = IRW1
           DO 110 I =1, N
             S(JJ) = RHS(I)
             JJ    = JJ +1
 110       CONTINUE
           JJ = IRW1
           DO 112 I=1,N
             JPERM = IS(I)
             IF (JPERM.EQ.I) GO TO 112
             RHS(JPERM) = S(JJ+I-1)
 112       CONTINUE
         ENDIF
       ENDIF
       IF ( (ERANAL).AND.(.NOT.PBIR) ) THEN
            DO 52 K=1,N
                  S(ISTW3+K-1) = RHS(K)
  52        CONTINUE
       ENDIF
       IF (LSCAL) THEN
         IF (ICNTL(9).EQ.1) THEN
           DO 55 K=1,N
            RHS(K)       = RHS(K)*ROWSCA(K)
  55       CONTINUE
         ELSE
           DO 56 K=1,N
            RHS(K)       = RHS(K)*COLSCA(K)
  56       CONTINUE
         ENDIF
       ENDIF
       CALL MA41C(N, S, LA, IS(ISTIW), LIW, S(IRW2),
     1            MAXFRT, RHS, IS(PTLUST), NSTEPS,
     2            S(IRW1), ICNTL(9), ICNTL)
        IF (LSCAL) THEN
         IF (ICNTL(9).EQ.1) THEN
           DO 213 I=1,N
            RHS(I) = RHS(I)*COLSCA(I)
  213      CONTINUE
         ELSE
           DO 214 I=1,N
            RHS(I) = RHS(I)*ROWSCA(I)
  214      CONTINUE
         ENDIF
        ENDIF
CC
       IF (ERANAL) THEN
        IF (PBIR) THEN
         INFO(1) = -12
         INFO(2) = MAX( INFO(2),MAXS+6*N-IAV1)
        ELSE
        ISTW1  = ISTIW - 2*N
        IF ((ICNTL(10).GT.0).AND.(ICNTL(11).GT.0)) THEN
         GIVSOL = .FALSE.
         IF (MPRINT.GE.0) WRITE(MPRINT,99991)
         CALL MA41Q(ICNTL(9),INFO(1),N,NZ,ASPK,IRN,ICN,
     *         RHS,S(ISTW3),S(IRW1),S(IRW2),GIVSOL,RSOL,
     *         RINFO(4), RINFO(5), RINFO(6), MPRINT, CNTL, ICNTL)
        ENDIF
        CALL MA41V (ICNTL(9), ASPK, NZ, N, IRN, ICN,
     *      S(ISTW3), RHS,
     *      S, LA, IS(ISTIW), LIW, S(IRW1), IS(ISTW1),
     *      IS(PTLUST), NSTEPS, MAXFRT, RINFO,
     *      LSCAL, COLSCA, ROWSCA, CNTL, ICNTL, INFO)
        ENDIF
       ENDIF
 440   IF ((ICNTL(9).EQ.1).AND.(KEEP(23).EQ.1)) THEN
        JJ = IRW1
        DO 210 I =1, N
         S(JJ) = RHS(I)
         JJ    = JJ +1
 210    CONTINUE
        JJ = IRW1
        DO 212 I=1,N
         JPERM = IS(I)
         IF (JPERM.EQ.I) GO TO 212
         RHS(JPERM) = S(JJ+I-1)
 212    CONTINUE
       ENDIF
        IF (PROK) WRITE(MPRINT,99987) INFO(1), INFO(2)
        IF ((INFO(1).GE.0).AND.(ICNTL(4).GE.3).AND.
     *      (ICNTL(2).GT.0)) THEN
         K = MIN0(10,N)
         IF (ICNTL(4).EQ.4) K = N
         WRITE (ICNTL(2),99992)
         IF (N.GT.0) WRITE (ICNTL(2),99993) (RHS(I),I=1,K)
        ENDIF
      ENDIF
      GO TO 500
C*************************
C*************************
 450  INFO(1) = -11
      INFO(2) = MAXS + N + MAXFRT - IAV1
      GO TO 500
 460  INFO(1) = -7
      INFO(2) = 2*NZ+12*N+1
      GO TO 500
 475  INFO(1) = -7
      INFO(2) = NZ+12*N+1
      GO TO 500
 480  INFO(1) = -7
      INFO(2) = 2*NZ+11*N+IKEEP
      GO TO 500
 485  INFO(1) = -9
      INFO(2) = KEEP(30)
      IF (LP.GT.0) THEN
        WRITE(LP,99997) MAXS, MAXIS, KEEP(30), KEEP(29),
     *      KEEP(13),KEEP(14), KEEP(25), KEEP(26), KEEP(27)
      ENDIF
      GO TO 500
 490  INFO(1) = -8
      INFO(2) = KEEP(29)
      IF (LP.GT.0) THEN
        WRITE(LP,99997) MAXS, MAXIS, KEEP(30), KEEP(29),
     *      KEEP(13),KEEP(14), KEEP(25), KEEP(26), KEEP(27)
       ENDIF
 500  RETURN
99999 FORMAT (/'Statistics prior to the solve phase     ...........'/
     1        'MAXS                                          =',I12/
     2        'ICNTL (9)                                     =',I12/
     3        ' --- (10)                                     =',I12/
     4        ' --- (11)                                     =',I12)
99998 FORMAT(//'****** Solve & check step ********'/)
99997 FORMAT(/'Statistics prior to the numerical factorization ...'/
     1        'MAXS                                 =',I12/
     2        'MAXIS                                =',I12/
     3        'Minimum (estimated) value of MAXS    =',I12/
     4        'Minimum (estimated) value of MAXIS   =',I12/
     5        'Real space for original matrix       =',I12/
     6        'Integer space for original matrix    =',I12/
     7        'Real space for factors               =',I12/
     8        'Integer space for factors            =',I12/
     1        'Maximum frontal size (estimated)     =',I12)
99993 FORMAT ('RHS'/(1X,1P,5E14.6))
99992 FORMAT (//'Leaving solver with solution' )
99991 FORMAT (//'Error analysis' )
99990 FORMAT(/'Size of working arrays for factorization:'/
     1        'MAXS                                 =',I12/
     2        'MAXIS                                =',I12)
99989 FORMAT('INFO (7) minimum value of MAXIS    (estimated) =',I12/
     1       'INFO (8) minimum value of MAXS     (estimated) =',I12)
99987 FORMAT(//'Leaving solver with:  INFO(1) ............ =',I12/
     1         '                      INFO(2) ............ =',I12)
      END
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA41F(N, NZ, IRN, ICN, IW, LIW, IKEEP, PTRAR,
     *            IORD, NFSIZ, FILS, FRERE,
     *            ICNTL, INFO, RINFO, KEEP)
C*******
      INTEGER N,NZ,LIW,IORD
      INTEGER PTRAR(N,4), NFSIZ(N), FILS(N), FRERE(N)
      INTEGER   IRN(NZ), ICN(NZ), IW(LIW), IKEEP(N,3)
      REAL RINFO(20)
      INTEGER INFO(20), ICNTL(20), KEEP(50)
      INTEGER K,I,L1,L2,IWFR,NCMPA,LLIW
      INTEGER NEMIN, MPRINT, LP, MP, LDIAG
      LOGICAL PROK
      EXTERNAL MC47B, MA41G, MA41J, MA41K, MA41L, MA41M, MA41N
      MPRINT= ICNTL(3)
      PROK  = (MPRINT.GE.0)
      LP    = ICNTL(1)
      MP    = ICNTL(2)
      LDIAG = ICNTL(4)
CC
CC   CHECK parameters
CC
      IF (KEEP(1).LT.1) KEEP(1) = 1
      NEMIN = KEEP(1)
      IF (LDIAG.LE.2 .OR. MP.LT.0) GO TO 10
      WRITE (MP,99999) N, NZ, LIW, INFO(1)
      K = MIN0(10,NZ)
      IF (LDIAG.EQ.4) K = NZ
      IF (K.GT.0) WRITE (MP,99998) (IRN(I),ICN(I),I=1,K)
      K = MIN0(10,N)
      IF (LDIAG.EQ.4) K = N
      IF (IORD.EQ.1 .AND. K.GT.0) THEN
        WRITE (MP,99997) (IKEEP(I,1),I=1,K)
      ENDIF
   10 LLIW = LIW - 2*N
      L1 = LLIW + 1
      L2 = L1 + N
      IF (IORD.NE.1) THEN
C*******************************
C*******************************
        CALL MA41G(N, NZ, IRN, ICN, IW, LLIW, PTRAR, PTRAR(1,2),
     *              IW(L2), IW(L1), IWFR, KEEP(13), KEEP(14),
     *              INFO(1), INFO(2), ICNTL)
        CALL MC47B(N, LLIW, PTRAR, IWFR, PTRAR(1,2), IW,
     *   IW(L1), IKEEP,
     *   IKEEP(1,2), NCMPA, FILS, IKEEP(1,3), IW(L2), PTRAR(1,3))
      ELSE
C***********************************
C***********************************
       DO 15 K=1,N
         IF ((IKEEP(K,1).LE.0).OR.(IKEEP(K,1).GT.N))
     *    GO TO 40
 15    CONTINUE
       CALL MA41J(N, NZ, IRN, ICN, IKEEP, IW, LLIW, PTRAR,
     *             PTRAR(1,2), IW(L1), IWFR,
     *             INFO(1),INFO(2), MP)
       CALL MA41K(N, PTRAR, IW, LLIW, IWFR, IKEEP,
     *    IKEEP(1,2), IW(L1),
     *    IW(L2), NCMPA)
      ENDIF
C********************************************
C********************************************
      CALL MA41L(N, PTRAR, IW(L1), IKEEP, IKEEP(1,2), IKEEP(1,3),
     *          NFSIZ, INFO(6), FILS, FRERE, PTRAR(1,3), NEMIN)
C*****************************************************************
C*****************************************************************
      CALL MA41M(N, NZ, IRN, ICN, IKEEP(1,3), IKEEP(1,2),
     * PTRAR(1,3), INFO(6), PTRAR(1,2), IORD,
     * KEEP(13), KEEP(14), INFO(3), INFO(4),
     * KEEP(16), INFO(8),KEEP(15),INFO(5), RINFO(1),
     * KEEP(12))
C*****************************************************************
C*****************************************************************
      CALL MA41N(N, NZ, IKEEP(1,1), FILS, FRERE,
     *      IKEEP(1,3), IKEEP(1,2),
     *      IRN, ICN, PTRAR(1,3), PTRAR(1,4), PTRAR)
C*************************
C*************************
      IF (PROK) THEN
        WRITE (MPRINT,99992) INFO(1), INFO(2), INFO(3), INFO(4),
     *      INFO(5), INFO(6), ICNTL(6), ICNTL(7), RINFO(1)
      ENDIF
      IF (LDIAG.GT.2 .AND. MP.GE.0) THEN
       IF (.NOT.PROK) THEN
        WRITE (MP,99992) INFO(1), INFO(2), INFO(3), INFO(4),
     *      INFO(5),INFO(6), ICNTL(6), ICNTL(7), RINFO(1)
       ENDIF
       K = MIN0(10,N)
       IF (LDIAG.EQ.4) K = N
       IF (K.GT.0) WRITE (MP,99997) (IKEEP(I,1),I=1,K)
       IF (K.GT.0) WRITE (MP,99991) (IKEEP(I,2),I=1,K)
       IF (K.GT.0) WRITE (MP,99990) (IKEEP(I,3),I=1,K)
       IF (K.GT.0) WRITE (MP,99986) (PTRAR(I,1),I=1,K)
       IF (K.GT.0) WRITE (MP,99985) (PTRAR(I,2),I=1,K)
       IF (K.GT.0) WRITE (MP,99984) (PTRAR(I,3),I=1,K)
       IF (K.GT.0) WRITE (MP,99983) (PTRAR(I,4),I=1,K)
       IF (K.GT.0) WRITE (MP,99987) (NFSIZ(I),I=1,K)
       IF (K.GT.0) WRITE (MP,99989) (FILS(I),I=1,K)
       IF (K.GT.0) WRITE (MP,99988) (FRERE(I),I=1,K)
      ENDIF
      GO TO 90
C*************************
C*************************
   40 INFO(1) = -4
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99996) INFO(1)
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1)) WRITE (LP,99982) INFO(2)
   90 RETURN
99999 FORMAT (/'Entering analysis phase (MA41F) with ...'/
     * '                N         NZ         LIW       INFO(1)'/,
     * 9X, I8, I11, I12, I14)
99998 FORMAT ('Matrix entries:    IRN()   ICN()'/
     * (I12, I7, I12, I7, I12, I7))
99997 FORMAT ('IKEEP(.,1)=', 10I6/(12X, 10I6))
99996 FORMAT (/'** Error return ** from MA41F *  INFO(1)=', I3)
99994 FORMAT ('Value of NZ out-of-range .. =', I10)
99993 FORMAT ('LIW too small, must be increased from', I10, ' to at ',
     * 'least', I10)
99992 FORMAT(/'Leaving analysis phase (MA41F) with  ...'/
     1       'INFO (1)                                       =',I12/
     2       'INFO (2)                                       =',I12/
     3       ' --  (3) Real space for factors    (estimated) =',I12/
     4       ' --  (4) Integer space for factors (estimated) =',I12/
     5       ' --  (5) Maximum frontal size      (estimated) =',I12/
     6       ' --  (6) Number of nodes in the tree           =',I12/
     7       'ICNTL(6) Maximum transversal option            =',I12/
     8       'ICNTL(7) Pivot order option                    =',I12/
     9   'RINFO(1) Operations during elimination (estim) =  ',1PE10.3)
99991 FORMAT ('IKEEP(.,2)=', 10I6/(12X, 10I6))
99990 FORMAT ('IKEEP(.,3)=', 10I6/(12X, 10I6))
99989 FORMAT ('FILS (.)  =', 10I6/(12X, 10I6))
99988 FORMAT ('FRERE(.)  =', 10I6/(12X, 10I6))
99987 FORMAT ('NFSIZ(.)  =', 10I6/(12X, 10I6))
99986 FORMAT ('PTRAR(.,1)=', 10I6/(12X, 10I6))
99985 FORMAT ('PTRAR(.,2)=', 10I6/(12X, 10I6))
99984 FORMAT ('PTRAR(.,3)=', 10I6/(12X, 10I6))
99983 FORMAT ('PTRAR(.,4)=', 10I6/(12X, 10I6))
99982 FORMAT ('Out-of-range value in KEEP        INFO(2)=', I3)
      END
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA41G(N,NZ, IRN, ICN, IW, LW, IPE, LEN,
     *  IQ, FLAG, IWFR,
     * NRORM, NIORM, IFLAG,IERROR, ICNTL)
      INTEGER N,NZ,LW,IFLAG,IERROR,NRORM,NIORM,IWFR
      INTEGER IRN(NZ), ICN(NZ), IW(LW), FLAG(N), LEN(N)
      INTEGER IPE(N), IQ(N), ICNTL(20)
      INTEGER MP
      INTEGER I,K,J,N1,LAST,NDUP,K1,K2,L
      INTEGER NBERR
      MP = ICNTL(2)
      NIORM  = 3*N
      IERROR = 0
      DO 10 I=1,N
        IPE(I) = 0
   10 CONTINUE
      DO 50 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     *                          .OR.(J.LT.1)) THEN
           IERROR = IERROR + 1
        ELSE
          IF (I.NE.J) THEN
           IPE(I) = IPE(I) + 1
           IPE(J) = IPE(J) + 1
           NIORM  = NIORM + 1
          ENDIF
        ENDIF
   50 CONTINUE
      IF (IERROR.GE.1) THEN
         NBERR  = 0
         IFLAG  = IFLAG+1
         IF ((MP.GE.0).AND.(ICNTL(4).GE.2))  THEN
          WRITE (MP,99999)
          DO 70 K=1,NZ
           I = IRN(K)
           J = ICN(K)
           IF ((I.GT.N).OR.(J.GT.N).OR.(I.LT.1)
     *                            .OR.(J.LT.1)) THEN
            NBERR = NBERR + 1
            IF (NBERR.LE.10)  THEN
               IF (MOD(K,10).GT.3 .OR. MOD(K,10).EQ.0 .OR.
     *             (10.LE.K .AND. K.LE.20)) THEN
                 WRITE (MP,'(I8,A,I8,A,I8,A)')
     *             K,'th entry (in row',I,' and column',J,') ignored'
               ELSE
                 IF (MOD(K,10).EQ.1) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'st entry (in row',I,' and column',J,') ignored'
                 IF (MOD(K,10).EQ.2) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'nd entry (in row',I,' and column',J,') ignored'
                 IF (MOD(K,10).EQ.3) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'rd entry (in row',I,' and column',J,') ignored'
               ENDIF
            ELSE
               GO TO 100
            ENDIF
           ENDIF
   70     CONTINUE
         ENDIF
      ENDIF
  100 NRORM = NIORM - 2*N
      IQ(1) = 1
      N1 = N - 1
      IF (N1.GT.0) THEN
        DO 110 I=1,N1
            IQ(I+1) = IPE(I) + IQ(I)
  110   CONTINUE
      ENDIF
      LAST = MAX(IPE(N)+IQ(N)-1,IQ(N))
      DO 115 I = 1,N
         FLAG(I) = 0
         IPE(I)  = IQ(I)
  115 CONTINUE
      DO 130 K=1,LAST
        IW(K) = 0
  130 CONTINUE
      IWFR = LAST + 1
      DO 200 K=1,NZ
         I = IRN(K)
         J = ICN(K)
         IF (I.NE.J) THEN
          IF (I.LT.J) THEN
            IF ((I.GE.1).AND.(J.LE.N)) THEN
             IW(IQ(I)) = -J
             IQ(I)     = IQ(I) + 1
            ENDIF
          ELSE
            IF ((J.GE.1).AND.(I.LE.N)) THEN
             IW(IQ(J)) = -I
             IQ(J)     = IQ(J) + 1
            ENDIF
          ENDIF
         ENDIF
  200 CONTINUE
      NDUP = 0
      DO 260 I=1,N
        K1 = IPE(I)
        K2 = IQ(I) -1
        IF (K1.GT.K2) THEN
         LEN(I) = 0
         IPE(I) = 0
         IQ(I)  = 0
        ELSE
         DO 240 K=K1,K2
           J     = -IW(K)
           IF (J.LE.0) GO TO 250
           L     = IQ(J)
           IQ(J) = L + 1
           IF (FLAG(J).EQ.I) THEN
            NDUP = NDUP + 1
            IW(L) = 0
            IW(K) = 0
           ELSE
            IW(L)   = I
            IW(K)   = J
            FLAG(J) = I
           ENDIF
  240    CONTINUE
  250    IQ(I) = IQ(I) - IPE(I)
         IF (NDUP.EQ.0) LEN(I) = IQ(I)
        ENDIF
  260 CONTINUE
      IF (NDUP.NE.0) THEN
       IWFR = 1
       DO 280 I=1,N
         K1 = IPE(I)
         IF (K1.EQ.0) GO TO 280
         K2 = K1 + IQ(I) - 1
         L = IWFR
         IPE(I) = IWFR
         DO 270 K=K1,K2
           IF (IW(K).NE.0) THEN
            IW(IWFR) = IW(K)
            IWFR     = IWFR + 1
           ENDIF
  270    CONTINUE
         LEN(I) = IWFR - L
  280  CONTINUE
      ENDIF
      RETURN
99999 FORMAT (/'*** Warning message from subroutine MA41F ***')
      END
CCC
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA41H(N, NZ, NZ1, ASPK, A, LA, IRN, ICN, IW,
     *   LIW ,PERM, IW4, PTRAR, PTRARW, PTRAIW,
     *   LSCAL,COLSCA,ROWSCA)
C***********************
      INTEGER N,NZ,LA,LIW
      REAL A(LA), ASPK(NZ), COLSCA(*), ROWSCA(*)
      INTEGER IW4(N,2), NZ1(2), PTRAR(N,2), PTRARW(N), PTRAIW(N)
      LOGICAL LSCAL
      INTEGER   IRN(NZ), ICN(NZ), IW(LIW), PERM(N), ISR, ISI
      REAL ZERO
      INTEGER IOLD,IA,IIW,I,I1,K,JOLD,INEW,JNEW
      INTEGER IS1,ISHIFT,IS,IAS
      PARAMETER (ZERO=0.0E0)
      DO 10 IOLD=1,N
        IW4(IOLD,1) = PTRAR(IOLD,1)
        IW4(IOLD,2) = PTRAR(IOLD,2)
   10 CONTINUE
      IA  = LA - NZ1(1) + 1
      IIW = LIW - NZ1(2) + 1
      ISR = IA - PTRARW(1)
      ISI = IIW - PTRAIW(1)
      IF ((ISR.NE.0).OR.(ISI.NE.0)) THEN
        DO 90 I=1,N
         PTRARW(I) = PTRARW(I) + ISR
         PTRAIW(I) = PTRAIW(I) + ISI
  90    CONTINUE
      ENDIF
      DO 100 I=1,N
       I1       = PTRAIW(I)
       IA       = PTRARW(I)
       IW(I1)   = IW4(I,1)
       IW(I1+1) = -IW4(I,2)
       IW(I1+2) = I
       A(IA)    = ZERO
  100 CONTINUE
CCC
CC store off diagonal values
CCC
      DO 120 K=1,NZ
        IOLD = IRN(K)
        JOLD = ICN(K)
        IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     *                 .OR.(JOLD.LT.1) ) GOTO 120
        IF (IOLD.EQ.JOLD) THEN
         IA = PTRARW(IOLD)
         IF (LSCAL) THEN
          A(IA) = A(IA)+ASPK(K)*ROWSCA(IOLD)*COLSCA(JOLD)
         ELSE
          A(IA) = A(IA)+ASPK(K)
         ENDIF
        ELSE
         INEW = PERM(IOLD)
         JNEW = PERM(JOLD)
         IF (INEW.LT.JNEW) THEN
CCC
CCC
           IS1         = PTRAIW(IOLD)
           ISHIFT      = IW(IS1) + IW4(IOLD,2)
           IW4(IOLD,2) = IW4(IOLD,2) - 1
           IIW         = IS1 + ISHIFT + 2
           IW(IIW)     = JOLD
           IS          = PTRARW(IOLD)
           IAS         = IS + ISHIFT
           IF (LSCAL) THEN
             A(IAS)      = ASPK(K)*ROWSCA(IOLD)*COLSCA(JOLD)
           ELSE
             A(IAS)      = ASPK(K)
           ENDIF
         ELSE
CCC
CCC
           ISHIFT      = PTRAIW(JOLD)+IW4(JOLD,1)+2
           IW(ISHIFT)  = IOLD
           IAS         = PTRARW(JOLD)+IW4(JOLD,1)
           IW4(JOLD,1) = IW4(JOLD,1) - 1
           IF (LSCAL) THEN
             A(IAS)    = ASPK(K)*ROWSCA(IOLD)*COLSCA(JOLD)
           ELSE
             A(IAS)    = ASPK(K)
           ENDIF
         ENDIF
        ENDIF
  120 CONTINUE
      RETURN
      END
CCC
      SUBROUTINE MA41J(N, NZ, IRN, ICN, PERM, IW, LW, IPE, IQ, FLAG,
     * IWFR, IFLAG, IERROR, MP)
      INTEGER N,NZ,LW,IWFR,IFLAG,IERROR
      INTEGER IRN(NZ), ICN(NZ), IW(LW), PERM(N), FLAG(N)
      INTEGER IPE(N), IQ(N)
      INTEGER MP
      INTEGER I,J,K,LBIG,L,ID,IN,LEN,JDUMMY,K1,K2
      IERROR = 0
      DO 10 I=1,N
        IQ(I) = 0
   10 CONTINUE
      DO 80 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IW(K) = -I
        IF (I-J.LT.0) THEN
          IF (I.GE.1 .AND. J.LE.N) GO TO 60
        ELSE IF (I-J.GT.0) THEN
          IF (J.GE.1 .AND. I.LE.N) GO TO 60
        ELSE
          IW(K) = 0
          IF (I.GE.1 .AND. I.LE.N) GO TO 80
        ENDIF
        IERROR = IERROR + 1
        IW(K) = 0
        IF (IERROR.LE.1 .AND. MP.GE.0) WRITE (MP,99999)
        IF (IERROR.LE.10 .AND. MP.GE.0) THEN
          IF (MOD(K,10).GT.3 .OR. MOD(K,10).EQ.0 .OR.
     *      (10.LE.K .AND. K.LE.20)) THEN
            WRITE (MP,'(I8,A,I8,A,I8,A)')
     *             K,'th entry (in row',I,' and column',J,') ignored'
          ELSE
            IF (MOD(K,10).EQ.1) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'st entry (in row',I,' and column',J,') ignored'
            IF (MOD(K,10).EQ.2) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'nd entry (in row',I,' and column',J,') ignored'
            IF (MOD(K,10).EQ.3) WRITE(MP,'(I8,A,I8,A,I8,A)')
     *             K,'rd entry (in row',I,' and column',J,') ignored'
          ENDIF
        ENDIF
        GO TO 80
   60   IF (PERM(J).GT.PERM(I)) GO TO 70
        IQ(J) = IQ(J) + 1
        GO TO 80
   70   IQ(I) = IQ(I) + 1
   80 CONTINUE
      IF (IERROR.GE.1) IFLAG = IFLAG+1
      IWFR = 1
      LBIG = 0
      DO 100 I=1,N
        L = IQ(I)
        LBIG = MAX0(L,LBIG)
        IWFR = IWFR + L
        IPE(I) = IWFR - 1
  100 CONTINUE
      DO 140 K=1,NZ
        I = -IW(K)
        IF (I.LE.0) GO TO 140
        L = K
        IW(K) = 0
        DO 130 ID=1,NZ
          J = ICN(L)
          IF (PERM(I).LT.PERM(J)) GO TO 110
          L = IPE(J)
          IPE(J) = L - 1
          IN = IW(L)
          IW(L) = I
          GO TO 120
  110     L = IPE(I)
          IPE(I) = L - 1
          IN = IW(L)
          IW(L) = J
  120     I = -IN
          IF (I.LE.0) GO TO 140
  130   CONTINUE
  140 CONTINUE
      K = IWFR - 1
      L = K + N
      IWFR = L + 1
      DO 170 I=1,N
        FLAG(I) = 0
        J = N + 1 - I
        LEN = IQ(J)
        IF (LEN.LE.0) GO TO 160
        DO 150 JDUMMY=1,LEN
          IW(L) = IW(K)
          K = K - 1
          L = L - 1
  150   CONTINUE
  160   IPE(J) = L
        L = L - 1
  170 CONTINUE
      IF (LBIG+N.LE.LBIG) GO TO 190
      DO 180 I=1,N
        K = IPE(I)
        IW(K) = IQ(I)
        IF (IQ(I).EQ.0) IPE(I) = 0
  180 CONTINUE
      GO TO 230
  190 IWFR = 1
      DO 220 I=1,N
        K1 = IPE(I) + 1
        K2 = IPE(I) + IQ(I)
        IF (K1.LE.K2) GO TO 200
        IPE(I) = 0
        GO TO 220
  200   IPE(I) = IWFR
        IWFR = IWFR + 1
        DO 210 K=K1,K2
          J = IW(K)
          IF (FLAG(J).EQ.I) GO TO 210
          IW(IWFR) = J
          IWFR = IWFR + 1
          FLAG(J) = I
  210   CONTINUE
        K = IPE(I)
        IW(K) = IWFR - K - 1
  220 CONTINUE
  230 RETURN
99999 FORMAT (/'*** Warning message from subroutine MA41J ***')
      END
C**END of MA41J***
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA41K(N, IPE, IW, LW, IWFR, IPS, IPV, NV, FLAG,
     *                  NCMPA)
      INTEGER N,LW,IWFR,NCMPA
      INTEGER IPE(N)
      INTEGER IW(LW), IPS(N), IPV(N), NV(N), FLAG(N)
      INTEGER I,J,ML,MS,ME,IP,MINJS,IE,KDUMMY,JP
      INTEGER LN,JP1,JS,LWFR,JP2,JE
      EXTERNAL MA41D
      DO 10 I=1,N
        FLAG(I) = 0
        NV(I) = 0
        J = IPS(I)
        IPV(J) = I
   10 CONTINUE
      NCMPA = 0
      DO 100 ML=1,N
        MS = IPV(ML)
        ME = MS
        FLAG(MS) = ME
        IP = IWFR
        MINJS = N
        IE = ME
        DO 70 KDUMMY=1,N
          JP = IPE(IE)
          LN = 0
          IF (JP.LE.0) GO TO 60
          LN = IW(JP)
          DO 50 JP1=1,LN
            JP = JP + 1
            JS = IW(JP)
            IF (FLAG(JS).EQ.ME) GO TO 50
            FLAG(JS) = ME
            IF (IWFR.LT.LW) GO TO 40
            IPE(IE) = JP
            IW(JP) = LN - JP1
            CALL MA41D(N, IPE, IW, IP-1, LWFR,NCMPA)
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 30
            DO 20 JP=IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   20       CONTINUE
   30       IP = LWFR
            JP = IPE(IE)
   40       IW(IWFR) = JS
            MINJS = MIN0(MINJS,IPS(JS)+0)
            IWFR = IWFR + 1
   50     CONTINUE
   60     IPE(IE) = -ME
          JE = NV(IE)
          NV(IE) = LN + 1
          IE = JE
          IF (IE.EQ.0) GO TO 80
   70   CONTINUE
   80   IF (IWFR.GT.IP) GO TO 90
        IPE(ME) = 0
        NV(ME) = 1
        GO TO 100
   90   MINJS = IPV(MINJS)
        NV(ME) = NV(MINJS)
        NV(MINJS) = ME
        IW(IWFR) = IW(IP)
        IW(IP) = IWFR - IP
        IPE(ME) = IP
        IWFR = IWFR + 1
  100 CONTINUE
      RETURN
      END
C** end of MA41K**
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA41L(N, IPE, NV, IPS, NE, NA, NFSIZ, NSTEPS,
     *                  FILS, FRERE,NDD,NEMIN)
      INTEGER N,NSTEPS
      INTEGER NDD(N)
      INTEGER IPE(N), FILS(N), FRERE(N)
      INTEGER NV(N), IPS(N), NE(N), NA(N), NFSIZ(N)
      INTEGER NEMIN
      INTEGER I,IF,IS,NR,NR1,INS,INL,INB,INF,INFS,INSW
      INTEGER K,L,ISON,IN,INP,IFSON,INC,INO
      INTEGER INOS,IB,IL,INT
      DO 10 I=1,N
        IPS(I) = 0
        NE(I) = 0
   10 CONTINUE
      DO 20 I=1,N
        IF (NV(I).GT.0) GO TO 20
        IF = -IPE(I)
        IS = -IPS(IF)
        IF (IS.GT.0) IPE(I) = IS
        IPS(IF) = -I
   20 CONTINUE
      NR = N + 1
      DO 50 I=1,N
        IF (NV(I).LE.0) GO TO 50
        IF = -IPE(I)
        IF (IF.NE.0) THEN
         IS = -IPS(IF)
         IF (IS.GT.0) IPE(I) = IS
         IPS(IF) = -I
        ELSE
         NR = NR - 1
         NE(NR) = I
        ENDIF
   50 CONTINUE
      DO 999 I=1,N
       FILS(I) = IPS(I)
 999  CONTINUE
CTEST GO TO 1151
      NR1 = NR
      INS = 0
 1000 IF (NR1.GT.N) GO TO 1151
C--
C---  PICK UP NEXT ROOT
C--
      INS = NE(NR1)
      NR1 = NR1 + 1
 1070 INL = FILS(INS)
      IF (INL.LT.0) THEN
       INS = -INL
       GO TO 1070
      ENDIF
 1080 IF (IPE(INS).LT.0) THEN
       INS       = -IPE(INS)
C--
C---- INS IS A FATHER OF REORGANIZED SONS SO WE CAN
C---- CLEAR THE POINTER TO THE SONS.
C--
       FILS(INS) = 0
       GO TO 1080
      ENDIF
      IF (IPE(INS).EQ.0) THEN
       INS = 0
       GO TO 1000
      ENDIF
      INB = IPE(INS)
CTEST
      IF (NV(INB).EQ.0) THEN
       INS = INB
       GO TO 1070
      ENDIF
      IF (NV(INB).GE.NV(INS)) THEN
       INS = INB
       GO TO 1070
      ENDIF
      INF = INB
 1090 INF = IPE(INF)
      IF (INF.GT.0) GO TO 1090
C--
C---- -INF IS THE FATHER
C--
      INF  = -INF
      INFS = -FILS(INF)
      IF (INFS.EQ.INS) THEN
C--
C----  PREVIOUS BROTHER OF INB: INS IS THE FIRST SON OF INF
C--
       FILS(INF) = -INB
       IPS(INF)  = -INB
       IPE(INS)  = IPE(INB)
       IPE(INB)  = INS
       INS       = INB
       GO TO 1070
      ENDIF
      INSW = INFS
 1100 INFS = IPE(INSW)
      IF (INFS.NE.INS) THEN
       INSW = INFS
       GO TO 1100
      ENDIF
      IPE(INS) = IPE(INB)
      IPE(INB) = INS
      IPE(INSW)= INB
      INS      =INB
      GO TO 1070
 1151 CONTINUE
      DO 51 I=1,N
       FRERE(I) = IPE(I)
       FILS(I) = IPS(I)
 51   CONTINUE
      IS = 1
      I = 0
      DO 160 K=1,N
        IF (I.GT.0) GO TO 60
        I = NE(NR)
        NE(NR) = 0
        NR = NR + 1
        IL = N
        NA(N) = 0
   60   DO 70 L=1,N
          IF (IPS(I).GE.0) GO TO 80
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
   80   IPS(I) = K
        NE(IS) = NE(IS) + 1
        IF (NV(I).GT.0) GO TO 89
      IN = I
 81   IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 81
      IF = -IN
      IN = IF
 82   INL = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 82
      IFSON = -IN
      FILS(INL) = I
      IN = I
 83   INP = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 83
      IF (IFSON .EQ. I) GO TO 86
      FILS(INP) = -IFSON
      IN = IFSON
 84   INC =IN
      IN = FRERE(IN)
      IF (IN.NE.I) GO TO 84
      FRERE(INC) = FRERE(I)
      GO TO 120
 86   IF (FRERE(I).LT.0) FILS(INP) = 0
      IF (FRERE(I).GT.0) FILS(INP) = -FRERE(I)
      GO TO 120
   89   IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        NDD(IS) = NV(I)
        NFSIZ(I) = NV(I)
        IF (NA(IS).LT.1) GO TO 110
        IF (NDD(IS-1)-NE(IS-1).EQ.NDD(IS)) GO TO 100
        IF (NE(IS).GE.NEMIN) GO TO 110
        IF (NE(IS-1).GE.NEMIN) GO TO 110
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        NDD(IS-1) = NDD(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
      IN=I
 101  INL = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 101
      IFSON = -IN
      IN = IFSON
 102  INO = IN
      IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 102
      FILS(INL) = INO
      NFSIZ(I) = NDD(IS-1)
      IN = INO
 103  INP = IN
      IN = FILS(IN)
      IF (IN.GT.0) GO TO 103
      INOS = -IN
      IF (IFSON.EQ.INO) GO TO 107
      IN = IFSON
      FILS(INP) = -IFSON
 105  INS = IN
      IN =  FRERE(IN)
      IF (IN.NE.INO) GO TO 105
      IF (INOS.EQ.0) FRERE(INS) = -I
      IF (INOS.NE.0) FRERE(INS) =  INOS
      IF (INOS.EQ.0) GO TO 109
 107  IN = INOS
      IF (IN.EQ.0) GO TO 109
 108  INT = IN
      IN =  FRERE(IN)
      IF (IN.GT.0) GO TO 108
      FRERE(INT) = -I
 109  CONTINUE
        GO TO 120
  110   IS = IS + 1
  120   IB = IPE(I)
        IF (IB.GE.0) THEN
          IF (IB.GT.0) NA(IL) = 0
          I = IB
        ELSE
          I = -IB
          IL = IL + 1
        ENDIF
  160 CONTINUE
      NSTEPS = IS - 1
      DO 170 I=1,N
        K = FILS(I)
        IF (K.GT.0) THEN
          FRERE(K)  = N + 1
          NFSIZ(K)  = 0
        ENDIF
 170  CONTINUE
      RETURN
      END
      SUBROUTINE MA41M(N, NZ, IRN, ICN, NA, NE, ND, NSTEPS,
     * LSTKR, IORD,
     * NRORM, NIORM, NRLADU, NIRADU, NRLNEC, NRLTOT, NIRTOT,
     * MAXFR, OPSA, PERLU)
      INTEGER  N,NZ,NSTEPS,IORD
      INTEGER  NRORM, NIORM, NRLADU, NIRADU, NRLNEC
      INTEGER  NRLTOT, NIRTOT, MAXFR, PERLU
      INTEGER  NA(NSTEPS), NE(NSTEPS), ND(NSTEPS),
     *         IRN(NZ), ICN(NZ)
      INTEGER  LSTKR(N)
      INTEGER  NRLBDU, NIRBDU
      INTEGER I,IOLD,JOLD,ISTKR,ITOP,ITREE,NELIM,NFR
      INTEGER K,LSTK,NSTK
      REAL OPSA
CLLTC OPSLLT: ESTIMATION OF NB OF FLOPS FOR LLT
CC
      INTRINSIC MIN,INT,REAL
      IF (IORD.EQ.1) THEN
        NIORM = N*3
        DO 40 I=1,NZ
         IOLD = IRN(I)
         JOLD = ICN(I)
         IF (MIN(IOLD,JOLD).LT.1 .OR. IOLD.GT.N
     *    .OR. JOLD.GT.N .OR. IOLD.EQ.JOLD) GO TO 40
         NIORM = NIORM + 1
   40   CONTINUE
        NRORM = NIORM - 2*N
      ENDIF
      ISTKR = 0
      OPSA = 0.0E0
      NRLADU = 0
      NIRADU = 1
      MAXFR  = 0
      NRLTOT = 0
      NRLNEC = 0
      ITOP = 0
      DO 90 ITREE=1,NSTEPS
        NELIM = NE(ITREE)
        NFR = ND(ITREE)
        IF (NFR.GT.MAXFR) MAXFR = NFR
        NSTK = NA(ITREE)
        IF (NSTK.GT.0) THEN
         DO 70 K=1,NSTK
          LSTK = LSTKR(ITOP)
          ISTKR = ISTKR - LSTK
          ITOP = ITOP - 1
   70    CONTINUE
        ENDIF
        NRLADU = NRLADU + NELIM*(2*NFR-NELIM)
        NRLNEC = MAX0(NRLNEC,NRLADU+ISTKR)
CCC We store the frontal matrices
CCC  directly in the LU area:  NIRADU = NIRADU + 3 + 2*NFR  is modified
CC   as follow:
        NIRADU = NIRADU + 5 + 2*NFR
        OPSA = OPSA + REAL(2*NFR*NELIM)*REAL(NFR-NELIM-1)+
     *         REAL(NELIM*(NELIM+1))*REAL(2*NELIM+1)/REAL(3)
        OPSA = OPSA + REAL(((2*NFR-NELIM-1)*NELIM)/2)
        IF (ITREE.EQ.NSTEPS) GO TO 90
        IF (NFR.EQ.NELIM) GO TO 90
        ITOP = ITOP + 1
        LSTKR(ITOP) = (NFR-NELIM)*(NFR-NELIM)+4
        ISTKR = ISTKR + LSTKR(ITOP)
        NRLNEC = MAX0(NRLNEC,NRLADU+ISTKR)
   90 CONTINUE
      NIRADU = MAX(NIRADU,2*N)
      NRLBDU = NRLNEC + PERLU*(NRLNEC/100)
      NIRBDU = NIRADU + PERLU*(NIRADU/100)
      NRLTOT = NRORM  + MAX(NRLBDU+MAXFR*MAXFR+3,N)
      NIRTOT = NIRBDU + NIORM
      RETURN
      END
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA41N( N, NZ, PERM, FILS, FRERE,
     * NSTK, NA, IRN, ICN,
     * PTRAIW, PTRARW, IW4)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER N,NZ
      INTEGER FILS(N), FRERE(N), PERM(N)
      INTEGER NSTK(N), NA(N)
      INTEGER PTRAIW(N), PTRARW(N), IRN(NZ), ICN(NZ), IW4(N,2)
      INTEGER NBLEAF
CCC
      INTEGER NBROOT,ILEAF,I,IN,IPERM,II1,JDEB,J,JTEMP,ISON
      INTEGER IOLD,K,JOLD,INEW,JNEW,ISHIFT,I1
      NA(N)   = 0
      IF (N.GT.1) NA(N-1) = 0
      NBROOT  = 0
      ILEAF   = 1
      DO 11 I=1,N
      NSTK(I) = 0
      NA(I)   = 0
      IF (FRERE(I).EQ. N+1) GO TO 11
      IF (FRERE(I).EQ.0) NBROOT = NBROOT + 1
      IN = I
 12   IN = FILS(IN)
      IF (IN.GT.0) GO TO 12
      IF (IN.LT.0) GO TO 14
      IPERM = PERM(I)
      IF (ILEAF.EQ.1) THEN
       NA(ILEAF) = I
       ILEAF        = ILEAF + 1
       GO TO 11
      ENDIF
      IF (IPERM.GE.PERM(NA(ILEAF-1))) THEN
       NA(ILEAF) = I
       ILEAF        = ILEAF + 1
      ELSE
       DO 20 I1=1,ILEAF-1
        II1 = NA(I1)
        IF (IPERM.LT.PERM(II1)) GO TO 21
 20    CONTINUE
 21    JDEB = I
       DO 25 J=I1,ILEAF-1
        JTEMP    = NA(J)
        NA(J)    = JDEB
        JDEB     = JTEMP
 25    CONTINUE
       NA(ILEAF)  = JDEB
       ILEAF      = ILEAF + 1
      ENDIF
      GO TO 11
 14   ISON = -IN
 13   NSTK(I) = NSTK(I) + 1
      ISON = FRERE(ISON)
      IF (ISON.GT.0) GO TO 13
 11   CONTINUE
      NBLEAF = ILEAF-1
      IF (N.GT.1) THEN
       IF (NBLEAF.GT.N-2) THEN
        IF (NBLEAF.EQ.N-1) THEN
         NA(N-1) = -NA(N-1)-1
         NA(N)   = NBROOT
        ELSE
         NA(N) = -NA(N)-1
        ENDIF
       ELSE
        NA(N-1) = NBLEAF
        NA(N)   = NBROOT
       ENDIF
      ENDIF
CSTAT C
CSTAT C Statistics concerning the three first levels of the
CSTAT C elimination tree (root is level 1)
CSTAT C
CSTAT if defined(perfanal)
CSTAT      IF (ICNTL(3).GT.0) THEN
CSTAT        CALL STATRE(ICNTL(3), N, FRERE,FILS,ND, KEEP(1), NBLEAF,
CSTAT     *     INFO, RINFO)
CSTAT      ENDIF
CSTAT endif
      DO 50 IOLD=1,N
        IW4(IOLD,1) = 0
        IW4(IOLD,2) = 0
   50 CONTINUE
      DO 70 K=1,NZ
        IOLD = IRN(K)
        JOLD = ICN(K)
        IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     *                 .OR.(JOLD.LT.1) ) GO TO 70
        IF (IOLD.NE.JOLD) THEN
         INEW = PERM(IOLD)
         JNEW = PERM(JOLD)
         IF (INEW.LT.JNEW) THEN
CCC
CC
           IW4(IOLD,2) = IW4(IOLD,2) + 1
         ELSE
CCC
CCC
           IW4(JOLD,1) = IW4(JOLD,1) + 1
         ENDIF
        ENDIF
CCCCCCC
CC we increment the nb of nonzeros not on diagonal
CCCCCCC
   70 CONTINUE
CCC
CC
CCC
CC we built the pointers PTRARW(I) & PTRAIW(I)
CC
      PTRARW(1) = 1
      PTRAIW(1) = 1
      IF (N.GT.1) THEN
      DO 90 I=1,N-1
       ISHIFT = IW4(I,1) + IW4(I,2)
       PTRARW(I+1) = PTRARW(I) + ISHIFT + 1
       PTRAIW(I+1) = PTRAIW(I) + ISHIFT + 3
  90  CONTINUE
      ENDIF
      RETURN
      END
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
C**********************************************
C**********************************************
      SUBROUTINE MA41O (N, NZ, MTRANS, PERM,
     *      IRN, ICN, IW, LIW, ICNTL, INFO)
CCCC*************************************************************
C**********************************
      INTEGER N, NZ, LIW, PERM(N), MTRANS
      INTEGER IRN(NZ),ICN(NZ), IW(LIW)
      INTEGER ICNTL(20), INFO(20)
C***** working variables
      INTEGER MPRINT,LP, MP
      INTEGER ISPW, JPERM
      INTEGER NZREAL, NUMNZ, I, J, JPOS, K
      INTEGER PLENR, IP, IRNW
      LOGICAL PROK, IDENT
      EXTERNAL MC21A
CC***************************************************************
      MPRINT = ICNTL(3)
      LP     = ICNTL(1)
      MP     = ICNTL(2)
      PROK = (MPRINT.GE.0)
      IF (PROK) WRITE(MPRINT,101)
 101    FORMAT(/'****** Preprocessing of original matrix '/)
CC
       IF (N.EQ.1) THEN
        MTRANS=0
        GO TO 500
       ENDIF
       PLENR = 1
       IP    = PLENR + N
       ISPW  = IP    + N
       IRNW  = ISPW  + 4*N
       IF (PROK) THEN
         WRITE(MPRINT,'(A)')
     *     'Compute column permutation'
         WRITE(MPRINT,'(A)') 'Permuted matrix has no zeros on diagonal'
       ENDIF
       DO 5 J=PLENR,PLENR+N-1
        IW(J) = 0
  5    CONTINUE
       NZREAL = 0
       DO 10 K=1,NZ
         I = IRN(K)
         J = ICN(K)
         IF ((J.LE.N).AND.(J.GE.1).AND.
     *           (I.LE.N).AND.(I.GE.1)) THEN
             IW(PLENR+J-1) = IW(PLENR+J-1) + 1
             NZREAL        = NZREAL+1
         ENDIF
  10   CONTINUE
       IW(IP)   = 1
       IW(ISPW) = 1
       DO 20 J=1,N-1
        IW(IP+J)   = IW(IP+J-1)+IW(PLENR+J-1)
        IW(ISPW+J) = IW(IP+J)
  20   CONTINUE
       DO 30 K=1,NZ
         I = IRN(K)
         J = ICN(K)
         IF ((J.LE.N).AND.(J.GE.1).AND.
     *           (I.LE.N).AND.(I.GE.1)) THEN
          JPOS            = IW(ISPW+J-1)
          IW(IRNW+JPOS-1) = I
          IW(ISPW+J-1)    = IW(ISPW+J-1) + 1
         ENDIF
  30   CONTINUE
       CALL MC21A(N,IW(IRNW),NZREAL,IW(IP),IW(PLENR),PERM,NUMNZ,
     *        IW(ISPW))
       IF (NUMNZ.LT.N) GO TO 400
       IDENT = .TRUE.
       DO 80 J=1,N
        JPERM = PERM(J)
        IW(ISPW+JPERM-1) = J
        IF (JPERM.NE.J) IDENT = .FALSE.
 80    CONTINUE
       IF (IDENT) MTRANS=0
       IF (MTRANS.EQ.1) THEN
        IF (LIW.LT.2*NZ+11*N+1) THEN
         INFO(1) = -7
         INFO(2) = 2*NZ+12*N+1
         GO TO 500
        ENDIF
        DO 100 K=1,NZ
         J = ICN(K)
         IF ((J.LE.0).OR.(J.GT.N)) GO TO 100
         ICN(K) = IW(ISPW+ICN(K)-1)
 100    CONTINUE
       ENDIF
       IF (PROK .AND. (MTRANS.NE.1)) THEN
        WRITE (MPRINT,'(A)')
     *  'Column permutation of original matrix is the identity matrix'
       ENDIF
      GO TO 500
C=========================
C** Error return message
C=========================
 400  IF ((LP.GE.0).AND.(ICNTL(4).GE.1))
     *   WRITE (LP,'(/A)') '** Error: Matrix is structurally singular'
      INFO(1) = -6
      INFO(2) = NUMNZ
 500  RETURN
      END
CCC
      SUBROUTINE MA41P(N, R, W)
      INTEGER N, I
      REAL  R(N), W(N)
      DO 100 I = 1, N
        R(I) = R(I) * W(I)
 100  CONTINUE
      RETURN
      END
CCC
      SUBROUTINE MA41Q(MTYPE,IFLAG,N,NZ,ASPK,IRN,ICN,LHS,WRHS,
     *           W,RHS,GIVSOL,SOL,
     *           ANORM, XNORM, SCLNRM, MPRINT, CNTL, ICNTL)
      INTEGER MTYPE,N,NZ,IFLAG,ICNTL(20)
      INTEGER IRN(NZ),ICN(NZ)
      REAL ASPK(NZ),RHS(N),LHS(N),CNTL(10)
      REAL WRHS(N),SOL(*)
      REAL ANORM,RESMAX,RESL2,
     *     XNORM,ZERO,
     *     ERMAX,ERL2,MAXSOL,EPSI,
     *     COMAX,ERREL,SCLNRM,W(N)
      LOGICAL GIVSOL,PROK
      INTEGER MPRINT, MP
      INTEGER I,J,K
      INTRINSIC ABS, MAX, SQRT
      PARAMETER (ZERO = 0.0E0)
      MP = ICNTL(2)
      PROK = (MPRINT.GE.0)
      EPSI = CNTL(2)
      DO 10 K=1,N
      W(K)    = ZERO
      RHS(K)  = WRHS(K)
   10 CONTINUE
      IF (MTYPE.EQ.1) THEN
        DO 20 K=1,NZ
        I     = IRN(K)
        J     = ICN(K)
        IF ((I.LE.0).OR.(I.GT.N).OR.
     *      (J.LE.0).OR.(J.GT.N)) GO TO 20
        RHS(I)= RHS(I)-ASPK(K)*LHS(J)
        W(I)  = W(I)+ABS(ASPK(K))
   20   CONTINUE
      ELSE
        DO 25 K=1,NZ
          I     = IRN(K)
          J     = ICN(K)
          IF (I.LE.0 .OR. I.GT.N .OR.
     *        J.LE.0 .OR. J.GT.N) GO TO 25
          RHS(J)= RHS(J)-ASPK(K)*LHS(I)
          W(J)  = W(J)+ABS(ASPK(K))
   25   CONTINUE
      ENDIF
      ANORM=ZERO
      RESMAX= ZERO
      RESL2 = ZERO
      DO 65 K=1,N
        RESMAX= MAX(RESMAX,ABS(RHS(K)))
        RESL2 = RESL2+RHS(K)*RHS(K)
        ANORM = MAX(ANORM,W(K))
   65 CONTINUE
      XNORM =ZERO
      DO 67 K=1,N
        XNORM  = MAX(XNORM,ABS(LHS(K)))
   67 CONTINUE
      IF (XNORM.GT.EPSI) THEN
        SCLNRM=RESMAX/(ANORM*XNORM)
      ELSE
       IFLAG = IFLAG+ 2
       IF (( MP.GE.0).AND.(ICNTL(4).GE.2))
     * WRITE(MP,'(/A)') 'Max-norm of computed solution is zero'
       SCLNRM = RESMAX/ANORM
      ENDIF
      RESL2=SQRT(RESL2)
      ERMAX=ZERO
      COMAX=ZERO
      ERL2=ZERO
      IF (.NOT.GIVSOL) THEN
        IF (PROK) WRITE(MPRINT,104) RESMAX,RESL2,ANORM,XNORM,SCLNRM
      ELSE
        MAXSOL=ZERO
        DO 68 K=1,N
          MAXSOL = MAX(MAXSOL,ABS(SOL(K)))
   68   CONTINUE
        DO 70 K=1,N
          ERL2=(LHS(K)-SOL(K))**2+ERL2
          ERMAX=MAX(ERMAX,ABS(LHS(K)-SOL(K)))
   70   CONTINUE
        DO 71 K=1,N
          IF (ABS(SOL(K)).GT.EPSI) THEN
            COMAX=MAX(COMAX,(ABS(LHS(K)-SOL(K))/
     *                      ABS(SOL(K))))
          ENDIF
   71   CONTINUE
        ERL2=SQRT(ERL2)
        IF (MAXSOL.GT.EPSI) THEN
          ERREL = ERMAX/MAXSOL
        ELSE
          IFLAG = IFLAG+2
          IF ((MP.GE.0).AND.(ICNTL(4).GE.2))
     *         WRITE(MP,'(/A)') 'Max-norm of exact solution is zero'
          ERREL = ERMAX
        ENDIF
        IF (PROK) WRITE(MPRINT,103)
     *                 ERMAX,ERL2,ERREL,COMAX,RESMAX,RESL2,ANORM,
     *                 XNORM,SCLNRM
      ENDIF
 103   FORMAT (/'Error is     ............ (max-norm)       =',1PE9.2/
     *       '             ............ (2-norm)         =',1PE9.2/
     *       'Relative error........... (max-norm)       =',1PE9.2/
     *       'Component wise error..... (max-norm)       =',1PE9.2/
     *       'and residual is ......... (max-norm)       =',1PE9.2/
     *       '                       .. (2-norm)         =',1PE9.2/
     *       'Norm of input  matrix ... (max-norm)       =',1PE9.2/
     *       'Norm of computed soln ... (max-norm)       =',1PE9.2/
     *       'Scaled residual ......... (max-norm)       =',1PE9.2)
 104   FORMAT (/
     *       'Residual is ............        (max-norm) =',1PE9.2/
     *       '                      ..        (2-norm)   =',1PE9.2/
     *       'RINFO(4):norm of input  matrix  (max-norm) =',1PE9.2/
     *       'RINFO(5):norm of computed soln  (max-norm) =',1PE9.2/
     *       'RINFO(6):scaled residual ...... (max-norm) =',1PE9.2)
      RETURN
      END
CCC
      SUBROUTINE MA41R(N, A, LA, IW, LIW, W, MAXFRT, RHS,
     * PTLUST, NBLK)
      INTEGER N,LA,LIW,MAXFRT,NBLK
      REAL A(LA), RHS(N), W(MAXFRT)
      INTEGER PTLUST(NBLK)
      INTEGER   IW(LIW)
      INTEGER APOS
      INTEGER IBLK,IPOS,LIELL,IFR,JJ
      INTEGER J1,J2,J3,J,K,APOS1
      INTEGER NPIV,NCB,IRHS,IPIV
      REAL W1
      REAL ALPHA,ONE
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0)
      EXTERNAL SGEMV, STRSV
      DO 70 IBLK=1,NBLK
        IPOS  = PTLUST(IBLK) + 2
        LIELL = IW(IPOS)
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
        IPOS = IPOS + 1
        APOS = IW(IPOS)
        J1   = IPOS + 1
        J2   = IPOS + LIELL
        J3   = IPOS + NPIV
        IF (NPIV.EQ.0) GO TO 70
        IF ((NPIV.GE.10).AND.(LIELL.GE.16)) THEN
          IFR = 0
          DO 10 JJ=J1,J2
            J = IW(JJ)
            IFR = IFR + 1
            W(IFR) = RHS(J)
   10     CONTINUE
          CALL STRSV('U','T','U', NPIV, A(APOS), LIELL,
     *           W(1),1)
          IF (NPIV.LT.LIELL) THEN
             APOS1 = APOS + NPIV*LIELL
             NCB   = LIELL - NPIV
             CALL SGEMV('T', NPIV, NCB, ALPHA, A(APOS1),NPIV,
     *           W(1), 1, ONE, W(NPIV+1), 1)
          ENDIF
          IFR = 0
          DO 20 JJ=J1,J2
            J = IW(JJ)
            IFR = IFR + 1
            RHS(J) = W(IFR)
   20     CONTINUE
        ELSE
          APOS = APOS + LIELL
          DO 40 IPIV = 1, NPIV
            IRHS = IW(J1)
            W1   = RHS(IRHS)
            J1   = J1 + 1
            K    = APOS
            IF (J1.LE.J3) THEN
             DO 30 J=J1,J3
              IRHS      = IW(J)
              RHS(IRHS) = RHS(IRHS) - A(K)*W1
              K         = K + LIELL
   30        CONTINUE
            ENDIF
            IF (J3.LT.J2) THEN
              DO 35 J=J3+1,J2
               IRHS      = IW(J)
               RHS(IRHS) = RHS(IRHS) - A(K)*W1
               K         = K + NPIV
   35         CONTINUE
            ENDIF
            APOS = APOS + LIELL + 1
   40     CONTINUE
        ENDIF
   70 CONTINUE
      RETURN
      END
CCC
      SUBROUTINE MA41S(N, A, LA, IW, LIW, W, MAXFRT, RHS,
     *    PTLUST, NBLK,W2)
      INTEGER N,LA,LIW,MAXFRT,NBLK
      REAL A(LA), RHS(N), W(MAXFRT), W2(N)
      INTEGER PTLUST(NBLK)
      INTEGER   IW(LIW)
      INTEGER APOS,APIV,NPIV,IPIV,K
      INTEGER IBLK,IPOS,LIELL,IFR,JJ
      INTEGER J1,J2,J,IST,NCB,JPOS,IPSPIV,IIRHS,IRHS
      REAL W1
      REAL ALPHA,ONE
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0)
      EXTERNAL SGEMV, STRSV
      DO 110 IBLK=NBLK,1,-1
        IPOS  = PTLUST(IBLK)+2
        LIELL = IW(IPOS)
        IPOS  = IPOS + 1
        NPIV  = IW(IPOS)
        IF (NPIV.EQ.0) GO TO 110
        IPOS = IPOS + 1
        APOS = IW(IPOS)
        IF ((NPIV.GE.10).AND.(LIELL.GE.16)) THEN
          J1 = IPOS + 1
          J2 = IPOS + NPIV
          IFR = 0
          DO 10 JJ=J1,J2
            J = IW(JJ)
            IFR = IFR + 1
            W(IFR) = RHS(J)
   10     CONTINUE
          IF (LIELL.GT.NPIV) THEN
           J1 = IPOS + LIELL + NPIV + 1
           J2 = IPOS + 2*LIELL
           DO 20 JJ=J1,J2
             J = IW(JJ)
             IFR = IFR + 1
             W(IFR) = W2(J)
   20      CONTINUE
           NCB   = LIELL - NPIV
           IST   = APOS + NPIV
C**
           CALL SGEMV('T',NCB,NPIV, ALPHA, A(IST),LIELL,
     *            W(NPIV+1),1,ONE, W(1), 1)
          ENDIF
C**
          CALL STRSV('L','T','N',NPIV,A(APOS),LIELL,
     *           W(1),1)
          IFR = 0
          J1 = IPOS + LIELL + 1
          J2 = IPOS + LIELL + NPIV
          DO 30 JJ=J1,J2
            J     = IW(JJ)
            IFR   = IFR + 1
            W2(J) = W(IFR)
   30     CONTINUE
        ELSE
          JPOS = IPOS + NPIV + LIELL
          J2   = IPOS + LIELL*2
          APIV = APOS + (LIELL+1)*(NPIV-1)
          DO 50 IPIV=NPIV,1,-1
           IPSPIV = JPOS - LIELL
           IIRHS  = IW(IPSPIV)
           W1     = RHS(IIRHS)
           J1     = JPOS + 1
           IF (J1.LE.J2) THEN
            K = APIV + 1
            DO 40 J=J1,J2
             IRHS = IW(J)
             W1   = W1 - A(K)*W2(IRHS)
             K    = K + 1
   40       CONTINUE
           ENDIF
           IIRHS     = IW(JPOS)
           W2(IIRHS) = W1/A(APIV)
           JPOS      = JPOS - 1
           APIV      = APIV - LIELL - 1
   50     CONTINUE
        ENDIF
  110 CONTINUE
  120 RETURN
      END
CCC
      SUBROUTINE MA41T(N, A, LA, IW, LIW, W, MAXFRT, RHS, PTLUST,
     *                  NBLK)
      INTEGER N,LA,LIW,MAXFRT,NBLK
      REAL A(LA), RHS(N), W(MAXFRT)
      INTEGER PTLUST(NBLK)
      INTEGER   IW(LIW)
      INTEGER APOS, APOS1, NCB
      INTEGER IBLK,IPOS,LIELL,IFR,JJ
      INTEGER J1,J2,J
      INTEGER NPIV
      REAL ALPHA
      REAL ONE
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0)
      EXTERNAL SGEMV, STRSV
      DO 70 IBLK=1,NBLK
        IPOS = PTLUST(IBLK) + 2
        LIELL= IW(IPOS)
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
        IPOS = IPOS + 1
        APOS = IW(IPOS)
        J1 = IPOS + LIELL + 1
        J2 = J1 + LIELL -1
        IF (NPIV.EQ.0) GO TO 70
        IFR = 0
        DO 10 JJ=J1,J2
          J = IW(JJ)
          IFR = IFR + 1
          W(IFR) = RHS(J)
   10   CONTINUE
        CALL STRSV('L','N','N', NPIV, A(APOS), LIELL,
     *         W(1),1)
        IF (NPIV.LT.LIELL) THEN
           APOS1 = APOS + NPIV
           NCB   = LIELL - NPIV
           CALL SGEMV('N', NCB, NPIV, ALPHA, A(APOS1), LIELL,
     *         W(1), 1, ONE, W(NPIV+1), 1)
        ENDIF
        IFR = 0
        DO 40 JJ=J1,J2
          J = IW(JJ)
          IFR = IFR + 1
          RHS(J) = W(IFR)
   40   CONTINUE
   70 CONTINUE
      RETURN
      END
CCC
      SUBROUTINE MA41U(N, A, LA, IW, LIW, W, MAXFRT, RHS,
     *        PTLUST, NBLK, W2)
      INTEGER N,LA,LIW,MAXFRT,NBLK
      REAL A(LA), RHS(N), W(MAXFRT), W2(N)
      INTEGER PTLUST(NBLK)
      INTEGER   IW(LIW)
      INTEGER APOS, NCB
      INTEGER IBLK,IPOS,LIELL,IFR,JJ
      INTEGER J1,J2,J,IST,NPIV
      REAL ALPHA
      REAL ONE
      PARAMETER (ONE=1.0E0, ALPHA=-1.0E0)
      EXTERNAL SGEMV, STRSV
      DO 110 IBLK=NBLK,1,-1
        IPOS = PTLUST(IBLK)+2
        LIELL = IW(IPOS)
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
        IF (NPIV.EQ.0) GO TO 110
        IPOS = IPOS + 1
        APOS = IW(IPOS)
        J2 = IPOS + LIELL
        J1 = IPOS + 1 + LIELL
        J2 = IPOS + NPIV + LIELL
        IFR = 0
        DO 10 JJ=J1,J2
          J = IW(JJ)
          IFR = IFR + 1
          W(IFR) = RHS(J)
   10   CONTINUE
        IF (LIELL.GT.NPIV) THEN
          J1 = IPOS + NPIV + 1
          J2 = IPOS + LIELL
          DO 20 JJ=J1,J2
            J = IW(JJ)
            IFR = IFR + 1
            W(IFR) = W2(J)
   20     CONTINUE
          NCB   = LIELL - NPIV
          IST   = APOS + NPIV*LIELL
C**
          CALL SGEMV('N',NPIV,NCB, ALPHA, A(IST),NPIV,
     *           W(NPIV+1),1,ONE, W(1), 1)
        ENDIF
C**
        CALL STRSV('U','N','U',NPIV,A(APOS),LIELL,
     *         W(1),1)
        IFR = 0
        J1 = IPOS + 1
        J2 = IPOS + NPIV
        DO 70 JJ=J1,J2
          J = IW(JJ)
          IFR = IFR + 1
          W2(J) = W(IFR)
   70   CONTINUE
  110 CONTINUE
      RETURN
      END
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
        SUBROUTINE MA41V (MTYPE, ASPK, NZ, N, IRN, ICN,
     *    RHS, RHSO, A, LA, IW, LIW, W1, IW1, PTLUST,
     *    NSTEPS, MAXFRT, RINFO,
     *    LSCAL, COLSCA, ROWSCA, CNTL, ICNTL, INFO)
      INTEGER INFO(20)
      INTEGER MTYPE, NZ, N, LA,LIW,NSTEPS,MAXFRT
      INTEGER PTLUST(NSTEPS), IW(LIW)
      INTEGER IRN(NZ), ICN(NZ), IW1(2*N), ICNTL(20)
      REAL A(LA), ASPK(NZ), RHS(N), CNTL(10),
     *     RHSO(N), W1(5*N),
     *     RINFO(20), COLSCA(*), ROWSCA(*), SOL(1)
      LOGICAL LSCAL
      LOGICAL GIVSOL, LPROK, LCOND(2)
      REAL COND(2), DXIMAX, DXMAX, OM1, OLDOMG(2)
      INTEGER JOB,NOITER,JUMP,IFLAG,BDKEEP(4)
      INTEGER MPRINT, LP, MP, NITREF
      INTEGER I,KASE,K,SOLVET
      REAL ONE
      PARAMETER (ONE=1.0E0)
      EXTERNAL MA41C, MA41Q, MA41W
      MPRINT = ICNTL(3)
      LP     = ICNTL(1)
      MP     = ICNTL(2)
      NITREF = ICNTL(10)
      JOB    = ICNTL(11)
      LPROK  = (LP.GT.0)
        IF ((MPRINT.GE.0).AND.(NITREF.GT.0)) THEN
              WRITE(MPRINT,99992)
              WRITE(MPRINT,99995)
     *    'Maximum number of steps                   =',NITREF
        ENDIF
        DO 10 I=1,N
          W1(N+I) = ONE
  10    CONTINUE
        KASE   = 0
 1012   CONTINUE
         IF (MTYPE.EQ.1) THEN
            CALL MA41W(ASPK, NZ, N, IRN, ICN, RHS, RHSO,
     *             W1(1), W1(N+1), W1(2*N+1), IW1, KASE,
     *             RINFO(7), RINFO(9), JOB,
     *             COND, NITREF, NOITER, LCOND(1), LCOND(2),
     *             JUMP, DXIMAX, DXMAX, OM1, OLDOMG, IFLAG,
     *             BDKEEP)
         ELSE
            CALL MA41W(ASPK, NZ, N, ICN, IRN, RHS, RHSO,
     *             W1(1), W1(N+1), W1(2*N+1), IW1, KASE,
     *             RINFO(7), RINFO(9), JOB,
     *             COND, NITREF, NOITER, LCOND(1), LCOND(2),
     *             JUMP, DXIMAX, DXMAX, OM1, OLDOMG, IFLAG,
     *             BDKEEP)
         ENDIF
         IF (KASE.GT.0) THEN
           IF (MTYPE.EQ.1) THEN
             SOLVET = KASE -1
           ELSE
             SOLVET = KASE
           ENDIF
           IF (SOLVET.EQ.1) THEN
            IF (LSCAL) THEN
             DO 555 K=1,N
              W1(K)     = W1(K)*ROWSCA(K)
 555         CONTINUE
            ENDIF
            CALL MA41C(N,A,LA,IW,LIW,W1(2*N+1),MAXFRT,
     1          W1(1),PTLUST,NSTEPS,W1(3*N+1),SOLVET,ICNTL)
            IF (LSCAL) THEN
             DO 556 K=1,N
              W1(K)     = W1(K)*COLSCA(K)
 556         CONTINUE
            ENDIF
           ELSE
            IF (LSCAL) THEN
             DO 557 K=1,N
              W1(K)     = W1(K)*COLSCA(K)
 557         CONTINUE
            ENDIF
            CALL MA41C(N,A,LA,IW,LIW,W1(2*N+1),MAXFRT,
     1         W1(1),PTLUST,NSTEPS,W1(3*N+1),SOLVET,ICNTL)
            IF (LSCAL) THEN
             DO 558 K=1,N
              W1(K)     = W1(K)*ROWSCA(K)
 558         CONTINUE
            ENDIF
           ENDIF
           GO TO 1012
         ELSEIF (KASE.LT.0) THEN
           INFO(1) = INFO(1) + 8
         ENDIF
C-----------------------------------------------
       INFO(15) = NOITER
       IF ((MPRINT.GE.0).AND.(NITREF.GT.0)) THEN
          WRITE(MPRINT,99994)
          WRITE(MPRINT,99995)
     *    'Number of steps of iterative refinement =',NOITER
          WRITE(MPRINT,99993)
       ENDIF
       IF (ICNTL(11).GT.0) THEN
C==========================
C==========================
        IF ((MPRINT.GE.0).AND.(NITREF.GT.0)) WRITE(MPRINT,99990)
        IF ((MPRINT.GE.0).AND.(NITREF.LE.0)) WRITE(MPRINT,99991)
        GIVSOL = .FALSE.
        CALL MA41Q(MTYPE,INFO(1),N,NZ,ASPK,IRN,ICN,RHSO,RHS,
     *     W1(1),W1(N+1),GIVSOL,SOL,RINFO(4),RINFO(5),RINFO(6),
     *     MPRINT, CNTL, ICNTL)
        IF (MPRINT.GE.0) THEN
          WRITE(MPRINT, 2100)
     *    'RINFO(7):componentwise scaled residual(w1) =',RINFO(7)
          WRITE(MPRINT, 2100)
     *    '-----(8):---------------------------- (w2) =',RINFO(8)
          WRITE(MPRINT, 2100)
     *    '-----(9):Upper bound error ................=',RINFO(9)
          WRITE(MPRINT, 2100)
     *    'Condition number (1) ......................=',COND(1)
          WRITE(MPRINT, 2100)
     *    'Condition number (2) ......................=',COND(2)
        ENDIF
C==========================
C==========================
       ENDIF
CCC
CCC
 500  RETURN
 2100 FORMAT(A44,1P,E9.2)
99990 FORMAT (//'Error analysis after iterative refinement')
99991 FORMAT (//'Error analysis')
99992 FORMAT (//'Begin iterative refinement')
99993 FORMAT (/'End   iterative refinement')
99994 FORMAT (/'Statistics after iterative refinement')
99995 FORMAT(A42,I4)
      END
CCC
       SUBROUTINE MA41W ( A, NZ, N, IRN, ICN, RHS, X,
     *                    Y, D, W, IW, KASE, OMEGA, ERX, JOB,
     *                    COND, MAXIT, NOITER, LCOND1, LCOND2,
     *                    JUMP,  DXIMAX, DXMAX, OM1, OLDOMG, IFLAG,
     *                    BDKEEP)
      INTEGER NZ, N, KASE, BDKEEP(4)
      INTEGER IW(N,2), JOB
      INTEGER IRN(NZ), ICN(NZ)
      REAL A(NZ), RHS(N)
      REAL X(N), Y(N), D(N), W(N,3)
      INTEGER MAXIT, NOITER
      REAL COND(2), OMEGA(2)
      REAL CGCE, CTAU
      PARAMETER (CTAU = 1.0E3, CGCE=0.50E0)
      LOGICAL LCOND1, LCOND2
      INTEGER IFLAG, JUMP, I, IMAX
      REAL ERX, DXMAX, DXIMAX
      REAL DD, OM1, OM2,
     *     TAU, ZERO, ONE
      PARAMETER (ZERO=0.0E0, ONE = 1.0E0)
      REAL OLDOMG(2)
      INTEGER ISAMAX
      EXTERNAL MA41X, MA41Y, MA41P, MC51B, ISAMAX
      INTRINSIC     ABS, MAX
      IF ( KASE .EQ. 0 ) THEN
        LCOND1  = .FALSE.
        LCOND2  = .FALSE.
        COND(1) = ONE
        COND(2) = ONE
        ERX    = ZERO
        OM1    = ZERO
        IFLAG  = 0
        NOITER = 0
        JUMP   = 1
        CALL MA41X ( A, NZ, N, IRN, W(1,2) )
      ENDIF
      GO TO ( 1000, 2000, 3000, 4000 ) JUMP
 2000   CONTINUE
        DO 110 I=1,N
          X(I) = X(I) + Y(I)
  110   CONTINUE
      IF ( NOITER .GT. MAXIT ) THEN
        IFLAG  = IFLAG+8
        GO TO 150
      ENDIF
 1000 CONTINUE
        IMAX = ISAMAX ( N, X, 1 )
        DXMAX = ABS(X(IMAX))
        CALL MA41Y ( A, NZ, N, IRN, ICN, RHS, X, Y, W )
        OMEGA(1) = ZERO
        OMEGA(2) = ZERO
        DO 90 I=1,N
          TAU = ( W(I,2)*DXMAX + ABS(RHS(I)) ) * N * CTAU
          DD  = W(I,1) + ABS(RHS(I))
          IF ( ( DD + TAU ) .GT. TAU ) THEN
            OMEGA(1) = MAX(OMEGA(1),ABS(Y(I)/DD))
            IW(I,1) = 1
          ELSE
            IF ( TAU .GT. ZERO ) THEN
              OMEGA(2) = MAX(OMEGA(2),ABS(Y(I)/(DD+W(I,2)*DXMAX)))
            ENDIF
            IW(I,1) = 2
          ENDIF
   90   CONTINUE
        OM2 = OMEGA(1) + OMEGA(2)
        IF ( ( OM2 + ONE ) .LE. ONE ) GO TO 150
        IF ( MAXIT .EQ. 0 ) GO TO 150
        IF ( NOITER .GT. 1 .AND. OM2 .GT. OM1*CGCE ) THEN
          IF ( OM2 .GT. OM1 ) THEN
            OMEGA(1) = OLDOMG(1)
            OMEGA(2) = OLDOMG(2)
            DO 130 I = 1, N
              X(I) = W(I,3)
  130       CONTINUE
          ENDIF
          GO TO 150
        ENDIF
        DO 120 I = 1, N
          W(I,3) = X(I)
  120   CONTINUE
        OLDOMG(1) = OMEGA(1)
        OLDOMG(2) = OMEGA(2)
        OM1 = OM2
        NOITER = NOITER + 1
        KASE   = 2
        JUMP   = 2
        RETURN
  150 KASE = 0
      IF ( JOB .LE. 0 ) GO TO 160
      DO 100 I = 1, N
        IF ( IW(I,1) .EQ. 1 ) THEN
          W(I,1) = W(I,1) + ABS( RHS(I) )
          W(I,2) = ZERO
          LCOND1 = .TRUE.
        ELSE
          W(I,2) = W(I,2) * DXMAX + W(I,1)
          W(I,1) = ZERO
          LCOND2 = .TRUE.
        ENDIF
  100 CONTINUE
      DO 200 I = 1, N
          W(I,3) = X(I) * D(I)
  200 CONTINUE
      IMAX = ISAMAX ( N, W(1,3), 1 )
      DXIMAX = ABS(W(IMAX,3))
      IF ( .NOT. LCOND1 ) GO TO 155
 300    CALL MC51B ( N, KASE, Y, COND(1), W(1,3), IW(1,2),
     *               BDKEEP(1),  BDKEEP(2), BDKEEP(3), BDKEEP(4) )
        IF (KASE .EQ. 0) GO TO 500
        IF (KASE .EQ. 1)  CALL MA41P ( N, Y, D )
        IF (KASE .EQ. 2)  CALL MA41P ( N, Y, W )
          JUMP = 3
          RETURN
 3000     CONTINUE
        IF (KASE .EQ. 1)  CALL MA41P ( N, Y, W )
        IF (KASE .EQ. 2)  CALL MA41P ( N, Y, D )
        GO TO 300
 500  IF ( DXIMAX .GT. ZERO )    COND(1) = COND(1) / DXIMAX
      ERX = OMEGA(1) * COND(1)
 155  IF ( .NOT. LCOND2 ) GO TO 160
      KASE = 0
 400    CALL MC51B ( N, KASE, Y, COND(2), W(1,3), IW(1,2),
     *               BDKEEP(1),  BDKEEP(2), BDKEEP(3), BDKEEP(4) )
        IF (KASE .EQ. 0) GO TO  800
        IF (KASE .EQ. 1)  CALL MA41P ( N, Y, D )
        IF (KASE .EQ. 2)  CALL MA41P ( N, Y, W(1,2) )
          JUMP = 4
          RETURN
 4000     CONTINUE
        IF (KASE .EQ. 1)  CALL MA41P ( N, Y, W(1,2) )
        IF (KASE .EQ. 2)  CALL MA41P ( N, Y, D )
        GO TO 400
  800 IF ( DXIMAX .GT. ZERO )    THEN
        COND(2) = COND(2) / DXIMAX
      ENDIF
      ERX = ERX + OMEGA(2) * COND(2)
  160 KASE = -IFLAG
      RETURN
      END
CCC
      SUBROUTINE MA41X(A, NZ, N, IRN, Z)
      INTEGER NZ, N, I, K
      INTEGER IRN(NZ)
      REAL A(NZ), Z(N), ZERO
      PARAMETER (ZERO = 0.0E0)
      INTRINSIC     ABS
      DO 100 I = 1, N
        Z(I) = ZERO
  100 CONTINUE
      DO 300 K = 1, NZ
        I = IRN(K)
        IF (I.GE.1 .AND. I.LE.N) Z(I) = Z(I) + ABS(A(K))
  300 CONTINUE
      RETURN
      END
CCC
      SUBROUTINE MA41Y(A, NZ, N, IRN, ICN, RHS, X, R, W)
      INTEGER NZ, N, I, K, J
      INTEGER IRN(NZ), ICN(NZ)
      REAL A(NZ), RHS(N), X(N), R(N), W(N), D, ZERO
      INTRINSIC ABS
      PARAMETER (ZERO = 0.0E0)
      DO 100 I=1,N
        R(I) = RHS(I)
        W(I) = ZERO
  100 CONTINUE
      DO 300 K = 1, NZ
        I = IRN(K)
        J = ICN(K)
        IF (I.GT.N .OR. J.GT.N .OR. I.LT.1 .OR. J.LT.1) GO TO 300
        D = A(K)*X(J)
        R(I) = R(I) - D
        W(I) = W(I) + ABS(D)
  300 CONTINUE
      RETURN
      END
CCC
C--------------------------------------------------------------------
C-             Copyright Rutherford Appleton Laboratory
C--------------------------------------------------------------------
      SUBROUTINE MA41Z(N, NZ, NSTEPS, A, LA, IW, LIW, NA, NE,
     * NSTK, MAXFRT, IFLAG, IERROR, ND, FILS, FRERE, PTRARW, PTRAIW,
     * PTRAST, PTRIST, PTLUST, TLKJ, PERMW,NPROCS, ITLOC,
     * NIRBDU, IPTA, IPOOL, LPOOL, IPOOLB,
     * UULOC, ICNTL, INFO, RINFO, KEEP)
C*******************************
      INTEGER NBUD
      PARAMETER(NBUD=29)
      INTEGER N, LIW, LA, IPTA(NBUD)
      INTEGER NZ, NSTEPS, LPOOL, KEEP(50)
      REAL A(LA)
      INTEGER IW(LIW)
      INTEGER NA(N), NE(N), NSTK(N), MAXFRT,
     *        IFLAG, IERROR, ND(N), FILS(N), FRERE(N),
     *        PTRARW(N), PTRAIW(N),PTLUST(NSTEPS),
     *        PTRAST(N), PTRIST(N), TLKJ(N), PERMW(N),
     *        NPROCS, ITLOC(N,NPROCS)
      INTEGER INFO(20), IPOOL(LPOOL), IPOOLB(N),
     *        NIRBDU
      REAL    RINFO(20)
      INTEGER ICNTL(20)
      REAL UULOC
      INTEGER NPROCW, ASTK, ISTEP, PERLU, NRLBDU
      INTEGER LENBUD,IPTEND
      INTEGER MPRINT,LP,MMP
      INTEGER I, ISP, LMAX, NBLEAF, MAXSZ, LFRONT, IPT
      INTEGER SIZFB, IPTOLD, IASS, LMIN, IPTINI
      LOGICAL ALOCLU,AGAIN
      REAL ZERO
      INTEGER LCKS(20)
      LOGICAL ENDTAG
      INTEGER NBACTI, NTOTPV, NBNODE
      INTEGER LRLU, IPTRLU, LRLUS,
     *    ISPON1, ISPON2, ISPON3, ISPON4,LENA,
     *    NBLEC, III, IIIB, LEAF, LEAFB, NBROOT
C***************
C***************************************
      INTRINSIC REAL, INT, LOG
      PARAMETER (ZERO=0.0E0)
      EXTERNAL MC51E
CCC
      MPRINT = ICNTL(3)
      LP     = ICNTL(1)
      ASTK = LA - NZ
      IF (KEEP(4).LE.0) KEEP(4)=32
      IF (KEEP(5).LE.0) KEEP(5)=16
      IF (KEEP(6).LE.0) KEEP(6)=24
      IF (KEEP(3).LE.KEEP(4)) KEEP(3)=KEEP(4)*2
      IF (NPROCS.EQ.1) THEN
         ISPON1 = 99999
         ISPON2 = 99999
         ISPON3 = 99999
         ISPON4 = 99999
      ELSE
         ISPON1 = KEEP(7)
         ISPON2 = KEEP(8)
         ISPON3 = KEEP(9)
         ISPON4 = KEEP(10)
         IF (ISPON2.LE.0) ISPON2 = 16
         IF (ISPON4.LE.0) ISPON4 = 16
         IF (ISPON1.LE.ISPON2) ISPON1 = ISPON2*2
         IF (ISPON3.LE.ISPON4) ISPON3 = ISPON4*2
      ENDIF
C*****************************
C*****************************
      ENDTAG = .FALSE.
      NBACTI = 0
      NBNODE = 0
      NPROCW = NPROCS
      INFO(9)   = 1
      INFO(10)  = 1
      MAXFRT = 0
      ALOCLU = .FALSE.
      AGAIN  = .FALSE.
      IF (NPROCS.EQ.1) THEN
C--------------------------------------------------------------------
C--------------------------------------------------------------------
        LRLU     = ASTK
        IPTRLU   = LRLU
        NRLBDU   = LRLU
        PERLU    = 0
        SIZFB  = 0
        LENA   = 0
        LENBUD = 0
      ELSE
C-----------------------------------------------------------------
C------------------------------------------------------------------
      PERLU  = KEEP(12)
      IF ((LA.GT.4*KEEP(30)).AND.(PERLU.NE.0)) THEN
          PERLU = 200
          IF (LA.GT.10*KEEP(30)) PERLU = 400
          PERLU = MAX(PERLU,KEEP(12))
      ENDIF
      NRLBDU   = KEEP(16) + PERLU*(KEEP(16)/100)
      LRLU   = NRLBDU
      IPTINI = NRLBDU
      IPTRLU = IPTINI
      LRLUS  = LRLU
      NBLEC  = 0
      IF (KEEP(2).LE.4) THEN
         LENBUD = -1
      ELSE
         LENBUD = KEEP(2)
      ENDIF
      DO 9 I=1,NBUD
        IPTA(I)=0
   9  CONTINUE
C**********************************
C**********************************
      MAXSZ  = KEEP(27)*KEEP(27)+3
      LENA   =  MAXSZ/16 +3
      ISP    = ASTK - NRLBDU
      NPROCW = MIN0(NPROCW,MAX(1,ISP/MAXSZ))
      IF (LENBUD.GT.4) THEN
        LMAX   = INT(LOG(REAL(LENBUD+1))/LOG(2.0))
        LENBUD = 2**LMAX
        ISP    = MIN0(ASTK - NRLBDU ,LENBUD)
      ENDIF
      IF (ISP.NE.LENBUD) THEN
        LMAX = INT(LOG(REAL(ISP))/LOG(2.0))
      ENDIF
       LMAX = MIN0(NBUD-1,LMAX)
       A(NRLBDU+1) = REAL(LMAX)
       IPTA(LMAX) = NRLBDU+1
       A(NRLBDU+2) = 0
       A(NRLBDU+3) = 0
       LENBUD = 2**LMAX
C********************************
C********************************
       LFRONT =  INT(LOG(REAL(MAXSZ+1))/LOG(2.0))  + 1
       LMIN   = LFRONT+INT(LOG(REAL(NPROCS+1))/LOG(2.0))
       IF ((LENBUD.LT.2*NRLBDU).OR.(LMAX.LT.LMIN+1)) ALOCLU=.TRUE.
       IF (LMAX.LT.LMIN)  THEN
         LENA   = MAX (LENA, (ASTK-NRLBDU-LENBUD)/NPROCS)
         LENA   = MIN0(MAXSZ+KEEP(12)*(MAXSZ/100),LENA)
       ELSE
         IF (LENA*(NSTEPS/2+1).LT.
     *                      (ASTK-NRLBDU-LENBUD))
     *     LENA = (ASTK-NRLBDU-LENBUD)/(NSTEPS/2+1)
       ENDIF
       LENA   = MIN0(ASTK-NRLBDU-LENBUD,LENA)
       IF (LMAX.LT.LMIN) THEN
        IF ((ICNTL(2).GE.0).AND.(ICNTL(4).GE.2))  THEN
         MMP = ICNTL(2)
         WRITE(MMP,'(/A)') 'Warning from factorization phase:'
         WRITE(MMP,'(A)')
     *     'The size of the real working array might be increased'
         WRITE(MMP,'(A)')
     *     'to improve performance'
        ENDIF
        IFLAG  = IFLAG+4
        IERROR = 2**LMIN-2**LMAX+1
        ALOCLU = .TRUE.
        IF (LMAX.LT.LFRONT+1) THEN
         LENBUD = 0
         IPTA(LMAX) = 0
         LMAX = 0
         AGAIN  = .TRUE.
        ENDIF
        IPTEND = NRLBDU + LENBUD + 1
        IPT   = IPTEND
        SIZFB = ASTK - IPT +1
        IF (LENBUD.EQ.0) THEN
         IF (SIZFB.GE.(MAXSZ+KEEP(12)*(MAXSZ/100)) ) THEN
            LENA = MAXSZ+KEEP(12)*(MAXSZ/100)
         ELSE
            LENA = SIZFB
         ENDIF
         LENA   = MAX(LENA,SIZFB/NSTEPS)
         LENA   = MAX(LENA,SIZFB/NPROCS)
         NPROCW = MAX(1,SIZFB/LENA)
         NPROCW = MIN0(NPROCS,NPROCW)
       ENDIF
      ENDIF
      IPTEND     = NRLBDU + LENBUD + 1
      SIZFB      = ASTK - IPTEND +1
      LENA       = MIN(LENA,SIZFB)
      IF (LENA.GT.0) THEN
       IPT   = IPTEND
       IF (IPT+LENA-1.LE.ASTK) THEN
        IPTA(NBUD) = IPT
        IPTOLD = 0
 237    A(IPT) = REAL(LENA)
        A(IPT+1) =REAL(IPTOLD)
        A(IPT+2) = REAL(IPT+LENA)
        IPTOLD = IPT
        IPT = IPT + LENA
        IF (IPT+LENA-1.LE.ASTK) GO TO 237
        A(IPTOLD+2) = 0
       ENDIF
      ENDIF
      IF (UULOC.EQ.ZERO) ALOCLU = .TRUE.
C---------------------------
C---------------------------
      ENDIF
      KEEP(18) = NRLBDU
      KEEP(19) = LENBUD
      KEEP(20) = SIZFB
      KEEP(21) = LENA
C********************************************
      ISTEP  = 0
      NTOTPV = 0
      DO 11 I=1,N
        NSTK(I) = NE(I)
 11   CONTINUE
      IF (N.EQ.1) THEN
        NBROOT = 1
        NBLEAF        = 1
        IPOOL(1)      = 1
      ELSEIF (NA(N).LT.0) THEN
        NBLEAF = N
        NBROOT = N
        IPOOL(N) = -NA(N)-1
        DO 12 I=1,NBLEAF-1
           IPOOL(I) = NA(I)
 12     CONTINUE
      ELSEIF (NA(N-1).LT.0) THEN
        NBLEAF = N-1
        NBROOT = NA(N)
        IPOOL(NBLEAF) = -NA(N-1)-1
        IF (NBLEAF-1.GT.0) THEN
         DO 13 I=1,NBLEAF-1
            IPOOL(I) = NA(I)
 13      CONTINUE
        ENDIF
      ELSE
        NBLEAF = NA(N-1)
        NBROOT = NA(N)
        DO 14 I = 1,NBLEAF
           IPOOL(I) = NA(I)
 14     CONTINUE
      ENDIF
      LEAF   = NBLEAF+1
      LEAFB  = 1
C******************************************
      III = 1
      IIIB= 1
      NPROCW = MIN0(NPROCW,MIN0(ISPON2,ISPON4))
      IF ((MPRINT.GE.0).AND.(NPROCS.GT.NPROCW))
     *                WRITE(MPRINT,99980) NPROCW
C******************************************
      DO 600 IASS=1,NPROCS
        CALL MC51E(NPROCS,NPROCW,N,IW,LIW,A,LA,NSTK,
     *         PERMW,IFLAG,ND,FILS,FRERE,MAXFRT,
     *         NTOTPV,PTRIST,PTRAST,PTRARW,PTRAIW,
     *         IPTA,TLKJ,ALOCLU,AGAIN,ITLOC(1,IASS),
     *         IERROR, IPOOL, LPOOL, IPOOLB, N,
     *         NIRBDU,NRLBDU,LENBUD,
     *         LENA, IPTEND, RINFO,LCKS,
     *         ISPON1,ISPON2,ISPON3,ISPON4,NBNODE,NBACTI,
     *         ENDTAG,INFO(9),INFO(10),LRLU,IPTRLU, IPTINI,
     *         LRLUS, NBLEC, III, IIIB, LEAF, LEAFB, NBROOT,
     *         UULOC,ICNTL,PTLUST,NSTEPS,ISTEP,INFO,KEEP)
 600  CONTINUE
      IF (NBNODE.LE.-1) THEN
         IF (MPRINT.GE.0) THEN
           WRITE(MPRINT,99981) RINFO(2)+RINFO(3)
         ENDIF
         NPROCW = 1
         ISPON1 = 9999
         ISPON2 = 9999
         ISPON3 = 9999
         ISPON4 = 9999
         NBNODE = 0
         CALL MC51E(NPROCS,NPROCW,N,IW,LIW,A,LA,NSTK,
     *         PERMW,IFLAG,ND,FILS,FRERE,MAXFRT,
     *         NTOTPV,PTRIST,PTRAST,PTRARW,PTRAIW,
     *         IPTA,TLKJ,ALOCLU,AGAIN,ITLOC(1,1),
     *         IERROR, IPOOL, LPOOL, IPOOLB, N,
     *         NIRBDU,NRLBDU,LENBUD,
     *         LENA, IPTEND, RINFO,LCKS,
     *         ISPON1,ISPON2,ISPON3,ISPON4,NBNODE,NBACTI,
     *         ENDTAG,INFO(9),INFO(10),LRLU,IPTRLU, IPTINI,
     *         LRLUS, NBLEC, III, IIIB, LEAF, LEAFB, NBROOT,
     *         UULOC,ICNTL,PTLUST,NSTEPS,ISTEP,INFO,KEEP)
      ENDIF
      INFO(9)  = INFO(9)  - 1 + KEEP(22)
      INFO(10) = INFO(10) - 1
      IF (IFLAG.EQ.-10) IERROR = NTOTPV
      IF (IFLAG.LT.0) GOTO 500
      IF (NTOTPV.NE.N) THEN
        IF ((LP.GE.0).AND.(ICNTL(4).GE.1)) THEN
         WRITE(LP,'(A)')
     *     'Error total number of pivots is not equal to'
         WRITE(LP,'(A)')
     *     'the order of the matrix '
        ENDIF
        IFLAG = -10
        IERROR = NTOTPV
      ENDIF
  500 RETURN
C===============
C===============
99980 FORMAT('Maximum number of concurrent nodes            =',I12)
99981 FORMAT('Number of operations before switch: ',1PE10.3)
      END
