C COPYRIGHT (c) 1995 Council for the Central Laboratory
*                    of the Research Councils
C######DATE 30 November 1995
C  April 2001: threadsafe version of MA51
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
       SUBROUTINE MC51A (N, NZ, NSCA, ASPK, IRN, ICN, COLSCA, ROWSCA,
     *                    S, MAXS, ICNTL, INFO)
CCCC*************************************************************
C INPUT (not modified)
C -----
C    IRN and ICN are integer arrays of size NZ
C    IRN(k) and ICN(k), k=1..NZ must be set on entry to hold
C        the row and column indices respectively.
C        Not altered in this routine
C    ASPK is a REAL array of length NZ.
C    The user must set ASPK(k) to the value
C    of the entry in row IRN(k) and column ICN(k).
C    It is not altered by the subroutine.
C
C INPUT/OUPUT
C -----------
C        ROWSCA and COLSCA are arrays or order N that
C                   must be set on entry to one. Contains on output
C                   respectively the row and the column scalings
C                   that will be applied later to the original matrix.
C WORKING ARRAY:
C  S Real working array of size MAXS that need not be
C                       set on input.
C**********************************
      INTEGER N, NZ, NSCA, MAXS
      INTEGER IRN(NZ), ICN(NZ)
      INTEGER ICNTL(20), INFO(20)
      REAL    ASPK(NZ), COLSCA(*), ROWSCA(*)
      REAL    S(MAXS)
C***** working variables
      INTEGER MPRINT,LP, MP
      INTEGER ISPW1, IWNOR
      INTEGER I, K, ITOT
      LOGICAL PROK
      REAL ONE
      PARAMETER (ONE=1.0E0)
      EXTERNAL MC51V, MC51W, MC51X, MC51Y
C***************************************************************
C ASPK is not modified during scaling
C
      LP     = ICNTL(1)
      MP     = ICNTL(2)
      MPRINT = ICNTL(3)
      PROK   = (MPRINT.GE.0)
      IF (PROK) WRITE(MPRINT,101)
 101    FORMAT(/'****** Scaling of original matrix '/)
CC
C  NSCA = 1, 2, 3, 4, 5, 6
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
C     initialize scaling arrays
        DO 10 I=1,N
            COLSCA(I) = ONE
            ROWSCA(I) = ONE
 10     CONTINUE
C
C********************************
C Compute scaling arrays
C
C
C check if enough space
C
        IF ((NSCA.EQ.5).OR.(NSCA.EQ.6)) THEN
          ITOT = 5*N + NZ
          IF (ITOT.GT.MAXS) GOTO 400
C
C         S(ISPW1+K-1) holds the matrix values so that
C         we do not modify ASPK in MC51A:
C
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
C
          IF (NSCA.EQ.1) THEN
            CALL MC51V(N,NZ,ASPK,IRN,ICN,S(IWNOR),
     *        COLSCA,ROWSCA,MPRINT)
          ELSEIF (NSCA.EQ.2) THEN
C
C MC51W calls HARWELL Library Subroutine MC29A
C
            CALL MC51W(N,NZ,ASPK,IRN,ICN,
     *      ROWSCA,COLSCA,S(IWNOR),MPRINT,MP,NSCA)
          ELSEIF (NSCA.EQ.3) THEN
            CALL MC51Y(N,NZ,ASPK,IRN,ICN,S(IWNOR),COLSCA,
     *      MPRINT)
          ELSEIF (NSCA.EQ.4) THEN
            CALL MC51X (N,NZ,IRN,ICN,ASPK,
     *      S(IWNOR),S(IWNOR+N),COLSCA,ROWSCA,MPRINT)
          ELSEIF (NSCA.EQ.5) THEN
            CALL MC51W(N,NZ,S(ISPW1),IRN,ICN,
     *      ROWSCA,COLSCA,S(IWNOR),MPRINT,MP,NSCA)
            CALL MC51Y(N,NZ,S(ISPW1),IRN,ICN,S(IWNOR),COLSCA,
     *          MPRINT)
          ELSEIF (NSCA.EQ.6) THEN
            CALL MC51W(N,NZ,S(ISPW1),IRN,ICN,
     *      ROWSCA,COLSCA,S(IWNOR),MPRINT,MP,NSCA)
            CALL MC51X (N,NZ,IRN,ICN,S(ISPW1),
     *      S(IWNOR),S(IWNOR+N),COLSCA,ROWSCA,MPRINT)
          ENDIF
      GOTO 500
C
C** ERROR return message
 400  INFO(1) = -5
      INFO(2) = ITOT
      IF ((LP.GE.0).AND.(ICNTL(4).GE.1))
     * WRITE(LP,'(/A)') '*** Error: Not enough space to scale matrix'
 500  RETURN
      END
CCC
      SUBROUTINE MC51B (N, KASE, X, EST, W, IW, ITER, J, JLAST, JUMP)
      INTEGER N, IW(N), KASE, ITER, J, JLAST, JUMP
      REAL W(N), X(N), EST
C
C      MC51B ESTIMATES THE 1-NORM OF A SQUARE, REAL
C      MATRIX A.
C      REVERSE COMMUNICATION IS USED FOR EVALUATING
C      MATRIX-VECTOR PRODUCTS.
C
C
C         N       INTEGER
C                 THE ORDER OF THE MATRIX.  N .GE. 1.
C
C         KASE    INTEGER
C                 SET INITIALLY TO ZERO . IF N .LE. 0 SET TO -1
C                 ON INTERMEDIATE RETURN
C                 = 1 OR 2.
C                ON FINAL RETURN
C                 =  0  ,IF SUCCESS
C                 = -1  ,IF N .LE.0
C
C         X       REAL ARRAY OF DIMENSION (N)
C                 IF 1-NORM IS REQUIRED
C                 MUST BE OVERWRITTEN BY
C
C                      A*X,             IF KASE=1,
C                      TRANSPOSE(A)*X,  IF KASE=2,
C
C                 AND MC41 MUST BE RE-CALLED, WITH ALL THE OTHER
C                 PARAMETERS UNCHANGED.
C                 IF INFINITY-NORM IS REQUIRED
C                 MUST BE OVERWRITTEN BY
C
C                      TRANSPOSE(A)*X,  IF KASE=1,
C                      A*X,             IF KASE=2,
C
C                 AND MC41 MUST BE RE-CALLED, WITH ALL THE OTHER
C                 PARAMETERS UNCHANGED.
C
C         EST     REAL
C                 CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C         W       REAL ARRAY OF DIMENSION (N)
C                 = A*V,   WHERE  EST = NORM(W)/NORM(V)
C                          (V  IS NOT RETURNED).
C         IW      INTEGER(N) USED AS WORKSPACE.
C
C
C      REFERENCE
C      N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C      THE ONE-NORM OF A
C      REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C      TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
C      UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.
C
C      SUBROUTINES AND FUNCTIONS
      INTRINSIC REAL, ABS, NINT, SIGN
      INTEGER ISAMAX
C
      INTEGER ITMAX
      PARAMETER (ITMAX = 5)
C
C      INTERNAL VARIABLES
C
      INTEGER I
      REAL ALTSGN, TEMP, ZERO, ONE
      PARAMETER (ZERO=0.0E0, ONE=1.0E0)
      EXTERNAL ISAMAX
C
C     IF (N .LE. 0) THEN
C       KASE = -1
C       RETURN
C     ENDIF
      IF (KASE .EQ. 0) THEN
         DO 10 I = 1,N
            X(I) = ONE/REAL(N)
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      ENDIF
C
      GOTO (100, 200, 300, 400, 500) JUMP
C
C      ................ ENTRY   (JUMP = 1)
C
  100 CONTINUE
      IF (N .EQ. 1) THEN
         W(1) = X(1)
         EST = ABS(W(1))
C         ... QUIT
         GOTO 510
      ENDIF
C
      DO 110 I = 1,N
         X(I) = SIGN(ONE,X(I))
         IW(I) = NINT(X(I))
  110 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
C
C      ................ ENTRY   (JUMP = 2)
C
  200 CONTINUE
      J = ISAMAX(N,X,1)
      ITER = 2
C
C      MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
  220 CONTINUE
      DO 230 I = 1,N
         X(I) = ZERO
  230 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      RETURN
C
C      ................ ENTRY   (JUMP = 3)
C
  300 CONTINUE
C
C      COPY X INTO W
C
      DO 310 I = 1,N
         W(I) = X(I)
  310 CONTINUE
      DO 320 I = 1,N
         IF ( NINT( SIGN(ONE,X(I)) ) .NE. IW(I) ) GOTO 330
  320 CONTINUE
C
C      REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GOTO 410
C
  330 CONTINUE
      DO 340 I = 1,N
         X(I) = SIGN(ONE,X(I))
         IW(I) = NINT(X(I))
  340 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
C
C      ................ ENTRY   (JUMP = 4)
C
  400 CONTINUE
      JLAST = J
      J = ISAMAX(N,X,1)
      IF (   (  ABS(X(JLAST)) .NE. ABS(X(J)) ) .AND.
     +      (ITER .LT. ITMAX)   ) THEN
         ITER = ITER + 1
         GOTO 220
      ENDIF
C
C      ITERATION COMPLETE.  FINAL STAGE.
C
  410 CONTINUE
      EST = ZERO
      DO 420 I = 1,N
         EST = EST + ABS(W(I))
  420 CONTINUE
C
      ALTSGN = ONE
      DO 430 I = 1,N
         X(I) = ALTSGN * (ONE+REAL(I-1)/REAL(N-1))
         ALTSGN = -ALTSGN
  430 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
C
C      ................ ENTRY   (JUMP = 5)
C
  500 CONTINUE
      TEMP = ZERO
      DO 520 I = 1,N
         TEMP = TEMP + ABS(X(I))
  520 CONTINUE
      TEMP = 2.0*TEMP/REAL(3*N)
      IF (TEMP. GT. EST) THEN
C
C      COPY X INTO W
C
         DO 530 I = 1,N
            W(I) = X(I)
  530    CONTINUE
         EST = TEMP
      ENDIF
C
  510 KASE = 0
      RETURN
C
      END
C
CCC
CC
CC
      SUBROUTINE MC51C
C Dummy subroutine ... only used in parallel version
      RETURN
      END
      SUBROUTINE MC51D
C Dummy subroutine ... only used in parallel version
      RETURN
      END
CCC
C--------------------------------------------------------------------
C            HARWELL SUBROUTINE LIBRARY Release 12 (1995)
C        --
C-             Copyright Rutherford Appleton Laboratory
C        --
C--------------------------------------------------------------------
      SUBROUTINE MC51E(NTASKS,NPROCW,N,IW,LIW,A,LA,
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
C
C  .. Purpose ..
C     =======
C Parallel routine doing node factorization.
C Based on a greedy method, we process all nodes of the
C elimination tree toward the root nodes. Two pool of tasks
C are used to give work to the routine.
C Each pool of task is associated  to one level of parallelism
C of the method.
C This routine extracts work from the pool, analyse the work to
C be done and perform the corresponding operations.
C If the task extracted
C corresponds to a node process then the operations involved are :
C assembly step to build the frontal matrix, factorization (may generate
C parallel tasks), stacking of the contribution block and factors, and
C eventually activate father node.
C
C Remark:
C ------
C  Circular pools of tasks are only necessary for parallel execution of
C  of the code.
C  Few other statements can only be used during parallel execution.
C
C
C  .. Parameters ..
C     ==========
      INTEGER NBUD
      PARAMETER (NBUD=29)
      INTEGER NTASKS,NPROCW,N,IFLAG,NTOTPV,MAXFRT,LA,LIW,
     *        IERROR,NIRBDU,IPTBEG,LENA, IPTEND,
     *        NSTEPS, ISTEP, INFO(20)
      INTEGER ISPON1,ISPON2,ISPON3,ISPON4,NBNODE,NBACTI
      REAL A(LA)
      INTEGER TLKJ(N),ITLOC(N), KEEP(50)
      INTEGER IW(LIW), PERM(N), NSTK(N), IPTA(NBUD)
      INTEGER PTRARW(N), PTRAIW(N), ND(N)
      INTEGER FILS(N),FRERE(N),PTRIST(N),PTRAST(N)
      INTEGER RC,POSELT,LENBUD,ICNTL(20),PTLUST(NSTEPS)
      INTEGER LPOOL, LPOOLB, IPOOL(LPOOL), IPOOLB(LPOOLB)
      INTEGER LCKS(20)
      REAL RINFO(20)
      LOGICAL ENDTAG
      INTEGER POSFAC,IWPOS,LRLU,
     *    IPTRLU, IPTINI, LRLUS,
     *    NBLEC, III, IIIB,
     *    LEAF, LEAFB, NBROOT
      REAL UU
      LOGICAL IOK,AGAIN,ALOCLU
      INTRINSIC MOD
C
C   .. Description of the parameters ..
C
C NTASKS corresponds to the number of tasks activated
C        if NTASKS=1 then the memory management is simplified,
C          all working space is given to LU area,
C          active matrices are stored at the top of the stack
C          of the factors.
C
C NPROCW indicates  how many frontal matrices can be stored
C        simultaneously in the active area. It is particularly
C        useful when our input working array for real is
C        such that we have to subdivided in fixed blocks of
C        the size of the largest frontal matrix. If NPROCW
C        is equal to one then we process the elimination with
C        a depth first search order and we get only the parallelism
C        comming from the node process.
C N is an INTEGER variable which must be set by the user to the
C   order n of the matrix A.  It is not altered by the subroutine.
C   Restriction: N greater or equal to 1.
C
C IW  Integer array of size LIW; must hold on entry integer informations
C     on the sorted arrows.
C     IW(1+I), I=1,IWPOS holds on output the LU factors.
C      The node factors are stacked in the ordered in which
C      they have been assembled. For each node the first entry is
C      the size of the contribution block, followed by the number of
C      uneliminated fully-summed variables, followed by the size of the
C      frontal matrix, followed by the number of eliminated pivots.
C      After the header follows the the row indices in the columns,
C      and the column indices in the rows.
C LIW is an INTEGER variable which must be set by the user to the
C     size of array IW. Not altered in this routine.
C
C A   real array of size LA; must hold on entry real informations
C      on the sorted arrows.
C      A(I), I=1,POSFAC holds on output the LU factors. The LU factors
C      are stacked (in subroutine STACK) in the ordered
C      in which their factorization has
C      been completed.
C      A(I), I=1,IPTBEG is used during factorization to hold the
C      LU factors and the contributions blocks. Both informations are
C      stacked on each side of the area of length IPTBEG.
C LA  is an INTEGER variable which must be set by the user to the
C     size of array A. Not altered in this routine.
C
C NSTK, ND, FILS, FRERE, PTRIST, PTRAST, PTRAIW, TLKJ, ITLOC
C     are integer arrays of size N described in MA41Z/ZD.
C
C PERM INTEGER ARRAY of SIZE N THAT NEED NOT BE SET ON ENTRY.
C     It is used in MC51E/ED as a working array to manage the
C     second level of parallelism. PERMW(INODE) holds the number of
C     parallel tasks of 2nd level that have still to be done.
C
C LPOOL,IPOOL(LPOOL),III,ILEAF are used to manage the pool of tasks
C      associated to the parallelism of the elimination tree (level 1)
C LPOOLB,IPOOLB(LPOOLB),IIIB,ILEAFB are used to manage the
C      pool of tasks
C      associated to the node level parallelism (level 2).
C      III (resp IIIB) points to the first task of level 1 (resp 2)
C      ready to be activated.
C      LEAF (resp LEAFB) points to the first free position in the
C      circular pool IPOOL (resp IPOOLB) to store a
C      new task of level 1 (resp 2).
C
C MAXFRT is the maximum frontal matrix size computed
C      during factorization.
C      It is a critical variable that is updated in a critical section.
C
C NTOTPV IS THE TOTAL NUMBER OF PIVOTS SELECTED. THIS IS USED
C      TO DETERMINE WHETHER THE MATRIX IS SINGULAR.
C      It is a critical variable that is updated in a critical section.
C
C .. Memory management ...
C LRLU  : Size of the real space free in the LU area
C LRLUS : Size of the real space that could be obtained
C         (after compress) in the LU area.
C IPTRLU: pointer to the first free position in Real LU area
C LENBUD holds the size of the space that will be manage
C     using a buddy system memory allocation scheme. The dynamic memory
C     allocation scheme is mainly used to allocate the frontal matrices
C     on which we will be working. LENBUD must be at least greater than
C     the 2**LMIN where LMIN is the smaller integer such that
C      2**LMIN > (order of the largest frontal matrix) **2 +3
C     If LENBUD is set equal to -1 then
C     the splitting of the working space is done automatically by the
C     numerical factorization routine MA41B/BD. The default value of
C     LENBUD is -1.
C LENA If only a part of the working space is allocated with a
C      buddy system memory allocation scheme then the remaining space
C      will be subdivided in equal pieces of size LENA.
C      The fixed block strategy is useful to handle efficiently the
C      small frontal matrices often near the leaves of
C      the elimination tree.
C      The default value is
C      the maximum size of the frontal matrices**2/16
C
C .. Blocking during frontal matrix factorization ..
C ICNTL(7-10) and
C ISPON1,ISPON2,ISPON3,ISPON4  are used to control the blocking
C       strategy used to factorize a frontal matrix.
C       A adapted block KJI scheme is used to factorize
C       each frontal matrix.
C       Let us now suppose
C       that we have to factorize a frontal matrix INODE,
C       where NASS is the number fully assembled var. of INODE
C       and NFRONT is the frontal matrix size.
C       If we denote by TLKJ(INODE) the block size
C       used in the KJI scheme
C       then we have:
C    -IF (NASS.GT.KEEP(3)) then TLKJ(INODE)= KEEP(6)
C                         ELSE TLKJ(INODE)= -KEEP(5)
C    (the signed is used to detect quickly that we don't
C     parallelize the elimination step).
C        Blocking is effective only if NASS > KEEP(4)
C
C    -Let us now suppose that we have factorized a block of size
C    IABS(TLKJ(INODE)).
C    Let NPIV be the number of already eliminated var THEN:
C    IF ( (NASS-NPIV).GT.ISPON1) and (TLKJ(INODE).GT.0) ) THEN
C    we update the following NASS-NPIV rows by blocks of ISPON2.
C
C    The updating of the rows containing the contribution block
C    can also be done by block of ISPON4 rows if the number of
C    pivots eliminated in the current front is bigger than KEEP(6)
C    and if (NFRONT-NASS) is greater than ISPON3.
C ISTEP is an integer variable set in a critical section to the
C    current step number of the factorization
C***************************************
C LOCAL ARRAYS (EACH TASK HAS ITS OWN COPY)
C       ITLOC IS USED IN THE ASSEMBLY
C***
C**************************************
C
C*********************
C LOCAL VARIABLES
C*********************
      REAL ZERO
      INTEGER I,NFRONT,NPIV,INODE,INOPV,IOLDPS
      INTEGER LKJIB,IFINB,NASS,NEL1,IROW,NPIVB
      INTEGER IN,IF,NPIVE
      EXTERNAL MC51F, MC51G, MC51H, MC51I, MC51M, MC51N,
     *         MC51O, MC51P, MC51Q, MC51R, MC51S, MC51T
C MAXFRW, NPVW, NOFFW, NELVAW,OPASSW, OPELIW, NCMPBW
C are variable that will hold information local to each task
C The updating of the global quantities is done at the
C end of the routine in a critical section. This avoids multiple
C critical sections during computation.
C Per process:
C  -OPASSW holds the number of operations performed during assembly
C  -OPELIW holds the number of operations performed during assembly
C  -NCMPBW holds the number of compresses required during factorization.
      INTEGER MAXFRW, NPVW, NOFFW, NELVAW, NCMPBW
      REAL    OPASSW, OPELIW
      LOGICAL LASTSON
      PARAMETER (ZERO=0.0E0)
CFPP$ NOAUTOEXPAND R
C *******
C *** INIT OF ITLOC WORKING ARRAY USEFULL TO
C *** SIMPLIFY THE INDIRECT ADDRESSING IN THE ASSEMBLY
C *** ITLOC IS LOCAL TO THE TASK MC51E
C *******
      DO 1234 I=1,N
       ITLOC(I)=0
 1234 CONTINUE
C
C initialization of local variables
C
      MAXFRW = 0
      NPVW   = 0
      NOFFW  = 0
      NELVAW = 0
      NCMPBW  = 0
      OPASSW = ZERO
      OPELIW = ZERO
C
C
C           |-   M A I N   L O O P   -|
C           V                         V
C     We extract work from the pools and do the
C     corresponding work,
C     IPOOL holds node number ready to be assembled
C     IPOOLB holds task identifiers corresponding
C            to node level parallelism.
C
 123  CONTINUE
C
C  IFLAG negative implies unrecoverable error detected during
C        factorization
C  NBNODE = -1 implies that we have reached a chain of
C        nodes that will be processed by switching off
C        the tree level parallelism
      IF (IFLAG.LT.0) GOTO 500
C
C***************************************************
CC
CC Suppress the following 5 lines  in a uniproc version of the code
       IF (NBNODE.LE.-1) THEN
c switch to automatic loop level parallelism
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
C         This statement cannot be reached
C         with a pseudo-parallel execution of the code.
C         NBNODE is the number of frontal matrices
C                being processed.
C         NPROCW is the maximum number of nodes
C                that can be processed in parallel.
            RC = -1
          ELSE
            RC = 0
            III = III + 1
            IF (III.GT.LPOOL) III = 1
C we increment the number of frontal matrices
C being currently processed
C
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
C
 124  IF (RC.EQ.1) THEN
        INODE=IPOOLB(I)
      ELSE
        INODE=IPOOL(I)
      ENDIF
      IF (INODE.GT.N) GO TO 55
      IF (INODE.LT.-N) GO TO 555
C
C  55 : WE PERFORM SPAWNED CONTRIBUTION
C  555: WE PERFORM SPAWNING DURING THE ELIMINATION PROCESS
C
C*********************************
C Assembly process of a new node
C*********************************
 234  CALL MC51F(NTASKS,N,INODE,IW,LIW,A,LA,
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
C
C INODE goes again to the pool
C (not enough space to allocate a frontal matrix
C  while another frontal matrix is currently being factored).
C
        IPOOL(LEAF) = -INODE
        LEAF = LEAF + 1
        IF (LEAF.GT.LPOOL) LEAF = 1
       GOTO 123
      ELSE
         ISTEP = ISTEP + 1
         PTLUST(ISTEP) = PTRIST(INODE)
      ENDIF
 20   CALL MC51I(N,INODE,IW,LIW,A,LA,INOPV,NOFFW,
     *                IFLAG,PTRIST,PTRAST,UU)
C JUMP IF NO MORE PIVOTS CAN BE FOUND IN ELEMENT INODE.
      IF (IFLAG.LT.0) GOTO 500
      IF (INOPV.EQ.1) GO TO 50
      IF (INOPV.EQ.2) THEN
CC NO MORE PIVOTS CAN BE SELECTED BUT UPDATING OPERATIONS HAVE
CC STILL TO BE PERFORMED. WE MEET AN UNSTABLE LINE FOR PIVOTING
CC WHICH IS NOT IN THE LAST BLOCK OF THE KJI FACTORIZATION.
CC SO WE HAVE TO UPDATE THE REMAINING BLOCKS WITH RESPECT TO THE
CC SELECTED VARIABLES OF THE CURRENT BLOCK.
CC
         CALL MC51S(N,INODE,IW,LIW,A,LA,
     *            PTRIST,PTRAST,TLKJ(INODE),KEEP(4))
         GOTO 20
      ENDIF
C***
C ACCUMULATE TOTAL NUMBER OF PIVOTS.
C
      NPVW = NPVW + 1
C**************************************
      IOLDPS = PTRIST(INODE)
      IF (IW(IOLDPS+2).LE.1) THEN
       CALL MC51O(N,INODE,IW,LIW,A,LA,
     *                 PTRIST,PTRAST)
       IW(IOLDPS+1) = IW(IOLDPS+1) + 1
       GO TO 61
      ENDIF
C*******
C INODE CONTAINS MORE THAN ONE PIVOT TO
C ELIMINATE
C*******
         LKJIB = IABS(TLKJ(INODE))
       CALL MC51M(N,INODE,IW,LIW,A,LA,
     *             PTRIST,PTRAST,IFINB,LKJIB,KEEP(4))
       IW(IOLDPS+1) = IW(IOLDPS+1) + 1
C*******
C IFINB = 0 Try to get another variable to eliminate.
C IFINB = 1 MEANS THAT WE REACH THE LAST PIVOT
C           OF A BLOCK OF SIZE LKJIB AND NASS >=  LKJIT
C           Perform updating of the remaining block of rows
C           in the fully summed rows.
C IFINB = -1 We have now to perform the updating of
C            the contribution block.
C*******
       IF (IFINB.EQ.0) GOTO 20
       IF (IFINB.EQ.(-1)) GOTO 50
C*******
C WE ALLOW SPAWNING OF THE UPDATING
C FROM THE PREVIOUS BLOCK OF PIVOTS
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
         CALL MC51Q(A,LA,
     *           NFRONT,NPIV,NASS,POSELT,LKJIB)
         GO TO 20
       ENDIF
CC
C we update the remaining part of the matrix
C with a blocking technique.
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
C
C JUMP BACK TO BEGINNING OF DO LOOP
      GO TO 123
 555  INODE  = -INODE
      IROW   = INODE/(N+1)
      INODE  = INODE - (N+1)*IROW
         LKJIB = IABS(TLKJ(INODE))
      CALL MC51R(N,INODE,IW,LIW,A,LA,IROW,
     *       PTRIST,PTRAST,LKJIB,ISPON2)
      IOK = .FALSE.
C**************************************
C***************************************
C PERM  IS SET TO NUMBER OF ROWS TO BE ELIMINATED.
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
C THE CONTRIBUTION BLOCK IS OF SIZE NEL1 BY NEL1
C----
      IF (NEL1.LE.0) GO TO 60
C 60: NO ROWS HAVE TO BE ELIMINATED, WE PERFORM A STACK
C
C SET PERM ARRAY.  NO LOCKON NEEDED AT THIS POINT SINCE NO
C     WORK FOR THIS NODE IS IN THE POOL
C WE UPDATE THE SUBDIAGONAL BLOCK USING A TRIANG. SOLVER
C***
C IF (NEL1.LT.ISPON3) DON'T DO ANY SPAWNING AND
C    update the block of rows corresponding to the
C    contribution block.
      IF ((NEL1.LT.ISPON3).OR.(NPIV.LT.KEEP(6))) THEN
        IROW = 1
        CALL MC51P(A,LA,IROW,NFRONT,
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
C
C JUMP BACK TO BEGINNING OF DO LOOP
      GO TO 123
 55   IROW   = INODE/(N+1)
      INODE  = INODE - (N+1)*IROW
      IOLDPS = PTRIST(INODE)
      POSELT = PTRAST(INODE)
      NFRONT = IW(IOLDPS)
      NPIV   = IW(IOLDPS+1)
      NASS   = IW(IOLDPS+2)
      CALL MC51P(A,LA,IROW,NFRONT,
     *     NPIV,NASS,POSELT,ISPON3,ISPON4,KEEP(6))
      IOK = .FALSE.
C**************************************
C***************************************
C PERM  IS SET TO NUMBER OF ROWS TO BE ELIMINATED.
      PERM(INODE) = PERM(INODE) - 1
      IF (PERM(INODE).EQ.0) IOK = .TRUE.
C***************************************
C**************************************
      IF (.NOT.IOK) GO TO 123
CCCCCC
C BEFORE STACKING THE LU FACTORS WE
C CHECK IF THERE IS A BLOCK OF UNELIMINATED VAR IN
C THE CURRENT FRONT AND TRY A COLUMN ORIENTED NON BLOCKED
C PIVOTING AND ELIMINATION
CCCCCC
 60   IOLDPS = PTRIST(INODE)
      NPIV   = IW(IOLDPS+1)
      NASS   = IW(IOLDPS+2)
CCC
C we store the number of variables eliminated for
C which updating of the remaining matrix has been done
CCC
      IW(IOLDPS+4) = NPIV
      IF (NASS.EQ.NPIV) GOTO 61
 62   CALL MC51H(N,INODE,IW,LIW,A,LA,INOPV,NOFFW,
     *                PTRIST,PTRAST,UU)
C JUMP IF NO MORE PIVOTS CAN BE FOUND IN ELEMENT INODE.
      IF (INOPV.NE.1) THEN
C      Accumulate total number of pivots
       NPVW = NPVW + 1
C***************************************
       CALL MC51N(N,INODE,IW,LIW,A,LA,
     *                 PTRIST,PTRAST,IFINB)
       IW(IOLDPS+1) = IW(IOLDPS+1) + 1
C
C      IFINB .EQ. 1 means that we reach the last pivot
C      of a block of size LKJIB and that NASS .GE. LKJIT
       IF (IFINB.EQ.0) GOTO 62
      ENDIF
C
C     Update the off-diag block row and the contribution block
      IOLDPS = PTRIST(INODE)
      POSELT = PTRAST(INODE)
      NFRONT = IW(IOLDPS)
      NPIV   = IW(IOLDPS+1)
      NASS   = IW(IOLDPS+2)
      NPIVB  = IW(IOLDPS+4)
      NPIVE  = NPIV - NPIVB
      NEL1   = NFRONT - NASS
      IF ((NPIVE.LE.0).OR.(NEL1.EQ.0)) GO TO 61
        CALL MC51T(A,LA,NPIVB,
     *                NFRONT,NPIV,NASS,POSELT)
C
C=====================
C Stack the LU FACTORS
C=====================
C
 61   CALL MC51G(NTASKS,N,INODE,A,LA,IW,LIW,A,
     *      IFLAG,IERROR,OPELIW,NELVAW,
     *       PTRIST,PTRAST,IPTA,ALOCLU,AGAIN,
     *       LENA,IPTEND,IPTBEG,LCKS,NBACTI,
     *       POSFAC,LRLU,IPTRLU,IPTINI,LRLUS,NBLEC,
     *       NCMPBW,KEEP(17),KEEP(22),ICNTL)
C
      IF (IFLAG.LT.0) GOTO 500
C============================
C TRY TO ACTIVATE FATHER NODE
C============================
      IN = INODE
 30   IN = FRERE(IN)
      IF (IN.GT.0) GO TO 30
      IF (IN.EQ.0) THEN
C      IN is a root node
        NBROOT = NBROOT - 1
        IF (NBROOT.EQ.0) ENDTAG = .TRUE.
C      We update the number of active frontal matrices
         NBACTI = NBACTI-1
         NBNODE = NBNODE-1
       GOTO 123
      ENDIF
      IF = -IN
C
C    We try to activate the father node IF
C    SUBTRACT ONE FROM COUNT OF UNELIMINATED SONS OF FATHER NODE.
C    IF NOW ZERO, activate FATHER NODE.
        NSTK(IF)=NSTK(IF)-1
        LASTSON = (NSTK(IF).EQ.0)
      IF (LASTSON) THEN
       IF ((III.EQ.LEAF).AND.(NBNODE.EQ.1).AND.
     *   (NPROCW.GT.1))  THEN
C        switch to automatic loop level parallelism
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
C JUMP BACK TO BEGINNING OF DO LOOP
C==================================
C
C     We update the number of active frontal matrices
C
        NBACTI = NBACTI-1
        NBNODE = NBNODE-1
      GO TO 123
C
C           ^                                  ^
C           |-   END  OF  M A I N   L O O P   -|
C
500   CONTINUE
C======================================
C We update in critical section
C the global counters with the
C local updated values of the counters
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
C end of driver tasks
C====================
      END
CCC
C--------------------------------------------------------------------
C            HARWELL SUBROUTINE LIBRARY Release 12 (1995)
C        --
C-             Copyright Rutherford Appleton Laboratory
C        --
C--------------------------------------------------------------------
      SUBROUTINE MC51F(NTASKS,N,INODE,IW,LIW,
     *       A,LA,IFLAG,IERROR,ND, FILS,FRERE,
     *       MAXFRW,OPASSW, IPTA,PTRIST,PTRAST,PTRARW,
     *       PTRAIW,ITLOC, LKJIB, NIRBDU,AGAIN,ALOCLU,
     *       LENBUD, LENA, IPTEND, IPTBEG,
     *       NBNODE, NBACTI, LCKS, LRLU, IPTRLU,
     *       IPTINI, IWPOS, POSFAC, LRLUS, NBLEC,
     *       ICNTL, KEEP)
C
C
C Purpose
C =======
C
C Subroutine performing the assembly
C of a new frontal matrix of node INODE.
C The contributions from the son nodes and the
C arrowhead corresponding to the variables to be
C eliminated at this stage are assembled in a full
C so called frontal matrix. Note that non eliminated
C variables at the son step might increase the
C predicted size of the frontal matrix estimated by the
C analysis.
C
C Remark:
C ------
C NTASKS=1 => sequential processing
C    of the tree. Memory allocation
C    strategy is simplified.
C==================================
C
C Parameters
C ==========
C
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
      REAL A(LA), OPASSW
      LOGICAL AGAIN,ALOCLU
C
C
C      .. Error Return ..
C         ============
C
C     On exit, a negative value of IFLAG corresponds to an error.
C     Possible values are:
C         IFLAG = -8 not enough INTEGER working space.
C             The INTEGER working
C             space must be increased by at least IERROR.
C         IFLAG = -9 not enough REAL working space. The real working
C             space must be increased by at least IERROR.
C
C
C      .. Local variables ..
C         ===============
C
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
      INTRINSIC REAL
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
      EXTERNAL MC51K, MC51U, MC51Z
C ZERO HAS THE VALUE 0.
C     SINCE THIS IS UPDATED BY EACH PROCESS, IT SHOULD BE PROTECTED.
C AMAX IS THE SIZE OF THE LARGEST COEFFICIENT IN THE FULLY SUMMED
C IDUMMY IS A DO INDEX NOT USED WITHIN THE LOOP.
C IELL IS A DO INDEX FOR PICKING UP STACKED ELEMENTS.
C IORG IS A DO INDEX FOR PICKING UP ORIGINAL ARROWS.
C J IS A TEMPORARY USED TO HOLD A ROW OR COLUMN NUMBER.
C JJ IS A DO INDEX FOR A LOOP OVER IW.
C JJJ IS ALSO A DO INDEX FOR A LOOP OVER IW.
C J1,J2,J3,J4 ARE BOUNDS FOR LOOPS OVER IW.
C K IS A TEMPORARY VARIABLE.
C LAELL=NFRONT**2 IS THE SIZE OF THE FRONT MATRIX.
C LSTK IS THE NUMBER OF VARIABLES IN THE STACK ELEMENT BEING ASSEMBLED.
C NASS WILL BE SET TO THE NUMBER OF FULLY ASSEMBLED VARIABLES IN THE
C     CURRENT NEWLY CREATED ELEMENT.
C NPIV IS THE NUMBER OF PIVOTS SO FAR SELECTED.
C NUMORG IS THE NUMBER OF ORIGINAL ARROWS TO BE ASSEMBLED IN THE
C     CURRENT STEP.
C NUMSTK IS THE NUMBER OF STACKED ELEMENTS TO BE ASSEMBLED IN THE
C     CURRENT STEP.
C NEWEL IS A POINTER INTO IW TO CONTROL OUTPUT OF INTEGER INFORMATION
C     FOR NEWLY CREATED ELEMENT.
C NFRONT HOLDS THE NUMBER OF VARIABLES IN THE CURRENT FRONT MATRIX.
C POSELT IS A POINTER TO THE POSITION IN A OF THE START OF THE CURRENT
C     FRONT ELEMENT.
C***************************************
        IN = INODE
C
C NUMORG = number of arrowheads to be assembled =number of estimated
C          variables to be eliminated
C IFSON  = first son of INODE
C NUMSTK = number of sons of INODE
C LAELL  = holds the size required to store the frontal matrix of
C          node INODE.
C
        NUMORG = 0
 3      NUMORG = NUMORG + 1
        IN = FILS(IN)
        IF (IN.GT.0) GO TO 3
C       -IN is the node number of the first SON of INODE
C        or zero if INODE is a leaf
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
C       Update the maximum frontal size and the number of variables
C       to be eliminated
        MAXFRW = MAX0(MAXFRW,NFRONT)
        NASS1 = NASS + NUMORG
C
C CHECK TO SEE
C     IF THERE IS SUFFICIENT SPACE.
C     THEN ZERO OUT FRONTAL MATRIX AS APPROPRIATE FIRST
        LAELL = NFRONT*NFRONT
C       POSELT is set to the position in A to store the frontal matrix
        IF (NTASKS.EQ.1) THEN
C-------------------------------
C        simple memory stack
C------------------------------
         IF (LRLU.LT.LAELL) GOTO 630
         LRLU     = LRLU - LAELL
         POSELT   = POSFAC
         POSFAC   = POSFAC + LAELL
        ELSE
C---------------------------
C        NTASKS > 1
C--------------------------
         LREQ = LAELL + 3
CC**********
CC Allocate space of length LREQ
C
         CALL MC51Z(A, LA, IPTA, LREQ, IPT, LSIZ, IRES,
     *    LENA,NBACTI,LCKS)
C OUTPUT:
C        IPT  is the pointer to the free space
C        2**LSIZ is the space really reserved if buddy system
C                was used otherwise LSIZ=0
C        IRES =
C                 1 If fixed block was used to allocate space
C                 2 If buddy system was used to allocate space
C                -1 IF failure in allocating space
C                   (Reachable only during parallel execution).
C                -2 If failure and if no nodes being currently
C                   factored (NBACTI=0)
CC*********
C
         IF (IRES.LT.0) THEN
C         we decrement the number of active frontal matrices
            NBNODE = NBNODE -1
         ENDIF
         IF (IRES.EQ.-2) GOTO 625
         IF (IRES.EQ.-1) GOTO 620
C        Space allocated in position IPT
C        First position to store frontal matrix value is POSELT
         POSELT = IPT+3
C----------------------------
C       end of test on NTASKS
C----------------------------
        ENDIF
C**************************************
C we increment the number of frontal matrices
C currently factored
C  POSELT hold the position of the frontal matrix in A.
C
             NBACTI = NBACTI+1
         POSEL1 = POSELT - NFRONT
         LAPOS2 = POSELT + LAELL - 1
         DO 230 K=POSELT,LAPOS2
           A(K) = ZERO
  230    CONTINUE
CC
C=================================================
C        LREQ = length of space to be allocated
         LREQ = 2*NFRONT+ 5
         IOK = .FALSE.
C        IOK if enough space for integer LU factors
        IF ((IWPOS +LREQ).LE.NIRBDU) IOK = .TRUE.
        IOLDPS = IWPOS
        IWPOS  = IWPOS + LREQ
        IF (.NOT.IOK) GO TO 610
C=================================================
C
        IOLDP2 = IOLDPS + 4
        NEWEL = IOLDP2 + NASS1
        NEWEL1= NASS1
        IW(IOLDPS) = NFRONT
C==========================================
C Build index list of new frontal matrix
C==========================================
C
C       FIRST ASSEMBLE NUMORG INCOMING PIVOT ROWS.
C
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
C
C       NOW ASSEMBLE NUMSTK STACK ELEMENTS
C
        IF (NUMSTK.NE.0) THEN
         NTOTFS = NUMORG
         ISON =  IFSON
         ICT11  = IOLDP2 + NFRONT
         DO 100 IELL=1,NUMSTK
C         ASSEMBLE ELEMENT IELL (SON ISON)
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
C
C         RUN THROUGH non eliminated variables of the son
C
          DO 74 JJ=J1,J3
            NTOTFS            = NTOTFS + 1
            JT1               = IW(JJ)
            IW(ICT11+NTOTFS)  = JT1
            ITLOC(JT1)        = NTOTFS
            IW(JJ)            = NTOTFS
            IW(IOLDP2+NTOTFS) = IW(JJ-NFSON)
  74      CONTINUE
C
C         RUN THROUGH COLUMN INDEX LIST OF  REST OF STACK ELEMENT IELL.
C
  75      J1 = J3 + 1
          IF (NASS1.NE.NFRONT) THEN
           DO 90 JJ=J1,J2
            J = IW(JJ)
            IF (ITLOC(J).EQ.0) THEN
C THE INDEX J WAS NOT FIND IN THE ARRAY
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
C
C NOW INCORPORATE ORIGINAL ARROWS.  NOTE THAT THE INDICES IN THESE
C     ARROWS NEED NOT BE IN ORDER.
        IBROT = INODE
        DO 180 IORG=1,NUMORG
          J1 = PTRAIW(IBROT)+2
          IBROT = FILS(IBROT)
CCC
CC with new MA41H/HD iw(j1-1) is always the number of non zeros
CC in the row
CCC
          J2 = J1 + IW(J1-2) - IW(J1-1)
          J1 = J1 +1
          IF (J1.LE.J2) THEN
           DO 170 JJ=J1,J2
            J   = IW(JJ)
            IF (ITLOC(J).EQ.0) THEN
C THE INDEX J WAS NOT FIND IN THE LIST
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
C
C     WE DEFINE THE REMAINING COLUMN INDICES W.R.T. THE ROW INDICES
C     OF THE NEW CURRENT FRONTAL
C
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
C   We take advantage of the symmetry of the
C   all fully summed variables except those
C   that have been transmitted from the sons
C*********************************************
        K1 = IOLDPS + 5 + NUMORG
        K2 = K1 + NFRONT - 1 + NASS
        DO 554 K=K1,K2
         I        = IW(K)
         ITLOC(I) = 0
  554   CONTINUE
C
C===================================
C ASSEMBLE REALS INTO FRONTAL MATRIX.
C===================================
C
C       JUMP IF THERE ARE NO STACK ELEMENTS TO ASSEMBLE.
        IF (NUMSTK.EQ.0) GO TO 290
C       PLACE REALS CORRESPONDING TO STACK ELEMENTS IN
C       CORRECT POSITIONS IN A.
        ISON   = IFSON
        DO 280 IELL=1,NUMSTK
          ISTCHK = PTRIST(ISON)
          LSTK   = IW(ISTCHK)
C         LSTK is always greater than 0
          SIZFR  = LSTK*LSTK
C         count number of operations
          OPASSW = OPASSW + REAL(SIZFR)
          NFS    = IW(ISTCHK+1)
          NELIM  = IW(ISTCHK+3)
          NFSON  = NELIM + LSTK
          J1     = ISTCHK + NFSON + 5 + NELIM
          J2     = J1 + LSTK - 1
C         NFSON1 is set to the leading dimension of
C         the contribution block
C         IACHK points to the first element of the contribution block
          IF (NTASKS.EQ.1) THEN
C-------------------
C         NTASKS = 1
C-------------------
             IACHK  = PTRAST(ISON)
             NFSON1 = LSTK
          ELSE
C-------------------
C         NTASKS > 1
C-------------------
            ISONLU = (PTRAST(ISON).LE.IPTBEG)
            IF (ISONLU) THEN
C            contribution block was stored in LU area
C            We increase the number of contributions block readers
                NBLEC = NBLEC +1
            ENDIF
C           if contribution block is available in active area then
C           all frontal matrix was left and we must extract only
C           the interesting block
            IASTK  = PTRAST(ISON)
            IF (ISONLU) THEN
             IACHK =  IASTK
             NFSON1 = LSTK
            ELSE
             IACHK  = IASTK + NFSON*NELIM + NELIM
             NFSON1 = NFSON
            ENDIF
C         end of test on value of NTASKS
          ENDIF
C----------------------------------------------------
C         Assemble real values from the son node ISON
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
C
C We put the original indices in the son
C LU factors set of indices. And we use the symmetry of the
C indices of the sons
C
          J3 = J1 + NFS
          DO 275 JJ=J3,J2
           IW(JJ) = IW(JJ-NFSON)
  275     CONTINUE
          IF (NFS.NE.0) THEN
C
C we cannot in this case use the symmetry of the indices
C of the son to retrieve its original value
C then we obtain it from the father frontal matrix.
C
            J3 = J3 -1
            DO 278 JJ=J1,J3
             JPOS = IW(JJ) + ICT11
             IW(JJ) = IW(JPOS)
  278       CONTINUE
          ENDIF
C==========================================================
C     FREE real space of son ISON
C==========================================================
      IF (NTASKS.EQ.1) THEN
C-------------------
C         NTASKS = 1
C-------------------
C        although the current son might not be at the top
C        of the stack of the Cont. Blocks, this update
C        is globally correct at the end of the loop on
C        all the sons of INODE.
         IPTRLU = IPTRLU + SIZFR
         LRLU   = LRLU + SIZFR
      ELSE
C-------------------
C         NTASKS > 1
C-------------------
        IF (ISONLU) THEN
C       the contribution block of node ISON was stacked in LU area
C       LRLUS is updated in MC51K.
          CALL MC51K(N, A, ISON, PTRAST, IPTBEG,
     *       LRLU, LRLUS, IPTRLU, IPTINI,LCKS)
C         Contribution block of the son was allocated in LU
C         We decrement the number of contributions block readers
              NBLEC = NBLEC -1
        ELSE
C         Contribution block of the son is in active area
          IPT = IASTK - 3
          CALL MC51U(A,LA,IPTA, IPT,
     *       IPTEND, IPTBEG,LCKS)
        ENDIF
C      end of test on NTASKS
       ENDIF
C
C Loop 280 we process next son of INODE
C
          ISON   = FRERE(ISON)
  280   CONTINUE
C====================================================
C END OF REAL assembly of the contribution
C blocks of the sons into the frontal matrix of INODE
C====================================================
C
C====================================================
C BEGIN OF INCORPORATE REALS FROM ORIGINAL ARROWS.
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
C         ASSEMBLE COLUMN
C         COMPUTATION OUTSIDE LOOP
          ICT12 = POSELT - NFRONT+ IJROW - 1
C*****************************
C         Because of potential duplicate entries
C         in arrowhead, we cannot vectorize the next loop
C*****************************
          DO 300 JJ=J1,J2
            APOS2 = ICT12 + IW(JJ)*NFRONT
            A(APOS2) = A(APOS2) + A(AINPUT)
            AINPUT = AINPUT + 1
  300     CONTINUE
          IF (J3.GT.J4) GO TO 320
C         Assemble the row
C         Computation outside loop
          ICT13 = POSELT + (IJROW-1)*NFRONT
          NBCOL = J4 - J3 + 1
C*****************************
C         Because of potential duplicate entries
C         in arrowhead, we cannot vectorize the next loop
C*****************************
          DO 310 JJ=1,NBCOL
           JJ1 = ICT13+IW(J3+JJ-1)-1
           A(JJ1) = A(JJ1) + A(AINPUT+JJ-1)
  310     CONTINUE
  320   CONTINUE
C====================================================
C END OF INCORPORATE REALS FROM ORIGINAL ARROWS.
C====================================================
C
C  Update description of assembled frontal matrix
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
C
C*************
C ERROR RETURN
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
C No space enough to allocate frontal matrix in active area
C But another frontal matrix is currently being factored.
C Reachable only during parallel execution.
C
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
C
C we start allocating space in LU if not done before:
C        ALOCLU = .TRUE.
C we force contribution blocks to be allocated in LU area
C and if not enough free space we compress LU area
C        AGAIN  = .TRUE.
C
      AGAIN  = .TRUE.
      ALOCLU = .TRUE.
      INODE = -INODE
      GOTO 640
C
  625 CONTINUE
C Not enough space to allocate a frontal matrix
C in the active area and
C no more active frontal matrices.
C Potential Deadlock due to lack of memory
C for storing the active frontal matrix
C
       IFLAG = -9
       IERROR = LREQ
       GOTO 640
C
  630  CONTINUE
C NTASKS=1 and
C Not enough space to allocate a frontal matrix
C
       IFLAG = -9
       IERROR = LAELL-LRLU
C
 640  RETURN
      END
CCC
C--------------------------------------------------------------------
C            HARWELL SUBROUTINE LIBRARY Release 12 (1995)
C        --
C-             Copyright Rutherford Appleton Laboratory
C        --
C--------------------------------------------------------------------
      SUBROUTINE MC51G(NTASKS, N, INODE,
     *         AFACT, LA, IW, LIW, A,
     *         IFLAG,IERROR,OPELIW,NELVAW,
     *         PTRIST,PTRAST,IPTA,ALOCLU,AGAIN,
     *         LENA,IPTEND,IPTBEG,LCKS,NBACTI,
     *         POSFAC,LRLU,IPTRLU,IPTINI,LRLUS,NBLEC,
     *         NCMPB,PTLUAC,LINAC,ICNTL)
C**********************************************
C
C  .. Purpose ..
C     =======
C Routine to store the LU factors and the
C contribution block of a factored frontal matrix.
C    The LU factors are stored by rows the
C    with the following format.
C    BY ROWS:
C      the Npiv first rows followed
C      by the NFRONT-NPIV off-diagonal block of the L factors
C Remark:
C ------
C NTASKS=1 => sequential processing
C    of the tree. Memory allocation
C    strategy is simplified.
C
C        PB1: Parallel version of the code
C        ---
C        A problem in the memory allocation scheme
C        has been detected. It is probably due to
C        scalar optimization done by the compiler:
C        an instruction in a critical section
C        has probably been partially executed out
C        of the critical section.
C        One should try to compile mc51g.f mc51h.f
C        mc51e.f mc51LU.f and mc51buddy.f without
C        scalar optimization.
C******************************************
C
C
C  .. Parameters ..
C     ==========
      INTEGER NBUD
      PARAMETER (NBUD=29)
      INTEGER NTASKS, N, LA, LIW, INODE,IFLAG,IERROR
      INTEGER LENA,IPTEND,IPTBEG,NCMPB, NBACTI
      INTEGER POSFAC, LRLU, IPTRLU, IPTINI, LRLUS, NBLEC
      INTEGER IW(LIW), IPTA(NBUD), PTLUAC, LINAC
      INTEGER LCKS(20)
      INTEGER ICNTL(20)
      INTEGER PTRIST(N), PTRAST(N), NELVAW
      REAL    A(LA), AFACT(LA), OPELIW
      LOGICAL ALOCLU,AGAIN
C
C LRLU  : Size of the real space free in the LU area
C IPTRLU: pointer to the first free position in Real LU area
C POSFAC IS POINTER FOR FACTORS IN A.
C     IS UPDATED BY ROUTINE AND MUST BE PROTECTED IN A CRITICAL
C     SECTION.
C PTLUAC give the position of the first LU block left in
C               active area
C
C
C     LOCAL VARIABLES
C     ===============
      INTEGER APOS, POSELT, OPSFAC, LSIZ, LP
      INTEGER IOLDPS,NFRONT,NPIV,NASS,LREQCB,LCONT
      INTEGER IPT,I,ISP,IRES,NPOS,J,J1,J2,LREQLU
      REAL    FLOP1
      LOGICAL IOK, FRESP, IOKLU
C
C***************************************
      INTRINSIC REAL
      EXTERNAL MC51J, MC51L, MC51U, MC51Z
      IOLDPS = PTRIST(INODE)
      NFRONT = IW(IOLDPS)
      NPIV = IW(IOLDPS+1)
      NASS = IW(IOLDPS+2)
      LCONT   = NFRONT-NPIV
      POSELT = PTRAST(INODE)
C     update the local number of uneliminated variables
      NELVAW = NELVAW + NASS - NPIV
C
C
C=================================
C stack the LU factors (real part)
C=================================
C     update flops performed during node elimination
      FLOP1  = REAL(2*NFRONT*NPIV)*REAL(NFRONT-NPIV-1)+
     *       REAL(NPIV*(NPIV+1))*REAL(2*NPIV+1)/REAL(3)
      FLOP1  = FLOP1 + REAL(((2*NFRONT-NPIV-1)*NPIV)/2)
      OPELIW = OPELIW + FLOP1
C------------------------------------------------------
C     ISP holds the space necessary for REAL LU factors
      ISP = 2*NPIV*NFRONT - NPIV*NPIV
C------------------------------------
C----------------------------------
C  NTASKS =1, simple stack is used
C----------------------------------
      IF (NTASKS.EQ.1) THEN
C      save position of real factors in header of integer factors
       IW(IOLDPS+4) = POSELT
       IF (LCONT.EQ.0) GOTO 610
C      save contribution block at the top of the stack
       LREQCB = LCONT*LCONT
       IF (LREQCB.GT.LRLU) GOTO 630
C      we do not update LRLU since space used by CB will
C      be freed from factors
       IPTRLU = IPTRLU - LREQCB
       POSFAC = POSFAC - LREQCB
       NPOS   = IPTRLU+1
       PTRAST(INODE) = NPOS
C      -----------------------
C      copy contribution block
C      -----------------------
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
C      ------------------------------------
C      move block left lower part of factors
C      ------------------------------------
       APOS = POSELT + (NPIV+1)*NFRONT
C         vectorization of inner loop is possible
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
C  NTASKS >1
C----------------------------------
C     Jump if no eliminations have been performed, and there are
C     no factors to store.
      IF (NPIV.EQ.0) GOTO 525
C
C     Check if we have enough space in LU area and reserve it
      IPT    = POSELT -3
      IOKLU = .FALSE.
         IF (LRLUS.GE.ISP) THEN
            LRLUS = LRLUS -ISP
            IOKLU = .TRUE.
         ENDIF
C
      IF (IOKLU) THEN
C       Check if enough space to store the LU factors
C       in the LU area otherwise we do a compress of the space
        IOK = .FALSE.
          IF (LRLU.GE.ISP) THEN
             LRLU  = LRLU - ISP
             IOK = .TRUE.
          ENDIF
        IF (.NOT.IOK) THEN
C        Reachable only during parallel execution of the code.
C        We compress the LU area
C        Space is logically available,
C        compress must be successful
         CALL MC51J(N,A,PTRAST,ISP,IRES,IPTBEG,
     *       LRLU,IPTRLU,IPTINI,NBLEC,LCKS)
         NCMPB = NCMPB + 1
C        IRES should always be positive.
C        See PB1 described in header.
         IF (IRES.EQ.-1) THEN
           IERROR = ISP
           GOTO 640
         ENDIF
        ENDIF
        IOK = .TRUE.
C       get pointer to free space obtained
          OPSFAC = POSFAC
          POSFAC = POSFAC + ISP
          IF ((POSFAC-1).GT.IPTRLU) IOK = .FALSE.
C         POSFAC -1 > IPTRLU means that the LU factors will
C         overlap the contribution block at the top the stack
C         IF ((POSFAC+LRLU-1).GT.IPTBEG) IOK = .FALSE.
C       IOK should always be true. See PB1 described in header.
        IF (.NOT.IOK) THEN
          IERROR = ISP
          GO TO 640
        ENDIF
      ELSE
C        LINAC holds the size of the factors left in active area
C        PTLUAC
C             =0  only root nodes may have been left in active area
C                 (if LINAC >0)
C             >0  PTLUAC points to head of the list of
C                 factor matrices left
        IF (NPIV.EQ.NFRONT) THEN
C        we have reached a root, we do not move it.
C        note that this block of LU factors will not
C        be in the linked list
           LINAC = LINAC + ISP
           IW (IOLDPS+4) = POSELT
         GOTO 610
        ENDIF
C       Check if enough space to store the LU factors
C       in the active area otherwise stop factorization.
C       we add to the record of the active area:
C         - the pointer to the next block left in the active area
C         - the effective size of the factors (ISP)
C         - the node number
        LREQLU = ISP + 3 + 3
        CALL MC51Z(A, LA, IPTA, LREQLU, OPSFAC, LSIZ, IRES,
     *    LENA,NBACTI,LCKS)
        IF (IRES.LT.0) GOTO 635
C        PTLUAC give the position of the first LU block left in
C               active area
         LINAC = LINAC + ISP
         A(OPSFAC+3) = REAL(PTLUAC)
         A(OPSFAC+4) = REAL(ISP)
         A(OPSFAC+5) = REAL(INODE)
         PTLUAC      = OPSFAC
        OPSFAC = OPSFAC+6
      ENDIF
C      save position of real factors in header of integer factors
      IW(IOLDPS+4) = OPSFAC
C------------------------------------------------
C     COPY FACTORS TO APPROPRIATE PART OF STORAGE.
C------------------------------------------------
      APOS = POSELT
C       Copy the first NPIV rows (by rows) of the factors
      J1 = APOS
      J2 = APOS + NPIV*NFRONT - 1
      DO 500 J=J1,J2
            AFACT(OPSFAC) = A(J)
            OPSFAC        = OPSFAC + 1
  500 CONTINUE
C     Copy (by rows) the off-diagonal block of L factors
C     of size LCONT=NFRONT-NPIV by NPIV
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
C
C==========================================
C Handling of the Contribution block
C==========================================
C
  525   LREQCB   = LCONT*LCONT + 4
C       LREQCB hold the space required to store the contribution block
C          in the LU area
        IPT    = POSELT -3
        IF (LCONT.NE.0) THEN
C********************************************
C        ALOCLU =.FALSE. if we don't use the LU area to store the
C                     contribution blocks.
C        ALOCLU is set in MA41F/FD and might be modified during facto.
C        if we cannot get space for frontal matrix.
C
          IF (.NOT.ALOCLU) GOTO 610
          IOK = .FALSE.
C*********************************************
CC*********************************************
C IF AGAIN=.TRUE. Then we had difficulties to
C obtain space in buddy or fixed block for a frontal matrix
C We compress the space each time we cannot allocate
C contribution blocks in LU area and we
C allocate them in LU area
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
C We check if after compress enough space
C will be available in LU area for contribution blocks.
C
            FRESP = .FALSE.
            IF (LRLUS.GE.LREQCB) THEN
             LRLUS = LRLUS -LREQCB
             FRESP = .TRUE.
            ENDIF
            IF (.NOT.FRESP) GOTO 610
C           Reachable only during parallel execution of the code.
C           We compress the LU area. The compress must be successful.
            CALL MC51J(N,A,PTRAST,LREQCB,IRES,IPTBEG,
     *       LRLU,IPTRLU,IPTINI,NBLEC,LCKS)
            NCMPB = NCMPB + 1
C           IRES should always be positive.
C           See PB1 described in header.
            IF (IRES.EQ.-1) THEN
             IERROR = LREQCB
             GOTO 640
            ENDIF
           ENDIF
C          we are sure that space is available in LU for cont block
           CALL MC51L(N,A,LREQCB,INODE,PTRAST,IPTBEG,
     *         IPTRLU,IPTINI,LCKS)
C**********************************
C          PTRAST(INODE) gives the position to store the cont. block
C          We increment the number of processes
C          accessing the contributions block in the LU area.
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
C
C We decrement the number of processes
C accessing the contributions block in the LU area.
C
            NBLEC = NBLEC -1
        ENDIF
C**********************************
CC
CC we free the working space in buddy system
CC
            CALL MC51U(A,LA,IPTA,IPT, IPTEND, IPTBEG,LCKS)
CC
  610 IW(IOLDPS+3) = NPIV
      IW(IOLDPS+2) = NFRONT
      IW(IOLDPS+1) = NASS - NPIV
      IW(IOLDPS)   = LCONT
C
C (LCONT = NFRONT- NPIV)
C
      GO TO 650
C
C **** ERROR RETURNS ****
C
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
C Error of the allocation scheme
C        See PB1 described in header.
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
      SUBROUTINE MC51H(N,INODE,IW,LIW,A,LA,
     *   INOPV,NOFFW,PTRIST,PTRAST,UU)
C PERFORM PIVOTING ON ASSEMBLED ELEMENT AT NODE INODE.
C Column oriented partial pivoting strategy with threshold.
C UU is the threshold value for pivoting
      INTEGER N,LIW,LA,INODE,INOPV
      REAL UU, A(LA)
      INTEGER IW(LIW)
C AMROW IS THE SIZE OF THE LARGEST COEFFICIENT IN THE FULLY SUMMED
C     PART OF THE CURRENT COLUMN.
C RMAX IS THE SIZE OF THE LARGEST COEFFICIENT IN THE CURRENT COLUMN.
      REAL ZERO, RMAX, AMROW
      PARAMETER (ZERO=0.0E0)
      REAL  SWOP
      INTEGER APOS, POSELT, PTRIST(N), PTRAST(N)
      INTEGER NFRONT,NOFFW,IOLDPS,NPIV,NASS,IPIV
      INTEGER NPIVP1,JMAX,J1,J3,JJ,J2,IDIAG,ISW,ISWPS1
      INTEGER ISWPS2,KSW
      INTEGER ISAMAX
      INTRINSIC MAX
      EXTERNAL ISAMAX
C NOFFW (local) holds the NUMBER OF off diag. pivots PIVOTS .
C
        INOPV   = 0
        POSELT  = PTRAST(INODE)
        IOLDPS  = PTRIST(INODE)
        NFRONT  = IW(IOLDPS)
        NPIV    = IW(IOLDPS+1)
        NPIVP1  = NPIV + 1
        NASS    = IW(IOLDPS+2)
C       NASS > NPIV
C***************
C EACH PASS THROUGH THIS LOOP TRIES TO CHOOSE ONE PIVOT.
          DO 460 IPIV=NPIVP1,NASS
C APOS IS THE FIRST ENTRY OF THE COLUMN TO BE SEARCHED FOR A PIVOT.
C           APOS = POSELT + NFRONT*(IPIV-1) + NPIV
            APOS = POSELT + NFRONT*NPIV + (IPIV-1)
C JMAX IS THE RELATIVE ROW INDEX OF THE PROSPECTIVE PIVOT.
            JMAX = 1
C
C UU is in this case always greater than zero, otherwise error=-5
C would have been detected in MC51I
C while looking for pivots in the row direction
C THRESHOLD PIVOTING IS BEING USED.
            AMROW = ZERO
C FIND LARGEST ENTRY IN THE FULLY SUMMED PART OF THE PROSPECTIVE PIVOT
C          COLUMN. ALSO RECORD ROW OF THIS LARGEST ENTRY.
            J1 = APOS
            J3    = NASS -NPIV
            JMAX  = ISAMAX(J3,A(J1),NFRONT)
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
C JUMP IF ALL THE ROW IS ZERO.
  370       IF (RMAX.EQ.ZERO) GO TO 460
C JUMP IF STABILITY TEST SATISFIED FOR DIAGONAL ENTRY.
            IDIAG = APOS + (IPIV - NPIVP1)*NFRONT
            IF (ABS(A(IDIAG)).GE.UU*RMAX) JMAX = IPIV - NPIV
            IF (ABS(A(IDIAG)).GE.UU*RMAX) GO TO 380
C CHECK LARGEST OFF-DIAGONAL IN FULLY-SUMMED BLOCK FOR STABILITY.
            IF (AMROW.LT.UU*RMAX) GO TO 460
C***************************************
C ACCUMULATE NUMBER OF OFF-DIAGONAL PIVOTS CHOSEN.
            NOFFW = NOFFW + 1
C**************************************
C
C PIVOT HAS BEEN CHOSEN.
C THE FOLLOWING LOOP INTERCHANGES PIVOT ROW AND COLUMN.
C
  380       IF (IPIV.EQ.NPIVP1) GO TO 400
C SWOP COLUMNS IPIV (FROM J3) AND NPIV+1 (J1 TO J2).
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
C SWOP ROWS    NPIV+1 AND NPIV+JMAX.
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
C A COMPLETE PASS WAS MADE WITHOUT FINDING A PIVOT.
C   1ST CASE: NASSW.EQ.NASS (ONLY POSSIBLE CASE FOR PIVOTN)
C             we send to the father the non eliminated
C             variables
       INOPV = 1
  420 RETURN
      END
CCC
C--------------------------------------------------------------------
C            HARWELL SUBROUTINE LIBRARY Release 12 (1995)
C        --
C-             Copyright Rutherford Appleton Laboratory
C        --
C--------------------------------------------------------------------
C********************************************************
C***********PIVOT********
      SUBROUTINE MC51I(N,INODE,IW,LIW,A,LA,
     *    INOPV,NOFFW,IFLAG,PTRIST,PTRAST,UU)
C
C Purpose
C =======
C
C PERFORM PIVOTING ON ASSEMBLED ELEMENT AT NODE INODE.
C ROW oriented partial pivoting strategy with threshold.
C UU is the threshold value for pivoting
C
C Parameters
C ==========
C
      INTEGER N,LA,LIW,INODE,IFLAG,INOPV,NOFFW
      REAL A(LA)
      REAL UU
      INTEGER IW(LIW), PTRIST(N), PTRAST(N)
C NOFFW :
C     NUMBER OF PIVOTS CHOSEN FROM OFF THE DIAGONAL OF THE
C     FRONTAL MATRIX.
C     NOFFW IS local to MC51E/ED and need not be protected in a
C     critical section.
C
C
C      .. Local variables ..
C         ===============
C
C AMROW IS THE SIZE OF THE LARGEST COEFFICIENT IN THE FULLY SUMMED
C     PART OF THE CURRENT ROW.
C RMAX IS THE SIZE OF THE LARGEST COEFFICIENT IN THE CURRENT ROW.
      REAL SWOP
      INTEGER APOS, POSELT
      REAL ZERO, RMAX, AMROW
      INTEGER NFRONT,IOLDPS,NPIV,NASS,NASSW,IPIV
      INTEGER NPIVP1,JMAX,J1,J3,JJ,J2,IDIAG,ISW,ISWPS1
      INTEGER ISWPS2,KSW
      INTEGER ISAMAX
      INTRINSIC MAX
      PARAMETER (ZERO=0.0E0)
      EXTERNAL ISAMAX
C
        INOPV   = 0
        POSELT  = PTRAST(INODE)
        IOLDPS  = PTRIST(INODE)
        NFRONT  = IW(IOLDPS)
        NPIV    = IW(IOLDPS+1)
        NPIVP1  = NPIV + 1
        NASS    = IW(IOLDPS+2)
C*****
C  NASSW is the last possible pivot of
C  the current block
C*****
        NASSW   = IABS(IW(IOLDPS+3))
C***************
C EACH PASS THROUGH THIS LOOP TRIES TO CHOOSE ONE PIVOT.
          DO 460 IPIV=NPIVP1,NASSW
C APOS IS THE FIRST ENTRY OF THE ROW TO BE SEARCHED FOR A PIVOT.
            APOS = POSELT + NFRONT*(IPIV-1) + NPIV
C IF THE USER HAS INDICATED THAT THE MATRIX IS DIAGONALLY DOMINANT, WE
C     DO NOT NEED TO TEST FOR STABILITY BUT WE DO CHECK TO SEE IF THE
C     PIVOT IS NON-ZERO.  IF IT IS ZERO, WE EXIT WITH AN ERROR.
C JMAX IS THE RELATIVE COLUMN INDEX OF THE PROSPECTIVE PIVOT.
            JMAX = 1
            IF (UU.GT.ZERO) GO TO 340
            IF (A(APOS).EQ.ZERO) GO TO 630
C JUMP IF DIAGONAL IS BEING USED AS PIVOT.
            GO TO 380
C THRESHOLD PIVOTING IS BEING USED.
  340       AMROW = ZERO
C FIND LARGEST ENTRY IN THE FULLY SUMMED PART OF THE PROSPECTIVE PIVOT
C          ROW. ALSO RECORD COLUMN OF THIS LARGEST ENTRY.
            J1 = APOS
            J2 = APOS - NPIV + NASS - 1
CCC
             J3    = NASS -NPIV
             JMAX  = ISAMAX(J3,A(J1),1)
             JJ    = JMAX + J1 - 1
             AMROW = ABS(A(JJ))
CCC
C           DO 350 JJ=J1,J2
C             IF (ABS(A(JJ)).LE.AMROW) GO TO 350
C             JMAX  = JJ - J1 + 1
C             AMROW = ABS(A(JJ))
C 350       CONTINUE
C DO SAME AS ABOVE FOR NON-FULLY-SUMMED PART ONLY HERE WE DO NOT NEED
C     TO RECORD COLUMN SO LOOP IS SIMPLER.
            RMAX = AMROW
            J1 = J2 + 1
            J2 = APOS - NPIV + NFRONT - 1
            IF (J2.LT.J1) GO TO 370
            DO 360 JJ=J1,J2
              RMAX = MAX(ABS(A(JJ)),RMAX)
  360       CONTINUE
C JUMP IF ALL THE ROW IS ZERO.
  370       IF (RMAX.EQ.ZERO) GO TO 460
C JUMP IF STABILITY TEST SATISFIED FOR DIAGONAL ENTRY.
            IDIAG = APOS + IPIV - NPIVP1
            IF (ABS(A(IDIAG)).GE.UU*RMAX) JMAX = IPIV - NPIV
            IF (ABS(A(IDIAG)).GE.UU*RMAX) GO TO 380
C CHECK LARGEST OFF-DIAGONAL IN FULLY-SUMMED BLOCK FOR STABILITY.
            IF (AMROW.LT.UU*RMAX) GO TO 460
C**************************************
C ACCUMULATE NUMBER OF OFF-DIAGONAL PIVOTS CHOSEN.
C
            NOFFW = NOFFW + 1
C***************************************
C
C PIVOT HAS BEEN CHOSEN.
C THE FOLLOWING LOOP INTERCHANGES PIVOT ROW AND COLUMN.
C
  380       IF (IPIV.EQ.NPIVP1) GO TO 400
C SWOP ROWS IPIV (FROM J3) AND NPIV+1 (J1 TO J2).
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
C SWOP COLUMNS NPIV+1 AND NPIV+JMAX.
C THis has to be in a critical section
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
C A COMPLETE PASS WAS MADE WITHOUT FINDING A PIVOT.
C   1st case: nassw.eq.nass
C             we send to the father the non eliminated
C             variables
C   2nd case: nassw .ne. nass
C             we have to update the following block w.r.t
C             the already eliminated variables.
C             2nd we swap rows npiv+1, nassw
C
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
      SUBROUTINE MC51J(N,A,PTRAST,ISPA,
     *   IRES,IPTBEG, LRLU,IPTRLU,IPTINI,NBLEC,LCKS)
C     ..
C
C  Purpose
C  =======
C
C Garbage collection routine :
C We compress the space allocated in LU
C
C
C  Parameters
C  ==========
C
C INPUT:
C       ISPA      : is the size of the memory space
C                   that we want to allocate in the LU area.
C       IPTBEG
C INPUT/OUTPUT:
C       A         : original matrix
C       PTRAST(I) : gives location of frontal matrix (real entries)
C                   for node  I
C       LRLU is updated
C            ISPA is subtracted from LRLU if IRES.NE.-1
C OUTPUT:
C       IRES      :
C                  1 if no problem during compress
C                 -1 if the free block after compress is smaller
C                    than ISPA (i.e. LRLU <ISPA)
C
C  Remark:
C  ------
C    This subroutine can be reached only during parallel execution
C    of the code.
C
      INTEGER N,LRLU,IPTRLU,IPTINI,NBLEC
      INTEGER IPTBEG,IPREV,INODE,LREQ,IPTNEW
      INTEGER I,II,ILAST,IPTCU,IOLD,LBLOCK,ICB
      REAL A(IPTBEG)
      INTEGER  IRES,PTRAST(N),ISPA
      LOGICAL IOK
      INTEGER LCKS(20)
      INTRINSIC REAL
C***
C Statement included in serial version to avoid compiler warning
C***
      LCKS(1) = 0
CCC
C LRLU  : Size of the real space free in the LU area
C IPTRLU: pointer to the first free position in Real LU area
CCC
      IRES = 1
      IOK  = .FALSE.
C*************************************************
C Wait for the end for processes currently
C accessing the LU area for assembling or stacking
C contribution blocks
 220   CONTINUE
      IF (NBLEC.GT.0) THEN
CC
C  We wait for the
C  number of processes accessing the
C  contributions blocks in LU area to be zero
C  We check is free space is enough without compress
CC
         GOTO 220
      ENDIF
C
C NBLEC =0
C get access to LU area (LCKS(6))
C*****
           IF (LRLU.GE.ISPA) THEN
            IOK = .TRUE.
            LRLU = LRLU -ISPA
           ENDIF
      IF (IOK) GOTO 410
C*************************
CCCCC
CCCCC
C  we compress the contribution blocks stored in the LU area
C    -we start from the bottom of the stack
CCC
C IPTINI point to the first busy block
C
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
C
C the first block can be compressed
C
       ILAST  = IPTINI + LREQ - 1
       II     = IPTBEG
       IPTINI = IPTNEW
       IF (2*LREQ.LT.LBLOCK) THEN
C
C copy is vectorizable
C
        DO 10 I=1,LREQ
         A(II) = A(ILAST)
         ILAST = ILAST -1
         II    = II - 1
 10     CONTINUE
       ELSE
C
C copy is not vectorizable
C
        DO 20 I=1,LREQ
         A(II) = A(ILAST)
         ILAST = ILAST -1
         II    = II - 1
 20     CONTINUE
       ENDIF
C      update position of node INODE
       PTRAST(INODE) = IPTNEW+4
      ENDIF
C
C     Dummy loop we have at maximum
C     N-1 contribution blocks to compress
      DO 50 ICB=1,N
C
       IF (IPREV.EQ.0) GOTO 400
       A(IPREV+1)    = REAL(IPTNEW)
       IPTCU = IPREV
       IOLD   = IPTNEW
       IPTNEW = IPTCU
       LBLOCK = IOLD - IPTCU
       IPREV  =  INT (A(IPTCU+2) + 0.5)
       INODE  =  INT (A(IPTCU+3) + 0.5)
       LREQ   =  INT (A(IPTCU) + 0.5)
       IF (LREQ.GE.LBLOCK) GOTO 50
C
C compress of the frontal matrix (IOLD) in the
C block of length LBLOCK is possible
C
       IPTNEW = IOLD - LREQ
       ILAST  = IPTCU + LREQ - 1
       II     = IPTNEW + LREQ -1
       IF (2*LREQ.LT.LBLOCK) THEN
C
C copy is vectorizable
C
        DO 110 I=1,LREQ
          A(II) = A(ILAST)
          ILAST = ILAST -1
          II    = II - 1
 110    CONTINUE
       ELSE
C
C copy is not vectorizable
C
        DO 120 I=1,LREQ
         A(II) = A(ILAST)
         ILAST = ILAST -1
         II    = II - 1
 120    CONTINUE
       ENDIF
C      update position of node INODE
       PTRAST(INODE) = IPTNEW+4
       A(IOLD+2) = REAL(IPTNEW)
  50  CONTINUE
C     this statement should never be reached
      IRES =-1
C
 400  CONTINUE
       LBLOCK = IPTNEW - IPTRLU - 1
       IPTRLU = IPTNEW - 1
        IF (LRLU+LBLOCK .GE.ISPA) THEN
           LRLU   = LRLU + LBLOCK - ISPA
         ELSE
           IRES=-1
         ENDIF
C
C free access to LU area
 410  CONTINUE
C***************************
C***************************
C***************************
C***************************
 500  RETURN
      END
CCC
      SUBROUTINE MC51K(N,A,INODE,
     * PTRAST,IPTBEG,LRLU,LRLUS,IPTRLU,IPTINI,LCKS)
C******************************************************
C
C  Purpose
C  =======
C
C We free space in LU and compress only the list of
C free blocks at the top of the stack
C
C  We free Real space at position IPT
C        IPT = PTRAST(INODE) - 4
C
C  Parameters
C  ==========
C
C INPUT:
C        INODE = node index of cont. block becoming free
C        PTRAST(inode) = position of node INODE in LU area
C        LRLUS : Space logically available in LU area.
C        LRLU  : Size of the real space free at the
C                top of the stack of contribution blocks.
C        IPTRLU: pointer to the first free position
C                in LU area
C OUTPUT:
C        IPTRLU, LRLU, LRLUS are updated.
C
C************************************************************
      INTEGER N
      INTEGER INODE,IPTBEG,IPREV,INEXT
      INTEGER LRLU, LRLUS, IPTRLU, IPTINI
      REAL A(IPTBEG)
      INTEGER  IPT, PTRAST(N),SIZFR
      REAL ZERO
      INTEGER LCKS(20)
      INTRINSIC REAL, INT
      PARAMETER (ZERO=0.0E0)
C***
C Statement included in serial version to avoid compiler warning
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
C        working space used in LU becomes completely free
          IPTINI = IPTBEG
          IPTRLU = IPTBEG
          LRLU   = LRLU + IPTBEG - IPT + 1
         ELSE
C         next block is now at the top
          A(INEXT+2) = ZERO
          IPTRLU     = INEXT - 1
          LRLU       = LRLU + INEXT - IPT
         ENDIF
       ELSE
C      block is not at the top
         IF (INEXT.EQ.0) THEN
          IPTINI = IPREV
          A(IPREV+1) = ZERO
         ELSE
          A(IPREV+1) = REAL(INEXT)
          A(INEXT+2) = REAL(IPREV)
         ENDIF
       ENDIF
      RETURN
      END
CCC
C--------------------------------------------------------------------
C            HARWELL SUBROUTINE LIBRARY Release 12 (1995)
C        --
C-             Copyright Rutherford Appleton Laboratory
C        --
C--------------------------------------------------------------------
C     ..
C
C  Purpose
C  =======
C  Memory Allocation routines dealing with the LU area
C  (the LU area holds both factors and contribution blocks)
C
C The following set of routines allows to allocate (MC51L),
C to free (MC51K) and to compress (MC51J) space in the
C space reserved initially for the LU factors.
C The organization of the records is a stack of chained records.
C
C Each record, REC(1:LREQ), allocated consists in :
C    REC(1)  = LREQ
C    REC(2)  = pointer to the previous record (or 0 bottom of stack)
C    REC(3)  = pointer to next record (or 0 if top of stack)
C    REC(4)  = Node number
C
C parameters used in by the procedures
C IPTBEG: SIZE reserved for the LU area
C LRLU  : Size of the real space free in the LU area
C IPTRLU: pointer to the first free position in Real LU area
C IPTINI: position of the first busy block
C LRLUS : indicates the amount of space available in the LU
C          (compress might be required to effectively get it)
C**
C
      SUBROUTINE MC51L(N,A,LREQ,INODE,
     *   PTRAST,IPTBEG, IPTRLU,IPTINI,LCKS)
C
C  Purpose
C  =======
C
C We allocate Real space of length LREQ In the LU factor AREA
C to store the contribution block of node INODE
C
C
C  Parameters
C  ==========
C
C INPUT:
C        LREQ  : size of th memory space to be allocated
C              (= space for contribution block + 4)
C        INODE : node number
C OUTPUT:
C        PTRAST (INODE) is updated
C
C************************************************************
      INTEGER LREQ,IPTBEG,N
      REAL A(IPTBEG)
      INTEGER IPT,INODE,PTRAST(N)
      REAL ZERO
      INTEGER LCKS(20)
      INTEGER IPTRLU, IPTINI
      INTRINSIC REAL
      PARAMETER (ZERO=0.0E0)
C***
C Statement included in serial version to avoid compiler warning
C***
      LCKS(1) = 0
C*********************************
        IF (IPTRLU.EQ.IPTBEG) THEN
         IPTRLU = IPTBEG - LREQ
         IPT    = IPTRLU + 1
         IPTINI = IPT
         A(IPT) = REAL(LREQ)
         A(IPT+1) = ZERO
         A(IPT+2) = ZERO
         A(IPT+3) = REAL(INODE)
        ELSE
         IPT = IPTRLU -LREQ+1
         A(IPT)   = REAL(LREQ)
         A(IPT+1) = REAL(IPTRLU +1)
         A(IPT+2) = ZERO
         A(IPT+3) = REAL(INODE)
         A(IPTRLU +3) = REAL(IPT)
         IPTRLU = IPT - 1
        ENDIF
         PTRAST(INODE) = IPT+4
C**********************************
      RETURN
      END
CCC
C ************** MC51M ****************
      SUBROUTINE MC51M(N,INODE,IW,LIW,A,LA,
     *     PTRIST, PTRAST,IFINB,LKJIB,LKJIT)
C
C  .. Purpose ..
C     =======
C
C BLAS2 updating of the remaining rows in the current
C block being factored.
C If we have reached the last variable of the current block then
C we update the first index and the last index of the
C new current block.
C****
C WE ONLY APPLY THE PREVIOUS ELIMINATION TO THE
C REMAINING PART OF THE BLOCK SO THAT WE CAN
C LATER PERFORM THE UPDATING OF THE REMAINING BLOCKS USING
C A blocked KJI SAXPY BLOCK SCHEME
C
C  .. Parameters ..
C     ==========
C  IFINB is an Integer variable that need not be set on entry.
C        On exit possible values are:
C         IFINB =0 other variables have still to be eliminated
C               in current block
C               =1 we have reached the last pivot in the current block
C                  and the current block is not the last block of the
C                  fully summed rows.
C               =-1 we have reached the last pivot in the current block
C                  and the current block is the last block of the
C                  fully summed variables rows.
      INTEGER N,LA,LIW,INODE,IFINB,LKJIB
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    ALPHA, VALPIV
      INTEGER APOS, POSELT, UUPOS, PTRIST(N), PTRAST(N)
      INTEGER LKJIT
      REAL ONE
      INTEGER IOLDPS,NFRONT,NPIV,NASS,JROW2
      INTEGER NEL2,NPIVP1,KROW,LPOS,NEL
      PARAMETER (ONE=1.0E0)
      EXTERNAL SGER
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
C         We update the first index and the last index of the
C         new current block.
          IW(IOLDPS+3) = MIN0(JROW2+LKJIB,NASS)
          IW(IOLDPS+4) = NPIVP1+1
         ENDIF
        ELSE
C        perform updating operations.
         APOS   = POSELT + NPIV*(NFRONT + 1)
C        WE DEFINE VALPIV= -1/ (THE VALUE OF THE PIVOT) = 1/A (APOS)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + NFRONT
         DO 541 KROW = 1,NEL2
C        CALCULATE MULTIPLIER.
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT
 541     CONTINUE
         LPOS   = APOS + NFRONT
         UUPOS  = APOS+1
C        we use BLAS2 routine to perform the update of
C        the remaining rows in the current block
         ALPHA = -1.0E0
         CALL SGER(NEL,NEL2,ALPHA,A(UUPOS),1,A(LPOS),NFRONT,
     *              A(LPOS+1),NFRONT)
        ENDIF
        RETURN
        END
CCC
      SUBROUTINE MC51N(N,INODE,IW,LIW,A,LA,
     *       PTRIST,PTRAST,IFINB)
      INTEGER N,LA,LIW,INODE,IFINB
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    ALPHA,VALPIV
      INTEGER APOS, POSELT,UUPOS,
     *        PTRIST(N), PTRAST(N)
      INTEGER IOLDPS,NFRONT,NPIV,NASS,KROW
      INTEGER NEL,LPOS,ICOL,NEL2,IRWPOS
      INTEGER NPIVP1
      REAL ONE
      PARAMETER (ONE=1.0E0)
      EXTERNAL SAXPY
C PERFORM THE ELIMINATION ON ROW IROW OF ELEMENT INODE.
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
C SEE IF ALL ELIMINATIONS ARE BEING PERFORMED IN THIS ONE CALL.
CPA WE DEFINE VALPIV= -1/ (THE VALUE OF THE PIVOT) = 1/A (APOS)
        VALPIV = ONE/A(APOS)
        LPOS   = APOS + NFRONT
        DO 541 KROW = 1,NEL
C CALCULATE MULTIPLIER.
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT
 541    CONTINUE
        LPOS   = APOS + NFRONT
        UUPOS  = APOS+1
        DO 440 ICOL = 1,NEL
C CALCULATE MULTIPLIER.
             IRWPOS  = LPOS + 1
             ALPHA   = -A(LPOS)
             CALL SAXPY(NEL2,ALPHA,A(UUPOS),1,A(IRWPOS),1)
             LPOS    = LPOS + NFRONT
  440   CONTINUE
        RETURN
        END
CCC************* MC51O***************
      SUBROUTINE MC51O(N,INODE,IW,LIW,A,LA,PTRIST,PTRAST)
      INTEGER N,INODE,LA,LIW
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    ALPHA,VALPIV
      INTEGER APOS, POSELT, UUPOS, PTRIST(N), PTRAST(N)
      INTEGER IOLDPS,NFRONT,NPIV,NEL
      INTEGER LPOS,JROW,IRWPOS
      REAL ONE
      PARAMETER (ONE=1.0E0)

      EXTERNAL SAXPY

        POSELT = PTRAST(INODE)
        IOLDPS = PTRIST(INODE)
        NFRONT = IW(IOLDPS)
        NPIV   = IW(IOLDPS+1)
        NEL    = NFRONT - NPIV - 1
        APOS   = POSELT + (NPIV)*NFRONT + NPIV
C JUMP IF ALL ELIMINATIONS COMPLETE.
        IF (NEL.EQ.0) GO TO 650
C SEE IF ALL ELIMINATIONS ARE BEING PERFORMED IN THIS ONE CALL.
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
C CALCULATE MULTIPLIER.
             IRWPOS  = LPOS + 1
             ALPHA   = -A(LPOS)
             CALL SAXPY(NEL,ALPHA,A(UUPOS),1,A(IRWPOS),1)
             LPOS    = LPOS + NFRONT
  440   CONTINUE
  650   RETURN
        END
CCC
C
C*******************************************
C   Level 3 BLAS update of frontal matrices
C*******************************************
C
      SUBROUTINE MC51P(A,LA,JROW1,NFRONT,
     *       NPIV,NASS,POSELT,ISPON3,ISPON4,LKJPAR)
      INTEGER LA,POSELT,LKJPAR
      REAL    A(LA)
      INTEGER JROW1, NFRONT, NPIV, NASS, ISPON3, ISPON4
C local scalars
      REAL ONE
      INTEGER NEL1,NEL11,LPOS2,LPOS1,LPOS
      REAL ALPHA
      INTEGER JROW2,JROW3
      PARAMETER (ONE=1.0E0)

      EXTERNAL STRSM, SGEMM
C************************************
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
C SEE IF ALL ELIMINATIONS ARE BEING PERFORMED IN THIS ONE CALL.
        JROW2 = MIN0(JROW1+ISPON4-1,NEL1)
        IF ((NEL1.LT.ISPON3).OR.(NPIV.LT.LKJPAR)) JROW2 = NEL1
CPA JROW3 = NB OF ROWS ON WHICH WE DO THE ELIMINATION
        JROW3  = JROW2-JROW1+1
        NEL1   = JROW3
CPA LPOS  = POSITION IN A OF THE FIRST ELEMENT OF HE BLOCK
        LPOS2  = POSELT + (NASS+JROW1-1)*NFRONT
C*** WE HAVE FIRST TO COMPUTE THE L FACTORS (BELOW DIAGONAL BLOCK)
C    SO WE APPLY A TRIANGULAR SOLVER OF SIZE NPIV WITH
C    JROW3 SECOND MEMBERS.
C BLAS-3 STRSM
C
        CALL STRSM('L','L','N','N',NPIV,JROW3,ONE,A(POSELT),NFRONT,
     *              A(LPOS2),NFRONT)
C******
C WE NOW UPDATE THE CORRESPONDING CONTRIBUTION BLOCK
C******
        LPOS   = LPOS2 + NPIV
        LPOS1  = POSELT + NPIV
        ALPHA  = -1.0E0
        CALL SGEMM('N','N',NEL11,NEL1,NPIV,ALPHA,A(LPOS1),
     *          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
  500   RETURN
        END
CCC
      SUBROUTINE MC51Q(A,LA,NFRONT,NPIV,NASS,POSELT,LKJIB)
      INTEGER LA, NFRONT, NPIV, NASS, LKJIB
      REAL    A(LA)
      INTEGER POSELT
      REAL    ALPHA
      INTEGER NEL1, NEL11, NPBEG, LPOS, LPOS1, LPOS2
      REAL ONE
      PARAMETER (ONE=1.0E0)

      EXTERNAL STRSM, SGEMM
CC in this case JROW1 = 1
CC and JROW2 = NEL1
        NEL1   = NASS - NPIV
        NPBEG  = NPIV - LKJIB + 1
        NEL11  = NFRONT - NPIV
        LPOS2  = POSELT + NPIV*NFRONT + NPBEG - 1
C*** WE HAVE FIRST TO COMPUTE THE L FACTORS (BELOW DIAGONAL BLOCK)
C    SO WE APPLY A TRIANGULAR SOLVER OF SIZE LKJIB WITH
C    (NFRONT - NASS) SECOND MEMBERS.
C***
        POSELT = POSELT + (NPBEG-1)*NFRONT + NPBEG - 1
        CALL STRSM('L','L','N','N',LKJIB,NEL1,ONE,A(POSELT),
     *               NFRONT,A(LPOS2),NFRONT)
        LPOS   = LPOS2 + LKJIB
        LPOS1  = POSELT + LKJIB
        ALPHA  = -1.0E0
        CALL SGEMM('N','N',NEL11,NEL1,LKJIB,ALPHA,A(LPOS1),
     *       NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
  500   RETURN
        END
CCC
      SUBROUTINE MC51R(N,INODE,IW,LIW,A,LA,
     *  JROW1, PTRIST,PTRAST,LKJIB,ISPON2)
      INTEGER N,LA,LIW
      INTEGER IOLDPS, NFRONT, NASS, NPIV, NPBEG,JROW2
      REAL    A(LA)
      INTEGER IW(LIW)
      REAL    ALPHA
      INTEGER INODE,POSELT,JROW1,
     *        PTRIST(N), PTRAST(N), LKJIB, ISPON2
      INTEGER NEL1, NEL11, JROW3, LPOS, LPOS1, LPOS2
      REAL ONE
      PARAMETER (ONE=1.0E0)

      EXTERNAL STRSM, SGEMM

C***************************
C if JROW1 .eq.0 we perform
C updating of lkjib lines
C plus elimination of then
C else
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
C SEE IF ALL ELIMINATIONS ARE BEING PERFORMED IN THIS ONE CALL.
CPA JROW3 = NB OF ROWS ON WHICH WE DO THE ELIMINATION
        JROW3  = JROW2-JROW1+1
        NEL1   = JROW3
CPA LPOS  = POSITION IN A OF THE FIRST ELEMENT OF HE BLOCK
        LPOS2  = POSELT + (NPIV+JROW1-1)*NFRONT + NPBEG - 1
C*** WE HAVE FIRST TO COMPUTE THE L FACTORS (BELOW DIAGONAL BLOCK)
C    SO WE APPLY A TRIANGULAR SOLVER OF SIZE LKJIB WITH
C    (NFRONT - NASS) SECOND MEMBERS.
C***
        POSELT = POSELT + (NPBEG-1)*NFRONT + NPBEG - 1
        CALL STRSM('L','L','N','N',LKJIB,JROW3,ONE,A(POSELT),
     *               NFRONT,A(LPOS2),NFRONT)
        LPOS   = LPOS2 + LKJIB
        LPOS1  = POSELT + LKJIB
        ALPHA  = -1.0E0
        CALL SGEMM('N','N',NEL11,NEL1,LKJIB,ALPHA,A(LPOS1),
     *          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
  500   RETURN
        END
CCC
C**** ROUTINE MC51S **************
C  WE CALL THIS SUB ONLY WHEN NO PIVOT WERE SELECTED
C  AND WHEN THE LAST SELECTED PIVOT IS BELONGING TO A BLOCK
C  FOR WHICH THE UPDATING (IN TERMS OF KJI ALGO HAS TO BE
C  PERFORMED.
C*************************************
      SUBROUTINE MC51S(N,INODE,IW,LIW,A,LA,
     *    PTRIST,PTRAST,LKJIB,LKJIT)
      INTEGER N,LA,LIW
      REAL    A(LA)
      INTEGER IW(LIW), LKJIB, INODE
      REAL    ALPHA, ONE
      INTEGER POSELT, PTRIST(N), PTRAST(N)
      INTEGER IOLDPS, NFRONT, NPIV, NASS, JROW2, NPBEG
      INTEGER NONEL, LKABS, LKJIW, NEL1, NEL11
      INTEGER LBP, IPOS, IPOSEL, KPOS, LPOS2
      INTEGER LPOS1,LPOS,LBPT,I1,K1,II,ISWOP,LBP1
      INTEGER LKJIT
      PARAMETER (ONE=1.0E0)

      EXTERNAL STRSM, SGEMM, SSWAP

        POSELT = PTRAST(INODE)
        IOLDPS = PTRIST(INODE)
        IPOSEL = POSELT
        NFRONT = IW(IOLDPS)
        NPIV   = IW(IOLDPS+1)
        NASS   = IW(IOLDPS+2)
        JROW2  = IABS(IW(IOLDPS+3))
        NPBEG  = IW(IOLDPS+4)
C*******
C we update the following values
C and pursue the elimination
C with an adapted  blocking strategy
C because of instability for the current
C node. We don't modify the strategy for
C the other nodes.
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
C JUMP IF ALL ELIMINATIONS COMPLETE.
        IF ((NEL1.EQ.0).OR.(LKJIW.EQ.0)) GO TO 500
        LPOS2  = POSELT + JROW2*NFRONT + NPBEG - 1
C*** WE HAVE FIRST TO COMPUTE THE L FACTORS (BELOW DIAGONAL BLOCK)
C    SO WE APPLY A TRIANGULAR SOLVER OF SIZE LKJIW WITH
C     NEL1 SECOND MEMBERS.
C***
         POSELT = POSELT + (NPBEG-1)*NFRONT + NPBEG - 1
         CALL STRSM('L','L','N','N',LKJIW,NEL1,ONE,A(POSELT),NFRONT,
     *               A(LPOS2),NFRONT)
C******
C WE NOW UPDATE THE CORRESPONDING CONTRIBUTION BLOCK
C******
        LPOS   = LPOS2 + LKJIW
        LPOS1  = POSELT + LKJIW
        ALPHA  = -1.0E0
        CALL SGEMM('N','N',NEL11,NEL1,LKJIW,ALPHA,A(LPOS1),
     *          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
  500   LBP = JROW2 - NPIV
C       LBP > 0
C***
C LBP is the length of the block pivot of non eliminated variable.
C We swap it only if necessary.
C***
         LBPT  = LBP + LBP
        IF  ((NEL1.GE.LBPT).AND.(NEL1.GE.LKABS)) THEN
C -- Index swapping (first):
         I1 = IOLDPS + 5 + NPIV
         K1 = IOLDPS + 5 + NASS - LBP
         DO 10 II=1,LBP
          ISWOP  = IW(I1)
          IW(I1) = IW(K1)
          IW(K1) = ISWOP
          I1     = I1 +1
          K1     = K1 + 1
  10     CONTINUE
C -- WE SWAPP THE REAL VALUE BLOCK
         IPOS = IPOSEL + NPIV*NFRONT
         KPOS = IPOSEL + (NASS-LBP)*NFRONT
         LBP1 = LBP * NFRONT
         CALL SSWAP(LBP1,A(IPOS),1,A(KPOS),1)
        ENDIF
  700   RETURN
        END
CCC
      SUBROUTINE MC51T(A,LA,NPIVB,NFRONT,
     *                             NPIV,NASS,POSELT)
      INTEGER NPIVB,NASS,LA
      REAL    A(LA)
      REAL    ALPHA, ONE
      INTEGER APOS, POSELT
      INTEGER NFRONT, NPIV, NASSL
      INTEGER LPOS, LPOS1, LPOS2, NEL1, NEL11, NPIVE
      PARAMETER (ONE=1.0E0)

      EXTERNAL STRSM, SGEMM

        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        NPIVE  = NPIV - NPIVB
        NASSL  = NASS - NPIVB
        APOS   = POSELT + NPIVB*NFRONT + NPIVB
        LPOS2  = APOS + NASSL
C*** WE HAVE FIRST TO COMPUTE THE U FACTORS
C    SO WE APPLY A TRIANGULAR SOLVER OF SIZE NPIVE WITH
C    NEL1 SECOND MEMBERS.
C BLAS-3 STRSM
C
        CALL STRSM('R','U','N','U',NEL1,NPIVE,ONE,A(APOS),NFRONT,
     *              A(LPOS2),NFRONT)
C******
C WE NOW UPDATE THE CORRESPONDING CONTRIBUTION BLOCK
C******
        LPOS   = LPOS2 + NFRONT*NPIVE
        LPOS1  = APOS + NFRONT*NPIVE
        ALPHA  = -1.0E0
        CALL SGEMM('N','N',NEL1,NEL11,NPIVE,ALPHA,A(LPOS2),
     *          NFRONT,A(LPOS1),NFRONT,ONE,A(LPOS),NFRONT)
  500   RETURN
        END
CCC
      SUBROUTINE MC51U(A, LA, IPTA, IPT,
     *       IPTEND, IPTBEG, LCKS)
C************************************************************
C We free Real space at position IPT
C
C INPUT:
C        IPT = position of real space to be freed
C
C************************************************************
      INTEGER NBUD
      PARAMETER (NBUD=29)
      INTEGER LA,IPTA(NBUD),IPTEND, IPTBEG
      INTEGER LCKS(20)
      REAL A(LA)
C****************
C LOCAL VARIABLES
C****************
      INTEGER IPT,IPTOLD,LSIZ,IFW,LSIZP1,IPT1,IPTBUD
      INTEGER IBW,IHEAD
      REAL ZERO
      INTRINSIC REAL, INT
      PARAMETER (ZERO=0.0E0)
C***
C Statement included in serial version to avoid compiler warning
C***
      LCKS(1) = 0
C****************************
C****************************
C FREE SPACE AT POINT IPT
C****************************
      IF (IPT.GE.IPTEND) THEN
C      STORAGE WAS FROM NON-BUDDY PART
       IPTOLD = IPTA(NBUD)
       IPTA(NBUD) = IPT
       A(IPT+1) = ZERO
       A(IPT+2) = REAL(IPTOLD)
       IF (IPTOLD.GT.0) A(IPTOLD+1) = REAL(IPT)
      ELSE
C ====================
C FREE SPACE FROM BUDDY
C =====================
  45  LSIZ = INT(A(IPT)+0.5)
C     free space in buddy of size lsiz
      IFW    = -1
      LSIZP1 = LSIZ + 1
      IPT1   = IPT - IPTBEG
      IF (((IPT1+2**LSIZP1-1)/2**LSIZP1)*2**LSIZP1.EQ.
     *      IPT1+2**LSIZP1-1) IFW = 1
      IPTBUD = IPT + IFW*2**LSIZ
C     search direction for amalgamation is IFW
C     CHECK TO SEE IF BUDDY IS FREE AND OF CORRECT SIZE
      IF (IPTBUD.LT.IPTEND .AND. IPTBUD.GT.IPTBEG .AND.
     *    INT(A(IPTBUD+1)+0.5).GE.0 .AND.
     *    INT(A(IPTBUD)+0.5).EQ.LSIZ) THEN
C       REMOVE BUDDY FROM ITS FREE LIST
        IFW = INT(A(IPTBUD+2)+0.5)
        IBW = INT(A(IPTBUD+1)+0.5)
C       We compress 2 success blocks of size LSIZ
        IF (IBW.EQ.0) IPTA(LSIZ) = IFW
        IF (IBW.GT.0) A(IBW+2)   = REAL(IFW)
        IF (IFW.GT.0) A(IFW+1)   = REAL(IBW)
C       AMALGAMATE WITH BUDDY
        IPT    = MIN0(IPT,IPTBUD)
        A(IPT) = REAL(LSIZ + 1)
        GO TO 45
      ELSE
C       PUT RESULTING FREED SPACE INTO POOL
        IHEAD = IPTA(LSIZ)
        IF (IHEAD.GT.0) A(IHEAD+1) = REAL(IPT)
        A(IPT+1)   = ZERO
        A(IPT+2)   = REAL(IHEAD)
        IPTA(LSIZ) = IPT
C       Finally we obtain a block of size
      ENDIF
      ENDIF
C****************************
C****************************
C
      RETURN
      END
CCC
C There follows the scaling routines
      SUBROUTINE MC51V(N,NZ,VAL,IRN,ICN,RNOR,
     *      COLSCA,ROWSCA,MPRINT)
C
C Compute diagonal scaling
C******************************
C  ORIGINAL MATRIX is in format
C   IRN(NZ), ICN(NZ), VAL(NZ) (see MA41A/AD for description)
C******************************
      INTEGER   N, NZ
      REAL      VAL(NZ),RNOR(N),COLSCA(N)
      REAL      ROWSCA(N)
      INTEGER   IRN(NZ),ICN(NZ)
      REAL      VDIAG
      INTEGER   MPRINT,I,J,K
      INTRINSIC SQRT
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
CC
C WE USE A DIAGONAL SCALING
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
      SUBROUTINE MC51W(N, NZ, VAL, ROWIND, COLIND,
     *                  RNOR, CNOR, WNOR, MPRINT, MP, NSCA)
C*****************
C MC29 based scaling: The scaling factors are chosen so that  the
C   scaled matrix has its nonzeros near unity in the sense that the sum
C   of the squares of the logarithms of the nonzeros is minimized.
C   MC29A computes RNOR(I), CNOR(I), I=1,N
C   and the scaled matrix B=(bij)
C   factors are given by the formula:
C     bij = exp(ri)*aij*exp(cj)
C   see Curtis and Reid, J. Inst. Maths. Applics. (1972), 10, pp 118-124
C   for a detailed description of the algorithm
C******************************
C INPUT ARRAYS:
C*************
C  ORIGINAL MATRIX FORMAT described in MA41A/AD
C    3 ARRAYS: ROWIND(NZ)  , COLIND(NZ) , VAL(NZ)
C IF NSCA=5 or NSCA =6 then the values of
C the input matrix (VAL) are scaled.
C***************
C OUTPUT ARRAYS:
C   CNOR(N) and RNOR(N) holding the scaling arrays
C***************
C WORKING ARRAYS
C  that need not be set on input
C   COLIND(NZ) will hold indices of the columns
C   WNOR(5*N)
C******************************
      INTEGER N, NZ
      REAL    VAL(NZ),RNOR(N),CNOR(N),
     *        WNOR(5*N)
      INTEGER COLIND(NZ),ROWIND(NZ)
      INTEGER J,I,K
      INTEGER MPRINT,MP,NSCA
      INTEGER IFAIL9
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
      EXTERNAL MC29A
CC
C WE USE SCALING BASED ON MC29
CC
      DO 15 I=1,N
       RNOR(I)   = ZERO
       CNOR(I)   = ZERO
  15  CONTINUE
      CALL MC29A(N,N,NZ,VAL,ROWIND,COLIND,
     *   RNOR,CNOR,WNOR, MP,IFAIL9)
C     IFAIL9 should always be zero since
C     N and NZ have already been checked.
C     IF (IFAIL9.NE.0) THEN
C       scaling was not performed
CCVD$ NODEPCHK
CCVD$ VECTOR
CCVD$ CONCUR
C       DO 20 I=1,N
C         RNOR(I)   = ONE
C         CNOR(I)   = ONE
C  20   CONTINUE
C       GOTO 500 --> (RETURN)
C      ENDIF
      DO 30 I=1,N
       CNOR(I) = EXP(CNOR(I))
       RNOR(I) = EXP(RNOR(I))
  30  CONTINUE
C*********************************************
C  NSCA = 5 or 6 then this scaling will
C       be followed by another one and we must
C       compute the scaled matrix
      IF ((NSCA.EQ.5).OR.(NSCA.EQ.6)) THEN
        DO 100 K=1,NZ
          I   = ROWIND(K)
          J   = COLIND(K)
          IF (MIN(I,J).LT.1 .OR. I.GT.N .OR. J.GT.N) GOTO 100
          VAL(K) = VAL(K) * CNOR(J) * RNOR(I)
 100    CONTINUE
      ENDIF
C
      IF (MPRINT.GE.0)
     *   WRITE(MPRINT,'(A)') 'End of scaling using MC29'
      RETURN
      END
CCC
      SUBROUTINE MC51X(N,NZ,IRN,ICN,VAL,
     *    RNOR,CNOR,COLSCA,ROWSCA,MPRINT)
C
C Compute Row and column scaling
C*******************************************************
C  ORIGINAL MATRIX is in format
C   IRN(NZ), ICN(NZ), VAL(NZ) (see MA41A/AD for description)
C*******************************************************
      INTEGER N, NZ
      REAL    VAL(NZ),RNOR(N),CNOR(N)
      REAL    COLSCA(N),ROWSCA(N)
      REAL    CMIN,CMAX,RMIN,ARNOR,ACNOR
      INTEGER IRN(NZ), ICN(NZ)
      REAL    VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
CC
C WE USE A ROW AND COLUMN SCALING
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
       WRITE(MPRINT,'(A,1PE12.4)') 'Maximum max-norm of columns:',CMAX
       WRITE(MPRINT,'(A,1PE12.4)') 'Minimum max-norm of columns:',CMIN
       WRITE(MPRINT,'(A,1PE12.4)') 'Minimum max-norm of rows   :',RMIN
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
      SUBROUTINE MC51Y(N,NZ,VAL,IRN,ICN,
     *       CNOR,COLSCA,MPRINT)
C
C Column scaling is computed
C*******************************************************
C  ORIGINAL MATRIX is in format
C   IRN(NZ), ICN(NZ), VAL(NZ) (see MA41A/AD for description)
C*******************************************************
C  WE COMPUTE THE MAXIMUM ELEMENT IN EACH COLUMN
C  STORED IN CNOR(N).
C  WE USE IT TO DIVIDE EACH COLUMN (J) OF THE ORIGINAL
C  MATRIX BY CNOR(J).
C  THIS IS EQUIVALENT TO MULTIPLY THE ORIGINAL MATRIX
C  BY A DIAGONAL MATRIX (ON THE RIGHT SIDE) CONTAINING
C  THE INVERSE OF THE MAXIMUM VALUE IN THE COLUMN.
C******************************
      INTEGER N,NZ
      REAL VAL(NZ),CNOR(N),COLSCA(N)
      INTEGER IRN(NZ), ICN(NZ)
      REAL VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
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
C            HARWELL SUBROUTINE LIBRARY Release 12 (1995)
C        --
C-             Copyright Rutherford Appleton Laboratory
C        --
C--------------------------------------------------------------------
C***********************************************************
C Subroutines managing the active area:
C   - allocation of space in active area: MC51Z
C   - liberation of space in active area : MC51U
C   - Move list of LU blocks left in active area to LU area:
C**********************************************************
      SUBROUTINE MC51Z(A, LA, IPTA, LREQ, IPT, LSIZ, IRES,
     *                  LENA,NBACTI,LCKS)
C*************************************************************
C We allocate Real space of length LREQ
C
C INPUT:
C        LREQ = length of space to be allocated
C OUTPUT:
C        IPT  is the pointer to the free space
C        2**LSIZ is the space really reserved if buddy system
C                was used otherwise LSIZ =0 on exit.
C        IRES =
C                 1 If fixed block was used to allocate space
C                 2 If buddy system was used to allocate space
C                -1 IF failure in allocating space
C                -2 If failure in allocating space and NBACTI=0
C
C************************************************************
      INTEGER NBUD
      PARAMETER (NBUD=29)
      INTEGER LA, IPTA(NBUD), IRES, NBACTI
      INTEGER LENA,LSIZ,LREQ,IPT
      INTEGER LCKS(20)
      INTEGER ISIZ,IFW,KSIZ,IPTBUD
      LOGICAL IOK
      REAL A(LA)
      REAL ZERO
      INTRINSIC REAL, INT
      PARAMETER (ZERO=0.0E0)
C***
C Statement included in serial version to avoid compiler warning
C***
      LCKS(1) = 0
      LSIZ = 0
      LSIZ = LOG(REAL(LREQ))/LOG(2.0) + 1
C***************************************
C DUM_LCK is in label 400
C**************************************
C OBTAIN SPACE OF SIZE LSIZ
      IOK = .FALSE.
      IF (LREQ.LE.LENA) THEN
C OBTAIN SPACE FROM NON-BUDDY SYSTEM
        IPT = IPTA(NBUD)
        IF (IPT.NE.0) THEN
          IOK  = .TRUE.
          IRES = 1
          IPTA(NBUD) = INT(A(IPT+2)+0.5)
          IF (IPTA(NBUD).NE.0) A(IPTA(NBUD)+1) = 0
        ENDIF
      ENDIF
      IF (IOK) GOTO 400
C USE BUDDY SYSTEM
      DO 300 ISIZ=LSIZ,NBUD-1
        IF (IPTA(ISIZ).NE.0) THEN
         IPT  = IPTA(ISIZ)
         IRES = 2
         IFW  = INT(A(IPT+2)+0.5)
         IPTA(ISIZ) = IFW
         IF (IFW.NE.0) A(IFW+1) = ZERO
         A(IPT+1) = REAL(-2)
         A(IPT)   = REAL(LSIZ)
         IF (ISIZ.EQ.LSIZ) GO TO 400
          DO 28 KSIZ = ISIZ-1,LSIZ,-1
            IPTBUD     = IPT + 2**KSIZ
C ADD BUDDY TO FREE LIST (WHICH IS EMPTY)
            IPTA(KSIZ) = IPTBUD
            A(IPTBUD)  = REAL(KSIZ)
            A(IPTBUD+1) = ZERO
            A(IPTBUD+2) = ZERO
 28       CONTINUE
          GOTO 400
        ENDIF
 300  CONTINUE
C     working space was not found in active area
C=================================================================
C No free space available
C Check if at least one frontal matrix is currently being factored
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
