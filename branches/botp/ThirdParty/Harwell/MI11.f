C COPYRIGHT (c) 1995 Council for the Central Laboratory
*                    of the Research Councils
C######DATE 30 November 1995
      SUBROUTINE MI11ID(ICNTL,CNTL)

C ICNTL(1)   Stream number for error messages. Default 6.
C ICNTL(2)   Stream number for warnings. Default 6.
C ICNTL(3)   If nonzero, no checking of input data. Default 0.
C ICNTL(4)   If nonzero, matrix is input by row indices
C            and column pointers. Otherwise, input is by column
C            indices and row pointers. Default 0.
C
C CNTL(1),CNTL(2)    Threshold tolerances.
C           CNTL(1) has default  0.0001.
C           CNTL(2) has default  sq. root machine precision.
C
C     .. Array Arguments ..
      DOUBLE PRECISION CNTL(2)
      INTEGER ICNTL(4)
C     ..
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = 0
      ICNTL(4) = 0
      CNTL(1) = 0.0001D0
      CNTL(2) = 1.0D0

      RETURN
      END

      SUBROUTINE MI11AD(N,NZ,INDEX,A,IP,JCN,FACTOR,IFPTR,IW,ICNTL,CNTL,
     +                  INFO)
C
C
C     This subroutine does incomplete Gaussian elimination
C     on a sparse matrix. The matrix is stored by row pointers
C     and column indices.
C     note: no pivoting is done
C           and the factorization is incomplete,
C           i.e., only where there exists a storage location
C           will the operation take place.
C
C Argument list (*) denotes argument changed by routine).
C
C  N     INTEGER variable. Must be set to the order of the matrix.
C
C  NZ    INTEGER variable. Must be set to the length of the arrays
C        A and INDEX. Must be at least IP(N+1)-1.
C
C  INDEX INTEGER array of length NZ. If ICNTL(4) = 0, the first
C        IP(N+1)-1 entries must be set by the user to hold the
C        column indices of the entries of A. The entries must
C        be ordered by rows, with the entries in each row
C        contiguous and those of row I preceding those of row
C        I+1 (I =1, 2,...,N). The ordering within each row
C        is unimportant. If ICNTL(4) is nonzero,
C        INDEX must hold row indices of the entries of A.
C
C  A     REAL (DOUBLE PRECISION) array of length NZ.
C        Must hold the entries of the matrix, corresponding to
C        the column indices set in INDEX.
C
C  IP    INTEGER array of length  N+1.
C        Must hold pointers to the start of the rows (columns
C        if ICNTL(4) is nonzero) in the arrays
C        A and INDEX.
C
C *JCN   INTEGER array length NZ.
C        On exit, JCN holds the column indices
C        of the matrix L+U (L strict lower triangle and U upper
C        triangle)  ordered by rows, with the entries
C        in each row ordered by columns.
C
C *FACTOR  REAL (DOUBLE PRECISION) array length NZ.
C        On exit,  contains the entries of the matrix L+U
C
C *IFPTR INTEGER array of length  N+1.
C        On exit, holds pointers to the start of the rows in the arrays
C        JCN and FACTOR.
C
C *IW   INTEGER array of length 5N+NZ. Used as workspace.
C
C ICNTL,CNTL control parameters (see MI11I/ID)
C
C *INFO INTEGER array of length 6. Used to hold information.
c       INFO(1) is error flag.
C
C     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER N,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NZ),CNTL(2),FACTOR(NZ)
      INTEGER ICNTL(4),IFPTR(N+1),INDEX(NZ),INFO(6),IP(N+1),IW(5*N+NZ),
     +        JCN(NZ)
C     ..
C     .. Local Scalars ..
      INTEGER I,IDUP,IOUT,IPERM,IROW,IWORK,JROW,K,KE,KS,LP,MP,NE,NZA,
     +        NZJ,NZNEW
C     ..
      DOUBLE PRECISION FD05AD
      EXTERNAL FD05AD

C     .. External Subroutines ..
      EXTERNAL MC38AD,MI11CD
C     ..
      LP = ICNTL(1)
      MP = ICNTL(2)
      DO 10 I = 1,6
         INFO(I) = 0
   10 CONTINUE
C
C Reset CNTL(2) if necessary
      IF (CNTL(2).LT.SQRT(FD05AD(1))) CNTL(2) = ONE

      IF (N.LT.1) THEN
         INFO(1) = -1
      ELSE IF (NZ.LT.IP(N+1)-1 .OR. NZ.LT.N) THEN
         INFO(1) = -2
      END IF
      IF (INFO(1).LT.0) THEN
         IF (LP.GT.0) WRITE (LP,FMT=9000) INFO(1)
         RETURN
      END IF
C
C If ICNTL(4) is nonzero, use MC38 to reorder by rows (column indices
C JCN and row pointers IFPTR).
C Otherwise, copy matrix into (JCN,FACTOR,IFPTR)
C
      IF (ICNTL(4).EQ.0) THEN
         DO 20 I = 1,IP(N+1) - 1
            JCN(I) = INDEX(I)
            FACTOR(I) = A(I)
   20    CONTINUE
         DO 30 I = 1,N + 1
            IFPTR(I) = IP(I)
   30    CONTINUE
      ELSE
         NE = IP(N+1) - 1
         CALL MC38AD(N,N,NE,IP,INDEX,.TRUE.,A,IFPTR,JCN,FACTOR)
      END IF

C Check the user-supplied data for errors.
C No further checking is done if ICNTL(3) is nonzero.
      IOUT = 0
      IDUP = 0
      NZA = 0
      IF (ICNTL(3).NE.0) THEN
         NZNEW = IFPTR(N+1) - 1
         INFO(5) = NZNEW
         GO TO 70
      END IF

C Check for out-of-range entries and duplicates (out-of-range
C entries are ignored; duplicates are summed).
C Initialise
      DO 40 I = 1,N
         IW(I) = 0
   40 CONTINUE

      NZNEW = 0
      KS = 1
      NZJ = 0
C Loop over the rows.
      DO 60 JROW = 1,N
         KE = IFPTR(JROW+1) - 1
         IFPTR(JROW+1) = IFPTR(JROW)
C Loop over the entries in row JROW
         DO 50 K = KS,KE
            I = JCN(K)
C Check index is within range.
            IF (I.GT.N .OR. I.LT.1) THEN
               IOUT = IOUT + 1
               GO TO 50
            END IF
C Check the value of the entry is nonzero
            IF (FACTOR(K).EQ.ZERO) THEN
               NZA = NZA + 1
               GO TO 50
            END IF
            IF (IW(I).LE.NZJ) THEN
               NZNEW = NZNEW + 1
               JCN(NZNEW) = I
               FACTOR(NZNEW) = FACTOR(K)
               IFPTR(JROW+1) = IFPTR(JROW+1) + 1
               IW(I) = NZNEW
            ELSE
C We have a duplicate in row J
               IDUP = IDUP + 1
               FACTOR(IW(I)) = FACTOR(IW(I)) + FACTOR(K)
            END IF
   50    CONTINUE
         KS = KE + 1
         NZJ = NZNEW
   60 CONTINUE

C NZNEW is the number of entries after the elimination of out-of-range
C entries and duplicates and entries with value zero.
      INFO(5) = NZNEW

      IF (IOUT.GT.0) THEN
         INFO(1) = 1
         INFO(2) = IOUT
C Immediate return if ALL entries out-of-range.
         IF (IOUT.EQ.NZ) THEN
            IF (MP.GT.0) THEN
               WRITE (MP,FMT=9010) INFO(1)
               WRITE (MP,FMT=9030) INFO(2)
            END IF
            RETURN
         END IF
      END IF

      IF (IDUP.GT.0) THEN
         INFO(1) = 1
         INFO(3) = IDUP
      END IF

      IF (NZA.GT.0) THEN
         INFO(1) = 1
         INFO(4) = NZA
C Immediate return if ALL entries have value zero.
         IF (NZA.EQ.NZ) THEN
            IF (MP.GT.0) THEN
               WRITE (MP,FMT=9010) INFO(1)
               WRITE (MP,FMT=9050) INFO(4)
            END IF
            RETURN
         END IF
      END IF
C
   70 CONTINUE
C
C Immediate return if N=1
      IF (N.EQ.1) THEN
         IW(1) = 1
         IW(2) = 1
         RETURN
      END IF
C
C Divide up workspace
      IPERM = 1
      IROW = IPERM + N
      IWORK = IROW + N
      CALL MI11CD(N,NZNEW,JCN,FACTOR,IFPTR,IW(IPERM),IW(IROW),IW(IWORK),

     +            ICNTL,CNTL,INFO)

      IF (INFO(1).EQ.-3) THEN
C Matrix is structurally singular
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9020) IW(1)
         END IF
         RETURN
      END IF

      IF (INFO(6).GT.0) INFO(1) = 1
C Printing of warnings
      IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
         WRITE (MP,FMT=9010) INFO(1)
         IF (IOUT.GT.0) WRITE (MP,FMT=9030) INFO(2)
         IF (IDUP.GT.0) WRITE (MP,FMT=9040) INFO(3)
         IF (NZA.GT.0) WRITE (MP,FMT=9050) INFO(4)
         IF (INFO(6).GT.0) WRITE (MP,FMT=9060) INFO(6)
      END IF

      RETURN
 9000 FORMAT (/' Error message from MI11A/AD. INFO(1) = ',I4)
 9010 FORMAT (/' Warning message from MI11A/AD. INFO(1) = ',I4)
 9020 FORMAT (/' The permuted matrix has ',I6,' zeros on the diagonal.')

 9030 FORMAT (/I6,' indices are out-of-range.')
 9040 FORMAT (/I6,' duplicated entries were entered.')
 9050 FORMAT (/I6,' matrix entries with value zero were entered.')
 9060 FORMAT (/I6,' diagonal entries were increased during the ',
     +       'elimination process.')
      END


      SUBROUTINE MI11BD(TRANS,N,NZ,JCN,FACTOR,IFPTR,IW,Y,Z,W)

C This routine performs preconditioning operations using the
C incomplete LU factorisation as the preconditioner.
C     P = (LU)(-1)Q
C
C i.e perform y= Pz  and y = P(T)z
C
C  TRANS    LOGICAL variable.  If TRANS=.TRUE. the
C           preconditioning operation y = P(T)z is performed
C           and if TRANS=.FALSE.  the preconditioning operation
C           y = Pz is performed.
C
C N,NZ,JCN,FACTOR,IFPTR must all be unchanged
C           since the call to MI11A/AD. These arguments
C           are not altered by the routine.
C
C  IW       first 2N entries must be unchanged since call to MI11A/AD.
C           IW(I), I=1,..,N holds permutation
C           IW(N+I), I=1,...,N points to the diagonal in row I.
C
C *Y        REAL (DOUBLE PRECISION) array of length N. On exit,
C           Y holds the preconditioned vector y.
C
C  Z        REAL (DOUBLE PRECISION) array of length N.
C           On entry, must hold the vector z which is to
C           be preconditioned.
C
C W         REAL (DOUBLE PRECISION) array of length N.
C           Used as workspace  (Not accessed if TRANS = FALSE)
C
C
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER N,NZ
      LOGICAL TRANS
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION FACTOR(NZ),W(*),Y(N),Z(N)
      INTEGER IFPTR(N+1),IW(2*N),JCN(NZ)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION T
      INTEGER I,IE,IR,IS,J
C     ..
      IF (TRANS) GO TO 50

C  We have to solve
C                   LUy = Qz
C
C U is unit diagonal

C First solve  Lx = Qz for x (store x in y)
      DO 20 I = 1,N
         T = ZERO
         IS = IFPTR(I)
         IE = IW(N+I) - 1
         DO 10 IR = IS,IE
            T = T + FACTOR(IR)*Y(JCN(IR))
   10    CONTINUE
         Y(I) = (Z(IW(I))-T)/FACTOR(IW(N+I))
   20 CONTINUE

C Solve Uy = x for y (overwrite x by y)
      DO 40 I = N - 1,1,-1
         T = ZERO
         IS = IW(N+I) + 1
         IE = IFPTR(I+1) - 1
         DO 30 IR = IS,IE
            T = T + FACTOR(IR)*Y(JCN(IR))
   30    CONTINUE
         Y(I) = Y(I) - T
   40 CONTINUE

      RETURN

   50 CONTINUE
C
C  We have to solve
C                   y = Q(T)LU(-T)z
C
C       i.e. LU(T)w = z, Qy = w

C  Solve U(T)v = z for v then L(T)v = w for w, then y = Q(T)w

C Copy Z into W
      DO 60 I = 1,N
         W(I) = Z(I)
   60 CONTINUE

      DO 80 I = 1,N - 1
         IS = IW(N+I) + 1
         IE = IFPTR(I+1) - 1
         DO 70 IR = IS,IE
            J = JCN(IR)
            W(J) = W(J) - W(I)*FACTOR(IR)
   70    CONTINUE
   80 CONTINUE

      DO 100 I = N,1,-1
         W(I) = W(I)/FACTOR(IW(N+I))
         IS = IFPTR(I)
         IE = IW(N+I) - 1
         DO 90 IR = IS,IE
            J = JCN(IR)
            W(J) = W(J) - W(I)*FACTOR(IR)
   90    CONTINUE
  100 CONTINUE

C Need  permutation y= Q(T)w
      DO 110 I = 1,N
         Y(IW(I)) = W(I)
  110 CONTINUE

      RETURN
      END

      SUBROUTINE MI11CD(N,NZ,JCN,FACTOR,IFPTR,IPERM,IROW,IW,ICNTL,CNTL,
     +                  INFO)

C *IPERM INTEGER array of length N. On exit, IPERM(I)
C        is the position in the original matrix A of row I in the
C        permuted matrix.
C
C *IROW  INTEGER array of length N. On exit, IROW(I)
C        is the position of the diagonal in row I.
C
C     .. Scalar Arguments ..
      INTEGER N,NZ
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CNTL(2),FACTOR(NZ)
      INTEGER ICNTL(4),IFPTR(N+1),INFO(6),IPERM(N),IROW(N),IW(3*N+NZ),
     +        JCN(NZ)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AMAXC,AMAXR,ATEMP
      INTEGER I,IE,IK,IS,ITEMP,J,JROW,K,KDIAG,KDUMMY,KE,KJ,KLO,KMAX,KOR,

     +        KS,L,N1,N2,NUMNZ,INFO6
C     ..
C     .. External Subroutines ..
      EXTERNAL MC21AD,MC22AD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MAX,MIN
C     ..
C     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C     ..
C
C Look to see if the matrix has any zeros on the diagonal.
C If a zero is found, use MC21 to find a row permutation that
C makes the matrix nonzeros on the diagonal.
      DO 30 K = 1,N
C Find the diagonal in row K
         KS = IFPTR(K)
         KE = IFPTR(K+1) - 1
         DO 10 L = KS,KE
            IF (JCN(L).EQ.K .AND. ICNTL(3).EQ.0) GO TO 30
            IF (JCN(L).EQ.K) GO TO 20
   10    CONTINUE
C No diagonal entry
         GO TO 50
C Diagonal entry found
   20    KDIAG = L
C Check size of diagonal entry (it can only be zero
C if the user did not allow input data to be checked for errors).
         IF (FACTOR(KDIAG).EQ.ZERO) GO TO 50
   30 CONTINUE

C Matrix has no zeros on the diagonal. No need for MC21/MC22.
      DO 40 I = 1,N
         IPERM(I) = I
   40 CONTINUE
      GO TO 90

   50 CONTINUE

C Set IROW(I) to the number of nonzeros in row I.

      DO 60 I = 1,N
         IROW(I) = IFPTR(I+1) - IFPTR(I)
   60 CONTINUE

      CALL MC21AD(N,JCN,NZ,IFPTR,IROW,IPERM,NUMNZ,IW)
      IF (NUMNZ.LT.N) THEN
C Matrix is structurally singular
         INFO(1) = -3
         IPERM(1) = NUMNZ
         RETURN
      END IF
C
C Permute the matrix to QA using MC22.
C Set the column permutation to the identity.
      DO 70 I = 1,N
         IW(I) = I
   70 CONTINUE

      N1 = N + 1
      N2 = N1 + 2*N
      CALL MC22AD(N,JCN,FACTOR,NZ,IROW,IPERM,IW,IW(N1),IW(N2))

C Set IFPTR for permuted matrix
      DO 80 I = 1,N
         IFPTR(I+1) = IFPTR(I) + IROW(I)
   80 CONTINUE

   90 CONTINUE
C
C Order the entries in each row of QA by columns
C (This bit of code is taken from MC20B/BD)
      KMAX = NZ
C Loop over rows
      DO 140 JROW = 1,N
         J = N + 1 - JROW
         KLO = IFPTR(J) + 1
         IF (KLO.GT.KMAX) GO TO 130
         KOR = KMAX
         DO 120 KDUMMY = KLO,KMAX
C Items KOR, KOR+1, .... ,KMAX are in order
            ATEMP = FACTOR(KOR-1)
            ITEMP = JCN(KOR-1)
            DO 100 K = KOR,KMAX
               IK = JCN(K)
               IF (ITEMP.LE.IK) GO TO 110
               JCN(K-1) = IK
               FACTOR(K-1) = FACTOR(K)
  100       CONTINUE
            K = KMAX + 1
  110       JCN(K-1) = ITEMP
            FACTOR(K-1) = ATEMP
            KOR = KOR - 1
  120    CONTINUE
C Next row
  130    KMAX = KLO - 2
  140 CONTINUE

C Find the diagonal in each row
      DO 160 K = 1,N
         KS = IFPTR(K)
         KE = IFPTR(K+1) - 1
         DO 150 L = KS,KE
            IF (JCN(L).EQ.K) THEN
               IROW(K) = L
               GO TO 160
            END IF
  150    CONTINUE
  160 CONTINUE

C Now perform LU factorisation  ... we now have a matrix
C QA which has nonzeros on the diagonal

      DO 260 K = 1,N
         KDIAG = IROW(K)
C Find the largest entry in row K (above diagonal)
C Jump if the diagonal is the last entry in row K
C (but increase it in size first, if necessary).
         IF (KDIAG.EQ.KE) THEN
            IF (ABS(FACTOR(KDIAG)).LT.CNTL(2)) THEN
               FACTOR(KDIAG) = CNTL(2)
               INFO(6) = INFO(6) + 1
            END IF
            GO TO 260
         END IF

         AMAXR = ZERO
         KS = IROW(K) + 1
         KE = IFPTR(K+1) - 1
         DO 170 L = KS,KE
            AMAXR = MAX(AMAXR,ABS(FACTOR(L)))
  170    CONTINUE

C Check size of diagonal entry. If it is too small compared
C with the maximum entry in the row, look to see if
C it is too small compared with the maximum entry in
C the column. If it is too small,
C increase it to the size of the largest entry in the row
C or the largest entry in the column (depending on which
C is the smaller of the two).
         INFO6 = INFO(6)
         IF (ABS(FACTOR(KDIAG)).LT.CNTL(1)*AMAXR) THEN
C Search below diagonal in column K
            AMAXC = ZERO
            DO 190 I = K + 1,N
C Look for entry in position row I, column K
               IS = IFPTR(I)
               IE = IROW(I) - 1
               DO 180 L = IS,IE
                  IF (K.EQ.JCN(L)) THEN
                     AMAXC = MAX(AMAXC,ABS(FACTOR(L)))
                     GO TO 190
                  ELSE IF (K.LT.JCN(L)) THEN
                     GO TO 190
                  END IF
  180          CONTINUE
  190       CONTINUE
            IF (ABS(FACTOR(KDIAG)).LT.CNTL(1)*AMAXC) THEN
               FACTOR(KDIAG) = MIN(AMAXC,AMAXR)
               INFO6 = INFO6 + 1
            END IF
         END IF

         IF (ABS(FACTOR(KDIAG)).LT.CNTL(1)) THEN
            FACTOR(KDIAG) = CNTL(2)
            IF (INFO6.EQ.INFO(6)) INFO6 = INFO6 + 1
         END IF
         INFO(6) = INFO6

         DO 200 L = KS,KE
            FACTOR(L) = FACTOR(L)/FACTOR(KDIAG)
  200    CONTINUE

C Loop over rows K+1 to N
         DO 250 I = K + 1,N
C Look for entry in position row I, column K
C (I>K so we are in lower triangular part)
            IS = IFPTR(I)
            IE = IROW(I) - 1
            DO 210 L = IS,IE
               IF (K.EQ.JCN(L)) THEN
C Entry located
                  IK = L
                  GO TO 220
               ELSE IF (K.LT.JCN(L)) THEN
                  GO TO 250
               END IF
  210       CONTINUE
            GO TO 250

  220       CONTINUE
            KJ = KDIAG + 1
C Update row I
            IE = IFPTR(I+1) - 1
            DO 240 J = IK + 1,IE
  230          CONTINUE
               IF (JCN(KJ).GT.JCN(J)) GO TO 240
               IF (JCN(KJ).LT.JCN(J)) THEN
                  KJ = KJ + 1
                  IF (KJ.GT.KE) GO TO 250
                  GO TO 230
               ELSE
                  FACTOR(J) = FACTOR(J) - FACTOR(IK)*FACTOR(KJ)
               END IF
  240       CONTINUE
  250    CONTINUE
  260 CONTINUE
      RETURN

      END
