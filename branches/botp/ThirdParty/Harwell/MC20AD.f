C[[[[[[DATE  09 MAR 1989     COPYRIGHT UKAEA, HARWELL.
C[[[[[[ALIAS MC20AD MC20BD
      SUBROUTINE MC20AD(NC,MAXA,A,INUM,JPTR,JNUM,JDISP)
C
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER JDISP,MAXA,NC
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*)
      INTEGER INUM(*),JNUM(*),JPTR(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,JA,JB,JCE,JCEP,K,KR,LOC,NULL
C     ..
      NULL   = -JDISP
C**      CLEAR JPTR
      DO 60 J = 1,NC
   60 JPTR(J) = 0
C**      COUNT THE NUMBER OF ELEMENTS IN EACH COLUMN.
      DO 120 K = 1,MAXA
         J      = JNUM(K) + JDISP
         JPTR(J) = JPTR(J) + 1
  120 CONTINUE
C**      SET THE JPTR ARRAY
      K      = 1
      DO 150 J = 1,NC
         KR     = K + JPTR(J)
         JPTR(J) = K
  150 K      = KR
C
C**      REORDER THE ELEMENTS INTO COLUMN ORDER.  THE ALGORITHM IS AN
C        IN-PLACE SORT AND IS OF ORDER MAXA.
      DO 230 I = 1,MAXA
C        ESTABLISH THE CURRENT ENTRY.
         JCE    = JNUM(I) + JDISP
         IF (JCE.EQ.0) GO TO 230
         ACE    = A(I)
         ICE    = INUM(I)
C        CLEAR THE LOCATION VACATED.
         JNUM(I) = NULL
C        CHAIN FROM CURRENT ENTRY TO STORE ITEMS.
         DO 200 J = 1,MAXA
C        CURRENT ENTRY NOT IN CORRECT POSITION.  DETERMINE CORRECT
C        POSITION TO STORE ENTRY.
            LOC    = JPTR(JCE)
            JPTR(JCE) = JPTR(JCE) + 1
C        SAVE CONTENTS OF THAT LOCATION.
            ACEP   = A(LOC)
            ICEP   = INUM(LOC)
            JCEP   = JNUM(LOC)
C        STORE CURRENT ENTRY.
            A(LOC) = ACE
            INUM(LOC) = ICE
            JNUM(LOC) = NULL
C        CHECK IF NEXT CURRENT ENTRY NEEDS TO BE PROCESSED.
            IF (JCEP.EQ.NULL) GO TO 230
C        IT DOES.  COPY INTO CURRENT ENTRY.
            ACE    = ACEP
            ICE    = ICEP
  200    JCE    = JCEP + JDISP
C
  230 CONTINUE
C
C**      RESET JPTR VECTOR.
      JA     = 1
      DO 250 J = 1,NC
         JB     = JPTR(J)
         JPTR(J) = JA
  250 JA     = JB
      RETURN

      END
      SUBROUTINE MC20BD(NC,MAXA,A,INUM,JPTR)
C
C     ******************************************************************
C
C     .. Scalar Arguments ..
      INTEGER MAXA,NC
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*)
      INTEGER INUM(*),JPTR(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS
C     ..
      KMAX   = MAXA
      DO 30 JJ = 1,NC
         J      = NC + 1 - JJ
         KLO    = JPTR(J) + 1
         IF (KLO.GT.KMAX) GO TO 30
         KOR    = KMAX
         DO 25 KDUMMY = KLO,KMAX
C ITEMS KOR, KOR+1, .... ,KMAX ARE IN ORDER
            ACE    = A(KOR-1)
            ICE    = INUM(KOR-1)
            DO 10 K = KOR,KMAX
               IK     = INUM(K)
               IF (IABS(ICE).LE.IABS(IK)) GO TO 20
               INUM(K-1) = IK
   10       A(K-1) = A(K)
            K      = KMAX + 1
   20       INUM(K-1) = ICE
            A(K-1) = ACE
   25    KOR    = KOR - 1
C        NEXT COLUMN
   30 KMAX   = KLO - 2
      RETURN

      END
