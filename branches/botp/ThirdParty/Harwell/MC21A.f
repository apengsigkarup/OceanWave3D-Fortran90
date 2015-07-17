C[[[[[[ALIAS   MC21A
C[[[[[[DATE   09 MAR 1989     COPYRIGHT UKAEA, HARWELL.
      SUBROUTINE MC21A(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
C     ..
C     .. External Subroutines ..
      EXTERNAL MC21B
C     ..
      CALL MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),IW(1,3),
     +           IW(1,4))
      RETURN

      END
      SUBROUTINE MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
C   PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.
C IT IS USED AS A WORK ARRAY IN THE SORTING ALGORITHM.
C   ELEMENTS (IPERM(I),I) I=1, ... N  ARE NON-ZERO AT THE END OF THE
C ALGORITHM UNLESS N ASSIGNMENTS HAVE NOT BEEN MADE.  IN WHICH CASE
C (IPERM(I),I) WILL BE ZERO FOR N-NUMNZ ENTRIES.
C   CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I
C WAS VISITED.
C   ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT.
C   OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE MAIN LOOP.
C
C   INITIALIZATION OF ARRAYS.
C     .. Scalar Arguments ..
      INTEGER LICN,N,NUMNZ
C     ..
C     .. Array Arguments ..
      INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
C     ..
      DO 10 I = 1,N
         ARP(I) = LENR(I) - 1
         CV(I)  = 0
   10 IPERM(I) = 0
      NUMNZ  = 0
C
C
C   MAIN LOOP.
C   EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT
C OR GIVES A ROW WITH NO ASSIGNMENT.
      DO 130 JORD = 1,N
         J      = JORD
         PR(J)  = -1
         DO 100 K = 1,JORD
C LOOK FOR A CHEAP ASSIGNMENT
            IN1    = ARP(J)
            IF (IN1.LT.0) GO TO 60
            IN2    = IP(J) + LENR(J) - 1
            IN1    = IN2 - IN1
            DO 50 II = IN1,IN2
               I      = ICN(II)
               IF (IPERM(I).EQ.0) GO TO 110
   50       CONTINUE
C   NO CHEAP ASSIGNMENT IN ROW.
            ARP(J) = -1
C   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J.
   60       OUT(J) = LENR(J) - 1
C INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS.
            DO 90 KK = 1,JORD
               IN1    = OUT(J)
               IF (IN1.LT.0) GO TO 80
               IN2    = IP(J) + LENR(J) - 1
               IN1    = IN2 - IN1
C FORWARD SCAN.
               DO 70 II = IN1,IN2
                  I      = ICN(II)
                  IF (CV(I).EQ.JORD) GO TO 70
C   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS.
                  J1     = J
                  J      = IPERM(I)
                  CV(I)  = JORD
                  PR(J)  = J1
                  OUT(J1) = IN2 - II - 1
                  GO TO 100

   70          CONTINUE
C
C   BACKTRACKING STEP.
   80          J      = PR(J)
               IF (J.EQ.-1) GO TO 130
   90       CONTINUE
C
  100    CONTINUE
C
C   NEW ASSIGNMENT IS MADE.
  110    IPERM(I) = J
         ARP(J) = IN2 - II - 1
         NUMNZ  = NUMNZ + 1
         DO 120 K = 1,JORD
            J      = PR(J)
            IF (J.EQ.-1) GO TO 130
            II     = IP(J) + LENR(J) - OUT(J) - 2
            I      = ICN(II)
            IPERM(I) = J
  120    CONTINUE
C
  130 CONTINUE
C
C   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE
C PERMUTATION IPERM.
      IF (NUMNZ.EQ.N) RETURN
      DO 140 I = 1,N
  140 ARP(I) = 0
      K      = 0
      DO 160 I = 1,N
         IF (IPERM(I).NE.0) GO TO 150
         K      = K + 1
         OUT(K) = I
         GO TO 160

  150    J      = IPERM(I)
         ARP(J) = I
  160 CONTINUE
      K      = 0
      DO 170 I = 1,N
         IF (ARP(I).NE.0) GO TO 170
         K      = K + 1
         IOUTK  = OUT(K)
         IPERM(IOUTK) = I
  170 CONTINUE
      RETURN

      END
