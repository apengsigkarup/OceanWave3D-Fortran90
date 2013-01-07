* COPYRIGHT (c) 1988 AEA Technology
*######DATE 21 Jan 1993
C       Toolpack tool decs employed.
C       SAVE statement added.
C 1/10/98 RC(3) not initialized to avoid SUN f90 failure
C 16 October 2001: STOP and WRITE statements removed.

      REAL FUNCTION FD05A( INUM )
      INTEGER INUM
      REAL RC( 5 )
C
C----------------------------------------------------------------
C  Real constants for: IEEE single precision (4-byte arithmetic)
C
C  Obtained from H.S.L. subroutine ZE02AM.
C  Nick Gould and Sid Marlow, Harwell Laboratory, April 1988.
C----------------------------------------------------------------
C
C  OBTAINED FROM H.S.L. SUBROUTINE ZE02AM.
C  NICK GOULD AND SID MARLOW, HARWELL, JULY 1988.
C
C  RC(1) THE 'SMALLEST' POSITIVE NUMBER: 1 + RC(1) > 1.
C  RC(2) THE 'SMALLEST' POSITIVE NUMBER: 1 - RC(2) < 1.
C  RC(3) THE SMALLEST NONZERO +VE REAL NUMBER.
C  RC(4) THE SMALLEST FULL PRECISION +VE REAL NUMBER.
C  RC(5) THE LARGEST FINITE +VE REAL NUMBER.
C
      SAVE RC
      DATA RC(1)/1.19210E-07/
      DATA RC(2)/5.96047E-08/
C     DATA RC(3)/1.40131E-45/
      DATA RC(4)/1.17550E-38/
      DATA RC(5)/3.40281E+38/

      IF ( INUM .LE. 0 ) THEN
         FD05A = RC( 1 )
      ELSE IF ( INUM .GE. 6 ) THEN
         FD05A = RC( 5 )
      ELSE IF ( INUM .EQ. 3 ) THEN
         FD05A = RC(4)/2.0**23
      ELSE
         FD05A = RC( INUM )
      ENDIF
      RETURN
      END
