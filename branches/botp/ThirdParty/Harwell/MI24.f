C *******************************************************************
C COPYRIGHT (c) 1995 Council for the Central Laboratory
*                    of the Research Councils
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
C######DATE 29 March 2001
C  March 2001: threadsafe version of MI04
      SUBROUTINE MI24ID( ICNTL, CNTL, ISAVE, RSAVE, LSAVE )
      INTEGER          ICNTL( 8 )
      DOUBLE PRECISION CNTL( 4 )
      INTEGER ISAVE(17)
      DOUBLE PRECISION RSAVE(9)
      LOGICAL LSAVE(4)
      INTEGER I
      DOUBLE PRECISION FD05AD
      EXTERNAL FD05AD
      INTRINSIC SQRT
      ICNTL( 1 ) = 6
      ICNTL( 2 ) = 6
      ICNTL( 3 ) = 0
      ICNTL( 4 ) = 0
      ICNTL( 5 ) = 0
      ICNTL( 6 ) = - 1
      ICNTL(7) = 0
      ICNTL(8) = 0
      CNTL( 1 ) = SQRT( FD05AD( 1 ) )
      CNTL( 2 ) = 0.0
      CNTL(3) = 0.0
      CNTL(4) = 0.0
      DO 10 I = 1, 17
      ISAVE(I) = 0
   10 CONTINUE
      DO 20 I = 1, 9
      RSAVE(I) = 0.0
   20 CONTINUE
      DO 30 I = 1, 4
      LSAVE(I) = .FALSE.
   30 CONTINUE
      RETURN
      END
      SUBROUTINE MI24AD( IACT, N, M, W, LDW, LOCY, LOCZ, H, LDH, RESID,
     *                   ICNTL, CNTL, INFO, ISAVE, RSAVE, LSAVE )
      DOUBLE PRECISION RESID
      INTEGER          IACT, N, M, LDW, LOCY, LOCZ, LDH
      DOUBLE PRECISION CNTL( 4 ), W( LDW, M + 7 ), H( LDH, M + 2 )
      INTEGER          ICNTL( 8 ), INFO( 4 )
      INTEGER ISAVE(17)
      DOUBLE PRECISION RSAVE(9)
      LOGICAL LSAVE(4)
      DOUBLE PRECISION    ZERO, ONE, POINT1
      PARAMETER         ( ZERO = 0.0D+0, ONE = 1.0D+0, POINT1 = 1.0D-1 )
      INTEGER             I, K, ITMAX, CS, SN, R, S, V, UU, Y, RES
      INTEGER             B, X, IPOS, U
      LOGICAL             LEFT, RIGHT
      DOUBLE PRECISION    AA, BB, BNRM2 , RNORM , DDOT, DNRM2, RSTOP
      DOUBLE PRECISION    PRESID, PRSTOP, FD05AD
      EXTERNAL            DAXPY, DCOPY, DDOT, DNRM2, DROT, DROTG, DSCAL
      EXTERNAL            DTRSV, DGEMV, FD05AD
      IPOS   = ISAVE(1)
      ITMAX  = ISAVE(2)
      B      = ISAVE(3)
      I      = ISAVE(4)
      K      = ISAVE(5)
      R      = ISAVE(6)
      X      = ISAVE(7)
      U      = ISAVE(8)
      V      = ISAVE(9)
      S      = ISAVE(10)
      Y      = ISAVE(11)
      CS     = ISAVE(12)
      SN     = ISAVE(13)
      UU     = ISAVE(14)
      RES    = ISAVE(15)
      BNRM2  = RSAVE(1)
      AA     = RSAVE(2)
      BB     = RSAVE(3)
      RNORM  = RSAVE(4)
      PRESID = RSAVE(5)
      RSTOP  = RSAVE(6)
      PRSTOP = RSAVE(7)
      LEFT   = LSAVE(1)
      RIGHT  = LSAVE(2)
      IF ( IACT .NE. 0 ) THEN
         IF ( IACT .LT. 0 ) GO TO 1000
         IF ( IACT .EQ. 1 .AND. ICNTL( 4 ) .EQ. 0 ) GO TO 1000
         IF ( IACT .EQ. 1 .AND. BNRM2 .EQ. ZERO ) GO TO 1000
         GO TO ( 40, 60, 70, 100, 110, 120, 160 ), IPOS
      END IF
      INFO( 1 ) = 0
      IF ( N .LT. 1 ) THEN
         INFO( 1 ) = - 1
      ELSE IF ( M .LT. 1 ) THEN
         INFO( 1 ) = - 2
      ELSE IF ( LDW .LT. MAX( 1, N ) ) THEN
         INFO( 1 ) = - 3
      ELSE IF ( LDH .LT. M + 1 ) THEN
         INFO( 1 ) = - 4
      ENDIF
      IF ( INFO( 1 ) .LT. 0 ) THEN
         IACT = - 1
         IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
         GO TO 1000
      END IF
      INFO( 2 ) = 0
      IF ( ICNTL( 6 ) .GT. 0 ) THEN
         ITMAX = ICNTL( 6 )
      ELSE
         ITMAX = 2 * N
      END IF
      RES = 1
      X   = 2
      S   = 3
      B   = 4
      UU  = 5
      Y   = 6
      V   = 7
      CS = M + 1
      SN = CS + 1
      BNRM2 = DNRM2( N, W( 1, RES ), 1 )
      IF ( BNRM2. EQ. ZERO ) THEN
         IACT = 1
         DO 10 I = 1, N
            W( I, X ) = ZERO
            W( I, RES ) = ZERO
   10    CONTINUE
         RESID = ZERO
         GO TO 1000
      END IF
      IF ( ICNTL( 4 ) .EQ. 0 ) THEN
         IF ( CNTL( 1 ) .LT. FD05AD( 1 ) .OR. CNTL( 1 ) .GT. ONE ) THEN
            INFO( 1 ) = 1
            IF (ICNTL( 2 ) .GT. 0 ) THEN
               WRITE( ICNTL( 2 ), 2010 ) INFO( 1 )
               WRITE( ICNTL( 2 ), 2020 )
            END IF
            CNTL( 1 ) = SQRT( FD05AD( 1 ) )
         END IF
      END IF
      LEFT  = ICNTL( 3 ) .EQ. 1 .OR. ICNTL( 3 ) .EQ. 3
      RIGHT = ICNTL( 3 ) .EQ. 2 .OR. ICNTL( 3 ) .EQ. 3
      IF ( ICNTL( 5 ) .EQ. 0 ) THEN
         DO 20 I = 1, N
            W( I, X ) = ZERO
   20    CONTINUE
      END IF
      CALL DCOPY( N, W( 1, RES ), 1, W( 1, B ), 1 )
      IF ( DNRM2( N, W( 1, X ), 1 ) .EQ. ZERO ) GO TO 50
   30 CONTINUE
         IPOS = 1
         IACT = 2
         LOCY = Y
         LOCZ = X
         GO TO 1000
   40    CONTINUE
         CALL DAXPY( N, - ONE, W( 1, Y ), 1, W( 1, RES ), 1 )
   50    CONTINUE
         RESID = DNRM2( N, W( 1, RES ), 1 )
         IF ( INFO( 2 ) .EQ. 0 ) THEN
            RSTOP  = MAX( RESID * CNTL( 1 ), CNTL( 2 ) )
            PRSTOP = RSTOP
         END IF
         IF ( INFO( 1 ) .LT. 0 ) THEN
            IACT = - 1
            IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
            GO TO 1000
         END IF
         IF ( LEFT ) THEN
            R = UU
            IPOS = 2
            IACT = 3
            LOCY = R
            LOCZ = RES
            GO TO 1000
         ELSE
            R = RES
         END IF
   60    CONTINUE
         IF ( ICNTL( 4 ) .NE. 0 .OR. ( ICNTL( 4 ) .EQ. 0 .AND.
     *        RESID .LE. RSTOP ) ) THEN
            IACT = 1
            IPOS = 3
            GO TO 1000
         END IF
   70    CONTINUE
         CALL DCOPY( N, W( 1, R ), 1, W( 1, V ), 1 )
         RNORM = DNRM2( N, W( 1, V ), 1 )
         CALL DSCAL( N, ONE / RNORM, W( 1, V ), 1 )
         W( 1, S ) = RNORM
         DO 80 K = 2, N
            W( K, S ) = ZERO
   80    CONTINUE
         I = 0
   90    CONTINUE
            I = I + 1
            INFO( 2 ) = INFO( 2 ) + 1
            IF ( INFO( 2 ) .GT. ITMAX ) THEN
               I = I - 1
               INFO( 1 ) = - 5
               IF ( ICNTL( 1 ) .GT. 0 ) WRITE( ICNTL( 1 ), 2030 ) ITMAX
               IF ( I .NE. 0 ) GO TO 150
               IACT = - 1
               WRITE( ICNTL( 1 ), 2000 ) INFO( 1 )
               GO TO 1000
            END IF
            IF ( RIGHT ) THEN
               IPOS = 4
               IACT = 4
               LOCY = Y
               LOCZ = V + I - 1
               GO TO 1000
            END IF
  100       CONTINUE
            IPOS = 5
            IACT = 2
            IF ( RIGHT ) THEN
               LOCY = RES
               LOCZ = Y
            ELSE
               LOCY = RES
               LOCZ = V + I - 1
            END IF
            GO TO 1000
  110       CONTINUE
            IF ( LEFT ) THEN
               U = UU
               IPOS = 6
               IACT = 3
               LOCY = UU
               LOCZ = RES
               GO TO 1000
            ELSE
               U = RES
            END IF
  120       CONTINUE
            DO 130 K = 1, I
               H( K, I ) = DDOT( N, W( 1, U ), 1,
     *                              W( 1, V + K - 1 ), 1 )
               CALL DAXPY( N, - H( K, I ), W( 1, V + K - 1 ), 1,
     *                                     W( 1, U ), 1 )
  130       CONTINUE
            H( I + 1, I ) = DNRM2( N, W( 1, U ), 1 )
            CALL DCOPY( N, W( 1, U ), 1, W( 1, V + I ), 1 )
            CALL DSCAL( N, ONE / H( I + 1, I ), W( 1, V + I ), 1 )
            DO 140 K = 1, I - 1
               CALL DROT( 1, H( K, I ), LDH, H( K + 1, I ), LDH,
     *                    H( K, CS ), H( K, SN ) )
  140       CONTINUE
            AA = H( I, I )
            BB = H( I + 1, I )
            CALL DROTG( AA, BB, H( I, CS ), H( I, SN ) )
            CALL DROT( 1, H( I, I ), LDH, H( I + 1, I ), LDH,
     *            H( I, CS ), H( I, SN ) )
            IF ( I .LT. N ) THEN
               CALL DROT( 1, W( I, S ), LDW, W( I + 1, S ),
     *                    LDW, H( I, CS ), H( I, SN ) )
               PRESID = ABS( W( I + 1, S ) )
               IF ( PRESID .LE. PRSTOP .AND. ICNTL( 4 ) .EQ. 0 ) THEN
                  PRSTOP = PRSTOP * POINT1
                  GO TO 150
               END IF
               IF ( I .LT. M )  GO TO 90
            END IF
  150    CONTINUE
         CALL DCOPY( I, W( 1, S ), 1, W( 1, Y ), 1 )
         CALL DTRSV( 'UPPER', 'NOTRANS', 'NONUNIT', I, H, LDH,
     *               W( 1, Y ), 1 )
         CALL DGEMV( 'NOTRANS', N, I, ONE, W( 1, V ), LDW, W( 1, Y ), 1,
     *               ZERO, W( 1, UU ), 1 )
         IF ( RIGHT ) THEN
            IPOS = 7
            IACT = 4
            LOCY = Y
            LOCZ = UU
            GO TO 1000
         END IF
  160    CONTINUE
         IF ( RIGHT ) THEN
            CALL DAXPY( N, ONE, W( 1, Y  ), 1, W( 1, X ), 1 )
         ELSE
            CALL DAXPY( N, ONE, W( 1, UU ), 1, W( 1, X ), 1 )
         END IF
         CALL DCOPY( N, W( 1, B ), 1, W( 1, RES ), 1 )
      GO TO 30
 1000 CONTINUE
      ISAVE(1) = IPOS
      ISAVE(2) = ITMAX
      ISAVE(3) = B
      ISAVE(4) = I
      ISAVE(5) = K
      ISAVE(6) = R
      ISAVE(7) = X
      ISAVE(8) = U
      ISAVE(9) = V
      ISAVE(10) = S
      ISAVE(11) = Y
      ISAVE(12) = CS
      ISAVE(13) = SN
      ISAVE(14) = UU
      ISAVE(15) = RES
      RSAVE(1) = BNRM2
      RSAVE(2) = AA
      RSAVE(3) = BB
      RSAVE(4) = RNORM
      RSAVE(5) = PRESID
      RSAVE(6) = RSTOP
      RSAVE(7) = PRSTOP
      LSAVE(1) = LEFT
      LSAVE(2) = RIGHT
      RETURN
 2000 FORMAT( / ' Error message from MI24A/AD. INFO(1) = ', I4 )
 2010 FORMAT( / ' Warning message from MI24A/AD. INFO(1) = ', I4 )
 2020 FORMAT( ' Convergence tolerance out of range.' )
 2030 FORMAT( /, ' # iterations required exceeds the maximum of ',
     *        I8, ' allowed by ICNTL(6)' )
      END
