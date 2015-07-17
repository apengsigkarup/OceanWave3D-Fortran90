      PROGRAM MAIN
C
C  Solve the linear system A x = b, where A is tridiagonal with
C  superdiagonal values 1, subdiagonals -1 and diagonals 2, and where
C  the square roots of the inverse of the diagonal of A are used to
C  precondition from both sides
C
      INTEGER N, LDW, M, LDH, I, IACT, LOCY, LOCZ
      PARAMETER ( N = 10, M = 5, LDW = N, LDH = M + 1 )
      DOUBLE PRECISION TWO, ONE, THREE, SQRTWO, RESID
      PARAMETER ( TWO = 2.0D+0, ONE = 1.0D+0, THREE = 3.0D+0 )
      INTEGER ICNTL( 8 ), INFO( 4 )
      DOUBLE PRECISION CNTL( 4 ), W( LDW, M + 7 ), H( LDH, M + 2 )
      INTEGER ISAVE(17)
      DOUBLE PRECISION RSAVE(9)
      LOGICAL LSAVE(4)
      EXTERNAL MI24AD, MI24ID
      SQRTWO = SQRT( TWO )
      CALL MI24ID( ICNTL, CNTL, ISAVE, RSAVE, LSAVE )
C
C  Precondition from both sides
C
      ICNTL( 3 ) = 3
      ICNTL( 6 ) = 100
C
C  Set right hand side, b
C
      W( 1, 1 ) = THREE
      DO 10 I = 2, N - 1
         W( I, 1 ) = TWO
   10 CONTINUE
      W( N, 1 ) = ONE
C
C  Perform an iteration of the GMRES(m) method
C
      IACT = 0
   20 CONTINUE
      CALL MI24AD( IACT, N, M, W, LDW, LOCY, LOCZ, H, LDH, RESID,
     *             ICNTL, CNTL, INFO, ISAVE, RSAVE, LSAVE )
      IF ( IACT .LT. 0 ) THEN
         WRITE( 6, 2020 ) INFO( 1 )
         GO TO 60
      END IF
C
C  Perform the matrix-vector product
C
      IF ( IACT .EQ. 2 ) THEN
         W( 1, LOCY ) = TWO * W( 1, LOCZ ) + W( 2, LOCZ )
         DO 30 I = 2, N - 1
            W( I, LOCY ) = - W( I - 1, LOCZ ) + TWO * W( I, LOCZ ) +
     *                       W( I + 1, LOCZ )
   30    CONTINUE
         W( N, LOCY ) = - W( N - 1, LOCZ ) + TWO * W( N, LOCZ )
         GO TO 20
      END IF
C
C  Perform the left preconditioning operation
C
      IF ( IACT .EQ. 3 ) THEN
         DO 40 I = 1, N
            W( I, LOCY ) = W( I, LOCZ ) / SQRTWO
   40    CONTINUE
         GO TO 20
      END IF
C
C  Perform the right preconditioning operation
C
      IF ( IACT .EQ. 4 ) THEN
         DO 50 I = 1, N
            W( I, LOCY ) = W( I, LOCZ ) / SQRTWO
   50    CONTINUE
         GO TO 20
      END IF
C
C  Solution found
C
      WRITE( 6, 2000 ) INFO( 2 ), ( W( I, 2 ), I = 1, N )
      IF ( INFO( 1 ) .GT. 0 ) WRITE( 6, 2010 ) INFO( 1 )
   60 CONTINUE
      STOP
 2000 FORMAT( I6,' iterations required by MI24 ', // ' Solution = ',
     *       / ( 1P, 5D10.2 ) )
 2010 FORMAT( ' Warning: INFO( 1  ) = ',I2,' on exit ' )
 2020 FORMAT( ' Error return: INFO( 1  ) = ',I2,' on exit ' )
      END
