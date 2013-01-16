SUBROUTINE StoreDataAscii(nx,ny,E,P,FineGrid,nr)
!
! ****************************************************************
!
!>
!! Store data in an ascii file.
!!
!!
!! By Allan P. Engsig-Karup.
!<
USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: nx, ny, i, j, k, nr, i0=100, io
REAL(KIND=long), DIMENSION(nx*ny) :: E,P
TYPE (Level_def) :: FineGrid
CHARACTER(len=30) :: filename,form
CHARACTER*3 CHOUT3
io=i0+nr
open(unit=io, file='fort.'//CHOUT3(INT(io)),status='unknown')

DO j=1,ny
   DO i=1,nx
      k=(j-1)*nx+i
      WRITE (io, 445) FineGrid%x(i,j), FineGrid%y(i,j), E(k), P(k)
   ENDDO
   WRITE(io,444)
ENDDO
CLOSE(io)
445 FORMAT(4e14.6)
446 FORMAT(100e14.6)
444 FORMAT()
END SUBROUTINE StoreDataAscii

    CHARACTER*3 FUNCTION CHOUT3(INTIN)
!------------------------------------------------------------------------
!
!     Description:  Utility to turn ints less than 1000 into chars.
!
!-----------------------------------------------------------------------
!
!     Parameters:
!
!     Name         Description
!     ------------------------------------------------------------------
!     INTIN        Input integer .LE. 999.
!     ------------------------------------------------------------------
!
!     Variables:
!
!     Name         Description
!     ------------------------------------------------------------------
!     CHOUT3       3 characters returned.
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER HOLDER(3), I, INTIN, INTTEM, FACTOR
      IF (INTIN .GT. 999) THEN
         PRINT *, 'TW-CHOUT3: Attempt to convert too large an '
         PRINT *, 'integer to character.  Only 3 places honored.'
      END IF
      INTTEM = INTIN
      DO  I = 3,1,-1
         FACTOR = 10**(I-1)
         HOLDER(I) = INTTEM / FACTOR 
         INTTEM = INTTEM - (HOLDER(I) * FACTOR)
      end DO
      CHOUT3 = CHAR(48+HOLDER(3))//CHAR(48+HOLDER(2))//CHAR(48+HOLDER(1))
      RETURN
end function
