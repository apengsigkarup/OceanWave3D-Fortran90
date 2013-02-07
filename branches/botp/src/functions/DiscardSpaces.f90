SUBROUTINE DiscardSpaces(str)
! By Allan P. Engsig-Karup.
IMPLICIT NONE
CHARACTER :: str*(*)
INTEGER :: i
DO i = 1, LEN_TRIM(str)
	IF (str(i:i)==' ') THEN
		str(i:i) = '0'
	ENDIF
END DO
END SUBROUTINE DiscardSpaces
