SUBROUTINE FILTERING_SWENSE(Nx,Ny,U,U_I,filterNP,filterALPHA,filtercoefficients,oddeven,GhostGridX,GhostGridY)

USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: Nx, Ny, i, j,k, filterALPHA, filterNP, oddeven, GhostGridX,GhostGridY
REAL(KIND=long) :: U(Nx,Ny), U_I(Nx,Ny), tmp(Nx,Ny), filtercoefficients(filterNP), idx(filterNP)
! test...
REAL(KIND=long) :: filtercoeff_bis(filterNP)


DO k = 1 , 2
  SELECT CASE (MOD(oddeven+k,2))
	CASE (0)
     IF (Nx>1) THEN
	  tmp = U+U_I
	  ! INTERIOR
	  DO j = 1, Ny
		DO i = 1+filterALPHA, Nx-filterALPHA
			U(i,j) = DOT_PRODUCT(tmp(i-filterALPHA:i+filterALPHA,j),filtercoefficients)-U_I(i,j)
		END DO
	  END DO
!$$$$$$       !
!$$$$$$       ! GD : modification, assuming symetry with respect to the boundary...
!$$$$$$       !
!$$$$$$       ! Left
!$$$$$$       DO j = 1, Ny
!$$$$$$         DO i = 1+GhostGridX,1+filterALPHA !keep the 1+filteralpha to prevent from using Ghost point?
!$$$$$$             !U(i,j) = DOT_PRODUCT(tmp(i-filterALPHA:i+filterALPHA,j),filtercoefficients)-U_I(i,j)
!$$$$$$             U(i,j) = DOT_PRODUCT(tmp(filterALPHA+i:2*i-GhostGridX:-1,j),filtercoefficients(1:filterALPHA-(i-GhostGridX)+1)) !points reproduced
!$$$$$$             U(i,j) = U(i,j) + DOT_PRODUCT(tmp(1+GhostGridX:i+filterALPHA,j),filtercoefficients(filterALPHA-(i-GhostGridX)+2:filterNP)) ! points in the domain
!$$$$$$             U(i,j) = U(i,j)-U_I(i,j)
!$$$$$$         END DO
!$$$$$$       END DO
!$$$$$$       ! Right
!$$$$$$       DO j = 1, Ny
!$$$$$$         DO i = Nx-filterALPHA-GhostGridX+1, Nx-GhostGridX !keep the Nx-filterALPHA-GhostGridX to prevent from using Ghost point?
!$$$$$$             U(i,j) = DOT_PRODUCT(tmp(i-filterALPHA:Nx-GhostGridX,j),filtercoefficients(1:Nx-GhostGridX-(i-filterALPHA)+1)) !points in the domain
!$$$$$$             U(i,j) = U(i,j) + DOT_PRODUCT(tmp(Nx-GhostGridX-1:Nx-GhostGridX-1-(filterALPHA-(Nx-i-GhostGridX)-1):-1,j),filtercoefficients(Nx-GhostGridX-(i-filterALPHA)+2:filterNP)) ! points reproduced
!$$$$$$             U(i,j) = U(i,j)-U_I(i,j)
!$$$$$$         END DO
!$$$$$$       END DO

!$$$$$$        !Left
!$$$$$$        DO j = 1+GhostGridY, Ny-GhostGridY
!$$$$$$          DO i = 1+GhostGridX,1+filterALPHA !keep the 1+filteralpha to prevent from using Ghost point?
!$$$$$$            filtercoeff_bis(1:filterNP) = filtercoefficients(1:filterNP)
!$$$$$$            filtercoeff_bis(filterNP:filterNP-(filterALPHA-(i-GhostGridX)+1-1):-1) = filtercoeff_bis(filterNP:filterNP-(filterALPHA-(i-GhostGridX)+1-1):-1) &
!$$$$$$                + filtercoefficients(1:filterALPHA-(i-GhostGridX)+1)
!$$$$$$            ! filtercoeff_bis starts at filterALPHA-(i-GhostGridX)+2
!$$$$$$            U(i,j) = DOT_PRODUCT(tmp(1+GhostGridX:i+filterALPHA,j),filtercoeff_bis(filterALPHA-(i-GhostGridX)+2:filterNP))-U_I(i,j)
!$$$$$$            print*,i,i+filterALPHA-(1+GhostGridX)+1,filterNP-(filterALPHA-(i-GhostGridX)+2)+1
!$$$$$$            print*,i,filterNP-(filterNP-(filterALPHA-(i-GhostGridX)+1-1))+1,filterALPHA-(i-GhostGridX)+1
!$$$$$$          END DO
!$$$$$$        END DO
!$$$$$$        pause
!$$$$$$        ! Right
!$$$$$$         DO j = 1+GhostGridY, Ny-GhostGridY
!$$$$$$           DO i = Nx-filterALPHA-GhostGridX+1, Nx-GhostGridX !keep the Nx-filterALPHA-GhostGridX to prevent from using Ghost point?
!$$$$$$             filtercoeff_bis(1:filterNP) = filtercoefficients(1:filterNP)
!$$$$$$             filtercoeff_bis(1:filterNP-(Nx-GhostGridX-(i-filterALPHA)+2)+1) = filtercoeff_bis(1:filterNP-(Nx-GhostGridX-(i-filterALPHA)+2)+1) &
!$$$$$$                 + filtercoefficients(filterNP:Nx-GhostGridX-(i-filterALPHA)+2:-1)
!$$$$$$             U(i,j) = DOT_PRODUCT(tmp(i-filterALPHA:Nx-GhostGridX,j),filtercoeff_bis(1:Nx-GhostGridX-(i-filterALPHA)+1))-U_I(i,j)
!$$$$$$           END DO
!$$$$$$         END DO
	 ENDIF
	CASE (1)
	  IF (Ny>1) THEN
		tmp = U+U_I
		! INTERIOR
		DO j = 1+filterALPHA, Ny-filterALPHA
			DO i = 1, Nx
				U(i,j) = DOT_PRODUCT(tmp(i,j-filterALPHA:j+filterALPHA),filtercoefficients)-U_I(i,j)
			END DO
		END DO
!$$$$$$         !
!$$$$$$         ! GD : modification, assuming symetry with respect to the boundary...
!$$$$$$         !
!$$$$$$         ! Left
!$$$$$$         DO j = 1+GhostGridY, 1+filterALPHA
!$$$$$$           DO i = 1,Nx !keep the 1+filteralpha to prevent from using Ghost point?
!$$$$$$               U(i,j) = DOT_PRODUCT(tmp(i,filterALPHA+j:2*j-GhostGridY:-1),filtercoefficients(1:filterALPHA-(j-GhostGridY)+1)) !points reproduced
!$$$$$$               U(i,j) = U(i,j) + DOT_PRODUCT(tmp(i,1+GhostGridY:j+filterALPHA),filtercoefficients(filterALPHA-(j-GhostGridY)+2:filterNP)) ! points in the domain
!$$$$$$               U(i,j) = U(i,j)-U_I(i,j)
!$$$$$$           END DO
!$$$$$$         END DO
!$$$$$$         ! Right
!$$$$$$         DO j = Ny-filterALPHA-GhostGridY+1, Ny-GhostGridY !keep the Nx-filterALPHA-GhostGridX to prevent from using Ghost point?
!$$$$$$           DO i = 1, Nx
!$$$$$$               U(i,j) = DOT_PRODUCT(tmp(i,j-filterALPHA:Ny-GhostGridY),filtercoefficients(1:Ny-GhostGridY-(j-filterALPHA)+1)) !points in the domain
!$$$$$$               U(i,j) = U(i,j) + DOT_PRODUCT(tmp(i,Ny-GhostGridY-1:Ny-GhostGridY-1-(filterALPHA-(Ny-j-GhostGridY)-1):-1),filtercoefficients(Ny-GhostGridY-(j-filterALPHA)+2:filterNP)) ! points reproduced
!$$$$$$               U(i,j) = U(i,j)-U_I(i,j)
!$$$$$$           END DO
!$$$$$$         END DO
	  ENDIF
  END SELECT
END DO
PRINT*,'  SG-FILTER applied.'
END SUBROUTINE
