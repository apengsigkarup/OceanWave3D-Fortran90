SUBROUTINE BeachGen3Dver2(h,h0,h1,lx,ly,x,y,Nx,Ny)
! By Allan P. Engsig-Karup.
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: Nx,Ny,i,j
REAL(KIND=long) :: h(Nx,Ny), h0, h1, l, x(Nx), y(Ny), xi, lx, ly

IF (lx /= ly .OR. ABS(lx-one)>1e-10) THEN
	PRINT*,'Error: Lx and Ly have to be equal to size 1 for the chosen test case.'
	STOP
ELSE
	l = lx
ENDIF

! Initialize:
h = zero
DO j = 1,Ny
	DO i = 1,Nx
        IF (SQRT(x(i)**2 + y(j)**2) > one) THEN
            xi = one
        ELSE
            xi = SQRT(x(i)**2 + y(j)**2)
        ENDIF
	    IF (xi > zero .AND. xi < l) THEN
	        h(i,j)   = h0 - (h0-h1)/two*(one+TANH(SIN(pi*(xi-l/two)/l)/(one-(two*(xi-l/two)/l)**2)))
		ELSE IF (xi>=l) THEN
			h(i,j) = h1
		ELSE
		    h(i,j)   = h0
		ENDIF
    END DO
END DO
END SUBROUTINE BeachGen3Dver2
