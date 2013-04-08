SUBROUTINE FILTERING_NEW(Nx,Ny,Wavefield,filterNP,filterALPHA,filtercoefficients,oddeven, swenseONOFF,&
     filtercoefficients2, GhostGridX, GhostGridY)

USE Precision
USE Constants
USE DataTypes
IMPLICIT NONE
INTEGER :: Nx, Ny, i, j,k, filterALPHA, filterNP, oddeven, swenseONOFF
REAL(KIND=long) :: tmp(Nx,Ny), filtercoefficients(filterNP), idx(filterNP)
TYPE(wavefield_FS) :: Wavefield
! GD: test...
REAL(long) :: filtercoefficients2(filterNP,filterALPHA)
INTEGER :: GhostGridX, GhostGridY

IF (swenseONOFF==0) THEN
! Filtering of E
  DO k = 1 , 2
    SELECT CASE (MOD(oddeven+k,2))
      CASE (0)
       IF (Nx>1) THEN
        tmp = Wavefield%E
        ! INTERIOR
        DO j = 1+GhostGridY, Ny-GhostGridY
          DO i = 1+filterALPHA+GhostGridX, Nx-filterALPHA-GhostGridX
              Wavefield%E(i,j) = DOT_PRODUCT(tmp(i-filterALPHA:i+filterALPHA,j),filtercoefficients)
          END DO
        END DO
        ! GD test - change 10-13-2009
        ! LEFT
        DO j = 1+GhostGridY, Ny-GhostGridY
            DO i = 1,filterALPHA!-GhostGridX !keep up to filterALPHA to replace eventually first interior points which uses ghost value...
              Wavefield%E(i+GhostGridX,j) = DOT_PRODUCT(tmp(1+GhostGridX:filterNP+GhostGridX,j),&
                filtercoefficients2(1:filterNP,i))
          END DO
        END DO
        ! RIGHT
        DO j = 1+GhostGridY, Ny-GhostGridY
            DO i = 1,filterALPHA !1+GhostGridX,filterALPHA
              Wavefield%E(Nx-GhostGridX-filterALPHA+i,j) = DOT_PRODUCT(tmp(Nx-GhostGridX-(filterNP)+1:Nx-GhostGridX,j),&
                   filtercoefficients2(filterNP:1:-1,filterALPHA+1-i))
          END DO
        END DO
      ENDIF
      CASE (1)
        IF (Ny>1) THEN
          tmp = Wavefield%E
          ! INTERIOR
             DO j = 1+filterALPHA, Ny-filterALPHA
                 DO i = 1+GhostGridX, Nx-GhostGridX
                     Wavefield%E(i,j) = DOT_PRODUCT(tmp(i,j-filterALPHA:j+filterALPHA),filtercoefficients)
                 END DO
             END DO
          !
          ! GD test - change 10-13-2009
          ! LEFT
          DO j = 1,filterALPHA!-GhostGridY !keep up to filterALPHA to replace eventually first interior points which uses ghost value...
              DO i = 1+GhostGridX, Nx-GhostGridX
                Wavefield%E(i,j+GhostGridY) = DOT_PRODUCT(tmp(i,1+GhostGridY:filterNP+GhostGridY),&
                    filtercoefficients2(1:filterNP,j))
            END DO
          END DO
          ! RIGHT
          DO j = 1,filterALPHA!1+GhostGridY,filterALPHA
              DO i = 1+GhostGridX, Nx-GhostGridX
                Wavefield%E(i,j+Ny-filterALPHA-GhostGridY) = DOT_PRODUCT(tmp(i,Ny-(filterNP)+1-1:Ny-1),&
                     filtercoefficients2(filterNP:1:-1,filterALPHA+1-j))
            END DO
          END DO
        ENDIF
    END SELECT
  END DO
  ! Filtering of P
  DO k = 1 , 2
    SELECT CASE (MOD(oddeven+k,2))
      CASE (0)
       IF (Nx>1) THEN
         tmp = Wavefield%P
         ! INTERIOR
         DO j = 1+GhostGridY, Ny-GhostGridY
           DO i = 1+filterALPHA+GhostGridX, Nx-filterALPHA-GhostGridX
               Wavefield%P(i,j) = DOT_PRODUCT(tmp(i-filterALPHA:i+filterALPHA,j),filtercoefficients)
           END DO
         END DO
         ! GD test - change 10-13-2009
         ! LEFT
         DO j = 1+GhostGridY, Ny-GhostGridY
             DO i = 1,filterALPHA!-GhostGridX
               Wavefield%P(i+GhostGridX,j) = DOT_PRODUCT(tmp(1+1:filterNP+1,j), &
                 filtercoefficients2(1:filterNP,i))
           END DO
         END DO
         ! RIGHT
         DO j = 1+GhostGridY, Ny-GhostGridY
             DO i = 1,filterALPHA!1+GhostGridX,filterALPHA
               Wavefield%P(i+Nx-filterALPHA-GhostGridX,j) = DOT_PRODUCT(tmp(Nx-(filterNP)+1-1:Nx-1,j),&
                    filtercoefficients2(filterNP:1:-1,filterALPHA+1-i))
           END DO
         END DO
       ENDIF
      CASE (1)
        IF (Ny>1) THEN
          tmp = Wavefield%P
          ! INTERIOR
          DO j = 1+filterALPHA, Ny-filterALPHA
              DO i = 1+GhostGridX, Nx-GhostGridX
                  Wavefield%P(i,j) = DOT_PRODUCT(tmp(i,j-filterALPHA:j+filterALPHA),filtercoefficients)
              END DO
          END DO
          ! GD test - change 10-13-2009
          ! LEFT
          DO j = 1,filterALPHA!-GhostGridY
              DO i = 1+GhostGridX, Nx-GhostGridX
                Wavefield%P(i,j+GhostGridY) = DOT_PRODUCT(tmp(i,1+1:filterNP+1),&
                    filtercoefficients2(1:filterNP,j))
            END DO
          END DO
          ! RIGHT
          DO j = 1,filterALPHA!1+GhostGridY,filterALPHA
              DO i = 1+GhostGridX, Nx-GhostGridX
                Wavefield%P(i,j+Ny-filterALPHA-GhostGridY) = DOT_PRODUCT(tmp(i,Ny-(filterNP)+1-1:Ny-1),&
                     filtercoefficients2(filterNP:1:-1,filterALPHA+1-j))
            END DO
          END DO
        ENDIF
    END SELECT
  END DO
ELSE
  ! Filtering of E + E_I
  ! FIXME: it does not seem to work (in the treatment of boundary points!)
  DO k = 1 , 2
    SELECT CASE (MOD(oddeven+k,2))
      CASE (0)
       IF (Nx>1) THEN
        tmp = Wavefield%E+Wavefield%E_I
        ! INTERIOR
        DO j = 1+GhostGridY, Ny-GhostGridY
          DO i = 1+filterALPHA, Nx-filterALPHA
              Wavefield%E(i,j) = DOT_PRODUCT(tmp(i-filterALPHA:i+filterALPHA,j),filtercoefficients)-Wavefield%E_I(i,j)
          END DO
        END DO
        ! GD test
        ! LEFT
        DO j = 1+GhostGridY, Ny-GhostGridY
            DO i = 1,filterALPHA!-GhostGridX !keep up to filterALPHA to replace eventually first interior points which uses ghost value...
              Wavefield%E(i+GhostGridX,j) = DOT_PRODUCT(tmp(1+GhostGridX:filterNP+GhostGridX,j),&
                filtercoefficients2(1:filterNP,i))-Wavefield%E_I(i+GhostGridX,j)
          END DO
        END DO
        ! RIGHT
        DO j = 1+GhostGridY, Ny-GhostGridY
            DO i = 1,filterALPHA !1+GhostGridX,filterALPHA
              Wavefield%E(i+Nx-filterALPHA-GhostGridX,j) = DOT_PRODUCT(tmp(Nx-(filterNP)+1-GhostGridX:Nx-GhostGridX,j),&
                   filtercoefficients2(filterNP:1:-1,filterALPHA+1-i))-Wavefield%E_I(i+Nx-filterALPHA-GhostGridX,j)
          END DO
        END DO
        ! END GD test
       ENDIF
      CASE (1)
        IF (Ny>1) THEN
          tmp = Wavefield%E+Wavefield%E_I
          ! INTERIOR
          DO j = 1+filterALPHA, Ny-filterALPHA
              DO i = 1+GhostGridX, Nx-GhostGridX
                  Wavefield%E(i,j) = DOT_PRODUCT(tmp(i,j-filterALPHA:j+filterALPHA),filtercoefficients)&
                       -Wavefield%E_I(i,j)
              END DO
          END DO
          ! GD test
          ! LEFT
          DO j = 1,filterALPHA!-GhostGridY !keep up to filterALPHA to replace eventually first interior points which uses ghost value...
              DO i = 1+GhostGridX, Nx-GhostGridX
                Wavefield%E(i,j+GhostGridY) = DOT_PRODUCT(tmp(i,1+GhostGridY:filterNP+GhostGridY),&
                    filtercoefficients2(1:filterNP,j))-Wavefield%E_I(i,j+GhostGridY)
            END DO
          END DO
          ! RIGHT
          DO j = 1,filterALPHA!1+GhostGridY,filterALPHA
              DO i = 1+GhostGridX, Nx-GhostGridX
                Wavefield%E(i,j+Ny-filterALPHA-GhostGridY) = DOT_PRODUCT(tmp(i,Ny-(filterNP)+1-GhostGridY:Ny-GhostGridY),&
                     filtercoefficients2(filterNP:1:-1,filterALPHA+1-j))-Wavefield%E_I(i,j+Ny-filterALPHA-GhostGridY)
            END DO
          END DO
          ! END GD test
        ENDIF
    END SELECT
  END DO
  ! Filtering of P + P_I
  DO k = 1 , 2
    SELECT CASE (MOD(oddeven+k,2))
      CASE (0)
       IF (Nx>1) THEN
        tmp = Wavefield%P+Wavefield%P_I_s
        ! INTERIOR
        DO j = 1+GhostGridY, Ny-GhostGridY
          DO i = 1+filterALPHA, Nx-filterALPHA
              Wavefield%P(i,j) = DOT_PRODUCT(tmp(i-filterALPHA:i+filterALPHA,j),filtercoefficients)&
                   -Wavefield%P_I_s(i,j)
          END DO
        END DO
        ! GD test
        ! LEFT
        DO j = 1+GhostGridY, Ny-GhostGridY
            DO i = 1,filterALPHA!-GhostGridX
              Wavefield%P(i+GhostGridX,j) = DOT_PRODUCT(tmp(1+GhostGridX:filterNP+GhostGridX,j), &
                filtercoefficients2(1:filterNP,i))-Wavefield%P_I_s(i+GhostGridX,j)
          END DO
        END DO
        ! RIGHT
        DO j = 1+GhostGridY, Ny-GhostGridY
            DO i = 1,filterALPHA!1+GhostGridX,filterALPHA
              Wavefield%P(i+Nx-filterALPHA-GhostGridX,j) = DOT_PRODUCT(tmp(Nx-(filterNP)+1-GhostGridX:Nx-GhostGridX,j),&
                   filtercoefficients2(filterNP:1:-1,filterALPHA+1-i))-Wavefield%P_I_s(i+Nx-filterALPHA-GhostGridX,j)
          END DO
        END DO
        ! END GD test
       ENDIF
      CASE (1)
        IF (Ny>1) THEN
          tmp = Wavefield%P+Wavefield%P_I_s
          ! INTERIOR
          DO j = 1+filterALPHA, Ny-filterALPHA
              DO i = 1+GhostGridX, Nx-GhostGridX
                  Wavefield%P(i,j) = DOT_PRODUCT(tmp(i,j-filterALPHA:j+filterALPHA),filtercoefficients)-&
                       Wavefield%P_I_s(i,j)
              END DO
          END DO
          ! GD test
          ! LEFT
          DO j = 1,filterALPHA!-GhostGridY
              DO i = 1+GhostGridX, Nx-GhostGridX
                Wavefield%P(i,j+GhostGridY) = DOT_PRODUCT(tmp(i,1+GhostGridY:filterNP+GhostGridY),&
                    filtercoefficients2(1:filterNP,j))-Wavefield%P_I_s(i,j+GhostGridY)
            END DO
          END DO
          ! RIGHT
          DO j = 1,filterALPHA!1+GhostGridY,filterALPHA
              DO i = 1+GhostGridX, Nx-GhostGridX
                Wavefield%P(i,j+Ny-filterALPHA-GhostGridY) = DOT_PRODUCT(tmp(i,Ny-(filterNP)+1-GhostGridY:Ny-GhostGridY),&
                     filtercoefficients2(filterNP:1:-1,filterALPHA+1-j))-Wavefield%P_I_s(i,j+Ny-filterALPHA-GhostGridY)
            END DO
          END DO
          ! END GD test
        ENDIF
    END SELECT
  END DO
ENDIF
!PRINT*,'  SG-FILTER applied.'
END SUBROUTINE FILTERING_NEW
