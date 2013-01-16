SUBROUTINE timeintegration_fsbc
! By Allan P. Engsig-Karup.
USE GlobalVariables
IMPLICIT NONE
!
EXTERNAL rhsFreeSurface3D, rhsLinearFreeSurface3D
EXTERNAL rhsFreeSurface3D_SWENSE, rhsLinearFreeSurface3D_SWENSE

IF (LinearONOFF==0) THEN
   ! LINEAR COMPUTATION
   IF(swenseONOFF/=0) THEN
      SELECT CASE (timemethod)
          CASE (1) ! Classical RK4
              CALL Runge_Kutta_4(rhsLinearFreeSurface3D_SWENSE)
          CASE (2) ! Carpenter & Kennedy, low-storage RK45
              CALL Runge_Kutta_45(rhsLinearFreeSurface3D_SWENSE)
      END SELECT
   ELSE
      SELECT CASE (timemethod)
          CASE (1) ! Classical RK4
              CALL Runge_Kutta_4(rhsLinearFreeSurface3D)
          CASE (2) ! Carpenter & Kennedy, low-storage RK45
              CALL Runge_Kutta_45(rhsLinearFreeSurface3D)
      END SELECT
   ENDIF
ELSE
  ! NONLINEAR COMPUTATION
  IF(swenseONOFF/=0) THEN
     SELECT CASE (timemethod)
         CASE (1) ! Classical RK4
             CALL Runge_Kutta_4(rhsFreeSurface3D_SWENSE)
         CASE (2) ! Carpenter & Kennedy, low-storage RK45
             CALL Runge_Kutta_45(rhsFreeSurface3D_SWENSE)
     END SELECT
  ELSE
     SELECT CASE (timemethod)
         CASE (1) ! Classical RK4
             CALL Runge_Kutta_4(rhsFreeSurface3D)
         CASE (2) ! Carpenter & Kennedy, low-storage RK45
             CALL Runge_Kutta_45(rhsFreeSurface3D)
     END SELECT
  ENDIF
ENDIF

END SUBROUTINE timeintegration_fsbc
